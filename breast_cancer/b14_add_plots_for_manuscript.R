
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged >=30

# Content: additional plots for manuscripts

# Date: 12/9/2025

# Last update: N/A

# Author: Yiyi Zhou

################################################################################


# crude annual rates of breast cancer by HIV status, age, and race/ethnicity----

library(polars)
library(exactci)
library(ggplot2)

data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/crude_ir/"


# function to return yearly data
year_strat <- function(data, year_lower, year_upper) {
  
  lower_date <- as.Date(paste0(year_lower, "-01-01"))
  upper_date <- as.Date(paste0(year_upper, "-01-01"))
  
  year.dat <- data$filter(
    pl$col("baseline")$dt$year() < year_upper &
      pl$col("first_event_date")$dt$year() >= year_lower
  )$with_columns(
    
    baseline0 = pl$col("baseline"),
    first_event_date0 = pl$col("first_event_date")
  )$with_columns(
    
    baseline = pl$when(pl$col("baseline")$gt(lower_date))$
      then(pl$col("baseline"))$
      otherwise(lower_date),
    
    index = (pl$col("first_event_date")$dt$year() >= year_upper)$cast(pl$Int8)
  )$with_columns(
    first_event_date = pl$when(pl$col("index") == 1)$
      then(upper_date)$
      otherwise(pl$col("first_event_date")),
    
    event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      cast(pl$Int8),
    
    compete_event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("compete_event_flag"))$
      cast(pl$Int8)
  )$with_columns(
    event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
  )$with_columns(
    event_time = pl$when(pl$col("event_time") == 0)$
      then(pl$col("event_time") + 0.000001)$
      otherwise(pl$col("event_time"))
  )$with_columns(
    event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil(),
    year_cat = pl$lit(paste0("[", year_lower, ", ", year_upper, ")"))
  )$drop("index")
  
  if (year.dat$filter(pl$col("event_time") <= 0)$height > 0) {
    cat("samp.year.", year_lower, ".", year_upper, "has event_time <= 0!")
  }
  
  # year.dat$write_ipc(paste0(
  #   data_dir, "samp.year.", year_lower, ".", year_upper, ".w",".arrow"
  # ))
  return(year.dat)
}

# function to put year-stratified data into a list
create_year_stratified_base <- function(data, year_cats) {
  year_strat_list <- list()
  
  for (year_cat in year_cats) {
    year_lower <- year_cat[1]
    year_upper <- year_cat[2]
    year_strat_list[[as.character(year_lower)]] <- year_strat(data, year_lower, year_upper)
  }
  
  return(year_strat_list)
}

# function to estimate crude incidence rates
calc_crude_ir <- function(ir_dat) {
  
  ir_dat <- ir_dat$with_columns(
    followup_years = ((pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()$
                        cast(pl$Float64) / 365.25)
  )$with_columns(
    followup_years = pl$when(pl$col("followup_years") == 0)$
      then(0.000001)$
      otherwise(pl$col("followup_years"))
  )$select(
    total_cases = pl$sum("event_flag")$cast(pl$Int64),
    total_person_years = pl$sum("followup_years")$cast(pl$Float64)
  )$with_columns(
    crude_incidence = (pl$col("total_cases") / pl$col("total_person_years"))$cast(pl$Float64)
  )
  
  return(ir_dat)
}

# function to calculate overall crude IR over time
calculate_overall_metrics <- function(year_stratified_data) {
  ir_year <- list()
  
  for (year in names(year_stratified_data)) {
    year.dat <- year_stratified_data[[year]]
    ir_result <- calc_crude_ir(year.dat)
    
    ir_result_r <- ir_result$to_data_frame()
    
    poisson_result <- poisson.exact(
      x = ir_result_r$total_cases,
      T = ir_result_r$total_person_years,
      r = 1,
      alternative = "two.sided",
      conf.level = 0.95,
      midp = TRUE
    )
    
    ir_year[[year]] <- ir_result$with_columns(
      crude_incidence = pl$lit(poisson_result$estimate),
      lower = pl$lit(poisson_result$conf.int[1]),
      upper = pl$lit(poisson_result$conf.int[2]),
      year = pl$lit(year)
    )
  }
  
  combined_df = pl$concat(ir_year)
  return(combined_df)
}

year_cats <- list(
  c(2001, 2002), c(2002, 2003), c(2003, 2004), c(2004, 2005), c(2005, 2006),
  c(2006, 2007), c(2007, 2008), c(2008, 2009), c(2009, 2010), c(2010, 2011),
  c(2011, 2012), c(2012, 2013), c(2013, 2014), c(2014, 2015), c(2015, 2016),
  c(2016, 2017), c(2017, 2018), c(2018, 2019), c(2019, 2020)
)

yearly_ir <- list()

dat_names <- c(
  "samp.dat",
  "samp.white",
  "samp.black",
  "samp.age.30.55",
  "samp.age.55.65"
)

for (dat_name in dat_names) {
  dat <- pl$read_ipc(paste0(data_dir, dat_name, ".w.arrow"))
  for (hiv_status in c("0", "1")) {
    cat("Analyzing", dat_name, "hiv", hiv_status, "\n")
    dat_hiv <- dat$with_columns(pl$col("hiv")$cast(pl$String))$
      filter(pl$col("hiv") == hiv_status)
    year_stratified_base <- create_year_stratified_base(dat_hiv, year_cats)
    yearly_ir[[paste0(dat_name,".hiv", hiv_status)]] <- calculate_overall_metrics(year_stratified_base)$
      with_columns(
        hiv = pl$lit(hiv_status),
        dat = pl$lit(dat_name)
      )
  }
}

pl$concat(yearly_ir, how = "vertical")$
  write_csv(paste0(table_dir, "breast_cx_crude_ir_yearly_for_add_plots.csv"))

plot_data <- pl$concat(yearly_ir, how = "vertical")$
  with_columns(
    year = pl$col("year")$cast(pl$Int32),
    crude_incidence = (pl$col("crude_incidence") * 10000),
    upper = (pl$col("upper") * 10000),
    lower = (pl$col("lower") * 10000)
  )$to_data_frame()



# CIF by hiv status, calendar period, age, and race/ethnicity-------------------

library(tidyverse)
library(patchwork)

table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/compete_risk/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"


csv_files <- list.files(path = table_dir, pattern = "^boot.*\\.csv$", full.names = TRUE)

csv_files_for_plot2 <- NULL

for (fp in csv_files) {
  if (grepl("2001.2012|2012.2020", fp) & !grepl("ex.ny", fp)) {
    csv_files_for_plot2 <- c(csv_files_for_plot2, fp)
  }
  
}


file_info2 <- tibble(full_path = csv_files_for_plot2) %>%
  mutate(filename = basename(full_path)) %>%
  mutate(
    data_name = str_extract(filename, "(?<=boot\\.).*(?=\\.(crude|wt1\\.trunc|wt2\\.trunc))"),
    analysis_type = str_extract(filename, "(crude|wt1\\.trunc|wt2\\.trunc)(?=\\.csv)")
  ) %>%
  mutate(
    col_label = case_when(
      grepl("black", data_name) ~ "Non-Hispanic Black",
      grepl("white", data_name) ~ "Non-Hispanic white",
      grepl("age.30.55", data_name) ~ "Age 30-54",
      grepl("age.55.65", data_name) ~ "Age 55-64",
      TRUE ~ "Overall"
    ),
    row_label = case_when(
      grepl("2001.2012", data_name) ~ "2001-2011",
      grepl("2012.2020", data_name) ~ "2012-2019",
      TRUE ~ NA_character_
    ),
    type_label = case_when(
      grepl("crude", analysis_type) ~ "Crude",
      grepl("wt1.trunc", analysis_type) ~ "Weight 1",
      grepl("wt2.trunc", analysis_type) ~ "Weight 2"
    )
  )


#  CIF -----

cif_data_long <- file_info2 %>%
  select(full_path, col_label, row_label, type_label) %>%
  pmap_dfr(function(full_path, col_label, row_label, type_label) {
    read_csv(full_path, col_types = cols()) %>%
      pivot_longer(
        cols = starts_with(c("risk_", "rd_", "rr_")),
        names_to = "full_measure",
        values_to = "value"
      ) %>%
      mutate(
        measure = str_remove(full_measure, "_ll$|_ul$"),
        ci_type = case_when(
          str_ends(full_measure, "_ll") ~ "lower",
          str_ends(full_measure, "_ul") ~ "upper",
          TRUE ~ "estimate"
        )
      ) %>%
      select(-full_measure) %>%
      pivot_wider(
        names_from = ci_type,
        values_from = value
      ) %>%
      mutate(
        outcome = if_else(str_detect(measure, "breast"), "Breast cancer", "Death"),
        hiv_status = if_else(str_detect(measure, "1$"), "With HIV", "Without HIV"),
        col_label = col_label,
        row_label = row_label,
        plot_type = type_label
      ) %>%
      mutate(
        row_label = factor(row_label, levels = c("2001-2011", "2012-2019")),
        col_label = factor(col_label, levels = c("Overall", "Non-Hispanic Black", "Non-Hispanic white", "Age 30-54", "Age 55-64")),
        hiv_status = factor(hiv_status, levels = c("With HIV", "Without HIV")),
        outcome = factor(outcome, levels = c("Breast cancer", "Death"))
      )
  })

cif_plot_data <- cif_data_long %>%
  filter(str_starts(measure, "risk_"))

create_cif_plot_wo_death <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = estimate, color = hiv_status)) +
    geom_step(linewidth = 0.5) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper, fill = hiv_status),
      alpha = 0.2, color = NA
    ) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    scale_color_manual(
      name = "", 
      values = c("Without HIV" = "darkblue", "With HIV" = "#68ACE5")
    ) + 
    scale_fill_manual(
      name = "",
      values = c("Without HIV" = "darkblue", "With HIV" = "#68ACE5")
    ) +
    labs(
      title = plot_title,
      x = "Age",
      y = "Cumulative Incidence"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, face = "plain"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 14)
    )
}


cif_plot_crude_2_wo_death <- cif_plot_data %>% 
  filter(plot_type == "Crude" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Crude")
cif_plot_wt1_2_wo_death <- cif_plot_data %>% 
  filter(plot_type == "Weight 1" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Weight 1")
cif_plot_wt2_2_wo_death <- cif_plot_data %>%
  filter(plot_type == "Weight 2" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Weight 2")



# RD-------

rd_plot_data <- cif_data_long %>% 
  filter(str_starts(measure, "rd_"))


create_rd_plot_wo_death <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    geom_step(linewidth = 0.5) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.2, color = NA
    ) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    labs(
      title = plot_title,
      x = "Age",
      y = "Risk Difference"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, face = "plain"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 14)
    )
  
}

rd_plot_crude_2_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Crude" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Crude")
rd_plot_wt1_2_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Weight 1" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Weight 1")
rd_plot_wt2_2_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Weight 2" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Weight 2")



plots_to_save <- list(
  cif_period_crude_wo_death_new = cif_plot_crude_2_wo_death,
  rd_period_crude_wo_death_new = rd_plot_crude_2_wo_death,
  cif_period_wt1_wo_death_new = cif_plot_wt1_2_wo_death,
  rd_period_wt1_wo_death_new = rd_plot_wt1_2_wo_death,
  cif_period_wt2_wo_death_new = cif_plot_wt2_2_wo_death,
  rd_period_wt2_wo_death_new = rd_plot_wt2_2_wo_death
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, "boot_", pn, ".jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 10, height = 5, units = "in", res = 300)
  print(p)
  dev.off()
}

cif_comb_crude_wt2 <- plots_to_save[["cif_period_crude_wo_death_new"]] +
  # plots_to_save[["cif_period_wt1_wo_death_new"]] +
  plots_to_save[["cif_period_wt2_wo_death_new"]] +
  plot_layout(ncol = 1, guides = "collect") & 
  theme(
    legend.position = "bottom", 
    plot.title = element_blank(),
    plot.margin = margin(t=20, r=20, b=20, l=20)
  )

rd_comb_crude_wt2 <- plots_to_save[["rd_period_crude_wo_death_new"]] +
  # plots_to_save[["rd_period_wt1_wo_death_new"]] +
  plots_to_save[["rd_period_wt2_wo_death_new"]] +
  plot_layout(ncol = 1, guides = "collect") & 
  theme(
    legend.position = "bottom", 
    plot.title = element_blank(),
    plot.margin = margin(t=20, r=20, b=20, l=20)
  )

plots_to_save_comb <- list(
  cif_comb_crude_wt2 = cif_comb_crude_wt2,
  rd_comb_crude_wt2 = rd_comb_crude_wt2
)

for (pn in names(plots_to_save_comb)) {
  fp <- paste0(figure_dir, "boot_", pn, ".jpeg")
  p <- plots_to_save_comb[[pn]]
  jpeg(fp, width = 11, height = 10, units = "in", res = 300)
  print(p)
  dev.off()
}
