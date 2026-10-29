
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged >=30

# Content: subset the wide data and estimate crude IR overall, by age and race, by age and calendar period

# Date: 7/7/2025

# Last update: 7/10/2025

# Author: Yiyi Zhou

################################################################################

library(polars)
library(exactci)
library(ggplot2)

data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/crude_ir/"

df <- pl$read_ipc(paste0(data_dir, "samp.dat.w.arrow"))

df$schema



# Functions --------------------------------------------------------------------

# function to subset data by race
race_strat <- function(data, data_name, race) {
  race_mapping <- c("white" = "1", "black" = "2", "hisp" = "3")
  race_label <- race_mapping[[race]]
  race.dat <- data$filter(pl$col("race") == race_label)
  race.dat$write_ipc(paste0(data_dir, data_name, ".", race, ".w.arrow"))
}

# function to subset data by time-updated age group
age_strat <- function(data, age_lower, age_upper, race = NULL, exclude_suffix = NULL) {
  
  age.dat <- data$filter(
    # remove those who entered study after age_upper
    pl$col("dob")$dt$offset_by(paste0(age_upper, "y")) > pl$col("baseline") &
      
      # remove those whose event occurred before age_lower
      pl$col("dob")$dt$offset_by(paste0(age_lower, "y")) <= pl$col("first_event_date")
  )$with_columns(
    
    # keep original baseline and first_event_date
    baseline0 = pl$col("baseline"),
    first_event_date0 = pl$col("first_event_date")
  )$with_columns( 
    
    # recode baseline to account for age_lower
    baseline = pl$when(pl$col("baseline") > pl$col("dob")$dt$offset_by(paste0(age_lower, "y")))$
      then(pl$col("baseline"))$
      otherwise(pl$col("dob")$dt$offset_by(paste0(age_lower, "y"))),
    
    # create index for events after age_upper
    index = (pl$col("first_event_date") >= pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      cast(pl$Int8)
  )$with_columns(
    
    # recode first_event_date
    first_event_date = pl$when(pl$col("index") == 1)$
      then(pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      otherwise(pl$col("first_event_date")),
    
    # update event_flag
    event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      cast(pl$Int8),
    
    # update compete_event_flag
    compete_event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("compete_event_flag"))$
      cast(pl$Int8)
  )$with_columns(
    
    # update event_time (in days) within the specific age range
    event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
  )$with_columns(
    
    # deal with a special case where the recoded baseline and the recoded 
    # first_event_date are the same day (for first_event_date occurred on the date of age_lower)
    event_time = pl$when(pl$col("event_time") == 0)$
      then(pl$col("event_time") + 0.000001)$
      otherwise(pl$col("event_time"))
  )$with_columns(
    
    # re-calculate number of time intervals within the age range
    event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil(),
    
    # add age category label
    age_cat = pl$lit(paste0("[", age_lower, ", ", age_upper, ")"))
  )$drop("index")
  
  if (age.dat$filter(pl$col("event_time") <= 0)$height > 0) {
    cat("samp.age.", "samp.age.", age_lower, ".", age_upper, 
        ".", race, exclude_suffix, "has event_time <= 0!", "\n")
  }
  
  if (is.null(race)) {
    age.dat$write_ipc(paste0(
      data_dir, "samp.age.", age_lower, ".", age_upper, exclude_suffix, ".w",".arrow"
    ))
  } else {
    age.dat$write_ipc(paste0(
      data_dir, "samp.age.", age_lower, ".", age_upper, 
      ".", race, exclude_suffix, ".w",".arrow"
    ))
  }
  # return(age.dat)
}

# function to subset data by year
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
  
  year.dat$write_ipc(paste0(
    data_dir, "samp.year.", year_lower, ".", year_upper, ".w",".arrow"
  ))

}

# function to further stratify year-specific data by time-updated age group
age_strat2 <- function(data, age_lower, age_upper, race = NULL) {
  
  age.dat <- data$filter(
    # remove those who entered study after age_upper
    pl$col("dob")$dt$offset_by(paste0(age_upper, "y")) > pl$col("baseline") &
      # remove those whose event occurred before age_lower
      pl$col("dob")$dt$offset_by(paste0(age_lower, "y")) <= pl$col("first_event_date")
  )$with_columns(
    
    # baseline0: original baseline, baseline1: recoded start date in year-stratified data
    baseline1 = pl$col("baseline"),
    
    # first_event_date0: original baseline, first_event_date1: recoded start date in year-stratified data
    first_event_date1 = pl$col("first_event_date")
  )$with_columns(
    
    # recode baseline to account for age_lower
    baseline = pl$when(pl$col("baseline") > pl$col("dob")$dt$offset_by(paste0(age_lower, "y")))$
      then(pl$col("baseline"))$
      otherwise(pl$col("dob")$dt$offset_by(paste0(age_lower, "y"))),
    
    # create index for events after age_upper
    index = (pl$col("first_event_date") >= pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      cast(pl$Int8)
  )$with_columns(
    
    # recode first_event_date
    first_event_date = pl$when(pl$col("index") == 1)$
      then(pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      otherwise(pl$col("first_event_date")),
    
    # update event_flag
    event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      cast(pl$Int8),
    
    # update compete_event_flag
    compete_event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("compete_event_flag"))$
      cast(pl$Int8)
  )$with_columns(
    
    # update event_time in days
    event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
  )$with_columns(
    # deal with a special case where the recoded baseline and the recoded 
    # first_event_date are the same day (for first_event_date occurred on the date of age_lower)
    event_time = pl$when(pl$col("event_time") == 0)$
      then(pl$col("event_time") + 0.000001)$
      otherwise(pl$col("event_time"))
  )$with_columns(
    # re-calculate number of time intervals
    event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil(),
    
    # add age category label
    age_cat = pl$lit(paste0("[", age_lower, ", ", age_upper, ")"))
  )$filter(pl$col("event_time") > 0)$
    drop("index")
  
  if (age.dat$filter(pl$col("event_time") <= 0)$height > 0) {
    cat("samp.age.", "samp.age.", age_lower, ".", age_upper, 
        ".", race, "has event_time <= 0!", "\n")
  }
  
  if (is.null(race)) {
    age.dat$write_ipc(paste0(
      data_dir, year_dat_name, ".age.", age_lower, ".", age_upper, ".w", ".arrow"
    ))
  } else {
    age.dat$write_ipc(paste0(
      data_dir, year_dat_name, ".age.", age_lower, ".", age_upper, ".", race, ".w", ".arrow"
    ))
  }

}



# Subset wide data -------------------------------------------------------------

# by race
races <- c("white", "black", "hisp")
for (race in races) {
  race_strat(df, "samp", race)
}


# by age
age_cats <- list(c(30, 55), c(55, 65))

for (age_cat in age_cats) {
  age_lower <- age_cat[1]
  age_upper <- age_cat[2]
  age_strat(df, age_lower, age_upper, race = NULL)
}


# by race and age
for (race in races) {
  race.dat <- pl$read_ipc(paste0(data_dir, "samp.", race, ".w.arrow"))
  for (age_cat in age_cats) {
    age_lower <- age_cat[1]
    age_upper <- age_cat[2]
    age_strat(race.dat, age_lower, age_upper, race)
  }
}

# by calendar period
year_cats <- list(c(2001, 2012), c(2012, 2020))

for (year_cat in year_cats) {
  year_lower <- year_cat[1]
  year_upper <- year_cat[2]
  year_strat(df, year_lower, year_upper)
}

# by calendar period and race
races <- c("white", "black")

file_names <- c(
  "samp.year.2001.2012",
  "samp.year.2012.2020"
)

for (fn in file_names) {
  fn.dat <- pl$read_ipc(paste0(data_dir, fn, ".w.arrow"))
  for (race in races) {
    race_strat(fn.dat, fn, race)
  }
}


# by age and calendar period
year_dat_names <- c("samp.year.2001.2012", "samp.year.2012.2020")

for (year_dat_name in year_dat_names) {
  dat <- pl$read_ipc(paste0(data_dir, year_dat_name,".w.arrow"))
  for (age_cat in age_cats) {
    age_lower <- age_cat[1]
    age_upper <- age_cat[2]
    age_strat2(dat, age_lower, age_upper, race = NULL)
  }
}



# Estimate crude incidence rates and CIs----------------------------------------

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

dat_names <- c(
  "samp.dat",
  "samp.white",
  "samp.black",
  "samp.hisp", # only 13 cases among Hispanic PWH
  
  "samp.age.30.55",
  "samp.age.30.55.white",
  "samp.age.30.55.black",
  # "samp.age.30.55.hisp",
  
  "samp.age.55.65",
  "samp.age.55.65.white",
  "samp.age.55.65.black",
  # "samp.age.55.65.hisp",
  
  "samp.year.2001.2012",
  "samp.year.2012.2020",
  "samp.year.2001.2012.age.30.55",
  "samp.year.2012.2020.age.30.55",
  "samp.year.2001.2012.age.55.65",
  "samp.year.2012.2020.age.55.65",
  
  "samp.year.2001.2012.white",
  "samp.year.2001.2012.black",
  "samp.year.2012.2020.white",
  "samp.year.2012.2020.black"
) 

ir_out <- list()

for (dat_name in dat_names) {
  
  dat <- pl$read_ipc(paste0(data_dir, dat_name, ".w.arrow"))$
    select(c("bene_id", "dob", "race", "state", "baseline", 
             "first_event_date", "event_flag", "hiv"))$
    with_columns(pl$col("hiv")$cast(pl$String))
  
  for (hiv_status in c("0", "1")) {
    cat("Analyzing", dat_name, "hiv", hiv_status, "\n")
    dat_hiv <- dat$filter(pl$col("hiv") == hiv_status)
    
    ir_temp <- calc_crude_ir(dat_hiv)$
      with_columns(
        data_name = pl$lit(dat_name),
        hiv = pl$lit(hiv_status)
      )
    
    ir_temp_r <- ir_temp$to_data_frame()
    
    ci_result <- poisson.exact(
      x = ir_temp_r$total_cases,
      T = ir_temp_r$total_person_years,
      r = 1,
      alternative = "two.sided",
      conf.level = 0.95,
      midp = TRUE
    )
    
    ir_temp <- ir_temp$with_columns(
      lower = pl$lit(ci_result$conf.int[1]),
      upper = pl$lit(ci_result$conf.int[2])
    )$select(c(
      "data_name", "hiv", "total_cases", "total_person_years", 
      "crude_incidence", "lower", "upper"
    ))
    
    ir_out[[paste0(dat_name, ".", "hiv", hiv_status)]] <- ir_temp
    
  }
  
}

ir_out_combined_df <- pl$concat(ir_out, how = "vertical")
ir_out_combined_df$filter(pl$col("total_cases") <= 10)$height # 0

ir_out_combined_df_formatted <- ir_out_combined_df$with_columns(
  total_person_years = pl$col("total_person_years")$round(0),
  crude_incidence = (pl$col("crude_incidence") * 10000)$round(2),
  lower = (pl$col("lower") * 10000)$round(2),
  upper = (pl$col("upper") * 10000)$round(2)
)$with_columns(
  pl$concat_str(
    pl$col("crude_incidence")$cast(pl$Utf8),
    pl$lit(" ("),
    pl$col("lower")$cast(pl$Utf8),
    pl$lit(", "),
    pl$col("upper")$cast(pl$Utf8),
    pl$lit(")")
  )$alias("crude_ir_ci")
)$select(c(
  "data_name", "hiv", "total_cases", "total_person_years", "crude_ir_ci"
))

ir_out_combined_df_formatted_r <- ir_out_combined_df_formatted$to_data_frame()

ir_out_combined_df$write_csv(paste0(table_dir, "breast_cx_crude_ir.csv"))
ir_out_combined_df_formatted$write_csv(paste0(table_dir, "breast_cx_crude_ir_formatted.csv"))
ir_out_combined_df_formatted[33:40,]$write_csv(paste0("/users/57285/c-yzhou172-57285/proposed/bx_update_ir_period_race/", "breast_cx_crude_ir_formatted_period_race.csv"))



# Visualize calendar time trends (yearly rates) by HIV status-------------------

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

for (hiv_status in c("0", "1")) {
  cat("Analyzing hiv", hiv_status, "\n")
  df_hiv <- df$with_columns(pl$col("hiv")$cast(pl$String))$
    filter(pl$col("hiv") == hiv_status)
  year_stratified_base <- create_year_stratified_base(df_hiv, year_cats)
  yearly_ir[[paste0("hiv", hiv_status)]] <- calculate_overall_metrics(year_stratified_base)$
    with_columns(hiv = pl$lit(hiv_status))
}

pl$concat(yearly_ir, how = "vertical")$
  write_csv(paste0(table_dir, "breast_cx_crude_ir_yearly_for_plot.csv"))

plot_data <- pl$concat(yearly_ir, how = "vertical")$
  with_columns(
    year = pl$col("year")$cast(pl$Int32),
    crude_incidence = (pl$col("crude_incidence") * 10000),
    upper = (pl$col("upper") * 10000),
    lower = (pl$col("lower") * 10000)
  )$to_data_frame()

nrow(plot_data[plot_data$total_cases <= 10,]) # 0

p <- ggplot(plot_data, aes(x = year, y = crude_incidence, color = hiv, group = hiv)) +
  geom_line(linewidth = 1) +
  geom_point() +
  # geom_smooth(method = "loess", se = FALSE, color = "gray50", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = hiv),
              alpha = 0.2, color = NA) +
  scale_x_continuous(breaks = seq(2001, 2019, by = 1)) +
  scale_color_manual(values = c("0" = "#1A85FF", "1" = "#d41159"),
                     labels = c("0" = "Without HIV", "1" = "With HIV")) +
  scale_fill_manual(values = c("0" = "#1A85FF", "1" = "#d41159"),
                    labels = c("0" = "Without HIV", "1" = "With HIV")) +
  labs(
    x = "Year",
    y = "Crude Incidence Rate (per 10,000 person-years)",
    color = "",
    fill = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

jpeg(paste0(figure_dir, "breast_cx_yearly_crude_ir", ".jpeg"), 
     width = 9, height = 5.5, units = "in", res = 300)
print(p)
dev.off()
