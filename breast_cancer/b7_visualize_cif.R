################################################################################
# 
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries >= 30
#
# Content: Visualize CIF
#
# Date: 7/18/2025
#
# Last update: 7/28/2025
#
# Author: Yiyi Zhou
#
################################################################################

library(tidyverse)
library(patchwork)

table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/compete_risk/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"



# point estimates only | focus: split by age------------------------------------

csv_files <- list.files(path = table_dir, pattern = "^samp.*\\.csv$", full.names = TRUE)

skip_csv_files <- NULL

for (fp in csv_files) {
  if (grepl(
    "samp.year.2001.2012.black|samp.year.2001.2012.white|samp.year.2012.2020.black|samp.year.2012.2020.white", 
    fp
  )) {
    skip_csv_files <- c(skip_csv_files, fp)
  }
}

csv_files_for_plot1 <- setdiff(csv_files, skip_csv_files)

file_info1 <- tibble(full_path = csv_files_for_plot1) %>%
  mutate(filename = basename(full_path)) %>%
  mutate(
    data_name = str_extract(filename, "^.*(?=_point_est_)"),
    analysis_type = str_extract(filename, "(?<=_point_est_).*(?=\\.csv$)")
  ) %>%
  mutate(
    col_label = case_when(
      grepl("black", data_name) ~ "Black",
      grepl("white", data_name) ~ "White",
      grepl("hisp", data_name) ~ "Hispanic",
      grepl("2001.2012", data_name) ~ "Pre-2012",
      grepl("2012.2020", data_name) ~ "Post-2012",
      TRUE ~ "Overall"),
    row_label = case_when(
      grepl("age.30.55", data_name) ~ "Under 55",
      grepl("age.55.65", data_name) ~ "55 and over",
      TRUE ~ "All ages"
    ),
    type_label = case_when(
      grepl("crude", analysis_type) ~ "Crude",
      grepl("wt1_trunc", analysis_type) ~ "Weight 1",
      grepl("wt2_trunc", analysis_type) ~ "Weight 2"
    )
  )



#  CIF -----

cif_data_long <- file_info1 %>%
  select(full_path, col_label, row_label, type_label) %>%
  pmap_dfr(function(full_path, col_label, row_label, type_label) {
    read_csv(full_path, col_types = cols()) %>%
      pivot_longer(
        cols = starts_with("risk_"),
        names_to = "measure",
        values_to = "cif_value"
      ) %>%
      mutate(
        outcome = if_else(str_detect(measure, "breast"), "Breast cancer", "Death"),
        hiv_status = if_else(str_detect(measure, "1$"), "With HIV", "Without HIV"),
        col_label = col_label,
        row_label = row_label,
        plot_type = type_label
      )
  })

cif_plot_data <- cif_data_long %>%
  mutate(
    row_label = factor(row_label, levels = c("All ages", "Under 55", "55 and over")),
    col_label = factor(col_label, levels = c("Overall", "White", "Black", "Hispanic", "Pre-2012", "Post-2012")),
    hiv_status = factor(hiv_status, levels = c("With HIV", "Without HIV")),
    outcome = factor(outcome, levels = c("Breast cancer", "Death"))
  )

create_cif_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = cif_value, color = hiv_status, linetype = outcome)) +
    geom_step(linewidth = 0.5) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    scale_color_manual(
      name = "", 
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
    ) + 
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
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


cif_plot_crude <- cif_plot_data %>% filter(plot_type == "Crude") %>% create_cif_plot("Crude")
cif_plot_wt1 <- cif_plot_data %>% filter(plot_type == "Weight 1") %>% create_cif_plot("Weight 1")
cif_plot_wt2 <- cif_plot_data %>% filter(plot_type == "Weight 2") %>% create_cif_plot("Weight 2")



# RD-------

rd_data_long <- file_info1 %>%
  select(full_path, col_label, row_label, type_label) %>%
  pmap_dfr(function(full_path, col_label, row_label, type_label) {
    read_csv(full_path, col_types = cols()) %>%
      pivot_longer(
        cols = starts_with("rd_"),
        names_to = "measure",
        values_to = "rd_value"
      ) %>%
      mutate(
        outcome = if_else(measure == "rd_breast", "Breast cancer", "Death"),
        hiv_status = if_else(str_detect(measure, "1$"), "With HIV", "Without HIV"),
        col_label = col_label,
        row_label = row_label,
        plot_type = type_label
      )
  })

rd_plot_data <- rd_data_long %>%
  mutate(
    row_label = factor(row_label, levels = c("All ages", "Under 55", "55 and over")),
    col_label = factor(col_label, levels = c("Overall", "White", "Black", "Hispanic", "Pre-2012", "Post-2012")),
    outcome = factor(outcome, levels = c("Breast cancer", "Death"))
  )

create_rd_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = rd_value, linetype = outcome)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    geom_step(linewidth = 0.5) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
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


rd_plot_crude <- rd_plot_data %>% filter(plot_type == "Crude") %>% create_rd_plot("Crude")
rd_plot_wt1 <- rd_plot_data %>% filter(plot_type == "Weight 1") %>% create_rd_plot("Weight 1")
rd_plot_wt2 <- rd_plot_data %>% filter(plot_type == "Weight 2") %>% create_rd_plot("Weight 2")


plots_to_save <- list(
  cif_age_crude = cif_plot_crude,
  rd_age_crude = rd_plot_crude,
  cif_age_wt1 = cif_plot_wt1,
  rd_age_wt1 = rd_plot_wt1,
  cif_age_wt2 = cif_plot_wt2,
  rd_age_wt2 = rd_plot_wt2
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, pn, "_point_est.jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 12, height = 7, units = "in", res = 300)
  print(p)
  dev.off()
}



# point estimates only | focus: split by calendar period------------------------


csv_files_for_plot2 <- NULL

for (fp in csv_files) {
  if (grepl("2001.2012|2012.2020", fp) & !grepl("age", fp)) {
    csv_files_for_plot2 <- c(csv_files_for_plot2, fp)
  }
  
}

file_info2 <- tibble(full_path = csv_files_for_plot2) %>%
  mutate(filename = basename(full_path)) %>%
  mutate(
    data_name = str_extract(filename, "^.*(?=_point_est_)"),
    analysis_type = str_extract(filename, "(?<=_point_est_).*(?=\\.csv$)")
  ) %>%
  mutate(
    col_label = case_when(
      grepl("black", data_name) ~ "Black",
      grepl("white", data_name) ~ "White",
      TRUE ~ "Overall"
    ),
    row_label = case_when(
      grepl("2001.2012", data_name) ~ "Pre-2012",
      grepl("2012.2020", data_name) ~ "Post-2012",
      TRUE ~ NA_character_
    ),
    type_label = case_when(
      grepl("crude", analysis_type) ~ "Crude",
      grepl("wt1_trunc", analysis_type) ~ "Weight 1",
      grepl("wt2_trunc", analysis_type) ~ "Weight 2"
    )
  )


#  CIF -----

cif_data_long <- file_info2 %>%
  select(full_path, col_label, row_label, type_label) %>%
  pmap_dfr(function(full_path, col_label, row_label, type_label) {
    read_csv(full_path, col_types = cols()) %>%
      pivot_longer(
        cols = starts_with("risk_"),
        names_to = "measure",
        values_to = "cif_value"
      ) %>%
      mutate(
        outcome = if_else(str_detect(measure, "breast"), "Breast cancer", "Death"),
        hiv_status = if_else(str_detect(measure, "1$"), "With HIV", "Without HIV"),
        col_label = col_label,
        row_label = row_label,
        plot_type = type_label
      )
  })

cif_plot_data <- cif_data_long %>%
  mutate(
    row_label = factor(row_label, levels = c("Pre-2012", "Post-2012")),
    col_label = factor(col_label, levels = c("Overall", "White", "Black")),
    hiv_status = factor(hiv_status, levels = c("With HIV", "Without HIV")),
    outcome = factor(outcome, levels = c("Breast cancer", "Death"))
  )

create_cif_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = cif_value, color = hiv_status, linetype = outcome)) +
    geom_step(linewidth = 0.5) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    scale_color_manual(
      name = "", 
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
    ) + 
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
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


cif_plot_crude_2 <- cif_plot_data %>% filter(plot_type == "Crude") %>% create_cif_plot("Crude")
cif_plot_wt1_2 <- cif_plot_data %>% filter(plot_type == "Weight 1") %>% create_cif_plot("Weight 1")
cif_plot_wt2_2 <- cif_plot_data %>% filter(plot_type == "Weight 2") %>% create_cif_plot("Weight 2")


# RD-------

rd_data_long <- file_info2 %>%
  select(full_path, col_label, row_label, type_label) %>%
  pmap_dfr(function(full_path, col_label, row_label, type_label) {
    read_csv(full_path, col_types = cols()) %>%
      pivot_longer(
        cols = starts_with("rd_"),
        names_to = "measure",
        values_to = "rd_value"
      ) %>%
      mutate(
        outcome = if_else(measure == "rd_breast", "Breast cancer", "Death"),
        hiv_status = if_else(str_detect(measure, "1$"), "With HIV", "Without HIV"),
        col_label = col_label,
        row_label = row_label,
        plot_type = type_label
      )
  })

rd_plot_data <- rd_data_long %>%
  mutate(
    row_label = factor(row_label, levels = c("Pre-2012", "Post-2012")),
    col_label = factor(col_label, levels = c("Overall", "White", "Black")),
    outcome = factor(outcome, levels = c("Breast cancer", "Death"))
  )

create_rd_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = rd_value, linetype = outcome)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
    geom_step(linewidth = 0.5) +
    facet_grid(
      rows = vars(row_label),
      cols = vars(col_label)
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
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


rd_plot_crude_2 <- rd_plot_data %>% filter(plot_type == "Crude") %>% create_rd_plot("Crude")
rd_plot_wt1_2 <- rd_plot_data %>% filter(plot_type == "Weight 1") %>% create_rd_plot("Weight 1")
rd_plot_wt2_2 <- rd_plot_data %>% filter(plot_type == "Weight 2") %>% create_rd_plot("Weight 2")


plots_to_save <- list(
  cif_period_crude = cif_plot_crude_2,
  rd_period_crude = rd_plot_crude_2,
  cif_period_wt1 = cif_plot_wt1_2,
  rd_period_wt1 = rd_plot_wt1_2,
  cif_period_wt2 = cif_plot_wt2_2,
  rd_period_wt2 = rd_plot_wt2_2
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, pn, "_point_est.jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 7, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
}




# with CIs | focus: split by age------------------------------------------------

csv_files <- list.files(path = table_dir, pattern = "^boot.*\\.csv$", full.names = TRUE)

skip_csv_files <- NULL

for (fp in csv_files) {
  if (grepl(
    "boot.samp.year.2001.2012.black|boot.samp.year.2001.2012.white|boot.samp.year.2012.2020.black|boot.samp.year.2012.2020.white", 
    fp
  )) {
    skip_csv_files <- c(skip_csv_files, fp)
  }
}

csv_files_for_plot1 <- setdiff(csv_files, skip_csv_files)

file_info1 <- tibble(full_path = csv_files_for_plot1) %>%
  mutate(filename = basename(full_path)) %>%
  mutate(
    data_name = str_extract(filename, "(?<=boot\\.).*(?=\\.(crude|wt1\\.trunc|wt2\\.trunc))"),
    analysis_type = str_extract(filename, "(crude|wt1\\.trunc|wt2\\.trunc)(?=\\.csv)")
  ) %>%
  mutate(
    col_label = case_when(
      grepl("black", data_name) ~ "Black",
      grepl("white", data_name) ~ "White",
      grepl("hisp", data_name) ~ "Hispanic",
      grepl("2001.2012", data_name) ~ "Pre-2012",
      grepl("2012.2020", data_name) ~ "Post-2012",
      TRUE ~ "Overall"),
    row_label = case_when(
      grepl("age.30.55", data_name) ~ "Under 55",
      grepl("age.55.65", data_name) ~ "55 and over",
      TRUE ~ "All ages"
    ),
    type_label = case_when(
      grepl("crude", analysis_type) ~ "Crude",
      grepl("wt1.trunc", analysis_type) ~ "Weight 1",
      grepl("wt2.trunc", analysis_type) ~ "Weight 2"
    )
  )



#  CIF -----

cif_data_long <- file_info1 %>%
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
        row_label = factor(row_label, levels = c("All ages", "Under 55", "55 and over")),
        col_label = factor(col_label, levels = c("Overall", "White", "Black", "Hispanic", "Pre-2012", "Post-2012")),
        hiv_status = factor(hiv_status, levels = c("With HIV", "Without HIV")),
        outcome = factor(outcome, levels = c("Breast cancer", "Death"))
      )
  })

cif_plot_data <- cif_data_long %>%
  filter(str_starts(measure, "risk_"))

create_cif_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = estimate, color = hiv_status, linetype = outcome)) +
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
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
    ) + 
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
    ) +
    scale_fill_manual(
      name = "",
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
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
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
    ) + 
    scale_fill_manual(
      name = "",
      values = c("Without HIV" = "#1A85FF", "With HIV" = "#d41159")
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

cif_plot_crude <- cif_plot_data %>% filter(plot_type == "Crude") %>% create_cif_plot("Crude")
cif_plot_wt1 <- cif_plot_data %>% filter(plot_type == "Weight 1") %>% create_cif_plot("Weight 1")
cif_plot_wt2 <- cif_plot_data %>% filter(plot_type == "Weight 2") %>% create_cif_plot("Weight 2")

cif_plot_crude_wo_death <- cif_plot_data %>% 
  filter(plot_type == "Crude" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Crude")
cif_plot_wt1_wo_death <- cif_plot_data %>% 
  filter(plot_type == "Weight 1" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Weight 1")
cif_plot_wt2_wo_death <- cif_plot_data %>% 
  filter(plot_type == "Weight 2" & outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Weight 2")


# RD-------

rd_plot_data <- cif_data_long %>%
  filter(str_starts(measure, "rd_"))

create_rd_plot <- function(data, plot_title) {
  ggplot(data, aes(x = time, y = estimate, linetype = outcome)) +
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
    scale_linetype_manual(
      name = "",
      values = c("Breast cancer" = "solid", "Death" = "dashed")
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

rd_plot_crude <- rd_plot_data %>% filter(plot_type == "Crude") %>% create_rd_plot("Crude")
rd_plot_wt1 <- rd_plot_data %>% filter(plot_type == "Weight 1") %>% create_rd_plot("Weight 1")
rd_plot_wt2 <- rd_plot_data %>% filter(plot_type == "Weight 2") %>% create_rd_plot("Weight 2")


rd_plot_crude_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Crude" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Crude")
rd_plot_wt1_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Weight 1" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Weight 1")
rd_plot_wt2_wo_death <- rd_plot_data %>% 
  filter(plot_type == "Weight 2" & outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Weight 2")


plots_to_save <- list(
  cif_age_crude = cif_plot_crude,
  rd_age_crude = rd_plot_crude,
  cif_age_wt1 = cif_plot_wt1,
  rd_age_wt1 = rd_plot_wt1,
  cif_age_wt2 = cif_plot_wt2,
  rd_age_wt2 = rd_plot_wt2,
  cif_age_crude_wo_death = cif_plot_crude_wo_death,
  rd_age_crude_wo_death = rd_plot_crude_wo_death,
  cif_age_wt1_wo_death = cif_plot_wt1_wo_death,
  rd_age_wt1_wo_death = rd_plot_wt1_wo_death,
  cif_age_wt2_wo_death = cif_plot_wt2_wo_death,
  rd_age_wt2_wo_death = rd_plot_wt2_wo_death
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, "boot_", pn, ".jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 12, height = 7, units = "in", res = 300)
  print(p)
  dev.off()
}



# with CIs | focus: split by calendar period------------------------


csv_files_for_plot2 <- NULL

for (fp in csv_files) {
  if (grepl("2001.2012|2012.2020", fp) & !grepl("age", fp)) {
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
      grepl("black", data_name) ~ "Black",
      grepl("white", data_name) ~ "White",
      TRUE ~ "Overall"
    ),
    row_label = case_when(
      grepl("2001.2012", data_name) ~ "Pre-2012",
      grepl("2012.2020", data_name) ~ "Post-2012",
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
        row_label = factor(row_label, levels = c("Pre-2012", "Post-2012")),
        col_label = factor(col_label, levels = c("Overall", "White", "Black")),
        hiv_status = factor(hiv_status, levels = c("With HIV", "Without HIV")),
        outcome = factor(outcome, levels = c("Breast cancer", "Death"))
      )
  })

cif_plot_data <- cif_data_long %>%
  filter(str_starts(measure, "risk_"))


cif_plot_crude_2 <- cif_plot_data %>% filter(plot_type == "Crude") %>% create_cif_plot("Crude")
cif_plot_wt1_2 <- cif_plot_data %>% filter(plot_type == "Weight 1") %>% create_cif_plot("Weight 1")
cif_plot_wt2_2 <- cif_plot_data %>% filter(plot_type == "Weight 2") %>% create_cif_plot("Weight 2")

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


rd_plot_crude_2 <- rd_plot_data %>% filter(plot_type == "Crude") %>% create_rd_plot("Crude")
rd_plot_wt1_2 <- rd_plot_data %>% filter(plot_type == "Weight 1") %>% create_rd_plot("Weight 1")
rd_plot_wt2_2 <- rd_plot_data %>% filter(plot_type == "Weight 2") %>% create_rd_plot("Weight 2")

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
  cif_period_crude = cif_plot_crude_2,
  rd_period_crude = rd_plot_crude_2,
  cif_period_wt1 = cif_plot_wt1_2,
  rd_period_wt1 = rd_plot_wt1_2,
  cif_period_wt2 = cif_plot_wt2_2,
  rd_period_wt2 = rd_plot_wt2_2,
  cif_period_crude_wo_death = cif_plot_crude_2_wo_death,
  rd_period_crude_wo_death = rd_plot_crude_2_wo_death,
  cif_period_wt1_wo_death = cif_plot_wt1_2_wo_death,
  rd_period_wt1_wo_death = rd_plot_wt1_2_wo_death,
  cif_period_wt2_wo_death = cif_plot_wt2_2_wo_death,
  rd_period_wt2_wo_death = rd_plot_wt2_2_wo_death
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, "boot_", pn, ".jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 7, height = 6, units = "in", res = 300)
  print(p)
  dev.off()
}
