################################################################################
# 
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries >= 30
#
# Content: Visualize CIF for sensitivity analysis (dropping NY 2012-2016)
#
# Date: 10/27/2025
#
# Last update: 
#
# Author: Yiyi Zhou
#
################################################################################

library(tidyverse)
library(patchwork)

table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/compete_risk/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"



# with CIs ----------------------------------------------------------------------

csv_files <- list.files(path = table_dir, pattern = "^boot\\..*ex\\.ny.*\\.csv$", full.names = TRUE)

file_info <- tibble(full_path = csv_files) %>%
  mutate(filename = basename(full_path)) %>%
  mutate(
    data_name = str_extract(filename, "(?<=boot\\.).*(?=\\.(crude|wt1\\.trunc|wt2\\.trunc))"),
    analysis_type = str_extract(filename, "(crude|wt1\\.trunc|wt2\\.trunc)(?=\\.csv)")
  ) %>%
  mutate(
    label = case_when(
      grepl("age.30.55", data_name) ~ "Under 55",
      grepl("age.55.65", data_name) ~ "55 and over",
      grepl("black", data_name) ~ "Black",
      grepl("white", data_name) ~ "White",
      TRUE ~ "Overall"
    ),
    type = case_when(
      grepl("crude", analysis_type) ~ "Crude",
      grepl("wt1.trunc", analysis_type) ~ "Weight 1",
      grepl("wt2.trunc", analysis_type) ~ "Weight 2"
    )
  )



#  CIF -----

cif_data_long <- file_info %>%
  select(full_path, label, type) %>%
  pmap_dfr(function(full_path, label, type) {
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
        label = label,
        plot_type = type
      ) %>%
      mutate(
        label = factor(label, levels = c(
          "Overall", "Under 55", "55 and over", "White", "Black"
        )),
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
      rows = vars(plot_type),
      cols = vars(label)
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
      rows = vars(plot_type),
      cols = vars(label)
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

cif_plot_ex_ny <- cif_plot_data %>% 
  create_cif_plot("Post-2012, excluding New York 2012-2016")


cif_plot_crude_wo_death_ex_ny <- cif_plot_data %>% 
  filter(outcome == "Breast cancer") %>% 
  create_cif_plot_wo_death("Post-2012, excluding New York 2012-2016")


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
      rows = vars(plot_type),
      cols = vars(label)
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
      rows = vars(plot_type),
      cols = vars(label)
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

rd_plot_ex_ny <- rd_plot_data %>% 
  create_rd_plot("Post-2012, excluding New York 2012-2016")


rd_plot_wo_death_ex_ny <- rd_plot_data %>% 
  filter(outcome == "Breast cancer") %>% 
  create_rd_plot_wo_death("Post-2012, excluding New York 2012-2016")



plots_to_save <- list(
  cif_plot_ex_ny = cif_plot_ex_ny,
  cif_plot_wo_death_ex_ny = cif_plot_wo_death_ex_ny,
  rd_plot_ex_ny = rd_plot_ex_ny,
  rd_plot_wo_death_ex_ny = rd_plot_wo_death_ex_ny
)

for (pn in names(plots_to_save)) {
  fp <- paste0(figure_dir, "boot_", pn, ".jpeg")
  p <- plots_to_save[[pn]]
  jpeg(fp, width = 10, height = 7, units = "in", res = 300)
  print(p)
  dev.off()
}



