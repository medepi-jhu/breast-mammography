
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged 30-64

# Content: estimate crude IR overall (CIs based on normal approximation), by age and race, by age and calendar period

# Date: 7/18/2025

# Author: Yiyi Zhou

################################################################################

library(tidyverse)

figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/crude_ir/"


df <- read_csv(paste0(table_dir, "breast_cx_crude_ir.csv")) %>%
  rename(
    lower_midp = lower,
    upper_midp = upper
  ) %>%
  mutate(
    lower_normal = crude_incidence - 1.96 * sqrt(total_cases / (total_person_years)^2),
    upper_normal = crude_incidence + 1.96 * sqrt(total_cases / (total_person_years)^2)
  )

df_formatted <- df %>%
  mutate(
    total_person_years = round(total_person_years, 0),
    across(c(crude_incidence, starts_with("lower"), starts_with("upper")), 
           ~round((. * 10000), 2))
  ) %>%
  mutate(
    crude_ir_ci_midp = paste0(crude_incidence, " (", lower_midp, ", ", upper_midp, ")"),
    crude_ir_ci_normal = paste0(crude_incidence, " (", lower_normal, ", ", upper_normal, ")")
  ) %>%
  select(
    data_name, hiv, total_cases, total_person_years, crude_ir_ci_midp, crude_ir_ci_normal
  )

yearly_df <- read_csv(paste0(table_dir, "breast_cx_crude_ir_yearly_for_plot.csv")) %>%
  rename(
    lower_midp = lower,
    upper_midp = upper
  ) %>%
  mutate(
    lower_normal = crude_incidence - 1.96 * sqrt(total_cases / (total_person_years)^2),
    upper_normal = crude_incidence + 1.96 * sqrt(total_cases / (total_person_years)^2)
  )

yearly_df_formatted <- yearly_df %>%
  mutate(
    total_person_years = round(total_person_years, 0),
    across(c(crude_incidence, starts_with("lower"), starts_with("upper")), 
           ~round((. * 10000), 2))
  ) %>%
  mutate(
    crude_ir_ci_midp = paste0(crude_incidence, " (", lower_midp, ", ", upper_midp, ")"),
    crude_ir_ci_normal = paste0(crude_incidence, " (", lower_normal, ", ", upper_normal, ")")
  ) %>%
  select(
    year, hiv, total_cases, total_person_years, crude_ir_ci_midp, crude_ir_ci_normal
  )


