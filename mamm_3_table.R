
################################################################################
#
# Project: MedEpi Mammography 
#
# Purpose: Build tables
#
# Last Update: 04 Nov 2025
#
################################################################################


library("tidyverse")


# Read in results ---------------------------------------------------------

crude <- read_csv("../results/crude_rates.csv")
wt <- read_csv("../results/wt_rates.csv")
prob <- read_csv("../results/prob2yr.csv")
detect <- read_csv("../results/screen_detect.csv")


# Table 2 -----------------------------------------------------------------

wt.hiv <- wt %>% 
  mutate(hiv = "HIV",
         rate_crude = round(rate_crude_hiv*100, 1), 
         irr_crude = round(exp(logirr_crude), 2), 
         irr_crude_ll = round(exp(logirr_crude_ll), 2), 
         irr_crude_ul = round(exp(logirr_crude_ul), 2), 
         irr_crude_ci = paste0("(", irr_crude_ll, ", ", irr_crude_ul, ")"),
         
         rate_wt1 = round(rate_wt1_hiv*100, 1),
         irr_wt1 = round(exp(logirr_wt1), 2), 
         irr_wt1_ll = round(exp(logirr_wt1_ll), 2), 
         irr_wt1_ul = round(exp(logirr_wt1_ul), 2), 
         irr_wt1_ci = paste0("(", irr_wt1_ll, ", ", irr_wt1_ul, ")"),
         
         rate_wt2 = round(rate_wt2_hiv*100, 1),
         irr_wt2 = round(exp(logirr_wt2), 2), 
         irr_wt2_ll = round(exp(logirr_wt2_ll), 2), 
         irr_wt2_ul = round(exp(logirr_wt2_ul), 2),
         irr_wt2_ci = paste0("(", irr_wt2_ll, ", ", irr_wt2_ul, ")")) %>% 
  select(period, group, hiv,
         rate_crude, irr_crude, irr_crude_ci, 
         rate_wt1, irr_wt1, irr_wt1_ci,
         rate_wt2, irr_wt2, irr_wt2_ci)

wt.nohiv <- wt %>% 
  mutate(hiv = "No HIV",
         rate_crude = round(rate_crude_nohiv*100, 1), 
         irr_crude = 1, 
         irr_crude_ci = NA,
         
         rate_wt1 = round(rate_wt1_nohiv*100, 1),
         irr_wt1 = 1, 
         irr_wt1_ci = NA,
         
         rate_wt2 = round(rate_wt2_nohiv*100, 1),
         irr_wt2 = 1, 
         irr_wt2_ci = NA) %>% 
  select(period, group, hiv,
         rate_crude, irr_crude, irr_crude_ci, 
         rate_wt1, irr_wt1, irr_wt1_ci,
         rate_wt2, irr_wt2, irr_wt2_ci)

wt.table <- bind_rows(wt.nohiv, wt.hiv) %>% 
  mutate(period = factor(period, levels=c("pre2012", "post2012")),
         group = factor(group, levels=c("Overall", "Black", "White"))) %>% 
  arrange(period, group, desc(hiv))

write_csv(wt.table, "../results/table2.csv")


# Table 4 -----------------------------------------------------------------

detect2 <- detect %>% 
  select(-n) %>% 
  filter(!(analysis %in% c("ipw1", "ipw2"))) %>% 
  mutate(hiv = ifelse(hiv==0, "No HIV", "HIV"),
         period = factor(period, levels=c("pre2012", "post2012")),
         group = factor(group, levels=c("Overall", "Black", "White")),
         across(!c(hiv, period, group, analysis), ~ round(.x, 1)),
         ci = paste0("(", rate_ll, ", ", rate_ul, ")")) %>% 
  pivot_wider(id_cols=c(period, group, hiv),
              names_from = analysis,
              values_from = c(rate, ci)) %>% 
  select(period, group, hiv, 
         rate_crude, ci_crude, 
         `rate_ipw1+ipsw`, `ci_ipw1+ipsw`,
         `rate_ipw2+ipsw`, `ci_ipw2+ipsw`)

write_csv(detect2, "../results/table4.csv")
