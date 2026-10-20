
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Summarize bootstraps
#
# Last Update: 30 Oct 2025
#
################################################################################


library("tidyverse")


# Function to read results ------------------------------------------------

read <- function(r, analysis, group) {

  if (analysis=="prob2yr") {
    if (file.exists(paste0("~/mammography/results/", analysis, "_", group, "_", r, ".csv"))) {
      res <- read_csv(paste0("~/mammography/results/", analysis, "_", group, "_", r, ".csv"))
    } else {
      res <- NULL
    }

  } else {

    if (file.exists(paste0("~/mammography/results/", analysis, "_", group, "_pre2012_", r, ".csv"))) {
      res_pre <- read_csv(paste0("~/mammography/results/", analysis, "_", group, "_pre2012_", r, ".csv")) %>%
        mutate(period="pre2012")
    } else {
      res_pre <- NULL
    }

    if (file.exists(paste0("~/mammography/results/", analysis, "_", group, "_post2012_", r, ".csv"))) {
      res_post <- read_csv(paste0("~/mammography/results/", analysis, "_", group, "_post2012_", r, ".csv")) %>%
        mutate(period="post2012")
    } else {
      res_post <- NULL
    }

    res <- bind_rows(res_pre, res_post)
  }

  return(res)

}


# Rates -------------------------------------------------------------------

res_overall <- lapply(c(0, seq(100, 1000, 100)), function(x){read(x, analysis="wt_rates", group="overall")})
res_overall <- bind_rows(res_overall) %>%
  mutate(group = "Overall") %>%
  group_by(period, rep) %>%
  filter(!duplicated(rep)) %>%  # I had to re-run boot=0 for the correct n values
  ungroup()

res_black <- lapply(c(0, seq(100, 1000, 100)), function(x){read(x, analysis="wt_rates", group="black")})
res_black <- bind_rows(res_black) %>%
  mutate(group = "Black") %>%
  group_by(period, rep) %>%
  filter(!duplicated(rep)) %>%  # I had to re-run boot=0 for the correct n values
  ungroup()

res_white <- lapply(c(0, seq(100, 1000, 100)), function(x){read(x, analysis="wt_rates", group="white")})
res_white <- bind_rows(res_white) %>%
  mutate(group = "White") %>%
  group_by(period, rep) %>%
  filter(!duplicated(rep)) %>%  # I had to re-run boot=0 for the correct n values
  ungroup()

res_rates <- bind_rows(res_overall, res_black, res_white)

res_rates_summ <- res_rates %>%
  group_by(group, period) %>%
  summarize(across(c(logirr_crude, logirr_wt1, logirr_wt2), ~quantile(.x, p=0.025), .names="{col}_ll"),
            across(c(logirr_crude, logirr_wt1, logirr_wt2), ~quantile(.x, p=0.975), .names="{col}_ul"))

res_rates <- filter(res_rates, rep==0) %>%
  left_join(res_rates_summ, by=c("group", "period")) %>%
  select(-rep)

sum(res_rates$n_hiv<11)

write_csv(res_rates, "~/mammography/results/wt_rates.csv")


# Probability of mammography in last 2 years ------------------------------

res_overall <- lapply(seq(100, 1000, 100), function(x){read(x, analysis="prob2yr", group="overall")})
res_overall <- bind_rows(res_overall) %>%
  mutate(group = "Overall")

res_black <- lapply(seq(100, 1000, 100), function(x){read(x, analysis="prob2yr", group="black")})
res_black <- bind_rows(res_black) %>%
  mutate(group = "Black")

res_white <- lapply(seq(100, 1000, 100), function(x){read(x, analysis="prob2yr", group="white")})
res_white <- bind_rows(res_white) %>%
  mutate(group = "White")

res_prob <- bind_rows(res_overall, res_black, res_white)

res_prob_summ <- res_prob %>%
  group_by(group, year, analysis, hiv) %>%
  summarize(p_screened_ll = quantile(p_screened, p=0.025),
            p_screened_ul = quantile(p_screened, p=0.975))

res_prob <- filter(res_prob, boot_rep==0) %>%
  left_join(res_prob_summ, by=c("group", "year", "analysis", "hiv")) %>%
  select(-boot_rep)

sum(res_prob$n<11)

write_csv(res_prob, "~/mammography/results/prob2yr.csv")


# Screen-detected breast cancers ------------------------------------------

res_overall <- lapply(c(1000), function(x){read(x, analysis="screen_detect", group="overall")})
res_overall <- bind_rows(res_overall) %>%
  mutate(group = "Overall")

res_black <- lapply(c(1000), function(x){read(x, analysis="screen_detect", group="black")})
res_black <- bind_rows(res_black) %>%
  mutate(group = "Black")

res_white <- lapply(c(1000), function(x){read(x, analysis="screen_detect", group="white")})
res_white <- bind_rows(res_white) %>%
  mutate(group = "White")

res_detect <- bind_rows(res_overall, res_black, res_white) %>%
  select(-rep)

sum(res_detect$n<11)

write_csv(res_detect, "~/mammography/results/screen_detect.csv")
