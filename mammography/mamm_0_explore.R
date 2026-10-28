
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Manage raw data and preliminary analyses
#
# Last Update: 05 May 2025
#
################################################################################


library("tidyverse")


# Read in data ------------------------------------------------------------

dat <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/bene.csv",
                col_types=paste0("cdcc", paste(rep("D", 22), collapse=""))) %>%
  rename(start_date = first_elig_period_start_date_65,
         end_date = first_elig_period_end_date_65) %>%
  select(bene_id, start_date, end_date, first_elig_state, dob, dod, breast_cx_type1_dt, breast_cx_type2_dt)

proc <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/porc.csv",
                 col_types=c("cDddddd")) %>%
  mutate(year = year(proc_date))


# Estimate rates by state and year ----------------------------------------

rates <- function(analysis) {

res <- NULL

for (state in c("AL", "AR", "AZ", "CA", "CO", "DC", "FL", "GA", "IL", "IN", "KY",
                "MA", "MD", "MI", "MO", "MS", "NC", "NJ", "NV", "NY", "OH", "OK",
                "PA", "TN", "TX", "SC", "WA")) {

  temp.res <- data.frame(analysis = analysis,
                         year = seq(2001, 2021, 1),
                         state = state)


# Manage individual data --------------------------------------------------
  # Person-time in those without breast cancer
  dat2 <- dat %>%
    filter(first_elig_state==state) %>%
    mutate(base = start_date + months(6),

           # Recode leap year birthdays
           dob = ifelse(month(dob)==2 & day(dob)==29, dob - days(1), dob),
           dob = as.Date(dob, origin=ymd("1970-01-01")),

           age40 = dob + years(40),
           age50 = dob + years(50))

  # Recode basline if subsetting by age
  if (analysis=="all") {
    dat2$baseline <- dat2$base

  } else if (analysis==">40") {
    dat2$baseline <- pmax(dat2$base, dat2$age40)

  } else if (analysis==">50") {
    dat2$baseline <- pmax(dat2$base, dat2$age50)

  }

  # Recode censoring date as needed
  dat2$event_date <- pmin(dat2$end_date, dat2$breast_cx_type2_dt, dat2$dod, na.rm=T)

  dat3 <- dat2 %>%
    mutate(start_year = year(baseline),
           end_year = year(event_date)) %>%
    filter(event_date>baseline) %>%
    filter(is.na(breast_cx_type1_dt) | breast_cx_type1_dt>baseline)

  # Make data long
  dat3$r <- dat3$end_year - dat3$start_year + 1
  long <- bind_rows(lapply(dat3, rep, dat3$r)) %>%
    group_by(bene_id) %>%
    mutate(add = row_number() - 1) %>%
    ungroup()

  long2 <- long %>%
    mutate(year = start_year + add,
           person_years = case_when(year==start_year & end_year==year ~ as.numeric(event_date - baseline)/365.25,
                                    year==start_year & end_year>year ~ as.numeric(ymd(paste0(year, "-12-31")) - baseline )/365.25,
                                    year>start_year & end_year==year ~ as.numeric(event_date - ymd(paste0(year, "-01-01")))/365.25,
                                    year>start_year & end_year>year ~ 1,
                                    T ~ NA_real_)) %>%
    filter(person_years>0) %>%
    select(-c(r, add))


# Estimate rates ---------------------------------------------------------

  time <- long2 %>%
    group_by(year) %>%
    summarize(tot_py = sum(person_years),
              n = n())

  # Rate of all mammographies
  all.mamm <- proc %>%
    left_join(select(dat3, bene_id, baseline, event_date), by="bene_id") %>%
    filter(bene_id %in% dat3$bene_id) %>%
    filter(proc_date>=baseline) %>%
    filter(proc_date<=event_date) %>%
    group_by(year) %>%
    summarize(n_proc = n())

  all.mamm2 <- left_join(time, all.mamm, by="year") %>%
    mutate(rate_total = n_proc/tot_py)

  # Rate of screening mammographies
  screen.mamm <- proc %>%
    filter(bene_id %in% dat3$bene_id) %>%
    left_join(select(dat2, bene_id, breast_cx_type2_dt, baseline, event_date), by="bene_id") %>%
    group_by(bene_id) %>%
    mutate(last_mamm = lag(proc_date)) %>%
    filter(proc_date>=baseline) %>%
    filter(proc_date<=event_date) %>%
    ungroup()

  screen.mamm2 <- screen.mamm %>%
    mutate(mamm_9month = case_when(is.na(last_mamm) ~ 0,
                                   ((last_mamm %--% proc_date)/months(1))>=9 ~ 0,
                                   ((last_mamm %--% proc_date)/months(1))<9 ~ 1,
                                   T ~ NA_real_),
           diag_1year = case_when(is.na(breast_cx_type2_dt) ~ 0,
                                  breast_cx_type2_dt>proc_date ~ 0,
                                  ((breast_cx_type2_dt %--% proc_date)/years(1))>=1 ~ 0,
                                  ((breast_cx_type2_dt %--% proc_date)/years(1))<1 ~ 1,
                                  T ~ NA_real_)) %>%
    filter(mammograph_typ==2 & mamm_9month==0 & diag_1year==0)

  screen.mamm3 <- screen.mamm2 %>%
    group_by(year) %>%
    summarize(n_proc = n())

  screen.mamm4 <- left_join(time, screen.mamm3, by="year") %>%
    mutate(rate_screen = n_proc/tot_py)

  temp.res <- temp.res %>%
    left_join(select(time, year, n), by="year") %>%
    left_join(select(all.mamm2, year, rate_total), by="year") %>%
    left_join(select(screen.mamm4, year, rate_screen), by="year")

  res <- bind_rows(res, temp.res)

}

return(res)
}

all.res <- lapply(c("all", ">40", ">50"), function(x){rates(x)})
all.res <- bind_rows(all.res)

write_csv(all.res, "./results/rates_state_yr.csv")
