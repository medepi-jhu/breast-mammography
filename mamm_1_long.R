
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Set up data for weighted analysis
#
# Last Update: 03 Jul 2025
#
################################################################################


library("tidyverse")

period <- "all"


# Read in data ------------------------------------------------------------

dat <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/bene.csv",
                col_types=paste0("cdcc", paste(rep("D", 22), collapse=""))) %>%
  rename(start_date = first_elig_period_start_date_65,
         end_date = first_elig_period_end_date_65,
         state = first_elig_state) %>%
  mutate(region = case_when(state %in% c("IA", "IL", "IN", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD",
                                         "WI") ~ "Midwest",
                            state %in% c("CT", "MA", "ME", "NH", "NJ", "NY", "PA", "RI", "VT") ~ "Northeast",
                            state %in% c("AL", "AR", "DC", "DE", "FL", "GA", "KY", "LA", "MD", "MS", "NC",
                                         "OK", "SC", "TN", "TX", "VA", "WV") ~ "South",
                            state %in% c("AK", "AZ", "CA", "CO", "HI", "ID", "MT", "NM", "OR", "NV", "UT",
                                         "WA", "WY") ~ "West"))

proc <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/porc.csv",
                 col_types=c("cDddddd")) %>%
  arrange(bene_id, proc_date)


# Manage individual data --------------------------------------------------

dat2 <- dat %>%
  mutate(base = start_date + months(6),

           # Recode leap year birthdays
           dob = ifelse(month(dob)==2 & day(dob)==29, dob - days(1), dob),
           dob = as.Date(dob, origin=ymd("1970-01-01")),

           age40 = dob + years(40),
           age50 = dob + years(50),
           age65 = dob + years(65)) %>%
    # Exclude MD (due to breast cx and mammography data quality) and all new states
    filter(state %in% c("AL", "CA", "CO", "FL", "GA", "IL", "MA", "NC", "NY", "OH", "PA", "TX", "WA"))

# Recode baseline and end date according to age ranges
dat2$baseline <- pmax(dat2$base, dat2$age50)
dat2$end_date <- pmin(dat2$end_date, dat2$age65)

# Make exclusions related to data quality and calendar time
if (period=="all") {
  dat2 <- dat2 %>%
    mutate(start_year = year(baseline),
           # Admin censor date could vary
           admin_censor = case_when(state=="MA" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="OH" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="GA" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="PA" & start_year<2013 ~ ymd("2012-12-31"),
                                    T ~ ymd("2019-12-31"))) %>%
    filter(!(state=="MA" & start_year==2011)) %>%
    filter(!(state=="OH" & start_year==2011)) %>%
    filter(!(state=="GA" & (start_year %in% c(2011, 2012)))) %>%
    filter(!(state=="PA" & start_year==2013))

} else if (period=="pre2012") {
  dat2 <- dat2 %>%
    mutate(start_year = year(baseline),
           # Admin censor date could vary
           admin_censor = case_when(state=="MA" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="OH" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="GA" & start_year<2011 ~ ymd("2010-12-31"),
                                    T ~ ymd("2011-12-31"))) %>%
    filter(!(state=="MA" & start_year==2011)) %>%
    filter(!(state=="OH" & start_year==2011)) %>%
    filter(!(state=="GA" & start_year==2011))

} else if (period=="post2012") {
  dat2 <- dat2 %>%
    mutate(start_year = year(baseline),
           # Admin censor date could vary
           admin_censor = case_when(state=="MA" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="OH" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="GA" & start_year<2011 ~ ymd("2010-12-31"),
                                    state=="PA" & start_year<2013 ~ ymd("2012-12-31"),
                                    T ~ ymd("2019-12-31"))) %>%
    filter(!(state=="MA" & start_year==2011)) %>%
    filter(!(state=="OH" & start_year==2011)) %>%
    filter(!(state=="GA" & (start_year %in% c(2011, 2012)))) %>%
    filter(!(state=="PA" & start_year==2013))

  # Allow people to enter in 2012 unless they are from one of the states with data quality issues at that time
  dat2$baseline <- pmax(dat2$baseline, ymd("2012-01-01"))

}

# Recode censoring date as needed (this differs from incidence analysis bc we don't censor at any cancer, just breast cx)
dat2$end_date <- pmin(dat2$end_date, dat2$breast_cx_type2_dt, dat2$dod,
                      dat2$admin_censor, na.rm=T)

# Make final exclusions and create variables necessary for long data
dat3 <- dat2 %>%
  mutate(end_year = year(end_date),
         hiv_base = ifelse(is.na(hiv_type2_dt), 0, ifelse(hiv_type2_dt<baseline, 1, 0)),
         age_base = floor((dob %--% baseline)/years(1)),
         yr_base = year(baseline)) %>%
  # Require event time
  filter(end_date>baseline) %>%
  # Exclude those with prevalent breast cancer
  filter(is.na(breast_cx_type1_dt) | breast_cx_type1_dt>baseline)

# Sample PWoH to improve computational efficiency
set.seed(123)
p1 <- 1
p0 <- 0.25

dat.hiv <- filter(dat3, hiv_base==1) %>%
  mutate(wt_s = 1/p1)
dat.nohiv <- filter(dat3, hiv_base==0)[rbinom(nrow(filter(dat3, hiv_base==0)), 1, p0)==1, ] %>%
  mutate(wt_s = 1/p0)

# Determine number of records per person
samp <- bind_rows(dat.hiv, dat.nohiv) %>%
  mutate(tot_pt = (baseline %--% end_date)/days(90),
         tot_int = ceiling(tot_pt))


# Manage procedures -------------------------------------------------------

screen.mamm <- proc %>%
  filter(bene_id %in% samp$bene_id) %>%
  left_join(select(samp, bene_id, breast_cx_type2_dt, baseline, end_date, hiv_base), by="bene_id") %>%
  group_by(bene_id) %>%
  mutate(last_mamm = lag(proc_date)) %>%
  filter(proc_date>=baseline) %>%
  filter(proc_date<=end_date) %>%
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

# Which interval did procedure occur during
screen.mamm3 <- screen.mamm2 %>%
  mutate(int = ceiling((baseline %--% proc_date)/days(90))) %>%
  select(bene_id, int, proc_date)


# Make data long ----------------------------------------------------------

long <- bind_rows(lapply(samp, rep, samp$tot_int)) %>%
    group_by(bene_id) %>%
    mutate(int = row_number()) %>%
    ungroup()

# Merge procedures and build variables
long2 <- left_join(long, screen.mamm3, by=c("bene_id", "int")) %>%
  mutate(int_start = baseline + days((int-1)*90),
         int_end = as.Date(ifelse(int==tot_int, end_date, int_start + days(89)), origin=ymd("1970-01-01")),

         # Code outcome
         had_mamm = as.numeric(!is.na(proc_date)),

         # Code drop out (follow-up ends prior to age 65 or admin censor date)
         drop = case_when(int<tot_int ~ 0,
                          end_date<age65 ~ 1,
                          end_date<admin_censor ~ 1,
                          T ~ 0),

         # Comorbidities
                  myocard_infarct = as.numeric(myocard_infarct_typ1_dt<=int_start &
                                        !is.na(myocard_infarct_typ1_dt)),
         cong_heart_fail = as.numeric(cong_heart_fail_typ1_dt<=int_start &
                                        !is.na(cong_heart_fail_typ1_dt)),
         periph_vas_dis = as.numeric(periph_vas_dis_typ1_dt<=int_start &
                                       !is.na(periph_vas_dis_typ1_dt)),
         cereb_vas_dis = as.numeric(cereb_vas_dis_typ1_dt<=int_start &
                                      !is.na(cereb_vas_dis_typ1_dt)),
         dementia = as.numeric(dementia_typ1_dt<=int_start & !is.na(dementia_typ1_dt)),
         chron_pulm_dis = as.numeric(chron_pulm_dis_typ1_dt<=int_start &
                                       !is.na(chron_pulm_dis_typ1_dt)),
         rheumatic_dis = as.numeric(rheumatic_dis_typ1_dt<=int_start &
                                      !is.na(rheumatic_dis_typ1_dt)),
         peptic_ulcer = as.numeric(peptic_ulcer_typ1_dt<=int_start &
                                     !is.na(peptic_ulcer_typ1_dt)),
         mild_liv_dis = as.numeric(mild_liv_dis_typ1_dt<=int_start &
                                     !is.na(mild_liv_dis_typ1_dt)),
         mod_sevr_liver_dis = as.numeric(mod_sevr_liver_dis_typ1_dt<=int_start &
                                           !is.na(mod_sevr_liver_dis_typ1_dt)),
         diabt_wo_chron_complict = as.numeric(diabt_wo_chron_complict_typ1_dt<=int_start &
                                                !is.na(diabt_wo_chron_complict_typ1_dt)),
         diabt_chron_complict = as.numeric(diabt_chron_complict_da_typ1_dt<=int_start &
                                             !is.na(diabt_chron_complict_da_typ1_dt)),
         hemiplegia = as.numeric(hemiplegia_paraplegia_d_typ1_dt<=int_start &
                                   !is.na(hemiplegia_paraplegia_d_typ1_dt)),
         renal_dis = as.numeric(renal_dis_typ1_dt<=int_start &
                                  !is.na(renal_dis_typ1_dt)),
         n_cond = myocard_infarct + cong_heart_fail + periph_vas_dis + cereb_vas_dis +
                  dementia + chron_pulm_dis + rheumatic_dis + peptic_ulcer + mild_liv_dis +
                  mod_sevr_liver_dis + diabt_wo_chron_complict + diabt_chron_complict +
                  hemiplegia + renal_dis,
         n_cond = case_when(n_cond==0 ~ 0,
                            n_cond==1 ~ 1,
                            T ~ 2)) %>%
  select(bene_id, wt_s, int, int_start, int_end, hiv_base, race, state, region, age_base,
         yr_base, n_cond, drop, had_mamm)

long3 <- long2 %>%
  group_by(bene_id) %>%
  mutate(n_cond_1 = lag(n_cond),
         n_cond_1 = ifelse(is.na(n_cond_1), n_cond, n_cond_1),
         n_cond_2 = lag(n_cond, n=2),
         n_cond_2 = ifelse(is.na(n_cond_2), n_cond, n_cond_2)) %>%
  ungroup()


# Output data -------------------------------------------------------------

if (period=="all") {
  write_csv(long3, "/cms01/data/dua/57285/users/c-jrudolp9-57285/mamm_long.csv")

} else {
  write_csv(long3, paste0("/cms01/data/dua/57285/users/c-jrudolp9-57285/mamm_long_", period, ".csv"))

}
 
