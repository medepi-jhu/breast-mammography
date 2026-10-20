
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Estimate breast cancer detection rate
#
# Last Update: 24 Jul 2025
#
################################################################################


library("tidyverse")
library("splines")
library("parallel")

group <- "white" # "overall", "black", "white"
period <- "post2012" # "all", "pre2012", "post2012"

nboot <- 1000
cores <- 10


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

mamm <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/porc.csv",
                 col_types=c("cDddddd")) %>%
  arrange(bene_id, proc_date) %>%
  mutate(year = year(proc_date))

surg <- read_csv("/cms01/data/dua/57285/medican/data_requests/jackie/20250714_jackie_mammography_v2/surgery.csv") %>%
  arrange(bene_id, surgery_dt) %>%
  group_by(bene_id) %>%
  filter(!duplicated(surgery_dt)) %>%
  ungroup()

diag <- read_csv("/cms01/data/dua/57285/medican/data_requests/jackie/20250714_jackie_mammography_v2/long_diag.csv") %>%
  arrange(bene_id, diag_dt) %>%
  group_by(bene_id) %>%
  filter(!duplicated(diag_dt)) %>%
  ungroup()


# Manage data -------------------------------------------------------------

dat2 <- dat %>%
  mutate(base = start_date + months(6),

         # Recode leap year birthdays
         dob = ifelse(month(dob)==2 & day(dob)==29, dob - days(1), dob),
         dob = as.Date(dob, origin=ymd("1970-01-01")),

         age50 = dob + years(50),
               age65 = dob + years(65)) %>%
  # Exclude MD (due to breast cx and mammography data quality) and all new states
  filter(state %in% c("AL", "CA", "CO", "FL", "GA", "IL", "MA", "NC", "NY", "OH", "PA", "TX", "WA"))

if (group=="overall") {
  dat2 <- dat2
} else if (group=="black") {
  dat2 <- filter(dat2, race==2)
} else if (group=="white") {
  dat2 <- filter(dat2, race==1)
}


# Recode baseline and censoring date according to age ranges
dat2$baseline <- pmax(dat2$base, dat2$age50)
dat2$censor_date <- pmin(dat2$end_date, dat2$age65)

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

# Recode censoring date as needed (this differs from other analyses bc no censoring by cancer)
dat2$censor_date <- pmin(dat2$censor_date, dat2$dod, dat2$admin_censor, na.rm=T)

# Make final exclusions
dat3 <- dat2 %>%
  # Require follow-up time
  filter(censor_date>baseline) %>%
  # Exclude those with prevalent breast cancer at enrollment baseline
  filter(is.na(breast_cx_type1_dt) | breast_cx_type1_dt>baseline)


# Identify screening mammography ------------------------------------------

# Pull eligible mammograms
screen.mamm <- mamm %>%
  filter(bene_id %in% dat3$bene_id) %>%
  left_join(dat3, by="bene_id") %>%
  group_by(bene_id) %>%
  mutate(last_mamm = lag(proc_date)) %>%
  filter(proc_date>=baseline) %>%
  filter(proc_date<=censor_date) %>%
  ungroup()

# Limit to screening mammograms
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


# Screen detected breast cancer -------------------------------------------

screen.mamm3 <- screen.mamm2 %>%
  mutate(age_screen = (dob %--% proc_date)/years(1)) %>%
  # For main analysis, pull first screening mammography
  filter(!duplicated(bene_id)) %>%
  # Exclude those with prevalent breast cancer at screening date
  filter(is.na(breast_cx_type1_dt) | breast_cx_type1_dt>proc_date) %>%
  # Need to have potential for 1 year of follow-up
  filter(age_screen<=64) %>%
  filter((proc_date %--% admin_censor)/years(1) >=1)

# Did they have a breast cx diagnosis within 123 days?
screen.mamm3 <- screen.mamm3 %>%
  mutate(# Only count if within 1st eligibility period
         breast_cx_type2_dt = ifelse(breast_cx_type2_dt>censor_date, NA, breast_cx_type2_dt),
         breast_cx_type2_dt = as.Date(breast_cx_type2_dt, origin=ymd("1970-01-01")),
         t_cancer = (proc_date %--% breast_cx_type2_dt)/days(1),
         breast_cx_123 = ifelse(is.na(t_cancer), 0, as.numeric(t_cancer<123)))

# Did they have a carcinoma in-situ within 286 days?
diag2 <- diag %>%
  filter(bene_id %in% screen.mamm3$bene_id) %>%
  left_join(select(dat3, bene_id, censor_date), by="bene_id") %>%
  filter(diag_dt<=censor_date) %>% # Only count if within 1st eligibility period
  left_join(select(screen.mamm3, bene_id, proc_date), by="bene_id") %>%
  filter(diag_dt>=proc_date) %>%
  filter(!duplicated(bene_id))

screen.mamm3 <- screen.mamm3 %>%
  left_join(select(diag2, bene_id, diag_dt), by="bene_id") %>%
  mutate(t_insitu = (proc_date %--% diag_dt)/days(1),
         insitu_286 = ifelse(is.na(t_insitu), 0, as.numeric(t_insitu<286)))

# Did they have a subsequent mammogram within 82 days?
mamm2 <- mamm %>%
  rename(mamm_dt = proc_date) %>%
  filter(bene_id %in% screen.mamm3$bene_id) %>%
  left_join(select(dat3, bene_id, censor_date), by="bene_id") %>%
  filter(mamm_dt<=censor_date) %>% # Only count if within 1st eligibility period
  left_join(select(screen.mamm3, bene_id, proc_date), by="bene_id") %>%
  filter(mamm_dt>proc_date) %>%
  filter(!duplicated(bene_id))

screen.mamm3 <- screen.mamm3 %>%
  left_join(select(mamm2, bene_id, mamm_dt), by="bene_id") %>%
  mutate(t_mamm = (proc_date %--% mamm_dt)/days(1),
         mamm_82 = ifelse(is.na(t_mamm), 0, as.numeric(t_mamm<82)))

# Did they have a breast-directed surgery within 365 days?
surg2 <- surg %>%
  filter(bene_id %in% screen.mamm3$bene_id) %>%
  left_join(select(dat3, bene_id, censor_date), by="bene_id") %>%
  filter(surgery_dt<=censor_date) %>% # Only count if within 1st eligibility period
  left_join(select(screen.mamm3, bene_id, proc_date), by="bene_id") %>%
  filter(surgery_dt>=proc_date) %>%
  filter(!duplicated(bene_id))

screen.mamm3 <- screen.mamm3 %>%
  left_join(select(surg2, bene_id, surgery_dt), by="bene_id") %>%
  mutate(t_surg = (proc_date %--% surgery_dt)/days(1),
         surg_365 = ifelse(is.na(t_surg), 0, as.numeric(t_surg<365)))

# Did they die or drop out before the end of the year?
screen.mamm3 <- screen.mamm3 %>%
  mutate(t_death = (proc_date %--% dod)/days(1),
         died_365 = ifelse(is.na(t_death), 0, as.numeric(t_death<365)),
         t_drop = (proc_date %--% end_date)/days(1),
         drop_365 = ifelse(is.na(t_drop), 0, as.numeric(t_drop<365)))

# Build variables necessary for analysis
screen.mamm4 <- screen.mamm3 %>%
  mutate(# Outcome
         delta = case_when(breast_cx_123==1 & surg_365==1 ~ "detected",
                           breast_cx_123==0 & insitu_286==1 & mamm_82==1 ~ "detected",
                           died_365==1 ~ "died",
                           drop_365==1 ~ "dropped",
                           T ~ "no event"),
         # HIV status
         hiv = as.numeric(!is.na(hiv_type2_dt) & hiv_type2_dt<=proc_date),

         # Comorbidities
         myocard_infarct = as.numeric(myocard_infarct_typ1_dt<=proc_date &
                                        !is.na(myocard_infarct_typ1_dt)),
         cong_heart_fail = as.numeric(cong_heart_fail_typ1_dt<=proc_date &
                                        !is.na(cong_heart_fail_typ1_dt)),
         periph_vas_dis = as.numeric(periph_vas_dis_typ1_dt<=proc_date &
                                       !is.na(periph_vas_dis_typ1_dt)),
         cereb_vas_dis = as.numeric(cereb_vas_dis_typ1_dt<=proc_date &
                                      !is.na(cereb_vas_dis_typ1_dt)),
         dementia = as.numeric(dementia_typ1_dt<=proc_date & !is.na(dementia_typ1_dt)),
         chron_pulm_dis = as.numeric(chron_pulm_dis_typ1_dt<=proc_date &
                                       !is.na(chron_pulm_dis_typ1_dt)),
         rheumatic_dis = as.numeric(rheumatic_dis_typ1_dt<=proc_date &
                                      !is.na(rheumatic_dis_typ1_dt)),
         peptic_ulcer = as.numeric(peptic_ulcer_typ1_dt<=proc_date &
                                     !is.na(peptic_ulcer_typ1_dt)),
         mild_liv_dis = as.numeric(mild_liv_dis_typ1_dt<=proc_date &
                                     !is.na(mild_liv_dis_typ1_dt)),
         mod_sevr_liver_dis = as.numeric(mod_sevr_liver_dis_typ1_dt<=proc_date &
                                           !is.na(mod_sevr_liver_dis_typ1_dt)),
         diabt_wo_chron_complict = as.numeric(diabt_wo_chron_complict_typ1_dt<=proc_date &
                                                !is.na(diabt_wo_chron_complict_typ1_dt)),
         diabt_chron_complict = as.numeric(diabt_chron_complict_da_typ1_dt<=proc_date &
                                             !is.na(diabt_chron_complict_da_typ1_dt)),
         hemiplegia = as.numeric(hemiplegia_paraplegia_d_typ1_dt<=proc_date &
                                   !is.na(hemiplegia_paraplegia_d_typ1_dt)),
         renal_dis = as.numeric(renal_dis_typ1_dt<=proc_date &
                                  !is.na(renal_dis_typ1_dt)),
         n_cond = myocard_infarct + cong_heart_fail + periph_vas_dis + cereb_vas_dis +
           dementia + chron_pulm_dis + rheumatic_dis + peptic_ulcer + mild_liv_dis +
           mod_sevr_liver_dis + diabt_wo_chron_complict + diabt_chron_complict +
           hemiplegia + renal_dis,
         n_cond = case_when(n_cond==0 ~ 0,
                            n_cond==1 ~ 1,
                            T ~ 2)) %>%
  select(bene_id, hiv, delta, race, state, age_screen, n_cond)


# Bootstrap resample ------------------------------------------------------

boot_rep <- function(r) {

  set.seed(r)

  samp <- table(screen.mamm4[sample(1:nrow(screen.mamm4), nrow(screen.mamm4), replace=T), (names(screen.mamm4)=="bene_id")])

  boot <- NULL
  if (r==0) {
    boot <- screen.mamm4 %>%
      rename(bid = bene_id)
  } else {
    for(zzz in 1:max(samp)){
      cc <- screen.mamm4[screen.mamm4$bene_id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$bene_id, zzz)
      boot <- rbind(boot, cc)
    }
    boot <- select(boot, -bene_id)
  }


# Breast cancer detection rate --------------------------------------------

# Drop those with missing outcome (we'll weight later to deal with them)
selected <- filter(boot, !(delta %in% c("dropped", "died"))) %>%
  mutate(detect = as.numeric(delta=="detected"))

# Crude
crude <- selected %>%
  group_by(hiv) %>%
  summarize(rate = (sum(detect)/n())*1000) %>%
  mutate(analysis = "crude")

# Weighted
if (group=="overall") {
  # Exposure weights
  p.h1 <- glm(hiv ~ as.factor(race) + state + bs(age_screen, df=3),
               data=selected, family=binomial(link="logit"))$fitted.values
  selected$wt_hiv1 <- ifelse(selected$hiv==1, 1, p.h1/(1-p.h1))

  p.h2 <- glm(hiv ~ as.factor(race) + state + bs(age_screen, df=3) + as.factor(n_cond),
              data=selected, family=binomial(link="logit"))$fitted.values
  selected$wt_hiv2 <- ifelse(selected$hiv==1, 1, p.h2/(1-p.h2))

  # Selection weights
  boot$included <- as.numeric(boot$bid %in% selected$bid)

  p.s <- glm(included ~ hiv + as.factor(race) + state + bs(age_screen, df=3) + as.factor(n_cond),
                 data=boot, family=binomial(link="logit"))$fitted.values
  selected$wt_s <- ifelse(boot$included==1, 1/p.s, 0)[boot$included==1]

  selected$wt_tot1 <- selected$wt_hiv1*selected$wt_s
  selected$wt_tot2 <- selected$wt_hiv2*selected$wt_s

} else if (group %in% c("black", "white")) {
  # Exposure weights
  p.h1 <- glm(hiv ~ state + bs(age_screen, df=3),
              data=selected, family=binomial(link="logit"))$fitted.values
  selected$wt_hiv1 <- ifelse(selected$hiv==1, 1, p.h1/(1-p.h1))

  p.h2 <- glm(hiv ~ state + bs(age_screen, df=3) + as.factor(n_cond),
              data=selected, family=binomial(link="logit"))$fitted.values
  selected$wt_hiv2 <- ifelse(selected$hiv==1, 1, p.h2/(1-p.h2))

  # Selection weights
  boot$included <- as.numeric(boot$bid %in% selected$bid)

  p.s <- glm(included ~ hiv + state + bs(age_screen, df=3) + as.factor(n_cond),
             data=boot, family=binomial(link="logit"))$fitted.values
  selected$wt_s <- ifelse(boot$included==1, 1/p.s, 0)[boot$included==1]

  selected$wt_tot1 <- selected$wt_hiv1*selected$wt_s
  selected$wt_tot2 <- selected$wt_hiv2*selected$wt_s

}

# Trim extreme weights
p1 <- quantile(selected$wt_hiv1, p=0.01)
p99 <- quantile(selected$wt_hiv1, p=0.99)

selected$wt_hiv1 <- ifelse(selected$wt_hiv1<p1, p1, selected$wt_hiv1)
selected$wt_hiv1 <- ifelse(selected$wt_hiv1>p99, p99, selected$wt_hiv1)

p1 <- quantile(selected$wt_hiv2, p=0.01)
p99 <- quantile(selected$wt_hiv2, p=0.99)

selected$wt_hiv2 <- ifelse(selected$wt_hiv2<p1, p1, selected$wt_hiv2)
selected$wt_hiv2 <- ifelse(selected$wt_hiv2>p99, p99, selected$wt_hiv2)

p1 <- quantile(selected$wt_tot1, p=0.01)
p99 <- quantile(selected$wt_tot1, p=0.99)

selected$wt_tot1 <- ifelse(selected$wt_tot1<p1, p1, selected$wt_tot1)
selected$wt_tot1 <- ifelse(selected$wt_tot1>p99, p99, selected$wt_tot1)

p1 <- quantile(selected$wt_tot2, p=0.01)
p99 <- quantile(selected$wt_tot2, p=0.99)

selected$wt_tot2 <- ifelse(selected$wt_tot2<p1, p1, selected$wt_tot2)
selected$wt_tot2 <- ifelse(selected$wt_tot2>p99, p99, selected$wt_tot2)

wt1 <- selected %>%
  group_by(hiv) %>%
  summarize(rate = (sum(detect*wt_hiv1)/sum(wt_hiv1))*1000) %>%
  mutate(analysis = "ipw1")

wt2 <- selected %>%
  group_by(hiv) %>%
  summarize(rate = (sum(detect*wt_hiv2)/sum(wt_hiv2))*1000) %>%
  mutate(analysis = "ipw2")

wt3 <- selected %>%
  group_by(hiv) %>%
  summarize(rate = (sum(detect*wt_tot1)/sum(wt_tot1))*1000) %>%
  mutate(analysis = "ipw1+ipsw")

wt4 <- selected %>%
  group_by(hiv) %>%
  summarize(rate = (sum(detect*wt_tot2)/sum(wt_tot2))*1000) %>%
  mutate(analysis = "ipw2+ipsw")


# Combine results ---------------------------------------------------------

est <- bind_rows(crude, wt1, wt2, wt3, wt4) %>%
  mutate(rep = r)

n.hiv <- boot %>%
  group_by(hiv) %>%
  summarize(n = n())

est <- left_join(est, n.hiv, by="hiv")

return(est)

}

#all.est <- lapply(0:nboot, function(x){boot_rep(x)})
all.est <- mclapply(0:nboot, function(x){boot_rep(x)}, mc.set.seed=F, mc.cores=cores)
all.est <- bind_rows(all.est)

res <- filter(all.est, rep==0)

res.summ <- all.est %>%
  group_by(analysis, hiv) %>%
  summarize(rate_ll = quantile(rate, p=0.025),
            rate_ul = quantile(rate, p=0.975))

res <- left_join(res, res.summ, by=c("analysis", "hiv"))


# Output results ----------------------------------------------------------

if (period=="all") {
  write_csv(res, paste0("./results/screen_detect_", group, ".csv"))

} else {
  write_csv(res, paste0("./results/screen_detect_", group, "_", period, "_", nboot, ".csv"))

}
