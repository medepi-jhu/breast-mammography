
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Create datasets for analyses requiring 1 record per person
#
# Last Update: 24 Jul 2025
#
################################################################################


library("tidyverse")
library("duckplyr")
library("splines")
library("parallel")

outcome <- "screen"
group <- "overall" # "overall", "black", "white"

boot_start <- 501
boot_end <- 1000
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

proc <- read_csv(file="/cms01/data/dua/57285/medican/data_requests/jackie/20250425_jackie_mammography_v1/porc.csv",
                 col_types=c("cDddddd")) %>%
  arrange(bene_id, proc_date) %>%
  mutate(year = year(proc_date))


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

# Recode baseline and end date according to age ranges and censoring date
dat2$baseline <- pmax(dat2$base, dat2$age50)

dat2$end_date <- pmin(dat2$end_date, dat2$age65, dat2$dod,
                      dat2$breast_cx_type2_dt, na.rm=T)

if (group=="overall") {
  dat2 <- dat2

} else if (group %in% c("black", "B_MW", "B_NE", "B_S", "B_W")) {
  dat2 <- filter(dat2, race==2)

  if (group=="B_MW") {
    dat2 <- filter(dat2, region=="Midwest")
  } else if (group=="B_NE") {
    dat2 <- filter(dat2, region=="Northeast")
  } else if (group=="B_S") {
    dat2 <- filter(dat2, region=="South")
  } else if (group=="B_W") {
    dat2 <- filter(dat2, region=="West")
  }

} else if (group %in% c("white", "W_MW", "W_NE", "W_S", "W_W")) {
  dat2 <- filter(dat2, race==1)

  if (group=="W_MW") {
    dat2 <- filter(dat2, region=="Midwest")
  } else if (group=="W_NE") {
    dat2 <- filter(dat2, region=="Northeast")
  } else if (group=="W_S") {
    dat2 <- filter(dat2, region=="South")
  } else if (group=="W_W") {
    dat2 <- filter(dat2, region=="West")
  }

} else if (group=="MW") {
  dat2 <- filter(dat2, region=="Midwest")
} else if (group=="NE") {
  dat2 <- filter(dat2, region=="Northeast")
} else if (group=="S") {
  dat2 <- filter(dat2, region=="South")
} else if (group=="W") {
  dat2 <- filter(dat2, region=="West")
}


# Set up annual data -----------------------------------------------------

est_prob <- function(data, yr) {

  # Make data quality exclusions for certain years
  if (yr %in% c(2011, 2012)) {
    data2 <- filter(data, !(state %in% c("MA", "OH", "GA")))
  } else if (yr==2013) {
    data2 <- filter(data, !(state %in% c("GA", "PA")))
  } else if (yr==2014) {
    data2 <- filter(data, state!="PA")
  } else {
    data2 <- data
  }

# Restrict to those eligible: enrolled on Dec 31 of year and no history of breast cancer
  # Build exposure and variables to control for
  data2 <- data2 %>%
    filter(baseline<=ymd(paste0(yr, "-12-31"))) %>%
    filter(end_date>=ymd(paste0(yr, "-12-31"))) %>%
    filter((is.na(breast_cx_type1_dt) | breast_cx_type1_dt>ymd(paste0(yr, "-12-31")))) %>%
    mutate(hiv = as.numeric(!is.na(hiv_type2_dt) & hiv_type2_dt<=ymd(paste0(yr, "-12-31"))),
           age = (dob %--% ymd(paste0(yr, "-12-31")))/years(1),
           # Comorbidities
           myocard_infarct = as.numeric(myocard_infarct_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                        !is.na(myocard_infarct_typ1_dt)),
           cong_heart_fail = as.numeric(cong_heart_fail_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                        !is.na(cong_heart_fail_typ1_dt)),
           periph_vas_dis = as.numeric(periph_vas_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                       !is.na(periph_vas_dis_typ1_dt)),
           cereb_vas_dis = as.numeric(cereb_vas_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                      !is.na(cereb_vas_dis_typ1_dt)),
           dementia = as.numeric(dementia_typ1_dt<=ymd(paste0(yr, "-12-31")) & !is.na(dementia_typ1_dt)),
           chron_pulm_dis = as.numeric(chron_pulm_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                       !is.na(chron_pulm_dis_typ1_dt)),
           rheumatic_dis = as.numeric(rheumatic_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                      !is.na(rheumatic_dis_typ1_dt)),
           peptic_ulcer = as.numeric(peptic_ulcer_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                     !is.na(peptic_ulcer_typ1_dt)),
           mild_liv_dis = as.numeric(mild_liv_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                     !is.na(mild_liv_dis_typ1_dt)),
           mod_sevr_liver_dis = as.numeric(mod_sevr_liver_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                           !is.na(mod_sevr_liver_dis_typ1_dt)),
           diabt_wo_chron_complict = as.numeric(diabt_wo_chron_complict_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                                !is.na(diabt_wo_chron_complict_typ1_dt)),
           diabt_chron_complict = as.numeric(diabt_chron_complict_da_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                             !is.na(diabt_chron_complict_da_typ1_dt)),
           hemiplegia = as.numeric(hemiplegia_paraplegia_d_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                   !is.na(hemiplegia_paraplegia_d_typ1_dt)),
           renal_dis = as.numeric(renal_dis_typ1_dt<=ymd(paste0(yr, "-12-31")) &
                                  !is.na(renal_dis_typ1_dt)),
           n_cond = myocard_infarct + cong_heart_fail + periph_vas_dis + cereb_vas_dis +
                    dementia + chron_pulm_dis + rheumatic_dis + peptic_ulcer + mild_liv_dis +
                    mod_sevr_liver_dis + diabt_wo_chron_complict + diabt_chron_complict +
                    hemiplegia + renal_dis,
           n_cond = case_when(n_cond==0 ~ 0,
                              n_cond==1 ~ 1,
                              T ~ 2)) %>%
    select(-c(myocard_infarct, cong_heart_fail, periph_vas_dis, cereb_vas_dis, dementia, chron_pulm_dis,
              rheumatic_dis, peptic_ulcer, mild_liv_dis, mod_sevr_liver_dis, diabt_wo_chron_complict,
              diabt_chron_complict, hemiplegia, renal_dis))


  # Bootstrap resample ------------------------------------------------------

boot_rep <- function(r) {

  set.seed(r)

  samp <- table(data2[sample(1:nrow(data2), nrow(data2), replace=T), (names(data2)=="bene_id")])

  boot <- NULL
  if (r==0) {
    boot <- data2 %>%
      mutate(bid = bene_id)
  } else {
    for(zzz in 1:max(samp)){
      cc <- data2[data2$bene_id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$bene_id, zzz)
      boot <- rbind(boot, cc)
    }

  }

  # Further require 2 prior years of follow-up, to be consistent with literature
  boot2 <- boot %>%
    filter(((baseline %--% ymd(paste0(yr, "-12-31")))/years(1)) >=2) %>%
    mutate(baseline2 = ymd(paste0(yr-2, "-12-31")))

  unique <- filter(boot2, !duplicated(bene_id))


# Merge procedures --------------------------------------------------------

  if (outcome=="screen") {

  screen.mamm <- proc %>%
    filter(bene_id %in% unique$bene_id) %>%
    left_join(select(unique, bene_id, breast_cx_type2_dt, baseline2), by="bene_id") %>%
    group_by(bene_id) %>%
    mutate(last_mamm = lag(proc_date)) %>%
    filter(proc_date>baseline2) %>%
    filter(proc_date<=ymd(paste0(yr, "-12-31"))) %>%
    ungroup()

  screen.mamm <- screen.mamm %>%
    mutate(mamm_9month = case_when(is.na(last_mamm) ~ 0,
                                   ((last_mamm %--% proc_date)/months(1))>=9 ~ 0,
                                   ((last_mamm %--% proc_date)/months(1))<9 ~ 1,
                                   T ~ NA_real_),
           diag_1year = case_when(is.na(breast_cx_type2_dt) ~ 0,
                                  breast_cx_type2_dt>proc_date ~ 0,
                                  ((breast_cx_type2_dt %--% proc_date)/years(1))>=1 ~ 0,
                                  ((breast_cx_type2_dt %--% proc_date)/years(1))<1 ~ 1,
                                  T ~ NA_real_)) %>%
    filter(mammograph_typ==2 & mamm_9month==0 & diag_1year==0) %>%
    filter(!duplicated(bene_id)) %>%
    rename(screen_date = proc_date)

  boot2 <- boot2 %>%
    left_join(select(screen.mamm, bene_id, screen_date), by="bene_id")

  } else {

  all.mamm <- proc %>%
    filter(bene_id %in% unique$bene_id) %>%
    left_join(select(unique, bene_id, baseline2), by="bene_id") %>%
    filter(proc_date>baseline2) %>%
    filter(proc_date<=ymd(paste0(yr, "-12-31"))) %>%
    filter(!duplicated(bid)) %>%
    rename(screen_date = proc_date)

  boot2 <- boot2 %>%
    left_join(select(all.mamm, bene_id, screen_date), by="bene_id")

  }


# Estimate probability of screening ---------------------------------------

  # Crude
  summ.crude <- boot2 %>%
    mutate(screen = !is.na(screen_date)) %>%
    group_by(hiv) %>%
    summarize(n = n(),
              p_screened = sum(screen)/n()) %>%
    mutate(analysis = "crude")

  # Weights to balance covariates by HIV status
  if (group=="overall") {

    p.h1 <- glm(hiv ~ as.factor(race) + state + bs(age, df=3),
                data=boot2, family=binomial(link="logit"))$fitted.values
    boot2$wt_hiv1 <- ifelse(boot2$hiv==1, 1, p.h1/(1-p.h1))

    p.h2 <- glm(hiv ~ as.factor(race) + state + bs(age, df=3) + as.factor(n_cond),
                data=boot2, family=binomial(link="logit"))$fitted.values
    boot2$wt_hiv2 <- ifelse(boot2$hiv==1, 1, p.h2/(1-p.h2))

  } else if (group %in% c("black", "white")) {

    p.h1 <- glm(hiv ~ state + bs(age, df=3),
                data=boot2, family=binomial(link="logit"))$fitted.values
    boot2$wt_hiv1 <- ifelse(boot2$hiv==1, 1, p.h1/(1-p.h1))

    p.h2 <- glm(hiv ~ state + bs(age, df=3) + as.factor(n_cond),
                data=boot2, family=binomial(link="logit"))$fitted.values
    boot2$wt_hiv2 <- ifelse(boot2$hiv==1, 1, p.h2/(1-p.h2))

  }

  # Trim extreme weights
  p1 <- quantile(boot2$wt_hiv1, p=0.01)
  p99 <- quantile(boot2$wt_hiv1, p=0.99)

  boot2$wt_hiv1 <- ifelse(boot2$wt_hiv1<p1, p1, boot2$wt_hiv1)
  boot2$wt_hiv1 <- ifelse(boot2$wt_hiv1>p99, p99, boot2$wt_hiv1)

  p1 <- quantile(boot2$wt_hiv2, p=0.01)
  p99 <- quantile(boot2$wt_hiv2, p=0.99)

  boot2$wt_hiv2 <- ifelse(boot2$wt_hiv2<p1, p1, boot2$wt_hiv2)
  boot2$wt_hiv2 <- ifelse(boot2$wt_hiv2>p99, p99, boot2$wt_hiv2)

  
  summ.hiv1 <- boot2 %>%
    mutate(screen = !is.na(screen_date)) %>%
    group_by(hiv) %>%
    summarize(n = n(),
              p_screened = sum(screen*wt_hiv1)/sum(wt_hiv1)) %>%
    mutate(analysis = "ipw1")

  summ.hiv2 <- boot2 %>%
    mutate(screen = !is.na(screen_date)) %>%
    group_by(hiv) %>%
    summarize(n = n(),
              p_screened = sum(screen*wt_hiv2)/sum(wt_hiv2)) %>%
    mutate(analysis = "ipw2")

  # Estimate weights to make those with 2 years of prior follow-up look like
  # all those eligible
  boot$included <- as.numeric(boot$bid %in% boot2$bid)

  if (group=="overall") {

    p.s <- glm(included ~ hiv + as.factor(race) + state + bs(age, df=3) + as.factor(n_cond),
               data=boot, family=binomial(link="logit"))$fitted.values
    boot2$wt_s <- ifelse(boot$included==1, 1/p.s, 0)[boot$included==1]

    boot2$wt1 <- boot2$wt_hiv1*boot2$wt_s
    boot2$wt2 <- boot2$wt_hiv2*boot2$wt_s

  } else if (group %in% c("black", "white")) {

    p.s <- glm(included ~ hiv + state + bs(age, df=3) + as.factor(n_cond),
               data=boot, family=binomial(link="logit"))$fitted.values
    boot2$wt_s <- ifelse(boot$included==1, 1/p.s, 0)[boot$included==1]

    boot2$wt1 <- boot2$wt_hiv1*boot2$wt_s
    boot2$wt2 <- boot2$wt_hiv2*boot2$wt_s

  }

# Trim extreme weights
  p1 <- quantile(boot2$wt1, p=0.01)
  p99 <- quantile(boot2$wt1, p=0.99)

  boot2$wt1 <- ifelse(boot2$wt1<p1, p1, boot2$wt1)
  boot2$wt1 <- ifelse(boot2$wt1>p99, p99, boot2$wt1)

  p1 <- quantile(boot2$wt2, p=0.01)
  p99 <- quantile(boot2$wt2, p=0.99)

  boot2$wt2 <- ifelse(boot2$wt2<p1, p1, boot2$wt2)
  boot2$wt2 <- ifelse(boot2$wt2>p99, p99, boot2$wt2)

  summ.select1 <- boot2 %>%
      mutate(screen = !is.na(screen_date)) %>%
      group_by(hiv) %>%
      summarize(n = n(),
                p_screened = sum(screen*wt1)/sum(wt1)) %>%
      mutate(analysis = "ipw1+ipsw")

  summ.select2 <- boot2 %>%
    mutate(screen = !is.na(screen_date)) %>%
    group_by(hiv) %>%
    summarize(n = n(),
              p_screened = sum(screen*wt2)/sum(wt2)) %>%
    mutate(analysis = "ipw2+ipsw")

  est <- bind_rows(summ.crude, summ.hiv1, summ.hiv2, summ.select1, summ.select2) %>%
    mutate(boot_rep = r)

  return(est)

}

#all.est <- lapply(boot_start:boot_end, function(x){boot_rep(x)})
all.est <- mclapply(boot_start:boot_end, function(x){boot_rep(x)}, mc.set.seed=F, mc.cores=cores)
all.est <- bind_rows(all.est) %>%
  mutate(year=yr)

#res <- filter(all.est, boot_rep==0)

#res.summ <- all.est %>%
#  group_by(analysis, hiv) %>%
#  summarize(p_screened_ll = quantile(p_screened, p=0.025),
#            p_screened_ul = quantile(p_screened, p=0.975))


#res <- left_join(res, res.summ, by=c("analysis", "hiv")) %>%
#  mutate(year = yr)

return(all.est)

}

res <- lapply(2003:2019, function(x){est_prob(data=dat2, yr=x)})
res <- bind_rows(res)

if (outcome=="screen") {
  write_csv(res, paste0("./results/prob2yr_", group, "_", boot_end, ".csv"))
} else {
  write_csv(res, paste0("./results/prob2yr_", group, "_", boot_end, "_total.csv"))
}
