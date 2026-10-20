
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Estimate crude annual rates of mammography
#
# Last Update: 27 Jun 2025
#
################################################################################


library("tidyverse")
library("tableone")

options(dplyr.summarise.inform = FALSE)

samp <- 0


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


# Estimate rates by state and year ----------------------------------------

rates <- function(analysis) {

res <- NULL

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
  if (analysis=="all") {
    dat2$baseline <- dat2$base
    dat2$end_date <- pmin(dat2$end_date, dat2$age65)

  } else if (analysis=="40-50") {
    dat2$baseline <- pmax(dat2$base, dat2$age40)
    dat2$end_date <- pmin(dat2$end_date, dat2$age50)

  } else if (analysis==">50") {
    dat2$baseline <- pmax(dat2$base, dat2$age50)
    dat2$end_date <- pmin(dat2$end_date, dat2$age65)

  }

# Make exclusions related to data quality
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

  # Recode censoring date as needed (this differs from incidence analysis bc we don't censor at any cancer, just breast cx)
  dat2$end_date <- pmin(dat2$end_date, dat2$breast_cx_type2_dt, dat2$dod, dat2$admin_censor, na.rm=T)

  dat3 <- dat2 %>%
    mutate(end_year = year(end_date),
           hiv_base = ifelse(is.na(hiv_type2_dt), 0, ifelse(hiv_type2_dt<baseline, 1, 0))) %>%
    # Require event time
    filter(end_date>baseline) %>%
    # Exclude those with prevalent breast cancer
    filter(is.na(breast_cx_type1_dt) | breast_cx_type1_dt>baseline)
  rm(dat2)

  # Get results for the sample
  if (samp==1) {
    set.seed(123)
    p1 <- 1
    p0 <- 0.25

    dat.hiv <- filter(dat3, hiv_base==1) %>%
      mutate(wt_s = 1/p1)
    dat.nohiv <- filter(dat3, hiv_base==0)[rbinom(nrow(filter(dat3, hiv_base==0)), 1, p0)==1, ] %>%
      mutate(wt_s = 1/p0)

    dat4 <- bind_rows(dat.hiv, dat.nohiv)

  } else {
    dat4 <- dat3
  }
  rm(dat3)

# Make data long
  dat4$r <- dat4$end_year - dat4$start_year + 1
  long <- bind_rows(lapply(dat4, rep, dat4$r)) %>%
    group_by(bene_id) %>%
    mutate(add = row_number() - 1) %>%
    ungroup()

  long2 <- long %>%
    mutate(year = start_year + add,
           person_years = case_when(year==start_year & end_year==year ~ as.numeric(end_date - baseline)/365.25,
                                    year==start_year & end_year>year ~ as.numeric(ymd(paste0(year, "-12-31")) - baseline )/365.25,
                                    year>start_year & end_year==year ~ as.numeric(end_date - ymd(paste0(year, "-01-01")))/365.25,
                                    year>start_year & end_year>year ~ 1,
                                    T ~ NA_real_)) %>%
    filter(person_years>0) %>%
    select(-c(r, add))
  rm(long)


# Estimate rates ---------------------------------------------------------

for (group in c("overall", "black", "white")) {
                #"Midwest", "Northeast", "South", "West"
                #"B_MW", "B_NE", "B_S", "B_W",
                #"W_MW", "W_NE", "W_S", "W_W"))

  if (group=="overall") {

  dat4 <- dat4 %>%
    mutate(age_base = (dob %--% baseline)/years(1),
                      # Comorbidities
           myocard_infarct = as.numeric(myocard_infarct_typ1_dt<=baseline &
                                        !is.na(myocard_infarct_typ1_dt)),
           cong_heart_fail = as.numeric(cong_heart_fail_typ1_dt<=baseline &
                                        !is.na(cong_heart_fail_typ1_dt)),
           periph_vas_dis = as.numeric(periph_vas_dis_typ1_dt<=baseline &
                                       !is.na(periph_vas_dis_typ1_dt)),
           cereb_vas_dis = as.numeric(cereb_vas_dis_typ1_dt<=baseline &
                                      !is.na(cereb_vas_dis_typ1_dt)),
           dementia = as.numeric(dementia_typ1_dt<=baseline & !is.na(dementia_typ1_dt)),
           chron_pulm_dis = as.numeric(chron_pulm_dis_typ1_dt<=baseline &
                                       !is.na(chron_pulm_dis_typ1_dt)),
           rheumatic_dis = as.numeric(rheumatic_dis_typ1_dt<=baseline &
                                      !is.na(rheumatic_dis_typ1_dt)),
           peptic_ulcer = as.numeric(peptic_ulcer_typ1_dt<=baseline &
                                     !is.na(peptic_ulcer_typ1_dt)),
           mild_liv_dis = as.numeric(mild_liv_dis_typ1_dt<=baseline &
                                     !is.na(mild_liv_dis_typ1_dt)),
           mod_sevr_liver_dis = as.numeric(mod_sevr_liver_dis_typ1_dt<=baseline &
                                           !is.na(mod_sevr_liver_dis_typ1_dt)),
           diabt_wo_chron_complict = as.numeric(diabt_wo_chron_complict_typ1_dt<=baseline &
                                                !is.na(diabt_wo_chron_complict_typ1_dt)),
           diabt_chron_complict = as.numeric(diabt_chron_complict_da_typ1_dt<=baseline &
                                             !is.na(diabt_chron_complict_da_typ1_dt)),
           hemiplegia = as.numeric(hemiplegia_paraplegia_d_typ1_dt<=baseline &
                                   !is.na(hemiplegia_paraplegia_d_typ1_dt)),
           renal_dis = as.numeric(renal_dis_typ1_dt<=baseline &
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

  print(CreateTableOne(data=dat4,
                       vars=c("age_base", "race", "state", "region", "start_year", "n_cond"),
                       factorVars=c("race", "start_year", "n_cond"),
                       strata="hiv_base", addOverall=T,
                       test=F),
        contDigits=1, nonnormal="age_base")


  }

  print(paste0("Starting analysis: ", group))

  # Person-time in those without breast cancer
  if (group=="overall") {
    long3 <- long2
    dat5 <- dat4
  } else if (group=="black") {
    long3 <- filter(long2, race==2)
    dat5 <- filter(dat4, race==2)
  } else if (group=="white") {
    long3 <- filter(long2, race==1)
    dat5 <- filter(dat4, race==1)
  } else if (group %in% c("Midwest", "Northeast", "South", "West")) {
    long3 <- filter(long2, region==group)
    dat5 <- filter(dat4, region==group)
  } else if (group=="B_MW") {
    long3 <- long2 %>%
      filter(region=="Midwest") %>%
      filter(race==2)
    dat5 <- dat4 %>%
      filter(region=="Midwest") %>%
      filter(race==2)
  } else if (group=="B_NE") {
    long3 <- long2 %>%
      filter(region=="Northeast") %>%
      filter(race==2)
    dat5 <- dat4 %>%
      filter(region=="Northeast") %>%
      filter(race==2)
   } else if (group=="B_S") {
    long3 <- long2 %>%
      filter(region=="South") %>%
      filter(race==2)
    dat5 <- dat4 %>%
      filter(region=="South") %>%
      filter(race==2)
  } else if (group=="B_W") {
    long3 <- long2 %>%
      filter(region=="West") %>%
      filter(race==2)
    dat5 <- dat4 %>%
      filter(region=="West") %>%
      filter(race==2)
  } else if (group=="W_MW") {
    long3 <- long2 %>%
      filter(region=="Midwest") %>%
      filter(race==1)
    dat5 <- dat4 %>%
      filter(region=="Midwest") %>%
      filter(race==1)
  } else if (group=="W_NE") {
    long3 <- long2 %>%
      filter(region=="Northeast") %>%
      filter(race==1)
    dat5 <- dat4 %>%
      filter(region=="Northeast") %>%
      filter(race==1)
  } else if (group=="W_S") {
    long3 <- long2 %>%
      filter(region=="South") %>%
      filter(race==1)
    dat5 <- dat4 %>%
      filter(region=="South") %>%
      filter(race==1)
  } else if (group=="W_W") {
    long3 <- long2 %>%
      filter(region=="West") %>%
      filter(race==1)
    dat5 <- dat4 %>%
      filter(region=="West") %>%
      filter(race==1)
  }

  # Person-time
  time <- long3 %>%
    group_by(year) %>%
    summarize(tot_py = sum(person_years),
              n = n())

  time.hiv <- long3 %>%
    group_by(hiv_base, year) %>%
    summarize(tot_py = sum(person_years),
              n = n())

  # Rate of all mammographies
  all.mamm <- proc %>%
    filter(bene_id %in% dat5$bene_id) %>%
    left_join(select(dat5, bene_id, baseline, end_date, hiv_base), by="bene_id") %>%
    filter(proc_date>=baseline) %>%
    filter(proc_date<=end_date)

  all.overall <- all.mamm %>%
    group_by(year) %>%
    summarize(n_total = n()) %>%
    mutate(hiv_base = 99)
  all.overall2 <- left_join(time, all.overall, by="year")

  all.hiv <- all.mamm %>%
    group_by(hiv_base, year) %>%
    summarize(n_total = n())
  all.hiv2 <- left_join(time.hiv, all.hiv, by=c("hiv_base", "year"))

  all.mamm2 <- bind_rows(all.overall2, all.hiv2) %>%
    mutate(rate_total = n_total/tot_py,
           se_total = sqrt(n_total/(tot_py**2))) %>%
  select(-tot_py)

  contrast <- data.frame(year = seq(2001, 2019, 1),
                         irr_total = all.mamm2$rate_total[all.mamm2$hiv_base==1]/all.mamm2$rate_total[all.mamm2$hiv_base==0],
                         se_lirr_total = sqrt((1/all.mamm2$n_total[all.mamm2$hiv_base==1]) + (1/all.mamm2$n_total[all.mamm2$hiv_base==0])))

  all.mamm2 <- left_join(all.mamm2, contrast, by="year")

  # Rate of screening mammographies
  screen.mamm <- proc %>%
    filter(bene_id %in% dat5$bene_id) %>%
    left_join(select(dat5, bene_id, breast_cx_type2_dt, baseline, end_date, hiv_base), by="bene_id") %>%
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
                                  T ~ NA_real_))

  if (group=="overall") {
    print("Total eligible procedures:")
    print(nrow(screen.mamm2))

    print("Total with screening codes:")
    print(sum(screen.mamm2$mammograph_typ==2))

    print("Above plus no prior mammography within 9 months:")
    print(sum(screen.mamm2$mammograph_typ==2 & screen.mamm2$mamm_9month==0))

    print("Above plus no breast cancer diagnosis within last year:")
    print(sum(screen.mamm2$mammograph_typ==2 & screen.mamm2$mamm_9month==0 & screen.mamm2$diag_1year==0))

  }

  screen.mamm2 <- screen.mamm2 %>%
    filter(mammograph_typ==2 & mamm_9month==0 & diag_1year==0)

  screen.overall <- screen.mamm2 %>%
    group_by(year) %>%
    summarize(n_screen = n()) %>%
    mutate(hiv_base = 99)
  screen.overall2 <- left_join(time, screen.overall, by="year") %>%
    select(-n) # don't need duplicate columns

  screen.hiv <- screen.mamm2 %>%
    group_by(hiv_base, year) %>%
    summarize(n_screen = n())
  screen.hiv2 <- left_join(time.hiv, screen.hiv, by=c("hiv_base", "year")) %>%
    select(-n)

  screen.mamm3 <- bind_rows(screen.overall2, screen.hiv2) %>%
    mutate(rate_screen = n_screen/tot_py,
           se_screen = sqrt(n_screen/(tot_py**2)))

  contrast <- data.frame(year = seq(2001, 2019, 1),
                         irr_screen = screen.mamm3$rate_screen[screen.mamm3$hiv_base==1]/screen.mamm3$rate_screen[screen.mamm3$hiv_base==0],
                         se_lirr_screen = sqrt((1/screen.mamm3$n_screen[screen.mamm3$hiv_base==1]) + (1/screen.mamm3$n_screen[screen.mamm3$hiv_base==0])))

  screen.mamm3 <- left_join(screen.mamm3, contrast, by="year")

  temp.res <- left_join(all.mamm2, screen.mamm3, by=c("hiv_base", "year")) %>%
    mutate(group = group,
           analysis = analysis)
  res <- bind_rows(res, temp.res)

}

return(res)
}


all.res <- lapply(c("40-50", ">50"), function(x){rates(x)}) # "all", "40-50", ">50"
all.res <- bind_rows(all.res) %>%
  mutate(n_screen = ifelse(n_screen<11 & !is.na(n_screen), NA, n_screen),
         n_total = ifelse(n_total<11 & !is.na(n_total), NA, n_total),
         n = ifelse(n<11, NA, n))

if (samp==1) {
  write_csv(all.res, "./results/crude_rates_samp.csv")

} else {
  write_csv(all.res, "./results/crude_rates.csv")

}
 
