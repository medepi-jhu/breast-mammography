
################################################################################
#
# Project: MedEpi Mammography
#
# Purpose: Set up data for weighted analysis
#
# Last Update: 29 Jul 2025
#
################################################################################

library("tidyverse")
library("duckplyr")
library("splines")
library("parallel")

boot_start <- 0
boot_end <- 0
cores <- 10

period <- "post2012" #all, pre2012, post2012
group <- "overall" #overall, black, white


# Read in data ------------------------------------------------------------

if (period=="all") {
  dat <- read_csv("/cms01/data/dua/57285/users/c-jrudolp9-57285/mamm_long.csv") %>%
    mutate(across(c(n_cond, n_cond_1, n_cond_2, race, int, yr_base), ~as.factor(.x)),
           py = (int_start %--% int_end)/years(1))

} else {
  dat <- read_csv(paste0("/cms01/data/dua/57285/users/c-jrudolp9-57285/mamm_long_", period, ".csv")) %>%
    mutate(across(c(n_cond, n_cond_1, n_cond_2, race, int, yr_base), ~as.factor(.x)),
           py = (int_start %--% int_end)/years(1))

}

if (group=="black") {
  dat <- filter(dat, race==2)

} else if (group=="white") {
  dat <- filter(dat, race==1)

}

base <- filter(dat, !duplicated(bene_id))


# Bootstrap resample ------------------------------------------------------

boot_rep <- function(r) {

  set.seed(r)

  samp <- table(base[sample(1:nrow(base), nrow(base), replace=T), (names(base)=="bene_id")])

  boot <- NULL
  if (r==0) {
    boot <- dat %>%
      rename(bid = bene_id)
  } else {
    for(zzz in 1:max(samp)){
      cc <- dat[dat$bene_id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid <- paste0(cc$bene_id, zzz)
      boot <- rbind(boot, cc)
    }
    boot <- select(boot, -bene_id)
  }

  boot.base <- filter(boot, !duplicated(bid))


# Estimate weights -------------------------------------------------------

  # Baseline odds weights
  if (group=="overall") {
  ps1 <- glm(hiv_base ~ bs(age_base, df=3) + race + state + yr_base,
                     family=binomial(link="logit"), weights=wt_s, data=boot.base)$fitted.values
  ps2 <- glm(hiv_base ~ bs(age_base, df=3) + race + state + yr_base + n_cond,
                     family=binomial(link="logit"), weights=wt_s, data=boot.base)$fitted.values
  } else {
  ps1 <- glm(hiv_base ~ bs(age_base, df=3) + state + yr_base,
                     family=binomial(link="logit"), weights=wt_s, data=boot.base)$fitted.values
  ps2 <- glm(hiv_base ~ bs(age_base, df=3) + state + yr_base + n_cond,
                     family=binomial(link="logit"), weights=wt_s, data=boot.base)$fitted.values

  }

  boot.base$wt_hiv1 <- ifelse(boot.base$hiv_base==1, 1, ps1/(1-ps1))
  boot.base$wt_hiv2 <- ifelse(boot.base$hiv_base==1, 1, ps2/(1-ps2))

  # Censoring weights in PWH
  if (group=="overall") {
  ps.c.den1 <- glm(I(drop==0) ~ bs(age_base, df=3) + race + state + yr_base + n_cond + n_cond_1 + n_cond_2 + int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==1))$fitted.values
  } else {
  ps.c.den1 <- glm(I(drop==0) ~ bs(age_base, df=3) + state + yr_base + n_cond + n_cond_1 + n_cond_2 + int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==1))$fitted.values
  }

  ps.c.num1 <- glm(I(drop==0) ~ int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==1))$fitted.values

  boot$wt_c[boot$hiv_base==1] <- ps.c.num1/ps.c.den1

  # Censoring weights in PWoH
  if (group=="overall") {
  ps.c.den0 <- glm(I(drop==0) ~ bs(age_base, df=3) + race + state + yr_base + n_cond + n_cond_1 + n_cond_2 + int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==0))$fitted.values
  } else {
  ps.c.den0 <- glm(I(drop==0) ~ bs(age_base, df=3) + state + yr_base + n_cond + n_cond_1 + n_cond_2 + int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==0))$fitted.values
  }

  ps.c.num0 <- glm(I(drop==0) ~ int,
                     family=binomial(link="logit"), weights=wt_s, data=filter(boot, hiv_base==0))$fitted.values

  boot$wt_c[boot$hiv_base==0] <- ps.c.num0/ps.c.den0

  if (r==1) {
  print(summary(boot.base$wt_hiv1))
  print(summary(boot.base$wt_hiv2))
  print(summary(boot$wt_c))
  }

  # Combine weights
  boot2 <- boot %>%
    left_join(select(boot.base, bid, wt_hiv1, wt_hiv2), by="bid") %>%
    group_by(bid) %>%
    mutate(cum_wt_c = cumprod(wt_c),
           wt1 = wt_s*wt_hiv1*cum_wt_c,
           wt2 = wt_s*wt_hiv2*cum_wt_c)

  # Truncate the weights (some censoring weights were extreme)
  wt1_p01 <- quantile(boot2$wt1, p=0.01)
  wt1_p99 <- quantile(boot2$wt1, p=0.99)
  boot2$wt1 <- ifelse(boot2$wt1<wt1_p01, wt1_p01, boot2$wt1)
  boot2$wt1 <- ifelse(boot2$wt1>wt1_p99, wt1_p99, boot2$wt1)

  wt2_p01 <- quantile(boot2$wt2, p=0.01)
  wt2_p99 <- quantile(boot2$wt2, p=0.99)
  boot2$wt2 <- ifelse(boot2$wt2<wt2_p01, wt2_p01, boot2$wt2)
  boot2$wt2 <- ifelse(boot2$wt2>wt2_p99, wt2_p99, boot2$wt2)


# Estimate rates ---------------------------------------------------------

  rates <- boot2 %>%
    group_by(hiv_base) %>%
    summarize(rate_crude = sum(had_mamm)/sum(py),
              rate_wt1 = sum(had_mamm*wt1)/sum(py*wt1),
              rate_wt2 = sum(had_mamm*wt2)/sum(py*wt2))

  first <- filter(boot2, !duplicated(bid)) %>%
    group_by(hiv_base) %>%
    summarize(n = n())

  est <- data.frame(rep = r,
                    n_hiv = first$n[first$hiv_base==1],
                    n_nohiv = first$n[first$hiv_base==0],
                                rate_crude_hiv = rates$rate_crude[rates$hiv_base==1],
                                rate_crude_nohiv = rates$rate_crude[rates$hiv_base==0],
                    logirr_crude = log(rates$rate_crude[rates$hiv_base==1]/rates$rate_crude[rates$hiv_base==0]),

                                rate_wt1_hiv = rates$rate_wt1[rates$hiv_base==1],
                    rate_wt1_nohiv = rates$rate_wt1[rates$hiv_base==0],
                                logirr_wt1 = log(rates$rate_wt1[rates$hiv_base==1]/rates$rate_wt1[rates$hiv_base==0]),

                                rate_wt2_hiv = rates$rate_wt2[rates$hiv_base==1],
                    rate_wt2_nohiv = rates$rate_wt2[rates$hiv_base==0],
                    logirr_wt2 = log(rates$rate_wt2[rates$hiv_base==1]/rates$rate_wt2[rates$hiv_base==0]))

  return(est)

}

#all.res <- lapply(0:nboot, function(x){boot_rep(x)})
all.res <- mclapply(boot_start:boot_end, function(x){boot_rep(x)}, mc.cores=cores, mc.set.seed=F)
all.res <- bind_rows(all.res)

# res <- filter(all.res, rep==0)
# res.summ <- all.res %>%
#   summarize(across(c(rate_crude_hiv, rate_crude_nohiv, logirr_crude, rate_wt1_hiv, rate_wt1_nohiv, logirr_wt1,
#                    rate_wt2_hiv, rate_wt2_nohiv, logirr_wt2), ~ quantile(.x, p=0.025), .names="{.col}_ll"),
#             across(c(rate_crude_hiv, rate_crude_nohiv, logirr_crude, rate_wt1_hiv, rate_wt1_nohiv, logirr_wt1,
#                      rate_wt2_hiv, rate_wt2_nohiv, logirr_wt2), ~ quantile(.x, p=0.975), .names="{.col}_ul"))
#
# res <- bind_cols(res, res.summ)

if (period=="all") {
  write_csv(all.res, paste0("./results/wt_rates_", group, "_", boot_end, ".csv"))

} else {
  write_csv(all.res, paste0("./results/wt_rates_", group, "_", period, "_", boot_end, ".csv"))

}
