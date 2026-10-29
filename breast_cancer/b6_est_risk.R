
# Project: Incidence of breast cancer by HIV status among Medicaid beneficiaries 

# Content: Estimate cumulative incidence with bootstrapped CIs for PWH and PWoH

# Date: 7/17/2025

# Last update: 7/21/2025

# Author: Yiyi Zhou

################################################################################

packages <- c(
  "polars", "tidyverse", "splines", "survival", "broom", "zoo", "tictoc", "parallel"
)
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
rm(packages)

file_name <- as.character(commandArgs(trailingOnly = T)[1])
analysis_type <- "compete"
data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/compete_risk/"



# read in data -----------------------------------------------------------------

dat <- pl$scan_ipc(paste0(data_dir, file_name, ".l.arrow"))$
  select(c(
    "bene_id","age_base", "year_base", "state", "race", "n_comorb_base_cat", 
    "n_comorb_cat", "n_comorb_cat_1", "n_comorb_cat_2", "hiv", "wt_s",
    "no_int", "event_int", "drop_compete", "drop_trad", "age_start", "age_end",
    "age_start_floor", "event_trad", "event_compete"
  ))$
  with_columns(
    pl$col("hiv")$cast(pl$Int8),
    pl$col("event_int")$cast(pl$Int64),
    pl$col("no_int")$cast(pl$Int64),
    age_end = pl$when(pl$col("age_end") == pl$col("age_start"))$
      then(pl$col("age_end") + 0.000001)$
      otherwise(pl$col("age_end"))
  )

# determine if the dataset is race-specific based on file name
is_race_specific <- any(grepl("white|black|hisp", file_name))


# a helper function to estimate weights

est_wt <- function(data, dat_name, type, is_race_specific, r) {
  
  data.base.r <- data$filter(pl$col("no_int") == 1)$collect()$to_data_frame()
  
  # --- build HIV weights ---
  formula_hiv <- "hiv ~ bs(age_base, df = 3) + as.factor(year_base) + state"
  
  if (!is_race_specific) {
    formula_hiv <- paste0(formula_hiv, " + race")
  }
  
  ps_hiv1 <- glm(formula_hiv, family = binomial(link = "logit"), data = data.base.r, weights = wt_s)$fitted.values
  
  formula_hiv <- paste0(formula_hiv, " + n_comorb_base_cat")
  
  ps_hiv2 <- glm(formula_hiv, family = binomial(link = "logit"), data = data.base.r, weights = wt_s)$fitted.values
  
  data.base.r$ps_hiv1 <- ps_hiv1
  data.base.r$ps_hiv2 <- ps_hiv2
  data.base.r$wt_hiv1 <- ifelse(data.base.r$hiv == 1, 1, ps_hiv1 / (1 - ps_hiv1))
  data.base.r$wt_hiv2 <- ifelse(data.base.r$hiv == 1, 1, ps_hiv2 / (1 - ps_hiv2))
  data.base <- as_polars_df(data.base.r)$lazy()$
    select(c("bid", "ps_hiv1", "ps_hiv2", "wt_hiv1", "wt_hiv2"))
  
  rm(data.base.r, formula_hiv, ps_hiv1, ps_hiv2)
  
  # --- build drop-out weights ---
  drop_ind <- ifelse(type == "compete", "I(drop_compete == 0)", "I(drop_trad == 0)")
  
  formula_drop <- paste0(drop_ind, " ~ bs(age_end, df = 3) + as.factor(year_base) + n_comorb_cat + n_comorb_cat_1 + n_comorb_cat_2 + state")
  
  if (!is_race_specific) {
    formula_drop <- paste0(formula_drop, " + race")
  }
  
  formula_drop_num <- paste0(drop_ind, " ~ bs(age_end, df = 3)")
  
  data <- data$with_columns(wt_drop = pl$lit(NA_real_))
  
  ## among those with hiv
  data.hiv.r <- data$filter(pl$col("hiv") == 1)$collect()$to_data_frame()
  ps_drop_den <- glm(formula_drop, family = binomial(link = "logit"), data = data.hiv.r, weights = wt_s)$fitted.values
  ps_drop_num <- glm(formula_drop_num, family = binomial(link = "logit"), data = data.hiv.r, weights = wt_s)$fitted.values
  data.hiv.r$ps_drop_den <- ps_drop_den
  data.hiv.r$ps_drop_num <- ps_drop_num
  data.hiv.r$wt_drop <- ps_drop_num / ps_drop_den
  data.hiv <- as_polars_df(data.hiv.r)$lazy()
  rm(data.hiv.r)
  
  ## among those without hiv
  data.nohiv.r <- data$filter(pl$col("hiv") == 0)$collect()$to_data_frame()
  ps_drop_den <- glm(formula_drop, family = binomial(link = "logit"), data = data.nohiv.r, weights = wt_s)$fitted.values
  ps_drop_num <- glm(formula_drop_num, family = binomial(link = "logit"), data = data.nohiv.r, weights = wt_s)$fitted.values
  data.nohiv.r$ps_drop_den <- ps_drop_den
  data.nohiv.r$ps_drop_num <- ps_drop_num
  data.nohiv.r$wt_drop <- ps_drop_num / ps_drop_den
  data.nohiv <- as_polars_df(data.nohiv.r)$lazy()
  
  rm(data.nohiv.r, ps_drop_den, ps_drop_num)
  
  # --- combine ---
  data <- pl$concat(list(data.hiv, data.nohiv), how = "vertical_relaxed")$
    join(data.base, on = "bid", how = "left")$
    sort(c("bid", "no_int"))$
    with_columns(
      cum_wt_drop = pl$col("wt_drop")$cum_prod()$over("bid")
    )$
    with_columns(
      wt1_final = pl$col("wt_hiv1") * pl$col("cum_wt_drop") * pl$col("wt_s"),
      wt2_final = pl$col("wt_hiv2") * pl$col("cum_wt_drop") * pl$col("wt_s")
    )$
    collect()
  
  # --- check distributions of PS, weights and truncate weights if needed ---
  
  get_weight_summary <- function(data_pl, wt_var, group_vars = NULL) {
    aggs <- list(
      n = pl$col(wt_var)$count(),
      min_wt = pl$col(wt_var)$min(),
      p1_wt = pl$col(wt_var)$quantile(0.01, "nearest"),
      median_wt = pl$col(wt_var)$median(),
      p99_wt = pl$col(wt_var)$quantile(0.99, "nearest"),
      max_wt = pl$col(wt_var)$max(),
      mean_wt = pl$col(wt_var)$mean(),
      sd_wt = pl$col(wt_var)$std()
    )
    if (is.null(group_vars)) {
      summary_df <- data_pl$select(aggs)
    } else {
      summary_df <- data_pl$group_by(group_vars)$agg(aggs)$sort(group_vars)
    }
    return(summary_df)
  }
  
  check_joint_weight_distribution <- function(data_pl, wt_final_var) {
    # joint weights
    wt_final_var_trunc <- paste0(wt_final_var, "_trunc")
    
    grouping_cols <- c("hiv", "age_start_floor")
    
    cat("--- Original combined weights:", wt_final_var, "' ---\n", sep = "")
    cat("\nOverall:\n")
    get_weight_summary(data_pl, wt_final_var)$print()
    
    # get_weight_summary(data_pl, "wt1_final_trunc_drop")$print()
    
    cat("\nBy HIV status and age:\n")
    get_weight_summary(data_pl, wt_final_var, grouping_cols)$print()
    
    cat("\n--- After truncation:'", wt_final_var_trunc, "' ---\n\n\n", sep = "")
    cat("\nOverall:\n")
    get_weight_summary(data_pl, wt_final_var_trunc)$print()
    cat("\nBy HIV status and age:\n")
    get_weight_summary(data_pl, wt_final_var_trunc, grouping_cols)$print()
    
    return(invisible(NULL))
  }
  
  # clip combined weights between their 1st and 99th percentiles
  truncate_weights <- function(data_pl) {
    data_pl$with_columns(
      wt1_final_trunc = pl$col("wt1_final")$clip(
        lower_bound = pl$col("wt1_final")$quantile(0.01, "nearest"),
        upper_bound = pl$col("wt1_final")$quantile(0.99, "nearest")
      ),
      wt2_final_trunc = pl$col("wt2_final")$clip(
        lower_bound = pl$col("wt2_final")$quantile(0.01, "nearest"),
        upper_bound = pl$col("wt2_final")$quantile(0.99, "nearest")
      )
    )
  }
  
  data <- truncate_weights(data)
  return(data)
}


# a helper function to estimate cumulative incidence (risk)

est_risk <- function(data, data_name, wt = NULL, type, r) {
  
  event_ind <- if (type == "compete") {
    data$event_compete
  } else if (type == "trad") {
    data$event_trad
  } else {
    print("check analysis type!")
  }
  
  wt_arg <- if (!is.null(wt)) {
    data[[wt]]
  } else {
    NULL
  }
  
  risk_fit <- survfit(
    Surv(age_start, age_end, event_ind) ~ hiv, 
    data = data, 
    weights = wt_arg, 
    id = bid, 
    se.fit = FALSE
  )
  
  
  res_long <- tibble(
    time = rep(risk_fit$time, times = length(risk_fit$states)),
    estimate = c(risk_fit$pstate),
    state = rep(risk_fit$states, each = nrow(risk_fit$pstate)),
    strata = rep(
      rep(names(risk_fit$strata), times = risk_fit$strata),
      times = length(risk_fit$states)
    )
  ) %>%
    filter(state != "" & state != "(s0)") %>%
    mutate(state = ifelse(state == 1, "breast", "death"))
  
  if (type == "trad") {
    return(res_long)
    
  } else if (type == "compete") {
    
    res_wide <- res_long %>%
      pivot_wider(
        names_from = c(state, strata),
        values_from = estimate,
        names_sep = "_"
      ) %>%
      bind_rows(tibble(time = min(data$age_start))) %>%
      arrange(time) %>%
      fill(everything(), .direction = "down") %>%
      mutate(across(everything(), ~replace_na(., 0))) %>%
      rename(
        risk_breast0 = `breast_hiv=0`,
        risk_death0 = `death_hiv=0`,
        risk_breast1 = `breast_hiv=1`,
        risk_death1 = `death_hiv=1`
      ) %>% 
      mutate(
        rd_breast = risk_breast1 - risk_breast0,
        rr_breast = risk_breast1 / risk_breast0,
        rd_death = risk_death1 - risk_death0,
        rr_death = risk_death1 / risk_death0,
        rep = r  
      ) %>%
      mutate(
        rr_breast = if_else(is.infinite(rr_breast), NA_real_, rr_breast),
        rr_death = if_else(is.infinite(rr_death), NA_real_, rr_death)
      )
    
    return(res_wide)
  }
}


# bootstrap when accounting for competing risk ---------------------------------

if (analysis_type == "compete") {
  
  boot_fun <- function(r) {
    
    # create bootstrap resampled dataset
    boot <- NULL
    if (r == 0) {
      boot <- dat$rename(bene_id = "bid") # assign bootstrap ID
    } else {
      # create a bootstrap sample of bene_id from the base data and count how many times
      # each bene_id was selected
      base <- dat$filter(pl$col("no_int") == 1)$collect()
      boot.samp <- base$
        sample(n = base$height, with_replacement = TRUE, seed = (r + 1))$
        select("bene_id")$
        with_columns(
          rep_n = pl$col("bene_id")$rank(method = "ordinal")$over("bene_id")
        )
      
      # if a bene_id was sampled k times, create k copies of that person's data from
      # the dat, giving each copy a unique bid
      boot <- dat$
        join(boot.samp$lazy(), on = "bene_id", how = "inner")$
        with_columns(
          bid = pl$concat_str(pl$col("bene_id"), pl$col("rep_n"))
        )$
        sort(c("bid", "no_int"))$
        drop(c("bene_id", "rep_n"))
    }
    
    # estimate weights
    new.boot <- est_wt(boot, file_name, analysis_type, is_race_specific, r)
    cat(file_name, "-", r, "- weights estimated\n")
    rm(boot)
    
    new.boot <- new.boot$to_data_frame() %>%  
      mutate(
        event_compete = factor(event_compete, levels = c(0, 1, 2)),
        event_trad = factor(event_trad, levels = c(0, 1)),
        hiv = factor(hiv, levels = c(0, 1))
      )
    
    # estimate crude risk
    res_crude <- est_risk(new.boot, file_name, wt = NULL, analysis_type, r)
    cat(file_name, "-", r, "- crude cif estimated\n")
    
    # estimate weighted risk
    res_wt1_trunc <- est_risk(new.boot, file_name, "wt1_final_trunc", analysis_type, r)
    cat(file_name, "-", r, "- wt1 cif estimated\n")
    
    res_wt2_trunc <- est_risk(new.boot, file_name, "wt2_final_trunc", analysis_type, r)
    cat(file_name, "-", r, "- wt2 cif estimated\n\n\n")
    
    rm(new.boot)
    
    toc()
    
    return(list(
      crude = res_crude, 
      wt1_trunc = res_wt1_trunc, 
      wt2_trunc = res_wt2_trunc
    ))
  }
  
  start_time <- Sys.time()
  r <- 1000
  # all.boot <- boot_fun(r) # r = 0
  all.boot <- mclapply(
    0:r, function(n) {boot_fun(n)},
    mc.cores = 5, mc.set.seed = FALSE
  )
  end_time <- Sys.time()
  print(paste0("Bootstrap success:", end_time - start_time))
  
  if (r == 0) {
    res_names <- c("crude", "wt1_trunc", "wt2_trunc")
    for (rn in res_names) {
      write_csv(all.boot[[rn]], paste0(table_dir, file_name, "_point_est_", rn, ".csv"))
    }
  } else {
    # summarize by time point
    boot_out <- function(boot_est) {
      boot_summ <- boot_est %>%
        group_by(time) %>%
        summarize(
          across(
            .cols = c(risk_breast0, risk_breast1, rd_breast, rr_breast,
                      risk_death0, risk_death1, rd_death, rr_death),
            .fns = list(
              ll = ~quantile(.x, p = 0.025, na.rm = TRUE),
              ul = ~quantile(.x, p = 0.975, na.rm = TRUE)
            ),
            .names = "{.col}_{.fn}"
          )
        )
      
      final_est <- boot_est %>% filter(rep == 0) %>% # get point estimates
        left_join(boot_summ, by = "time") %>%
        select(-rep)
      
      boot_est_name <- deparse(substitute(boot_est)) # extract boot summary data name as a string
      suffix <- sub(".*(\\.wt[1-2]\\.trunc|\\.crude)", "\\1", boot_est_name)
      rm(boot_est)
      
      write_csv(final_est, paste0(table_dir, "boot.", file_name, suffix, ".csv"))
    }
    
    all.boot.crude <- do.call(rbind, lapply(all.boot, `[[`, "crude"))
    write_csv(all.boot.crude, paste0(table_dir, "all_rep_res/", file_name, "_boot_cif_all_rep_crude.csv"))
    boot_out(all.boot.crude)
    rm(all.boot.crude)
    
    all.boot.wt1.trunc <- do.call(rbind, lapply(all.boot, `[[`, "wt1_trunc"))
    write_csv(all.boot.wt1.trunc, paste0(table_dir, "all_rep_res/", file_name, "_boot_cif_all_rep_wt1_trunc.csv"))
    boot_out(all.boot.wt1.trunc)
    rm(all.boot.wt1.trunc)
    
    all.boot.wt2.trunc <- do.call(rbind, lapply(all.boot, `[[`, "wt2_trunc"))
    write_csv(all.boot.wt2.trunc, paste0(table_dir, "all_rep_res/", file_name, "_boot_cif_all_rep_wt2_trunc.csv"))
    boot_out(all.boot.wt2.trunc) 
  }
}


