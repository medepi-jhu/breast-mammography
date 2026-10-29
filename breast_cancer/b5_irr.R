
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged >=30

# Content: Estimate crude and weighted IRR (PWH / PWoH) and CIs (sandwich estimator)

# Date: 7/14/2025

# Last update: 7/15/2025

# Author: Yiyi Zhou

################################################################################

library(polars)
library(splines)
library(tictoc)
library(parallel)


# directories ------------------------------------------------------------------

data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/irr/"
wt_check_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/irr/combined_wt/"



# function to estimate weights -------------------------------------------------

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



# function to estiamte crude and weighted incidence rates by HIV status---------

calc_ir <- function(ir_dat) {
  
  ir_dat <- ir_dat$lazy()$
    with_columns(
      followup_years = pl$col("age_end") - pl$col("age_start")
    )$
    with_columns(
      followup_years_wt1 = pl$col("followup_years") * pl$col("wt1_final_trunc"),
      event_trad_wt1 = pl$col("event_trad") * pl$col("wt1_final_trunc"),
      followup_years_wt2 = pl$col("followup_years") * pl$col("wt2_final_trunc"),
      event_trad_wt2 = pl$col("event_trad") * pl$col("wt2_final_trunc")
    )$
    group_by("hiv")$
    agg(
      total_cases = pl$col("event_trad")$sum(),
      total_person_years = pl$col("followup_years")$sum(),
      total_cases_wt1 = pl$col("event_trad_wt1")$sum(),
      total_person_years_wt1 = pl$col("followup_years_wt1")$sum(),
      total_cases_wt2 = pl$col("event_trad_wt2")$sum(),
      total_person_years_wt2 = pl$col("followup_years_wt2")$sum()
    )$
    with_columns(
      ir_crude = pl$col("total_cases") / pl$col("total_person_years"),
      ir_wt1 = pl$col("total_cases_wt1") / pl$col("total_person_years_wt1"),
      ir_wt2 = pl$col("total_cases_wt2") / pl$col("total_person_years_wt2")
    )$
    collect()
  
  return(ir_dat)
  
}



# function to bootstrap --------------------------------------------------------

boot_fun <- function(dat, dat_name, r) {
  
  is_race_specific <- any(grepl("white|black|hisp", dat_name))
  
  tic(paste0("bootstrap ", r))
  
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
  new.boot <- est_wt(boot, dat_name, analysis_type, is_race_specific, r)
  cat(dat_name, "-", r, "- weights estimated\n")
  rm(boot)
  
  # estimate crude and weighted 
  res <- calc_ir(new.boot)
  
  # estimate crude and weighted incidence rate ratio
  value_cols <- setdiff(res$columns, "hiv")
  res <- res$
    with_columns(data_name = pl$lit(dat_name))$
    pivot(
      values = value_cols,
      index = "data_name",
      on = "hiv"
    )$
    with_columns(
      irr_crude = pl$col("ir_crude_1") / pl$col("ir_crude_0"),
      irr_wt1 = pl$col("ir_wt1_1") / pl$col("ir_wt1_0"),
      irr_wt2 = pl$col("ir_wt2_1") / pl$col("ir_wt2_0"),
      rep = pl$lit(r)$cast(pl$Int64)
    )$
    select(c(
      "data_name", "total_cases_1", "total_cases_0", "total_person_years_1", 
      "total_person_years_0", "total_cases_wt1_1", "total_cases_wt1_0", 
      "total_person_years_wt1_1", "total_person_years_wt1_0", "total_cases_wt2_1", 
      "total_cases_wt2_0", "total_person_years_wt2_1", "total_person_years_wt2_0", 
      "ir_crude_1", "ir_crude_0", "ir_wt1_1", "ir_wt1_0", "ir_wt2_1", 
      "ir_wt2_0", "irr_crude", "irr_wt1", "irr_wt2", "rep"
    ))
  
  toc()
  
  return(res$to_data_frame())
  
}



# function to aggregate bootstrapped results -----------------------------------

boot_out <- function(boot_est) {
  
  boot_est <- as_polars_df(boot_est)
  
  estimates_to_summarize <- c(
    "total_cases_1", "total_cases_0", "total_person_years_1", 
    "total_person_years_0", "total_cases_wt1_1", "total_cases_wt1_0", 
    "total_person_years_wt1_1", "total_person_years_wt1_0", "total_cases_wt2_1", 
    "total_cases_wt2_0", "total_person_years_wt2_1", "total_person_years_wt2_0", 
    "ir_crude_1", "ir_crude_0", "ir_wt1_1", "ir_wt1_0", "ir_wt2_1", 
    "ir_wt2_0", "irr_crude", "irr_wt1", "irr_wt2"
  )
  
  boot_ci_expr <- unlist(lapply(estimates_to_summarize, function(col_name) {
    list(
      pl$col(col_name)$quantile(0.025)$alias(paste0(col_name, "_ll")),
      pl$col(col_name)$quantile(0.975)$alias(paste0(col_name, "_ul"))
    )
  }))
  
  boot_summ <- boot_est$select(boot_ci_expr)
  
  point_est <- boot_est$filter(pl$col("rep") == 0)
  
  final_est <- point_est$
    join(boot_summ, how = "cross")$
    drop("rep")
  
  format_expr <- function(col_name) {
    pl$concat_str(
      pl$col(col_name), pl$lit(" ("), 
      pl$col(paste0(col_name, "_ll")), pl$lit(", "),
      pl$col(paste0(col_name, "_ul")), pl$lit(")"),
      separator = ""
    )$alias(col_name) # overwrite using the formatted columns
  }
  
  final_est_formatted <- final_est$
    with_columns(
      pl$col("^total_.*$")$round(2),
      (pl$col("^ir_.*$") * 10000)$round(2),
      pl$col("^irr_.*$")$round(2)
    )$
    with_columns(
      lapply(estimates_to_summarize, format_expr)
    )$
    select(c(
      "data_name", 
      estimates_to_summarize
    ))
  
  return(list(final_est = final_est, final_est_formatted = final_est_formatted))
  
}



# a main function to process each dataset --------------------------------------

process_each_dataset <- function(dat_name, r) {
  
  dat <- pl$scan_ipc(paste0(data_dir, dat_name, ".l.arrow"))$
    select(c(
      "bene_id","age_base", "year_base", "state", "race", "n_comorb_base_cat", 
      "n_comorb_cat", "n_comorb_cat_1", "n_comorb_cat_2", "hiv", "wt_s",
      "no_int", "event_int", "drop_compete", "drop_trad", "age_start", "age_end",
      "age_start_floor", "event_trad", "event_compete"
    ))$
    # rename(bene_id = "bid")$
    with_columns(
      pl$col("hiv")$cast(pl$Int8),
      pl$col("event_int")$cast(pl$Int64),
      pl$col("no_int")$cast(pl$Int64),
      age_end = pl$when(pl$col("age_end") == pl$col("age_start"))$
        then(pl$col("age_end") + 0.000001)$
        otherwise(pl$col("age_end"))
    )
  
  # bootstrap
  # all.boot <- boot_fun(dat, dat_name, 0) # point estimates
  # return(all.boot)
  
  all.boot <- mclapply(0:r, function(n) {
    boot_fun(dat, dat_name, n)
  }, mc.cores = 5, mc.set.seed = FALSE)
  
  rm(dat)
  
  all.boot.res <- do.call(rbind, all.boot)
  
  readr::write_csv(all.boot.res, paste0(table_dir, "all_rep_res/", dat_name, "_boot_irr_all_rep.csv"))
  
  all.boot.res.out <- boot_out(all.boot.res)
  
  all.boot.res.out$final_est$write_csv(paste0(table_dir, dat_name, "_boot_irr.csv"))
  all.boot.res.out$final_est_formatted$write_csv(paste0(table_dir, dat_name, "_boot_irr_formatted.csv"))
  
  return(invisible(NULL))
  
}



# # analysis (point estimates only) ----------------------------------------------
# 
# analysis_type <- "compete"
# # test <- process_each_dataset("samp.dat", 3)
# # print(test)
# 
# dat_names <- c(
#   "samp.dat",
#   "samp.dat.white",
#   "samp.dat.black",
#   "samp.dat.hisp", # only 13 cases among Hispanic PWH
#   
#   "samp.age.30.55",
#   "samp.age.30.55.white",
#   "samp.age.30.55.black",
#   # "samp.age.30.55.hisp",
#   
#   "samp.age.55.65",
#   "samp.age.55.65.white",
#   "samp.age.55.65.black",
#   # "samp.age.55.65.hisp",
#   
#   "samp.year.2001.2012",
#   "samp.year.2012.2020",
#   "samp.year.2001.2012.age.30.55",
#   "samp.year.2012.2020.age.30.55",
#   "samp.year.2001.2012.age.55.65",
#   "samp.year.2012.2020.age.55.65",
#   
#   "samp.year.2001.2012.white",
#   "samp.year.2001.2012.black",
#   "samp.year.2012.2020.white",
#   "samp.year.2012.2020.black"
# )
# 
# # point estimates
# out_df <- data.frame(matrix(ncol = 23, nrow = 1))
# colnames(out_df) <- c(
#   "data_name", "total_cases_1", "total_cases_0", "total_person_years_1", 
#   "total_person_years_0", "total_cases_wt1_1", "total_cases_wt1_0", 
#   "total_person_years_wt1_1", "total_person_years_wt1_0", "total_cases_wt2_1", 
#   "total_cases_wt2_0", "total_person_years_wt2_1", "total_person_years_wt2_0", 
#   "ir_crude_1", "ir_crude_0", "ir_wt1_1", "ir_wt1_0", "ir_wt2_1", 
#   "ir_wt2_0", "irr_crude", "irr_wt1", "irr_wt2", "rep"
# )
# out_df <- as_polars_df(out_df)
# 
# for (dat_name in dat_names) {
#   cat("Analyzing:", dat_name, "\n")
#   start_time <- Sys.time()
#   out_res <- process_each_dataset(dat_name, 0)
#   end_time <- Sys.time()
#   cat(dat_name, ":point estimates done, ", end_time - start_time, "\n")
#   out_df <- pl$concat(list(out_df, out_res), how = "vertical_relaxed")
# }
# 
# out_df <- out_df$filter(pl$col("data_name")$is_not_null())
# out_df$write_csv(paste0(table_dir, "irr_point_est.csv"))
# 
# crude_ir_df <- pl$read_csv("/users/57285/c-yzhou172-57285/breast_cancer/tables/crude_ir/breast_cx_crude_ir_formatted.csv")
# 
# cols_not_to_unpivot <- c("data_name", "irr_crude", "irr_wt1", "irr_wt2", "rep")
# 
# out_df_long <- out_df$
#   unpivot(
#     index = cols_not_to_unpivot,
#     on = setdiff(out_df$columns, cols_not_to_unpivot),
#     variable_name = "variable",
#     value_name = "value"
#   )$
#   with_columns(
#     hiv = pl$col("variable")$str$slice(-1)$cast(pl$Int8),
#     base_name = pl$col("variable")$str$split("_")$list$slice(0, pl$col("variable")$str$split("_")$list$len() - 1)
#   )$
#   with_columns(
#     base_name = pl$col("base_name")$list$join("_")
#   )$
#   drop("variable")$
#   pivot(
#     values = "value",
#     index = c("hiv", cols_not_to_unpivot),
#     on = "base_name"
#   )$
#   select(c(
#     "data_name", "hiv", "total_cases", "total_person_years", "total_cases_wt1",
#     "total_person_years_wt1", "total_cases_wt2", "total_person_years_wt2", 
#     "ir_crude", "ir_wt1", "ir_wt2", "irr_crude", "irr_wt1", "irr_wt2", "rep"
#   ))$
#   rename(
#     total_cases = "total_cases_crude",
#     total_person_years = "total_person_years_crude"
#   )
# 
# out_df_long$write_csv(paste0(table_dir, "irr_point_est_long.csv"))
# 
# combined_df <- crude_ir_df$
#   with_columns(
#     pl$col("hiv")$cast(pl$Int8),
#     data_name = pl$when(pl$col("data_name") == "samp.white")$
#       then(pl$lit("samp.dat.white"))$
#       when(pl$col("data_name") == "samp.black")$
#       then(pl$lit("samp.dat.black"))$
#       when(pl$col("data_name") == "samp.hisp")$
#       then(pl$lit("samp.dat.hisp"))$
#       otherwise(pl$col("data_name"))
#   )$
#   join(
#     out_df_long, 
#     on = c("data_name", "hiv"),
#     how = "left"
#   )
# 
# combined_df$
#   with_columns(
#     pl$col("total_cases")$cast(pl$Int64), 
#     pl$col("total_person_years")$cast(pl$Int64), 
#     pl$col("total_cases_crude")$cast(pl$Int64), 
#     pl$col("total_person_years_crude")$cast(pl$Int64)
#   )$
#   filter(
#     (pl$col("total_cases") != pl$col("total_cases_crude")) | 
#       ((pl$col("total_person_years") - pl$col("total_person_years_crude"))$abs() > 1)
#   )$select(c(
#     "data_name", "hiv", "total_cases", "total_person_years", 
#     "total_cases_crude", "total_person_years_crude"
#   ))$height # 0
# 
# combined_df <- combined_df$
#   with_columns(
#     ir_wt1 = (pl$col("ir_wt1") * 10000)$round(2),
#     ir_wt2 = (pl$col("ir_wt2") * 10000)$round(2),
#     irr_crude = pl$col("irr_crude")$round(2),
#     irr_wt1 = pl$col("irr_wt1")$round(2),
#     irr_wt2 = pl$col("irr_wt2")$round(2)
#   )$
#   with_columns(
#     irr_crude = pl$when(pl$col("hiv") == 0)$
#       then(pl$lit("ref"))$
#       otherwise(pl$col("irr_crude")),
#     irr_wt1 = pl$when(pl$col("hiv") == 0)$
#       then(pl$lit("ref"))$
#       otherwise(pl$col("irr_wt1")),
#     irr_wt2 = pl$when(pl$col("hiv") == 0)$
#       then(pl$lit("ref"))$
#       otherwise(pl$col("irr_wt2"))
#   )$
#   select(c(
#     crude_ir_df$columns, "ir_wt1", "ir_wt2", "irr_crude", "irr_wt1", "irr_wt2"
#   ))
# 
# combined_df$write_csv(paste0(table_dir, "ir_irr_combined.csv"))



# analysis (with bootstrap CIs)-------------------------------------------------

analysis_type <- "compete"
dat_name <- as.character(commandArgs(trailingOnly = T)[1])
process_each_dataset(dat_name, 1000)

# organize results
csv_files <- list.files(path = table_dir, pattern = "^samp.*formatted\\.csv$", full.names = TRUE)

file_info <- as_polars_df(
  tibble(full_path = csv_files, filename = basename(full_path))
)$with_columns(
  data_name = pl$col("filename")$str$replace(
    pattern = "_boot_irr_formatted.csv",
    value = ""
  ))

list_of_irr_dfs <- list()

for (i in 1:file_info$height) {
  params <- file_info[i, ]
  irr_df <-  pl$read_csv(params$get_column("full_path"))$select(c("data_name", "irr_crude", "irr_wt1", "irr_wt2"))
  list_of_irr_dfs[[i]] <- irr_df
}

all_irr_dfs <- pl$concat(list_of_irr_dfs, how = "vertical_relaxed")

crude_ir_df <- pl$read_csv(paste0(table_dir, "ir_irr_combined.csv"))$
  select(c("data_name", "hiv", "total_cases", "total_person_years", "crude_ir_ci"))

combined_df <- crude_ir_df$
  join(
    all_irr_dfs,
    on = "data_name",
    how = "left"
  )$
  with_columns(
    irr_crude = pl$when(pl$col("hiv") == 0)$
      then(pl$lit("Reference"))$
      otherwise(pl$col("irr_crude")),
    irr_wt1 = pl$when(pl$col("hiv") == 0)$
      then(pl$lit("Reference"))$
      otherwise(pl$col("irr_wt1")),
    irr_wt2 = pl$when(pl$col("hiv") == 0)$
      then(pl$lit("Reference"))$
      otherwise(pl$col("irr_wt2"))
  )

combined_df$write_csv(paste0(table_dir, "ir_irr_combined_boot_ci.csv"))
