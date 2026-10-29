################################################################################
# 
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged 30-64
#
# Content: Stratify long data by time-fixed variable, e.g., race, or excluding a specific race group
#
# Date: 7/11/2025
#
# Author: Yiyi Zhou
#
################################################################################

library(polars)

data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/long_data_check/"

# by time-updated age or period
subset_long_data_tv <- function(file_name) {
  
  # --- make data long---
  
  # pull in age- or period- stratified wide data
  dat <- pl$scan_ipc(paste0(data_dir, file_name, ".w.arrow"))
  
  # determine age_upper/year_upper
  if (grepl("year", file_name) & grepl("age", file_name)) {
    stratification <- "by_year_age"
    year_upper <- sub('.*year\\.\\d{4}\\.(\\d{4}).*', '\\1', file_name)
    age_upper <- sub('.*age\\.\\d+\\.(\\d+).*', '\\1', file_name)
    age_lower <- sub('.*age\\.([0-9]+)\\.[0-9]+.*', '\\1', file_name)
  } else if (grepl("age", file_name)) {
    stratification <- "by_age"
    age_upper <- sub('.*age\\.\\d+\\.(\\d+).*', '\\1', file_name)
    age_lower <- sub('.*age\\.([0-9]+)\\.[0-9]+.*', '\\1', file_name)
  } else if (grepl("year", file_name)) {
    stratification <- "by_year"
    year_upper <- sub('.*year\\.\\d{4}\\.(\\d{4}).*', '\\1', file_name)
  } else {
    print("check the file name!")
  }
  
  dat1 <- dat$
    # repeat rows based on number of time intervals from baseline to end of follow-up
    with_columns(
      dummy = pl$lit(1)$repeat_by("event_int") # event_int: number of 90-day intervals for each id
    )$
    explode("dummy")$
    drop("dummy")$
    # group by 'bene_id' and create new variables
    with_columns(
      no_int = pl$col("bene_id")$cum_count()$over("bene_id"),
      drop_trad = pl$lit(0L),
      drop_compete = pl$lit(0L)
    )$with_columns(
      start_days_expr = ((pl$col("no_int") - 1) * 90)$cast(pl$Int64),
      end_days_expr = (pl$col("no_int") * 90)$cast(pl$Int64)
    )$with_columns(
      # start date of the record-specific time interval
      int_start = pl$col("baseline")$dt$offset_by(
        pl$concat_str(pl$col("start_days_expr"), pl$lit("d"), separator = "")
      ),
      # end date of the record_specific time interval
      int_end = pl$min_horizontal(
        pl$col("baseline")$dt$offset_by(
          pl$concat_str(pl$col("end_days_expr"), pl$lit("d"), separator = "")
        ),
        pl$col("first_event_date")
      )
    )$
    drop(c("start_days_expr", "end_days_expr"))$
    with_columns(
      event_trad = pl$when(pl$col("no_int") != pl$col("event_int"))$
        then(pl$lit(0))$
        otherwise(pl$col("event_flag"))$
        cast(pl$Int8),
      event_compete = pl$when(pl$col("no_int") != pl$col("event_int"))$
        then(pl$lit(0))$
        otherwise(pl$col("compete_event_flag"))$
        cast(pl$Int8)
    )$
    with_columns(
      # age on the start date of the record-specific interval
      age_start = (pl$col("int_start") - pl$col("dob"))$dt$total_days()$cast(pl$Float64) / 365.25,
      # age on the end date of the record-specific interval
      age_end = (pl$col("int_end") - pl$col("dob"))$dt$total_days()$cast(pl$Float64) / 365.25
    )$
    with_columns(
      # deal with rare cases where age_start==age_end in year- or age-stratified data
      # (first_event_date happening on the first day of the year period/age range)
      age_end = pl$when(pl$col("age_end") == pl$col("age_start"))$
        then(pl$col("age_end") + 0.000001)$
        otherwise(pl$col("age_end")),
      
      age_start_floor = pl$col("age_start")$floor()
    )
  
  # flag and count comorbidities (current and lagged)
  comorb_cols <- c(
    "myocard_infarct_typ1_dt", "cong_heart_fail_typ1_dt", "periph_vas_dis_typ1_dt", 
    "cereb_vas_dis_typ1_dt", "dementia_typ1_dt", "chron_pulm_dis_typ1_dt", 
    "rheumatic_dis_typ1_dt", "peptic_ulcer_typ1_dt", "mild_liv_dis_typ1_dt", 
    "diabt_wo_chron_complict_typ1_dt", "diabt_chron_complict_da_typ1_dt", 
    "hemiplegia_paraplegia_d_typ1_dt", "renal_dis_typ1_dt", "mod_sevr_liver_dis_typ1_dt"
  )
  
  tv_comorb_cols <- paste0(comorb_cols, "_tv")
  
  flag_comorb_logic <- function(comorb_col_name, starting_date) {
    pl$when(pl$col(comorb_col_name)$is_null())$
      then(pl$lit(0L))$
      otherwise(pl$col(comorb_col_name) <= pl$col(starting_date))$
      cast(pl$Int8)
  }
  
  dat2 <- dat1$
    with_columns(
      lapply(comorb_cols, function(x) {
        flag_comorb_logic(x, "int_start")$alias(paste0(x, "_tv"))
      })
    )$
    with_columns(
      n_comorb = pl$sum_horizontal(tv_comorb_cols)
    )$
    with_columns(
      n_comorb_cat = pl$when(pl$col("n_comorb") >= 2)$
        then(pl$lit("2+"))$
        otherwise(pl$col("n_comorb")$cast(pl$String))
    )$with_columns(
      # within each id's records, for the first row of the lagged variables, fill null with the original comorb cat
      n_comorb_cat_1 = pl$col("n_comorb_cat")$shift(1, fill_value = pl$col("n_comorb_cat"))$over("bene_id")
    )$with_columns(
      n_comorb_cat_2 = pl$col("n_comorb_cat")$shift(2, fill_value = pl$col("n_comorb_cat_1"))$over("bene_id")
    )$
    drop(tv_comorb_cols)$
    collect()
  
  # build drop-out indicators for those who were censored due to non-administrative censoring
  # e.g., disenrollment, death, other cancers
  censor_compete <- dat2$
    filter(
      (pl$col("int_end") == pl$col("first_event_date0")) & 
        (pl$col("compete_event_flag") == 0) &
        (pl$col("int_end") < pl$col("age_65_date")) & 
        (pl$col("int_end") < pl$col("admin_censor_date"))
    )
  
  censor_trad <- dat2$
    filter(
      (pl$col("int_end") == pl$col("first_event_date0")) & 
        (pl$col("event_flag") == 0) &
        (pl$col("int_end") < pl$col("age_65_date")) & 
        (pl$col("int_end") < pl$col("admin_censor_date"))
    )
  
  if (stratification == "by_age") {
    censor_compete <- censor_compete$
      filter(
        pl$col("int_end") < pl$col("dob")$dt$offset_by(
          pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
        )
      )$
      select("bene_id")
    
    censor_trad <- censor_trad$
      filter(
        pl$col("int_end") < pl$col("dob")$dt$offset_by(
          pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
        )
      )$
      select("bene_id")
    
  } else if (stratification == "by_year") {
    censor_compete <- censor_compete$
      filter(pl$col("int_end") < pl$lit(as.Date(paste0(year_upper, "-01-01"))))$
      select("bene_id")
    
    censor_trad <- censor_trad$
      filter(pl$col("int_end") < pl$lit(as.Date(paste0(year_upper, "-01-01"))))$
      select("bene_id")
    
  } else if (stratification == "by_year_age") {
    censor_compete <- censor_compete$
      filter(
        (pl$col("int_end") < pl$col("dob")$dt$offset_by(
          pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
        )) & 
          (pl$col("int_end") < pl$lit(as.Date(paste0(year_upper, "-01-01"))))
      )$
      select("bene_id")
    
    censor_trad <- censor_trad$
      filter(
        (pl$col("int_end") < pl$col("dob")$dt$offset_by(
          pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
        )) & 
          (pl$col("int_end") < pl$lit(as.Date(paste0(year_upper, "-01-01"))))
      )$
      select("bene_id")
  }
  
  dat3 <- dat2$
    with_columns(
      drop_compete = pl$when(
        (pl$col("bene_id")$is_in(censor_compete$get_column("bene_id"))) &
          (pl$col("int_end") == pl$col("first_event_date0"))
      )$then(pl$lit(1L))$
        otherwise(pl$col("drop_compete"))$
        cast(pl$Int8),
      
      drop_trad = pl$when(
        (pl$col("bene_id")$is_in(censor_trad$get_column("bene_id"))) &
          (pl$col("int_end") == pl$col("first_event_date0"))
      )$then(pl$lit(1L))$
        otherwise(pl$col("drop_trad"))$
        cast(pl$Int8)
    )$sort(c("bene_id", "no_int"))
  
  # --- export ---
  est_size <- round((dat3$estimated_size() / (1024^3)), 2)
  dat3$write_ipc(paste0(data_dir, file_name, ".l.arrow"))
  
  
  # --- data check ---
  sink(paste0(table_dir, file_name, ".txt"))
  cat("Checking wide data: \n")
  
  dat <- dat$collect()
  print(dat$describe()$to_data_frame())
  print("distribution of event_time and event_int:")
  print(dat$select(c("event_time", "event_int"))$describe())
  print("total records:")
  print(dat$height)
  print("any missingness for important variables:")
  cols_to_check <- c(
    "race", "dob", "start_date", "end_date", "baseline", "first_event_date", 
    "first_event", "event_flag", "compete_event_flag", "event_time", "event_int",
    "hiv", "age_base", "year_base", "state", "n_comorb_base_cat"
  )
  print(dat$select(pl$col(cols_to_check)$is_null()$sum())$to_data_frame())
  print("first_event:")
  print(dat$group_by("first_event")$agg(pl$col("bene_id")$count()))
  print("event_flag:")
  print(dat$group_by("event_flag")$agg(pl$col("bene_id")$count()))
  print("compete_event_flag:")
  print(dat$group_by("compete_event_flag")$agg(pl$col("bene_id")$count()))
  
  print("censored but not administratively censored (trad):")
  if (stratification == "by_age") {
    print(
      dat$filter(
        (pl$col("event_flag") == 0) &
          (pl$col("first_event_date") < pl$col("dob")$dt$offset_by(
            pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
          )) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height
    )
  } else if (stratification == "by_year") {
    print(
      dat$filter(
        (pl$col("event_flag") == 0) &
          (pl$col("first_event_date") < pl$lit(as.Date(paste0(year_upper, "-01-01")))) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height 
    )
  } else if (stratification == "by_year_age") {
    print(
      dat$filter(
        (pl$col("event_flag") == 0) &
          (pl$col("first_event_date") < pl$col("dob")$dt$offset_by(
            pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
          )) &
          (pl$col("first_event_date") < pl$lit(as.Date(paste0(year_upper, "-01-01")))) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height
    )
  }
  
  print("censored but not administratively censored (compete):")
  if (stratification == "by_age") {
    print(
      dat$filter(
        (pl$col("compete_event_flag") == 0) &
          (pl$col("first_event_date") < pl$col("dob")$dt$offset_by(
            pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
          )) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height 
    )
  } else if (stratification == "by_year") {
    print(
      dat$filter(
        (pl$col("compete_event_flag") == 0) &
          (pl$col("first_event_date") < pl$lit(as.Date(paste0(year_upper, "-01-01")))) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height
    ) 
  } else if (stratification == "by_year_age") {
    print(
      dat$filter(
        (pl$col("compete_event_flag") == 0) &
          (pl$col("first_event_date") < pl$col("dob")$dt$offset_by(
            pl$concat_str(pl$lit(age_upper), pl$lit("y"), separator = "")
          )) &
          (pl$col("first_event_date") < pl$lit(as.Date(paste0(year_upper, "-01-01")))) &
          (pl$col("first_event_date") < pl$col("age_65_date")) & 
          (pl$col("first_event_date") < pl$col("admin_censor_date"))
      )$height
    )
  }
  rm(dat)
  
  cat("\n\n\n")
  cat("Checking long data: \n")
  
  print("total records:")
  print(dat3$height)
  print("unique ids:")
  print(dat3$select("bene_id")$unique()$height)
  print("any missingness for important variables:")
  cols_to_check <- c(
    "race", "dob", "start_date", "end_date", "baseline", "first_event_date", 
    "first_event", "event_flag", "compete_event_flag", "event_time", "event_int",
    "hiv", "age_base", "year_base", "state", "n_comorb_base_cat", "no_int", 
    "int_start", "int_end", "age_start", "age_end", "n_comorb_cat", "n_comorb_cat_1",
    "n_comorb_cat_2", "event_trad", "event_compete", "drop_trad", "drop_compete"
  )
  print(dat3$select(pl$col(cols_to_check)$is_null()$sum())$to_data_frame())
  print("any no_int > event_int:")
  print(dat3$filter(pl$col("no_int") > pl$col("event_int"))$height)
  print("any int_start > int_end:")
  print(dat3$filter(pl$col("int_start") > pl$col("int_end"))$height)
  print("any int_start == int_end:")
  print(dat3$filter(pl$col("int_start") == pl$col("int_end"))$height)
  print("any age_start >= age_end:")
  print(dat3$filter(pl$col("age_start") >= pl$col("age_end"))$height)
  print("first_event, baseline records:")
  print(
    dat3$filter(pl$col("no_int") == 1)$
      group_by("first_event")$agg(pl$col("bene_id")$count())
  )
  print("event_trad, all records:")
  print(dat3$group_by("event_trad")$agg(pl$col("bene_id")$count()))
  print("event_trad, last records:")
  print(
    dat3$filter(pl$col("no_int") == pl$col("event_int"))$
      group_by("event_trad")$agg(pl$col("bene_id")$count())
  )
  print("event_compete, all records:")
  print(dat3$group_by("event_compete")$agg(pl$col("bene_id")$count()))
  print("event_compete, last records:")
  print(
    dat3$filter(pl$col("no_int") == pl$col("event_int"))$
      group_by("event_compete")$agg(pl$col("bene_id")$count())
  )
  print("drop_trad, all records:")
  print(dat3$group_by("drop_trad")$agg(pl$col("bene_id")$count()))
  print("drop_trad, last records:")
  print(
    dat3$filter(pl$col("no_int") == pl$col("event_int"))$
      group_by("drop_trad")$agg(pl$col("bene_id")$count())
  )
  print("drop_compete, all records:")
  print(dat3$group_by("drop_compete")$agg(pl$col("bene_id")$count()))
  print("drop_compete, last records:")
  print(
    dat3$filter(pl$col("no_int") == pl$col("event_int"))$
      group_by("drop_compete")$agg(pl$col("bene_id")$count())
  )
  print(dat3$select(cols_to_check)$describe()$to_data_frame())
  
  sink()
  return(cat(file_name, ": done, estimated size:", est_size, "GB\n"))
}

file_names <- c(
  "samp.age.30.55",
  "samp.age.55.65",
  "samp.year.2001.2012",
  "samp.year.2012.2020",
  "samp.year.2001.2012.age.30.55",
  "samp.year.2012.2020.age.30.55",
  "samp.year.2001.2012.age.55.65",
  "samp.year.2012.2020.age.55.65"
)

for (fn in file_names){
  subset_long_data_tv(fn)
}


# by race only, by race and age together, by race and calendar period together
race_strat_long <- function(data, data_name, race) {
  race_mapping <- c("white" = "1", "black" = "2", "hisp" = "3")
  race_label <- race_mapping[[race]]
  race.dat <- data$filter(pl$col("race") == race_label)
  race.dat$write_ipc(paste0(data_dir, data_name, ".", race, ".l.arrow"))
}

file_names <- c(
  "samp.dat",
  "samp.age.30.55",
  "samp.age.55.65",
  "samp.year.2001.2012",
  "samp.year.2012.2020"
)

for (fn in file_names) {
  fn.dat <- pl$read_ipc(paste0(data_dir, fn, ".l.arrow"))
  if (fn == "samp.dat") {
    races <- c("white", "black", "hisp")
  } else {
    races <- c("white", "black")
  }
  for (race in races) {
    race_strat_long(fn.dat, fn, race)
  }
}

