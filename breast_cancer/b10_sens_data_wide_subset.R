
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged >=30

# Content: subset the wide data for sensitivity analysis (removing NY 2012-2016)

# Date: 10/20/2025

# Last update: 

# Author: Yiyi Zhou

################################################################################

library(polars)
# library(exactci)
# library(ggplot2)

data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/crude_ir/"

df <- pl$read_ipc(paste0(data_dir, "samp.dat.ex.ny.w.arrow"))

# df$schema



# Functions --------------------------------------------------------------------

# function to subset data by year
year_strat <- function(data, year_lower, year_upper) {
  
  lower_date <- as.Date(paste0(year_lower, "-01-01"))
  upper_date <- as.Date(paste0(year_upper, "-01-01"))
  
  year.dat <- data$filter(
    pl$col("baseline")$dt$year() < year_upper &
      pl$col("first_event_date")$dt$year() >= year_lower
  )$with_columns(
    
    baseline0 = pl$col("baseline"),
    first_event_date0 = pl$col("first_event_date")
  )$with_columns(
    
    baseline = pl$when(pl$col("baseline")$gt(lower_date))$
      then(pl$col("baseline"))$
      otherwise(lower_date),
    
    index = (pl$col("first_event_date")$dt$year() >= year_upper)$cast(pl$Int8)
  )$with_columns(
    first_event_date = pl$when(pl$col("index") == 1)$
      then(upper_date)$
      otherwise(pl$col("first_event_date")),
    
    event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      cast(pl$Int8),
    
    compete_event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("compete_event_flag"))$
      cast(pl$Int8)
  )$with_columns(
    event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
  )$with_columns(
    event_time = pl$when(pl$col("event_time") == 0)$
      then(pl$col("event_time") + 0.000001)$
      otherwise(pl$col("event_time"))
  )$with_columns(
    event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil(),
    year_cat = pl$lit(paste0("[", year_lower, ", ", year_upper, ")"))
  )$drop("index")
  
  if (year.dat$filter(pl$col("event_time") <= 0)$height > 0) {
    cat("samp.year.", year_lower, ".", year_upper, "has event_time <= 0!")
  }
  
  year.dat$write_ipc(paste0(
    data_dir, "samp.year.", year_lower, ".", year_upper, ".ex.ny.w",".arrow"
  ))
  # return(year.dat)
}

# function to further stratify year-specific data by time-updated age group
age_strat2 <- function(data, age_lower, age_upper, race = NULL) {
  
  age.dat <- data$filter(
    # remove those who entered study after age_upper
    pl$col("dob")$dt$offset_by(paste0(age_upper, "y")) > pl$col("baseline") &
      # remove those whose event occurred before age_lower
      pl$col("dob")$dt$offset_by(paste0(age_lower, "y")) <= pl$col("first_event_date")
  )$with_columns(
    
    # baseline0: original baseline, baseline1: recoded start date in year-stratified data
    baseline1 = pl$col("baseline"),
    
    # first_event_date0: original baseline, first_event_date1: recoded start date in year-stratified data
    first_event_date1 = pl$col("first_event_date")
  )$with_columns(
    
    # recode baseline to account for age_lower
    baseline = pl$when(pl$col("baseline") > pl$col("dob")$dt$offset_by(paste0(age_lower, "y")))$
      then(pl$col("baseline"))$
      otherwise(pl$col("dob")$dt$offset_by(paste0(age_lower, "y"))),
    
    # create index for events after age_upper
    index = (pl$col("first_event_date") >= pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      cast(pl$Int8)
  )$with_columns(
    
    # recode first_event_date
    first_event_date = pl$when(pl$col("index") == 1)$
      then(pl$col("dob")$dt$offset_by(paste0(age_upper, "y")))$
      otherwise(pl$col("first_event_date")),
    
    # update event_flag
    event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      cast(pl$Int8),
    
    # update compete_event_flag
    compete_event_flag = pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("compete_event_flag"))$
      cast(pl$Int8)
  )$with_columns(
    
    # update event_time in days
    event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
  )$with_columns(
    # deal with a special case where the recoded baseline and the recoded 
    # first_event_date are the same day (for first_event_date occurred on the date of age_lower)
    event_time = pl$when(pl$col("event_time") == 0)$
      then(pl$col("event_time") + 0.000001)$
      otherwise(pl$col("event_time"))
  )$with_columns(
    # re-calculate number of time intervals
    event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil(),
    
    # add age category label
    age_cat = pl$lit(paste0("[", age_lower, ", ", age_upper, ")"))
  )$filter(pl$col("event_time") > 0)$
    drop("index")
  
  if (age.dat$filter(pl$col("event_time") <= 0)$height > 0) {
    cat("samp.age.", "samp.age.", age_lower, ".", age_upper, 
        ".", race, "has event_time <= 0!", "\n")
  }
  
  if (is.null(race)) {
    age.dat$write_ipc(paste0(
      data_dir, year_dat_name, ".age.", age_lower, ".", age_upper, ".w", ".arrow"
    ))
  } else {
    age.dat$write_ipc(paste0(
      data_dir, year_dat_name, ".age.", age_lower, ".", age_upper, ".", race, ".w", ".arrow"
    ))
  }
  # return(age.dat)
}



# Subset wide data -------------------------------------------------------------


# by calendar period
year_strat(df, 2012, 2020)

# by calendar period and race
races <- c("white", "black")

for (race in races) {
  race_strat(
    pl$read_ipc(paste0(data_dir, "samp.year.2012.2020.ex.ny", ".w.arrow")),
    "samp.year.2012.2020.ex.ny",
    race
  )
}

# by age and calendar period

age_cats <- list(c(30, 55), c(55, 65))
year_dat_name <- "samp.year.2012.2020.ex.ny"

for (age_cat in age_cats) {
  age_lower <- age_cat[1]
  age_upper <- age_cat[2]
  age_strat2(
    pl$read_ipc(paste0(data_dir, year_dat_name,".w.arrow")), 
    age_lower, 
    age_upper, 
    race = NULL
  )
}
