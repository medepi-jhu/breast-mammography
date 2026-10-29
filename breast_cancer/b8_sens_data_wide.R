
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged >=30

# Content: Manage data for sensitivity analysis removing NY from 2012-2016

# Date: 10/20/2025

# Last update: 

# Author: Yiyi Zhou

################################################################################

library(polars)
library(tableone)

raw_data_dir <- "/cms01/data/dua/57285/medican/data_requests/karine/20250422_cancer_incidence_v11/"
data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
table_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/tables/"



# read in new data ------------------------------------------------------------

# fetch only first row and get column names
df_temp <- pl$scan_csv(paste0(raw_data_dir, "incid.csv"))$limit(1)$collect()
col_names <- df_temp$columns

# define column types
dtypes_list <- lapply(col_names, function(col) {
  col_lower <- tolower(col)  
  if (grepl("id$", col_lower)) {
    pl$String
  } else if (grepl("flg|flag", col_lower)) {
    pl$Int8
  } else if (grepl("date|dt|typ1|typ2|dod|dob", col_lower)) {
    pl$Date
  } else if (grepl("sex|race|state|zip", col_lower)) {
    pl$String
  } else {
    pl$String
  }
})
names(dtypes_list) <- col_names

# read data as a lazy dataframe, restricting to women and the 14 original states
df <- pl$scan_csv(paste0(raw_data_dir, "incid.csv"), dtypes = dtypes_list)$
  rename(
    first_elig_state = "state",
    first_elig_period_start_date_65 = "start_date", 
    first_elig_period_end_date_65 = "end_date",
    first_elig_zip = "zip"
  )$
  filter(pl$col("sex") == "F")$
  filter(pl$col("state")$is_in(c(
    "AL", "CA", "CO", "FL", "GA", "IL", "MA", "MD", "NC", "NY", "OH", "PA", "TX", "WA"
  )))$
  select(c(
    "bene_id", "race", "state", "dob", "dod", "zip", "start_date", "end_date", 
    "hiv_cms_typ2_dt", "any_cx_typ1_dt", "any_cx_typ2_dt", "breast_cx_typ2_dt",
    "myocard_infarct_typ1_dt", "cong_heart_fail_typ1_dt", "periph_vas_dis_typ1_dt", 
    "cereb_vas_dis_typ1_dt", "dementia_typ1_dt", "chron_pulm_dis_typ1_dt", 
    "rheumatic_dis_typ1_dt", "peptic_ulcer_typ1_dt", "mild_liv_dis_typ1_dt", 
    "diabt_wo_chron_complict_typ1_dt", "diabt_chron_complict_da_typ1_dt", 
    "hemiplegia_paraplegia_d_typ1_dt", "renal_dis_typ1_dt", "mod_sevr_liver_dis_typ1_dt"
    # "metastatic_typ1_dt", "any_malignancy_typ1_dt"
  ))



# manage data ------------------------------------------------------------------

# define baseline and key variables at baseline
df1 <- df$
  with_columns(
    # 6 calendar months after start date
    baseline = pl$col("start_date")$dt$offset_by("6mo"),
    # calculate age 30 & age 65 dates
    age_30_date = pl$col("dob")$dt$offset_by("30y"),
    age_65_date = pl$col("dob")$dt$offset_by("65y")
  )$
  with_columns(
    # move baseline to >= age 30
    baseline = pl$max_horizontal(pl$col("age_30_date"), pl$col("baseline"))
    # administrative censoring on Dec 31, 2019
    # end_date = pl$min_horizontal(pl$col("end_date"), as.Date("2019-12-31"))
  )$
  # define date of administrative censoring
  with_columns(
    admin_censor_date = pl$when(
      (pl$col("state") == "MA") & (pl$col("baseline")$dt$year() < 2011)
    )$
      then(pl$lit(as.Date("2010-12-31")))$
      when((pl$col("state") == "OH") & (pl$col("baseline")$dt$year() < 2011))$
      then(pl$lit(as.Date("2010-12-31")))$
      when((pl$col("state") == "GA") & (pl$col("baseline")$dt$year() < 2011))$
      then(pl$lit(as.Date("2010-12-31")))$
      when((pl$col("state") == "PA") & (pl$col("baseline")$dt$year() < 2013))$
      then(pl$lit(as.Date("2012-12-31")))$
      when((pl$col("state") == "MD") & (pl$col("baseline")$dt$year() < 2016))$
      then(pl$lit(as.Date("2015-12-31")))$
      #  NY 2012-2016 due to HIV data quality issues
      when((pl$col("state") == "NY") & (pl$col("baseline")$dt$year() < 2012))$
      then(pl$lit(as.Date("2011-12-31")))$
      otherwise(pl$lit(as.Date("2019-12-31")))
  )$
  with_columns(
    # calculate age at baseline
    age_base = (pl$col("baseline") - pl$col("dob"))$dt$total_days()$cast(pl$Float64) / 365.25,
    # get baseline year
    year_base = pl$col("baseline")$dt$year(),
    # define HIV status at baseline
    hiv = pl$when(pl$col("hiv_cms_typ2_dt")$is_null())$
      then(0L)$
      otherwise(pl$col("hiv_cms_typ2_dt") <= pl$col("baseline"))$
      cast(pl$Int8)
  )$
  # with_columns(
  #   # ensure the baseline age falls within the range after approximate calculation
  #   age_base = pl$when(pl$col("age_base") < 30 & pl$col("age_base")$round(0) == 30)$
  #     then(pl$lit(30))$
  #     when(pl$col("age_base") >= 65 & pl$col("age_base")$round(0) == 65)$
  #     then(pl$lit(64.999999))$
  #     otherwise(pl$col("age_base"))
  # )$
  # make exclusions related to data quality
  filter(
    !((pl$col("state") == "MA") & (pl$col("year_base") == 2011)) &
      !((pl$col("state") == "OH") & (pl$col("year_base") == 2011)) &
      !((pl$col("state") == "GA") & (pl$col("year_base")$is_in(c(2011, 2012)))) &
      !((pl$col("state") == "PA") & (pl$col("year_base") == 2013)) &
      !((pl$col("state") == "MD") & (pl$col("year_base") == 2016)) &
      # removing NY 2012-2016
      !((pl$col("state") == "NY") & (pl$col("year_base")$is_in(2012:2016)))
  )$
  # remove anyone with a prior cancer (any type)
  filter(
    !(
      pl$col("any_cx_typ1_dt")$is_not_null() & 
        (pl$col("any_cx_typ1_dt") <= pl$col("baseline"))
    )
  )$
  # only keep people whose death (if any), end of eligibility, age 65, or 
  # date of administrative censoring occurred after baseline
  filter(
    (pl$col("dod")$is_null() | (pl$col("dod") > pl$col("baseline"))) &
      (pl$col("end_date") > pl$col("baseline")) &
      (pl$col("age_65_date") > pl$col("baseline")) &
      (pl$col("admin_censor_date") > pl$col("baseline"))
  )


# define end of follow-up: follow eligible individuals until any cancer 
# (type 2 date) or any censoring criteria are met
df2 <- df1$with_columns(
  first_event_date = pl$min_horizontal(
    pl$col("any_cx_typ2_dt"),
    pl$col("dod"),
    pl$col("end_date"),
    pl$col("age_65_date"),
    pl$col("admin_censor_date")
  )
)$with_columns( # define status at the end of follow-up
  first_event = pl$when(pl$col("first_event_date") == pl$col("age_65_date"))$
    then(pl$lit("age65"))$
    when(
      pl$col("any_cx_typ2_dt")$is_not_null() & 
        pl$col("breast_cx_typ2_dt")$is_not_null() &
        (pl$col("first_event_date") == pl$col("any_cx_typ2_dt")) &
        (pl$col("first_event_date") == pl$col("breast_cx_typ2_dt"))
    )$
    then(pl$lit("breast_cx"))$
    when(
      pl$col("dod")$is_not_null() & 
        pl$col("first_event_date") == pl$col("dod")
    )$
    then(pl$lit("death"))$
    when(pl$col("first_event_date") == pl$col("end_date"))$
    then(pl$lit("end_period"))$
    when(
      pl$col("any_cx_typ2_dt")$is_not_null() & 
        (
          (pl$col("breast_cx_typ2_dt")$is_not_null() &
             (pl$col("first_event_date") != pl$col("breast_cx_typ2_dt"))) |
            (pl$col("breast_cx_typ2_dt")$is_null())
        )
    )$
    then(pl$lit("other_cx"))$
    when((pl$col("first_event_date") == pl$col("admin_censor_date")))$
    then(pl$lit("admin_censor"))$
    otherwise(pl$lit("check_your_data"))
)$with_columns(
  event_flag = pl$when(
    pl$col("first_event")$is_not_null() & (pl$col("first_event") == "breast_cx")
  )$then(1L)$
    otherwise(0L)$
    cast(pl$Int8),
  
  compete_event_flag = pl$when(
    pl$col("first_event")$is_not_null() & (pl$col("first_event") == "breast_cx")
  )$then(1L)$
    when(
      pl$col("first_event")$is_not_null() &
        (pl$col("first_event") == "death")
    )$then(2L)$
    otherwise(0L)$
    cast(pl$Int8),
  
  # event time in days
  event_time = (pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()
)$with_columns(
  # number of event intervals
  event_int = (pl$col("event_time")$cast(pl$Float64) / 90)$ceil()
)


# define baseline comorbidities
comorb_cols <- c(
  "myocard_infarct_typ1_dt", "cong_heart_fail_typ1_dt", "periph_vas_dis_typ1_dt", 
  "cereb_vas_dis_typ1_dt", "dementia_typ1_dt", "chron_pulm_dis_typ1_dt", 
  "rheumatic_dis_typ1_dt", "peptic_ulcer_typ1_dt", "mild_liv_dis_typ1_dt", 
  "diabt_wo_chron_complict_typ1_dt", "diabt_chron_complict_da_typ1_dt", 
  "hemiplegia_paraplegia_d_typ1_dt", "renal_dis_typ1_dt", "mod_sevr_liver_dis_typ1_dt"
)

base_comorb_cols <- paste0(comorb_cols, "_base")

flag_comorb_logic <- function(comorb_col_name, starting_date) {
  pl$when(pl$col(comorb_col_name)$is_null())$
    then(pl$lit(0L))$
    otherwise(pl$col(comorb_col_name) <= pl$col(starting_date))$
    cast(pl$Int8)
}

df3 <- df2$with_columns(
  lapply(comorb_cols, function(x) {
    flag_comorb_logic(x, "baseline")$alias(paste0(x, "_base"))
  })
)$with_columns(
  n_comorb_base = pl$sum_horizontal(base_comorb_cols) 
)$with_columns(
  n_comorb_base_cat = pl$when(pl$col("n_comorb_base") >= 2)$
    then(pl$lit("2+"))$
    otherwise(pl$col("n_comorb_base")$cast(pl$String))
)$collect(engine = "streaming")


# df3$select(c(
#   "baseline", "first_event_date", "event_time", "event_int", 
#   "first_event", "event_flag", "hiv"
# ))$describe()



# sample data ------------------------------------------------------------------

samp.hiv <- df3$filter(pl$col("hiv") == 1) # 100% PWH

samp.nohiv <- df3$filter(pl$col("hiv") == 0)$ # 25% PWoH
  sample(fraction = 0.25, with_replacement = FALSE, seed = 123)

samp.df3 <- pl$concat(list(samp.hiv, samp.nohiv), how = "vertical_relaxed")$
  with_columns(
    wt_s = pl$when(pl$col("hiv") == 1)$then(1)$otherwise(1/0.25)
  )$sort("bene_id")

samp.df3$group_by("event_flag")$agg(pl$col("bene_id")$count())

samp.df3$group_by("hiv")$agg(pl$col("event_flag")$filter(pl$col("event_flag") == 1)$count())



# export -----------------------------------------------------------------------

df3$sort("bene_id")$write_ipc(paste0(data_dir, "dat.ex.ny.w.arrow"))
samp.df3$write_ipc(paste0(data_dir, "samp.dat.ex.ny.w.arrow"))


