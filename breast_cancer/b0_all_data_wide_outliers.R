
# Project: Incidence of breast cancer by HIV status among female Medicaid beneficiaries aged 30-64

# Content: Preliminary data management (overall data); identification of outlier states and years

# Date: 7/2/2025

# Author: Yiyi Zhou

################################################################################

library(polars)
library(ggplot2)
library(exactci) # not necessary if we don't estimate CIs for the crude rates (I'm just too lazy to change the code)

raw_data_dir <- "/cms01/data/dua/57285/medican/data_requests/karine/20250422_cancer_incidence_v11/"
data_dir <- "/cms01/data/dua/57285/users/c-yzhou172-57285/breast_cancer/interim_data/"
figure_dir <- "/users/57285/c-yzhou172-57285/breast_cancer/figures/"



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

# read data as a lazy dataframe, restricting to women
df <- pl$scan_csv(paste0(raw_data_dir, "incid.csv"), dtypes = dtypes_list)$
  rename(
    first_elig_state = "state",
    first_elig_period_start_date_65 = "start_date", 
    first_elig_period_end_date_65 = "end_date",
    first_elig_zip = "zip"
  )$
  filter(pl$col("sex") == "F")$
  select(c(
    "bene_id", "sex", "race", "state", "dob", "dod", "zip",
    "start_date", "end_date", "hiv_cms_typ2_dt", 
    "any_cx_typ1_dt", "any_cx_typ2_dt", "breast_cx_typ2_dt",
    "myocard_infarct_typ1_dt", "cong_heart_fail_typ1_dt", 
    "periph_vas_dis_typ1_dt", "cereb_vas_dis_typ1_dt", 
    "dementia_typ1_dt", "chron_pulm_dis_typ1_dt", 
    "rheumatic_dis_typ1_dt", "peptic_ulcer_typ1_dt", 
    "mild_liv_dis_typ1_dt", "diabt_wo_chron_complict_typ1_dt", 
    "diabt_chron_complict_da_typ1_dt", "hemiplegia_paraplegia_d_typ1_dt", 
    "renal_dis_typ1_dt", "mod_sevr_liver_dis_typ1_dt"
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
    baseline = pl$max_horizontal(pl$col("age_30_date"), pl$col("baseline")),
    # administrative censoring on Dec 31, 2019
    end_date = pl$min_horizontal(pl$col("end_date"), as.Date("2019-12-31"))
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
  # remove anyone with a prior cancer (any type)
  filter(
    !(
      pl$col("any_cx_typ1_dt")$is_not_null() & 
        (pl$col("any_cx_typ1_dt") <= pl$col("baseline"))
    )
  )$
  # only keep people whose death (if any), end of eligibility, age 65 occurred after baseline
  filter(
    (pl$col("dod")$is_null() | (pl$col("dod") > pl$col("baseline"))) &
      (pl$col("end_date") > pl$col("baseline")) &
      (pl$col("age_65_date") > pl$col("baseline"))
  )


# define end of follow-up: follow eligible individuals until any cancer 
# (type 2 date) or any censoring criteria are met
df2 <- df1$with_columns(
  first_event_date = pl$min_horizontal(
    pl$col("any_cx_typ2_dt"),
    pl$col("dod"),
    pl$col("end_date"),
    pl$col("age_65_date")
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
    cast(pl$Int8)
)$collect()


# estimate crude incidence rates of breast cancer by year and state ------------

# function to stratify data by year
year_strat <- function(data, year_lower, year_upper) {
  
  lower_date <- as.Date(paste0(year_lower, "-01-01"))
  upper_date <- as.Date(paste0(year_upper, "-01-01"))
  
  year.dat <- data$filter(
    pl$col("baseline")$dt$year() < year_upper &
      pl$col("first_event_date")$dt$year() >= year_lower
  )$with_columns(
    pl$col("baseline")$alias("baseline0"),
    pl$col("first_event_date")$alias("first_event_date0"),
    
    pl$when(pl$col("baseline")$gt(lower_date))$
      then(pl$col("baseline"))$
      otherwise(lower_date)$
      alias("baseline"),
    
    (pl$col("first_event_date")$dt$year() >= year_upper)$cast(pl$Int8)$alias("index")
  )$with_columns(
    pl$when(pl$col("index") == 1)$
      then(upper_date)$
      otherwise(pl$col("first_event_date"))$
      alias("first_event_date"),
    
    pl$when(pl$col("index") == 1)$
      then(0)$
      otherwise(pl$col("event_flag"))$
      alias("event_flag"),
    
    pl$lit(paste0("[", year_lower, ", ", year_upper, ")"))$alias("year_cat")
  )$drop("index")
  
  return(year.dat)
}

# function to put year-stratified data into a list
create_year_stratified_base <- function(data, year_cats) {
  year_strat_list <- list()
  
  for (year_cat in year_cats) {
    year_lower <- year_cat[1]
    year_upper <- year_cat[2]
    year_strat_list[[as.character(year_lower)]] <- year_strat(data, year_lower, year_upper)
  }
  
  return(year_strat_list)
}

# function to calculate crude IR
calc_crude_ir <- function(ir_dat) {
  
  ir_dat <- ir_dat$with_columns(
    ((pl$col("first_event_date") - pl$col("baseline"))$dt$total_days()$
       cast(pl$Float64) / 365.25)$alias("followup_years")
  )$with_columns(
    pl$when(pl$col("followup_years") == 0)$
      then(0.0000001)$
      otherwise(pl$col("followup_years"))$
      alias("followup_years")
  )$select(
    pl$sum("event_flag")$cast(pl$Int64)$alias("total_cases"),
    pl$sum("followup_years")$cast(pl$Float64)$alias("total_person_years")
  )$with_columns(
    (pl$col("total_cases") / pl$col("total_person_years"))$cast(pl$Float64)$alias("crude_incidence")
  )
  
  return(ir_dat)
}

# function to calculate crude IR over time stratified by time-fixed variables
# (sex, race, expansion status, birth cohort, and state)
calculate_stratified_metrics <- function(year_stratified_data, strata_var = NULL, 
                                         strata_levels = NULL, expansion_analysis = FALSE,
                                         time_updated_expansion = FALSE,
                                         years_sc_expansion_analysis = FALSE,
                                         ex_hiv_liver_dis = FALSE,
                                         restriction_analysis = FALSE,
                                         relaxed_restriction = FALSE,
                                         fibrosis_or_sobriety = NULL) {
  ir_year_strat <- list()
  ir_year_strat_formatted <- list()
  
  for (year in names(year_stratified_data)) {
    year.dat <- year_stratified_data[[year]]
    if (ex_hiv_liver_dis) {
      year.dat <- year.dat$filter(
        pl$col("hiv") == "0" & pl$col("mld") == "0" & pl$col("msld") == "0"
      )
    }
    year_lower <- as.numeric(year)
    
    if (expansion_analysis) {
      year.dat <- year.dat$with_columns(
        pl$lit(year_lower)$cast(pl$Int32)$alias("year")
      )
      
      expansion_map <- list(
        "2014" = c("AZ", "AR", "CA", "CO", "IL", "KY", "MD", "MA", "MI", 
                   "NV", "NJ", "NY", "OH", "WA", "DC"),
        "2015" = c("IN", "PA"),
        "2021" = c("MO", "OK"),
        "2023" = c("GA", "NC"),
        "2025" = c("MS"),
        "NA" = c("AL", "FL", "SC", "TN", "TX")
      )
      
      if (time_updated_expansion) {
        # time-updated expansion status
        year.dat <- year.dat$with_columns(
          pl$when(pl$col("state")$is_in(expansion_map[["NA"]]))$
            then(pl$lit("not expanded"))$
            when(pl$col("state")$is_in(expansion_map[["2014"]]) & pl$col("year") >= 2014)$
            then(pl$lit("expanded"))$
            when(pl$col("state")$is_in(expansion_map[["2015"]]) & pl$col("year") >= 2015)$
            then(pl$lit("expanded"))$
            when(pl$col("state")$is_in(expansion_map[["2021"]]) & pl$col("year") >= 2021)$
            then(pl$lit("expanded"))$
            when(pl$col("state")$is_in(expansion_map[["2023"]]) & pl$col("year") >= 2023)$
            then(pl$lit("expanded"))$
            when(pl$col("state")$is_in(expansion_map[["2025"]]) & pl$col("year") >= 2025)$
            then(pl$lit("expanded"))$
            otherwise(pl$lit("not expanded"))$
            alias("expansion")
        )
        
        strata_var <- "expansion"
        strata_levels <- c("expanded", "not expanded")
        
      } else {
        # time-fixed expansion status (ever expanded during 2001-2019 vs. never expanded during the period)
        year.dat <- year.dat$with_columns(
          pl$when(pl$col("state")$is_in(expansion_map[["NA"]]))$
            then(NA_real_)$
            when(pl$col("state")$is_in(expansion_map[["2014"]]))$
            then(2014)$
            when(pl$col("state")$is_in(expansion_map[["2015"]]))$
            then(2015)$
            when(pl$col("state")$is_in(expansion_map[["2021"]]))$
            then(2021)$
            when(pl$col("state")$is_in(expansion_map[["2023"]]))$
            then(2023)$
            when(pl$col("state")$is_in(expansion_map[["2025"]]))$
            then(2025)$
            otherwise(NA_real_)$
            alias("year_of_expansion")
        )$with_columns(
          pl$when(pl$col("year_of_expansion")$is_not_null() &
                    pl$col("year_of_expansion") <= 2019)$
            then(pl$lit("ever expanded"))$
            otherwise(pl$lit("never expanded"))$
            alias("ever_expansion")
        )$with_columns(
          pl$when(pl$col("ever_expansion") == "ever expanded")$
            then(pl$col("year") - pl$col("year_of_expansion"))$
            otherwise(NA_real_)$
            alias("years_sc_expansion")
        )
        
        if (years_sc_expansion_analysis) {
          strata_var <- "years_sc_expansion"
          strata_levels <- seq(-14, 5, 1)
          year.dat <- year.dat$filter(pl$col("ever_expansion") == "ever expanded")
        } else {
          strata_var <- "ever_expansion"
          strata_levels <- c("ever expanded", "never expanded")
        }
      }
      
    } else if (restriction_analysis) {
      year.dat <- year.dat$with_columns(
        pl$lit(year_lower)$cast(pl$Int32)$alias("year")
      )
      if (!relaxed_restriction) { # conservative approach
        df_restrict_conservative <- pl$read_csv("/users/57285/c-yzhou172-57285/hepc/simplified_hepc_restrictions_cleaned_conservative_long.csv")$
          with_columns(pl$col("year")$cast(pl$Int32))
        
        year.dat <- year.dat$join(
          other = df_restrict_conservative,
          on = c("state", "year"),
          how = "left"
        )
        
      } else { # relaxed
        df_restrict_relaxed <- pl$read_csv("/users/57285/c-yzhou172-57285/hepc/simplified_hepc_restrictions_cleaned_relaxed_long.csv")$
          with_columns(pl$col("year")$cast(pl$Int32))
        
        year.dat <- year.dat$join(
          other = df_restrict_relaxed,
          on = c("state", "year"),
          how = "left"
        )
      }
      
      if (is.null(fibrosis_or_sobriety)) {
        print("Check! Both fibrosis and sobriety are null!")
      } else {
        strata_var <- fibrosis_or_sobriety
        strata_levels <- c("Restricted", "No Restrictions")
      } 
    }
    
    if (!is.null(strata_var)) {
      for (level in strata_levels) {
        strat.dat <- year.dat$filter(pl$col(strata_var) == level)
        
        if (strat.dat$height > 0) {
          ir_result <- calc_crude_ir(strat.dat)
          ir_result_r <- ir_result$to_data_frame()
          
          poisson_result <- poisson.exact(
            x = ir_result_r$total_cases,
            T = ir_result_r$total_person_years,
            r = 1,
            alternative = "two.sided",
            conf.level = 0.95,
            midp = TRUE
          )
          
          ir_year_strat[[paste0(year, "_", level)]] <- ir_result$with_columns(
            pl$lit(poisson_result$estimate)$alias("est"),
            pl$lit(poisson_result$conf.int[1])$alias("lower"),
            pl$lit(poisson_result$conf.int[2])$alias("upper"),
            pl$lit(year)$alias("year"),
            pl$lit(level)$alias("stratum")
          )
          
          ir_year_strat_formatted[[paste0(year, "_", level)]] <- ir_result$with_columns(
            pl$col("total_person_years")$round(0)$alias("total_person_years"),
            (pl$lit(poisson_result$estimate) * 10000)$round(2)$alias("est"),
            (pl$lit(poisson_result$conf.int[1]) * 10000)$round(2)$alias("lower"),
            (pl$lit(poisson_result$conf.int[2]) * 10000)$round(2)$alias("upper"),
            pl$lit(year)$alias("year"),
            pl$lit(level)$alias("stratum")
          )$with_columns(
            pl$concat_str(
              pl$col("est")$cast(pl$Utf8),
              pl$lit(" ("),
              pl$col("lower")$cast(pl$Utf8),
              pl$lit(", "),
              pl$col("upper")$cast(pl$Utf8),
              pl$lit(")")
            )$alias("est_ci")
          )$select(
            "year",
            "stratum",
            "total_cases",
            "total_person_years",
            "est_ci"
          )
        }
      }
    }
  }
  
  list(
    combined_df = pl$concat(ir_year_strat),
    formatted_combined_df = pl$concat(ir_year_strat_formatted)
  )
}


year_cats <- list(
  c(2001, 2002), c(2002, 2003), c(2003, 2004), c(2004, 2005), c(2005, 2006),
  c(2006, 2007), c(2007, 2008), c(2008, 2009), c(2009, 2010), c(2010, 2011),
  c(2011, 2012), c(2012, 2013), c(2013, 2014), c(2014, 2015), c(2015, 2016),
  c(2016, 2017), c(2017, 2018), c(2018, 2019), c(2019, 2020)
)

year_stratified_base <- create_year_stratified_base(df2, year_cats)

ir_year_state = calculate_stratified_metrics(
  year_stratified_base,
  "state",
  c("AL", "AR", "AZ", "CA", "CO", 
    "DC", "FL", "GA", "IL", "IN",
    "KY", "MA", "MD", "MI", "MO",
    "MS", "NC", "NJ", "NV", "NY",
    "OH", "OK", "PA", "SC", "TN",
    "TX", "WA")
)

ir_year_state$combined_df$filter(pl$col("total_cases") <= 10)$height #0

# function to identify outliers, calculating thresholds across all years and states
identify_outliers <- function(data) {
  
  # Calculate log rates and thresholds
  data_log <- data$with_columns(
    (pl$col("est") + 1)$log()$alias("log_rate")
  )$with_columns(
    pl$col("log_rate")$mean()$alias("mean_log"),
    pl$col("log_rate")$std()$alias("sd_log")
  )$with_columns(
    (pl$col("log_rate") > (pl$col("mean_log") + 3 * pl$col("sd_log")) | 
       pl$col("log_rate") < (pl$col("mean_log") - 3 * pl$col("sd_log")))
    $alias("is_outlier")
  )
  
  # Extract outliers
  outliers <- data_log$filter(pl$col("is_outlier"))$sort(c("year", "state"))
  
  outliers <- outliers$with_columns(
    (pl$col("est") * 10000)$alias("crude_rate"),
    ((pl$col("mean_log")$exp() - 1) * 10000)$alias("mean_rate"),
    (((pl$col("mean_log") - 3 * pl$col("sd_log"))$exp() - 1) * 10000)$alias("lower_threshold"),
    (((pl$col("mean_log") + 3 * pl$col("sd_log"))$exp() - 1) * 10000)$alias("upper_threshold")
  )$select(
    "year",
    "state",
    "crude_rate",
    "total_cases",
    "total_person_years",
    "mean_rate",
    "lower_threshold",
    "upper_threshold"
  )
  
  return(list(
    outliers = outliers,
    data_with_flags = data_log
  ))
}

state_results <- ir_year_state$combined_df$rename(stratum = "state")
outlier_results <- identify_outliers(state_results)
outlier_results$outliers

# create a plot
plot_data <- ir_year_state$combined_df$
  with_columns(
    year = pl$col("year")$cast(pl$Int32),
    est = (pl$col("est") * 10000),
    upper = (pl$col("upper") * 10000),
    lower = (pl$col("lower") * 10000),
    state = pl$col("stratum")
  )$to_data_frame()


p <- ggplot(plot_data, aes(x = year, y = est, color = state, group = state)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(breaks = seq(2001, 2019, by = 1)) +
  scale_color_viridis_d(option = "turbo") + 
  labs(
    x = "Year",
    y = "Crude Incidence Rate (per 10,000 person-years)",
    color = "State"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  guides(color = guide_legend(nrow = 2))

jpeg(filename = paste0(figure_dir, "outlier_states_years.jpeg"),
     width = 9, height = 5.5, units = "in", res = 300)
print(p)
dev.off()
