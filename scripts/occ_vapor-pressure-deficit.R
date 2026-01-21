## Calculating and Adding Vapor Pressure Deficit w/ Master Data
# SSS - Fall 2025

# load in libraries
library(dplyr)

# bring in data set
gm_df <- read.csv("data/data_gmocc_all-ut.csv")

glimpse(gm_df)

# saturation vapor pressure (kPa) using FAO-56/Tetens
# defining the standard FAO-56/Tetens formula to compute saturation vapor pressure (es) from temperature
sat_vp_kPa <- function(Tc) 0.6108 * exp((17.27 * Tc) / (Tc + 237.3))

# calculate VPD at each col_point
gm_df <- gm_df %>%
  mutate(
    rh_input = rel_hum,                                 # rename RH column 
    rh_frac = case_when(                                # Normalize RH to a fraction
      is.na(rh_input) ~ NA_real_,
      rh_input > 1  ~ pmin(pmax(rh_input / 100, 0), 1), # clamping values to 0/1 for outliers
      TRUE ~ pmin(pmax(rh_input,0), 1)),
    es_kPa  = sat_vp_kPa(air_temp),                     # saturation vapor pressure at the measured air temperature
    vpd_kPa = pmax(es_kPa * (1 - rh_frac), 0))          # vapor pressure deficit at each row: VPD = es * (1 − RH)
                                                        # pmax prevents negative values from rounding/measurement noise

# summarize VPD by transect
vpd_by_transect <- gm_df %>%
  group_by(tran_id) %>%
  mutate(
    mean_vpd_kPa   = round(mean(vpd_kPa,   na.rm = TRUE), 3),
    median_vpd_kPa = round(median(vpd_kPa, na.rm = TRUE), 3))
  
# create a VPD subset dataframe
vpd_df <- vpd_by_transect %>%  
  select(date, site_id, tran_id, col_point, air_temp, rel_hum, rh_input, rh_frac, es_kPa, vpd_kPa, mean_vpd_kPa, median_vpd_kPa)

# export as subset csv
write.csv(vpd_df, "data/data_gmocc_vpd-ut.csv")

# select/rename columns to add to master datasheet
gm_vpd <- vpd_by_transect %>%
  select(-es_kPa, -rh_input, -rh_frac, -mean_vpd_kPa) %>%
  rename(vpd = vpd_kPa,
         med_vpd = median_vpd_kPa)

# update master data sheet
write.csv(gm_vpd, "data/data_gmocc_all-ut.csv")

