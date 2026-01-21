## Model selection Assignment 
# SSS 6 OCTOBER 2025

# ---- Set up ----
# libraries
library(ggplot2)
library(patchwork)
library(MuMIn)
library(easystats)
library(arm)
library(ggfortify)
library(grid)
library(readr)
library(dplyr)
library(AICcmodavg)
library(performance)
library(GGally)
library(olsrr)
library(lme4)
library(modelsummary)
library(broom.mixed)
library(lubridate)
library(psych)
library(tidyverse)
library(unmarked)

# bring in my dataset
gm <- read.csv("data/data_gmocc_all-ut.csv")

# ---- Data Prep ----
# converting date to julian
gm$date <- mdy(gm$date)
gm$julian_day <- yday(gm$date)

# calculating vars means per transect
gm_df <- gm %>%
  group_by(tran_id) %>%
  select(julian_day, site_id, tran_id, col_point, avg_wind_spd, vpd, rock_ind, shr_m_density, pp_total, elev, gr_temp, air_temp, detection) %>%
  mutate(m_air_temp = mean(air_temp, na.rm = TRUE),
         m_gr_temp = mean(gr_temp, na.rm = TRUE),
         m_rock = mean(rock_ind, na.rm = TRUE),
         m_elev = mean(elev, na.rm = TRUE),
         m_vpd = mean(vpd, na.rm = TRUE),
         m_wind = mean(avg_wind_spd, na.rm = TRUE),
         detection = as.integer(sum(detection, na.rm = TRUE) > 0)) %>%
  ungroup()

gm_df <- gm_df %>%
  group_by(tran_id) %>%
  slice_head(n = 1) %>%
  select(julian_day, site_id, tran_id, m_wind, m_vpd, m_rock, shr_m_density, pp_total, m_elev, m_gr_temp, m_air_temp, detection) %>%
  ungroup()

gm_df <- gm_df %>%
  rename(
    rockiness = m_rock,
    elevation = m_elev,
    air_temp = m_air_temp,
    gr_temp = m_gr_temp,
    prey_presence = pp_total,
    shrub_density = shr_m_density,
    transect_id = tran_id,
    vpd = m_vpd,
    wind_spd = m_wind) %>%
  arrange(julian_day)

glimpse(gm_df)

# ---- Detection Histories df ----
# rename transect IDs to "replicates" per site
tran_rename <- gm_df %>%
  group_by(site_id) %>%
  arrange(transect_id, .by_group = TRUE) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

# set max number of transects per site (should be 89)
Jmax <- tran_rename %>% 
  count(site_id) %>% 
  pull(n) %>% 
  max()

# wide df of detction histories (rows = site, cols = transect replicates)
y_wide <- tran_rename %>%
  select(site_id, rep_id, detection) %>%
  pivot_wider(names_from = rep_id, values_from = detection, names_prefix = "rep_")

# pad sites without max number transects with NAs
missing_cols <- setdiff(paste0("rep_", 1:Jmax), names(y_wide))

if (length(missing_cols)) y_wide[missing_cols] <- NA_integer_

# final wide df
df_det.history <- y_wide %>% 
  select(site_id, paste0("rep_", 1:Jmax))

# make an unmarked y matrix
mat.y <- as.matrix(y_wide[, -1, drop = FALSE])  

# ---- Covariate dfs ----
# air temp
atemp_wide <- tran_rename %>%
  select(site_id, rep_id, air_temp) %>%
  pivot_wider(names_from = rep_id, values_from = air_temp, names_prefix = "atemp_")
# make a matrix for unmarked
atemp <- as.matrix(atemp_wide[, -1, drop = FALSE]) 

# ground temp
gtemp_wide <- tran_rename %>%
  select(site_id, rep_id, gr_temp) %>%
  pivot_wider(names_from = rep_id, values_from = gr_temp, names_prefix = "gtemp_")
# make a matrix for unmarked
mat.gtemp <- as.matrix(gtemp_wide[, -1, drop = FALSE]) 

# vpd
vpd_wide <- tran_rename %>%
  select(site_id, rep_id, vpd) %>%
  pivot_wider(names_from = rep_id, values_from = vpd, names_prefix = "vpd_")
# make a matrix for unmarked
mat.vpd <- as.matrix(vpd_wide[, -1, drop = FALSE]) 

# shrub density
shrubden_wide <- tran_rename %>%
  select(site_id, rep_id, shrub_density) %>%
  pivot_wider(names_from = rep_id, values_from = shrub_density, names_prefix = "shrubden_")
# make a matrix for unmarked
mat.shrubden <- as.matrix(shrubden_wide[, -1, drop = FALSE]) 

# prey presence
prey_wide <- tran_rename %>%
  select(site_id, rep_id, prey_presence) %>%
  pivot_wider(names_from = rep_id, values_from = prey_presence, names_prefix = "pp_")
# make a matrix for unmarked
mat.prey <- as.matrix(prey_wide[, -1, drop = FALSE]) 

# rockiness
rock_wide <- tran_rename %>%
  select(site_id, rep_id, rockiness) %>%
  pivot_wider(names_from = rep_id, values_from = rockiness, names_prefix = "rock_")
# make a matrix for unmarked
mat.rock <- as.matrix(rock_wide[, -1, drop = FALSE]) 

# elevation
elev_wide <- tran_rename %>%
  select(site_id, rep_id, elevation) %>%
  pivot_wider(names_from = rep_id, values_from = elevation, names_prefix = "elev_")
# make a matrix for unmarked
at.elev <- as.matrix(elev_wide[, -1, drop = FALSE]) 

# wind speed
wind_wide <- tran_rename %>%
  select(site_id, rep_id, wind_spd) %>%
  pivot_wider(names_from = rep_id, values_from = wind_spd, names_prefix = "wind_")
# make a matrix for unmarked
mat.wind <- as.matrix(wind_wide[, -1, drop = FALSE]) 

# date
julian_wide <- tran_rename %>%
  select(site_id, rep_id, julian_day) %>%
  pivot_wider(names_from = rep_id, values_from = julian_day, names_prefix = "julian_")
# make a matrix for unmarked
mat.julian <- as.matrix(julian_wide[, -1, drop = FALSE]) 


# ---- Combined dfs ----
list_sitecovs <- list(shrubden_wide, prey_wide, rock_wide, elev_wide)
list_obscovs <- list(atemp_wide, gtemp_wide, vpd_wide, julian_wide, wind_wide)

df_sitecovs <- list_sitecovs %>%
  reduce(full_join, by = "site_id")

df_obscovs <- list_obscovs %>%
  reduce(full_join, by = "site_id")

df_covs <- full_join(df_sitecovs,df_obscovs, by = "site_id")

# combine det histories with covariate df
df_gmocc <- full_join(df_det.history, df_covs, by = "site_id")


# ---- Creating an unmarked frame ----
site <- y_wide[, "site_id"]

umf <- unmarkedFrameOccu(
  y = mat.y, siteCovs = list(df_sitecovs, obsCovs = df_obscovs)
