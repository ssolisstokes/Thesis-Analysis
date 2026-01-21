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
library(sjPlot)

# bring in my dataset
gm <- read.csv("data/data_gmocc_all-ut.csv")

# converting date to julian
gm$date <- mdy(gm$date)
gm$julian_day <- yday(gm$date)

# create a df to work with
  # selecting predictor variables: rockiness, pp, shrub density, air temp, and elevation (5)
  # create a transect-level df
    # median rock_ind, pp_total, shr_den, mean air_temp, median elev, "sum" of detection across col_points
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

gm_df$site_id <- as.factor(gm_df$site_id)
gm_df$transect_id <- as.factor(gm_df$transect_id)
glimpse(gm_df)

# ---- Checking for colinearity ---
#check for co-linearity among the five predictor variables
pairs(gm_df[,4:11], lower.panel = NULL) 

ggpairs(gm_df, columns = 4:11) + theme_minimal()
## Elevation is highly correlated with shrub density and air temperature. Biologically this makes sense to me, with higher elevations being associated with more shrub density
# (vegetative diversity/communities), and cooler temperatures.

# look at coefficients of the Global Model 
model_global01 <- glm(detection ~ rockiness + elevation + air_temp + prey_presence + shrub_density + vpd + wind_spd, family = binomial, data = gm_df)
model_global01

model_global02 <- glm(detection ~ rockiness + elevation + air_temp + prey_presence + shrub_density + vpd + wind_spd + site_id, family = binomial, data = gm_df)
model_global02
# Variance Inflation Factor
check_collinearity(model_global01)
check_collinearity(model_global02)

options(na.action = "na.fail")   # ensure common rows

# ---- Writing Models ----

# old glm models
model_01 <- glm(detection ~ rockiness + elevation, family = binomial, data = gm_df)
model_02 <- glm(detection ~ rockiness + air_temp, family = binomial, data = gm_df)
model_03 <- glm(detection ~ rockiness + shrub_density, family = binomial, data = gm_df)
model_04 <- glm(detection ~ rockiness + prey_presence, family = binomial, data = gm_df)
model_05 <- glm(detection ~ rockiness + shrub_density + prey_presence, family = binomial, data = gm_df)
model_06 <- glm(detection ~ shrub_density + prey_presence, family = binomial, data = gm_df)
model_07 <- glm(detection ~ rockiness + elevation + air_temp, family = binomial, data = gm_df)

# re-writing old glms to include site_id
model_08 <- glm(detection ~ rockiness + elevation + site_id, family = binomial, data = gm_df)
model_09 <- glm(detection ~ rockiness + air_temp + site_id, family = binomial, data = gm_df)
model_10 <- glm(detection ~ rockiness + shrub_density + site_id, family = binomial, data = gm_df)
model_11 <- glm(detection ~ rockiness + prey_presence + site_id, family = binomial, data = gm_df)
model_12 <- glm(detection ~ rockiness + shrub_density + prey_presence + site_id, family = binomial, data = gm_df)
model_13 <- glm(detection ~ shrub_density + prey_presence + site_id, family = binomial, data = gm_df)
model_14 <- glm(detection ~ rockiness + elevation + air_temp + site_id, family = binomial, data = gm_df)

# experimenting with site id
model_15 <- glm(detection ~ site_id, family = binomial, data = gm_df)
model_16 <- glm(detection ~ site_id + rockiness, family = binomial, data = gm_df)
model_17 <- glm(detection ~ site_id + shrub_density, family = binomial, data = gm_df)
model_18 <- glm(detection ~ site_id + prey_presence, family = binomial, data = gm_df)

# air temp quadratic
model_19 <- glm(detection ~ air_temp, family = binomial, data = gm_df)
model_20 <- glm(detection ~ I(air_temp^2), family = binomial, data = gm_df)

# old glmm models
model_21 <- glmer(detection ~ rockiness + elevation + (1|site_id/transect_id), family = binomial, data = gm_df)
model_22 <- glmer(detection ~ rockiness + air_temp + (1|site_id/transect_id), family = binomial, data = gm_df)
model_23 <- glmer(detection ~ rockiness + shrub_density + (1|site_id/transect_id), family = binomial, data = gm_df)
model_24 <- glmer(detection ~ rockiness + prey_presence + (1|site_id/transect_id), family = binomial, data = gm_df)
model_25 <- glmer(detection ~ rockiness + shrub_density + prey_presence + (1|site_id/transect_id), family = binomial, data = gm_df)
model_26 <- glmer(detection ~ shrub_density + prey_presence + (1|site_id/transect_id), family = binomial, data = gm_df)
model_27 <- glmer(detection ~ rockiness + elevation + air_temp + (1|site_id/transect_id), family = binomial, data = gm_df)

# rewriting old glms with new random effect
model_28 <- glmer(detection ~ rockiness + elevation + site_id + (1|transect_id), family = binomial, data = gm_df)
model_29 <- glmer(detection ~ rockiness + air_temp + site_id + (1|transect_id), family = binomial, data = gm_df)
model_30 <- glmer(detection ~ rockiness + shrub_density + site_id + (1|transect_id), family = binomial, data = gm_df)
model_31 <- glmer(detection ~ rockiness + prey_presence + site_id + (1|transect_id), family = binomial, data = gm_df)
model_32 <- glmer(detection ~ rockiness + shrub_density + site_id + (1|transect_id), family = binomial, data = gm_df)
model_33 <- glmer(detection ~ shrub_density + prey_presence + site_id + (1|transect_id), family = binomial, data = gm_df)
model_34 <- glmer(detection ~ rockiness + elevation + air_temp + site_id + (1|transect_id), family = binomial, data = gm_df)

# airtemp GLMMs
model_35 <- glmer(detection ~ air_temp + site_id + (1|transect_id), family = binomial, data = gm_df)
model_36 <- glmer(detection ~ I(air_temp^2) + site_id + (1|transect_id), family = binomial, data = gm_df)


# ---- Model Selection ----
# Compare models
model_comparison <- model.sel(model_01, model_02, model_03, model_04, model_05, model_06, model_07,
                              model_08, model_09, model_10, model_11, model_12, model_13, model_14,
                              model_15, model_16, model_17, model_18, model_19, model_20, model_21,
                              model_21, model_22, model_23, model_24, model_25, model_26, model_27, 
                              model_28, model_29, model_30, model_31, model_32, model_33, model_34,
                              model_35, model_36)
model_comparison

model_comparison2 <- model.sel(model_01, model_02, model_03, model_04, model_05, model_06, model_07,
                               model_08, model_09, model_10, model_11, model_12, model_13, model_14,
                               model_15, model_16, model_17, model_18, model_19, model_20, model_21,
                               model_21, model_22, model_23, model_24, model_25, model_27, 
                               model_28, model_29, model_30, model_31, model_32, model_33, model_34,
                               model_35, model_36)
model_comparison2
# model_02 (rock + air temp) best describes the data; followed by model_07 (rock, air temp, elevation)
# which variables are most influential in these models? 

topmodels <- model.sel(model_09, model_14, model_29, model_34)
sw(topmodels)

model.avg(topmodels, revised.var = TRUE)

# ---- Interpretation ----
# Can now report the model averaged coefficients for the predictor variables individual effects on detection probability.
# Used to average regression coefficients across multiple models with the ultimate goal of capturing a variable’s overall “effect.”
summary(model.avg(model_comparison2, subset = delta < 4)) 


model_09 <- glm(detection ~ rockiness + air_temp + site_id, family = binomial, data = gm_df)
# ---- Plots ----
p1 <- plot_model(model_09, type = "pred", terms = "rockiness")
p2 <- plot_model(model_09, type = "pred", terms = "air_temp")
p3 <- plot_model(model_09, type = "pred", terms = "site_id")

p1 + p2 + p3 + plot_layout(ncol = 3)


p_raw1 <- ggplot(gm_df, aes(rockiness, detection)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Raw data: Rockiness vs Detection")

p_raw2 <- ggplot(gm_df, aes(air_temp, detection)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Raw data: Air Temp vs Detection")

p_raw3 <- ggplot(gm_df, aes(site_id, detection)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Raw data: Site vs Detection")

p_raw1 + p_raw2 + p_raw3 + plot_layout(ncol = 3)

