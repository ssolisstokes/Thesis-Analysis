library(dplyr)
library(tidyr)
library(spOccupancy)
library(coda)
library(stars)
library(MCMCvis)
library(unmarked)

list.of.packages <- c("ggplot2", "corrplot", "gridExtra", "dplyr", "unmarked", "lubridate", "tibble", "sf", "MuMIn", "AICcmodavg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ---- Data Wrangling ----
# bring in master data
gm <- read.csv("data/data_gmocc_all-ut.csv", stringsAsFactors = FALSE)

# selecting only Paradise Canyon
gm <- gm %>%
  filter(site_id == "paradise_canyon")

# remove off transect detections 
gm$tran_id <- trimws(gm$tran_id) # fixing any spacing issues in tran_id
gm <- gm[!grepl("OTD", gm$tran_id, ignore.case = TRUE), ]

# bring in the transect length data
t.length <- read.csv("data/data_gmocc_tran-distance.csv", stringsAsFactors = FALSE)
t.length <- t.length %>%
  filter(site_id == "paradise_canyon")

# remove off transect detections
t.length$tran_id <- trimws(t.length$tran_id)
t.length <- t.length[!grepl("OTD", t.length$tran_id, ignore.case = TRUE), ]

# QAQC to remove duplicates if needed
#t.length <- t.length[!duplicated(t.length$tran_id), ]

# calculating time since midnight to use time as a variable
time <- strsplit(ifelse(is.na(gm$time), "NA:NA:NA", gm$time), ":", fixed = TRUE)
hh <- as.numeric(sapply(time, `[`, 1))
mm <- as.numeric(sapply(time, `[`, 2))
ss <- as.numeric(sapply(time, `[`, 3))
gm$time_sec <- hh*3600 + mm*60 + ss
gm$time_sec[!is.finite(gm$time_sec)] <- NA_real_

# dealing with transects/extra collection points for detections
transects  <- subset(gm, col_point %in% 1:5)  # getting only the 1-5 collection points per transect

# adding the transect length data
transects  <- merge(transects, t.length[, c("tran_id","path_km")], by = "tran_id", all.x = TRUE)
# sorted transect ids
all_tran  <- sort(unique(transects$tran_id))

# getting transects with detections
det_transects <- subset(gm, col_point == 6 & detection == 1)

# wide table of time at each collection point per transect
pts_wide <- transects %>%
  select(tran_id, col_point, time_sec) %>%
  distinct(tran_id, col_point, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_point, values_from = time_sec, names_prefix = "time_") %>%
  arrange(tran_id)

# assigning detection to the nearest time slot stop & creating a detection map 
det_map <- det_transects %>%
  select(site_id, tran_id, time_sec) %>%
  left_join(pts_wide, by = "tran_id") %>%
  rowwise() %>%
  mutate(
    rep_from_6 =
      if (!is.na(time_sec) && !is.na(time_1) && !is.na(time_2) &&
          ((time_1 <= time_sec && time_sec <= time_2) || (time_2 <= time_sec && time_sec <= time_1))) 1 else
            if (!is.na(time_sec) && !is.na(time_2) && !is.na(time_3) &&
                ((time_2 <= time_sec && time_sec <= time_3) || (time_3 <= time_sec && time_sec <= time_2))) 2 else
                  if (!is.na(time_sec) && !is.na(time_3) && !is.na(time_4) &&
                      ((time_3 <= time_sec && time_sec <= time_4) || (time_4 <= time_sec && time_sec <= time_3))) 3 else
                        if (!is.na(time_sec) && !is.na(time_4) && !is.na(time_5) &&
                            ((time_4 <= time_sec && time_sec <= time_5) || (time_5 <= time_sec && time_sec <= time_4))) 4) %>%
  ungroup() %>%
  mutate(det_rep = pmin(rep_from_6, 5L)) %>%
  select(tran_id, det_rep) %>%
  distinct()


# creating a transect map
tran_map <- transects %>%
  mutate(det_point = as.integer(detection == 1)) %>%
  select(tran_id, col_point, det_point) %>%
  group_by(tran_id, col_point) %>%
  summarise(det_point = as.integer(any(det_point == 1)), .groups = "drop")


# merging maps to make long df of transects
# make col_points replicates
y_long <- expand.grid(tran_id = all_tran, rep_idx = 1:5, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(tran_map %>% 
              rename(rep_idx = col_point), by = c("tran_id","rep_idx")) %>%
  left_join(det_map, by = "tran_id") %>%
  mutate(
    detection = ifelse(!is.na(det_rep) & det_rep == rep_idx, 1L,
                       ifelse(is.na(det_point), 0L, det_point))) %>%
  select(tran_id, rep_idx, detection) %>%
  arrange(tran_id, rep_idx)

# create wide df (transects as rows/replicates as columns)
y_wide <- y_long %>%
  pivot_wider(names_from = rep_idx, values_from = detection, names_prefix = "rep_") %>%
  arrange(tran_id)

# create full detection matrix; transects as sites, col_points as replicates
y <- as.matrix(y_wide %>% select(starts_with("rep_")))
rownames(y) <- y_wide$tran_id


# ---- Site Covariates ----
# all site covariates (occupancy)
site_covs <- transects %>%
  group_by(tran_id) %>%
  summarise(
    rockiness_mean = mean(rock_ind, na.rm = TRUE),  
    shrub_density_mean = mean(shr_pt_density, na.rm = TRUE),  
    elevation_mean = mean(elev, na.rm = TRUE),  
    prey_total = sum(pp_total, na.rm = TRUE),   
    #col_points = n_distinct(col_point),        # total col_points; future effort analysis?
    path_km = first(path_km),                      
    tran_utm_x = mean(utm_x, na.rm = TRUE),     # mean of utms to create a per-transect centroid for the spOcc model
    tran_utm_y = mean(utm_y, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(tran_id)

# align with detection matrix (y) rows
site_covs <- site_covs %>% 
  slice(match(rownames(y), tran_id))

# prey as rate per km 
#site_covs$prey_rate_km <- site_covs$prey_total / site_covs$path_km

# scale covariates
site_covs$rockiness_z <- as.numeric(scale(site_covs$rockiness_mean))
site_covs$shrubs_z <- as.numeric(scale(site_covs$shrub_density_mean))
site_covs$elev_z <- as.numeric(scale(site_covs$elevation_mean))
site_covs$prey_z <- as.numeric(scale(site_covs$prey_total))


# ---- Detection covariates ----
# all observation covariates (detection)
obs_covs <- expand.grid(tran_id = all_tran, col_point = 1:5, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(
    transects %>%
      select(tran_id, col_point, time_sec, air_temp, gr_temp, avg_wind_spd, vpd, rock_ind) %>%
      group_by(tran_id, col_point) %>%
      summarise(
        time_sec = mean(time_sec, na.rm = TRUE),
        air_temp = mean(air_temp, na.rm = TRUE),
        gr_temp = mean(gr_temp, na.rm = TRUE),
        wind = mean(avg_wind_spd, na.rm = TRUE),
        vpd = mean(vpd, na.rm = TRUE),
        rockiness = mean(rock_ind, na.rm = TRUE),
        .groups = "drop"),
    by = c("tran_id","col_point")) %>%
  arrange(tran_id, col_point)

# wide dfs per observation (detection) covariate
obs_time <- obs_covs %>% select(tran_id, col_point, time_sec) %>% pivot_wider(names_from = col_point, values_from = time_sec, names_prefix = "rep_") %>% arrange(tran_id)
obs_air <- obs_covs %>% select(tran_id, col_point, air_temp) %>% pivot_wider(names_from = col_point, values_from = air_temp, names_prefix = "rep_") %>% arrange(tran_id)
obs_gr <- obs_covs %>% select(tran_id, col_point, gr_temp) %>% pivot_wider(names_from = col_point, values_from = gr_temp, names_prefix = "rep_") %>% arrange(tran_id)
obs_wind <- obs_covs %>% select(tran_id, col_point, wind) %>% pivot_wider(names_from = col_point, values_from = wind, names_prefix = "rep_") %>% arrange(tran_id)
obs_vpd <- obs_covs %>% select(tran_id, col_point, vpd) %>% pivot_wider(names_from = col_point, values_from = vpd, names_prefix = "rep_") %>% arrange(tran_id)
obs_rock <- obs_covs %>% select(tran_id, col_point, rockiness)%>% pivot_wider(names_from = col_point, values_from = rockiness, names_prefix = "rep_") %>% arrange(tran_id)

# matrices 
time_mat <- as.matrix(obs_time %>% select(starts_with("rep_")))[match(rownames(y), obs_time$tran_id), ]
air_mat <- as.matrix(obs_air %>% select(starts_with("rep_")))[match(rownames(y), obs_air$tran_id), ]
gr_mat <- as.matrix(obs_gr %>% select(starts_with("rep_")))[match(rownames(y), obs_gr$tran_id), ]
wind_mat <- as.matrix(obs_wind %>% select(starts_with("rep_")))[match(rownames(y), obs_wind$tran_id), ]
vpd_mat <- as.matrix(obs_vpd %>% select(starts_with("rep_")))[match(rownames(y), obs_vpd$tran_id), ]
rock_mat <- as.matrix(obs_rock %>% select(starts_with("rep_")))[match(rownames(y), obs_rock$tran_id), ]

# z-score scale each detection matrix
time_mat[] <- (time_mat - mean(time_mat, na.rm = TRUE)) / sd(time_mat, na.rm = TRUE)
air_mat[] <- (air_mat - mean(air_mat, na.rm = TRUE)) / sd(air_mat, na.rm = TRUE)
gr_mat[] <- (gr_mat - mean(gr_mat, na.rm = TRUE)) / sd(gr_mat, na.rm = TRUE)
wind_mat[] <- (wind_mat - mean(wind_mat, na.rm = TRUE)) / sd(wind_mat, na.rm = TRUE)
vpd_mat[] <- (vpd_mat - mean(vpd_mat, na.rm = TRUE)) / sd(vpd_mat, na.rm = TRUE)
rock_mat[] <- (rock_mat - mean(rock_mat, na.rm = TRUE)) / sd(rock_mat, na.rm = TRUE)

# ---- Prep Covariates for spPGOcc ----
# store site covariates as a dataframe
occ_df <- data.frame(
  rockiness = site_covs$rockiness_z,
  shrubs = site_covs$shrubs_z,
  elev = site_covs$elev_z,
  prey = site_covs$prey_z)

# detection covariates: a named list of N x J matrices; names must match det.formula terms
det_df <- list(
  time = time_mat,
  vpd = vpd_mat,
  rockiness = rock_mat,
  t_air = air_mat,
  t_gr = gr_mat,
  wind = wind_mat)


umf <- unmarkedFrameOccu(y = y, siteCovs = occ_df, obsCovs = det_df)
summary(umf)

naive_occ <- sum(ifelse(rowSums(y, na.rm = TRUE)>0,1,0)/nrow(y))
naive_occ

null.model <- occu(formula = 
                     ~ 1 # detection
                   ~1, # occupancy
                   data = umf)

summary(null.model)
# null model occupancy = 6.04
# null model detection = -3.96

p<-backTransform(null.model, type = "det")
p
# the probability of detecting a Gila is ~ 0.02 (0.0186) (SE = 0.007) along a transect in paradise canyon

psi<- backTransform(null.model, type = "state")
psi
# there is roughly a 0.99 (99%?) chance that a Gila occupies a transect in paradise canyon 

global.model.um <- occu(formula = 
                          ~vpd + rockiness + I(t_air^2) + I(t_gr^2) + wind ~ 
                          ~ rockiness + shrubs + elev + prey,
                        data = umf)

summary(global.model.um)

mod.sel <- dredge(global.model.um)
mod.sel[1:10]

# top models results:
  # top model includes rockiness and avg. wind speed for det. covs
  # top model includes elevation, prey presence and rockiness for site (occupancy) covs
  # none of the top 10 models include shrub density
  # none of the top 10 models include vpd, or temperature data
  # 2nd best model: det vars: rockines/ occ vars: elev, prey presence, rockiness

occu.mb <- occu(formula = ~1 ~ rockiness + elev + prey, data = umf)
summary(occu.mb)

occu.psiPr <- predict(occu.mb, type = "state", appendData = TRUE)
head(occu.psiPr)
occu.psi <- occu(formula = ~rockiness + wind  ~rockiness + elev + prey, data = umf)
summary(occu.psi)

occu.psiPr <- predict(occu.psi,type = "state", appendData = TRUE)
head(occu.psiPr) 


fit <- occu(
  ~ rockiness + time + a.temp + g.temp + wind + effort + vpd + prey
  ~ shrub_density_site_mean + elev_site_mean,
  data = umf)

fit1 <- occu(~ rockiness ~ elev + shrubs, data = umf)
fit2 <- occu(~ 1 ~ elev, data = umf)
fit3 <- occu(~ rockiness + I(t_air^2) ~ prey, data = umf)


print(summary(fit1))
print(summary(fit2))
print(summary(fit3))




