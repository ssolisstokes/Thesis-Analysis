# running through occupancy model options
# option 2: spatial occupancy (point level)

library(dplyr)
library(tidyr)
library(spOccupancy)

# ---- Data Wrangling ----
# bring in master data
gm <- read.csv("data/data_gmocc_all-ut.csv", stringsAsFactors = FALSE)

# remove off transect detections 
gm$tran_id <- trimws(gm$tran_id)
gm <- gm[!grepl("OTD", gm$tran_id, ignore.case = TRUE), ]

# bring in the transect length data
t.length <- read.csv("data/data_gmocc_tran-distance.csv", stringsAsFactors = FALSE)

# remove off transect detections
t.length$tran_id <- trimws(t.length$tran_id)
t.length <- t.length[!grepl("OTD", t.length$tran_id, ignore.case = TRUE), ]
t.length <- t.length[!duplicated(t.length$tran_id), ]

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
# getting transects with detections
det_transects <- subset(gm, col_point == 6 & detection == 1)

# attaching the detection to the nearest collection point

# wide table of time at each collection point per transect
pts_wide <- transects %>%
  select(tran_id, col_point, time_sec) %>%
  distinct(tran_id, col_point, .keep_all = TRUE) %>%
  pivot_wider(names_from = col_point, values_from = time_sec, names_prefix = "time_") %>%
  arrange(tran_id)

# assigning detection to the nearest time slot stop
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
  mutate(encounter_rep = pmin(rep_from_6, 5L)) %>%
  select(tran_id, encounter_rep) %>%
  distinct()

all_tran <- sort(unique(transects$tran_id))

# all collection 
nondet_tran <- transects %>%
  mutate(det_point = as.integer(detection == 1)) %>%
  select(tran_id, col_point, det_point) %>%
  group_by(tran_id, col_point) %>%
  summarise(det_point = as.integer(any(det_point == 1)), .groups = "drop")

# full grid and merge detection/nondetection transects
# making col_points replicates
y_long <- expand.grid(tran_id = all_tran, rep_idx = 1:5, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(nondet_tran %>% rename(rep_idx = col_point), by = c("tran_id","rep_idx")) %>%
  left_join(det_map, by = "tran_id") %>%
  mutate(
    det_point = ifelse(!is.na(encounter_rep) & encounter_rep == rep_idx, 1L,
                       ifelse(is.na(det_point), 0L, det_point))) %>%
  select(tran_id, rep_idx, det_point) %>%
  arrange(tran_id, rep_idx)

# create wide df
y_wide <- y_long %>%
  pivot_wider(names_from = rep_idx, values_from = det_point, names_prefix = "rep_") %>%
  arrange(tran_id)
# create full detection matrix; transects (500) as sites, col_points as replicates (5)
y <- as.matrix(y_wide %>% select(starts_with("rep_")))
rownames(y) <- y_wide$tran_id

# ---- Site Covariates ----
site_covs <- transects %>%
  group_by(site_id, tran_id) %>%
  summarise(
    rockiness_mean = mean(rock_ind, na.rm = TRUE),  
    shrub_density_mean = mean(shr_pt_density, na.rm = TRUE),  
    elevation_mean = mean(elev, na.rm = TRUE),  
    prey_total = sum(pp_total,na.rm = TRUE),   
    n_points_tran = n_distinct(col_point),
    path_km = first(path_km),                      
    tran_utm_x = mean(utm_x, na.rm = TRUE),     # mean of utms to create a per-transect centroid for the spOcc model
    tran_utm_y = mean(utm_y, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(tran_id)

# align with detection matrix rows
site_covs <- site_covs %>% slice(match(rownames(y), tran_id))

# prey as rate per km 
#site_covs$prey_rate_km <- site_covs$prey_total / site_covs$path_km

# scale site (occupancy) covariates
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
        vpd_obs = mean(vpd, na.rm = TRUE),
        rock_obs = mean(rock_ind, na.rm = TRUE),
        .groups = "drop"),
    by = c("tran_id","col_point")) %>%
  arrange(tran_id, col_point)

# wide matrices per observation (detection) covariate
obs_time <- obs_covs %>% select(tran_id, col_point, time_sec) %>% pivot_wider(names_from = col_point, values_from = time_sec, names_prefix = "rep_") %>% arrange(tran_id)
obs_air <- obs_covs %>% select(tran_id, col_point, air_temp) %>% pivot_wider(names_from = col_point, values_from = air_temp, names_prefix = "rep_") %>% arrange(tran_id)
obs_gr <- obs_covs %>% select(tran_id, col_point, gr_temp) %>% pivot_wider(names_from = col_point, values_from = gr_temp, names_prefix = "rep_") %>% arrange(tran_id)
obs_wind <- obs_covs %>% select(tran_id, col_point, wind) %>% pivot_wider(names_from = col_point, values_from = wind, names_prefix = "rep_") %>% arrange(tran_id)
obs_vpd <- obs_covs %>% select(tran_id, col_point, vpd_obs) %>% pivot_wider(names_from = col_point, values_from = vpd_obs, names_prefix = "rep_") %>% arrange(tran_id)
obs_rock <- obs_covs %>% select(tran_id, col_point, rock_obs)%>% pivot_wider(names_from = col_point, values_from = rock_obs, names_prefix = "rep_") %>% arrange(tran_id)

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
X.p <- list(
  time = time_mat,
  vpd = vpd_mat,
  rock = rock_mat,
  t_air = air_mat,
  t_gr = gr_mat,
  wind = wind_mat)

# ---- Set Coordinates (UTM meters) ----
coords <- as.matrix(site_covs[, c("tran_utm_x", "tran_utm_y")])
rownames(coords) <- site_covs$tran_id

# ---- Fit Spatial Occupancy (NNGP GP) ----
# setting priors
# sigma.sq (spatial variance parameter)
# phi (spatial range parameter)
# w (the latent spatial random effects at each site)
# and nu (spatial smoothness parameter): only specified if adopting a Matern covariance function


priors <- list(
  beta.normal  = list(mean = 0,  var = 10),   # <-- scalar applies to all βs (incl. intercept)
  alpha.normal = list(mean = 0,  var = 10),   # <-- scalar applies to all αs (incl. intercept)
  sigma.sq.ig  = c(2, 2),
  phi.unif     = c(3/20000, 3/500))

set.seed(123)
fit <- spPGOcc(
  occ.formula = ~ rockiness + shrubs + elev + prey,                  # ψ (intercept included by default)
  det.formula = ~ time + vpd + rock + t_air + t_gr + wind,           # p
  data        = list(
    y        = y,
    occ.covs = occ_df,   # <-- key fix
    det.covs = X.p,      # <-- key fix
    coords   = coords
  ),
  priors      = priors,
  cov.model   = "exponential",    #spOccupancy supports four spatial covariance models (exponential, spherical, gaussian, and matern)
  NNGP        = TRUE, n.neighbors = 5,
  n.batch     = 2000, batch.length = 25,   # total iters = 50k
  n.burn      = 10000, n.thin = 10, n.report = 50,
  verbose     = TRUE
)

summary(fit)

