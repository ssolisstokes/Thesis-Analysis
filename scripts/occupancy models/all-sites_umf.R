# All Sites - transects as replicates unmarked occ. model
# Nov 2025

# loading in the libraries
library(dplyr)
library(tidyr)
library(spOccupancy)
library(coda)
library(stars)
library(MCMCvis)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(unmarked)
library(lubridate)
library(tibble)
library(sf)
library(MuMIn)
library(AICcmodavg)

# ---- Data Wrangling ----
# bring in master data
gm <- read.csv("data/data_gmocc_all-ut.csv", stringsAsFactors = FALSE)

# remove off transect detections 
gm$tran_id <- trimws(gm$tran_id) # fixing any spacing issues in tran_id
gm <- gm[!grepl("OTD", gm$tran_id, ignore.case = TRUE), ]

# bring in the transect length data
t.length <- read.csv("data/data_gmocc_tran-distance.csv", stringsAsFactors = FALSE)

# remove off transect detections
t.length$tran_id <- trimws(t.length$tran_id)
t.length <- t.length[!grepl("OTD", t.length$tran_id, ignore.case = TRUE), ]


# calculating time since midnight to use time as a variable
time <- strsplit(ifelse(is.na(gm$time), "NA:NA:NA", gm$time), ":", fixed = TRUE)
hh <- as.numeric(sapply(time, `[`, 1))
mm <- as.numeric(sapply(time, `[`, 2))
ss <- as.numeric(sapply(time, `[`, 3))
gm$time_sec <- hh*3600 + mm*60 + ss
gm$time_sec[!is.finite(gm$time_sec)] <- NA_real_


# just replacing the detection covariates
replacement_cols <- c("utm_x", "utm_y", "time", "detection","avg_wind_spd", "air_temp", "rel_hum", "gr_temp", "rock_ind", "shr_nw", "shr_ne", "shr_se", "shr_sw", "shr_pt_density", "vpd", "time_sec")

gm <- gm %>%
  group_by(site_id, tran_id) %>%
  dplyr::group_modify(~ {
    df <- .x
    # indices for candidate rows (col_point 1–5) and col_point 6 rows
    idx_1_5 <- which(df$col_point %in% 1:5 & !is.na(df$time_sec))
    idx_6   <- which(df$col_point == 6   & !is.na(df$time_sec))
    # if no candidates or no col_point 6, nothing to do
    if (length(idx_1_5) == 0 || length(idx_6) == 0) return(df)
    # for each col_point 6 row, find nearest-in-time col_point 1–5 and copy covariates
    for (i in idx_6) {
      diffs   <- abs(df$time_sec[idx_1_5] - df$time_sec[i])
      nearest <- idx_1_5[which.min(diffs)]
      df[nearest, replacement_cols] <- df[i, replacement_cols]
    }
    df
  }) %>%
  ungroup()

# keep only on-transect collection points 1–5 for covariates/effort summaries
transects <- subset(gm, col_point %in% 1:5)

transects %>% filter(detection == 1)
# add the transect length (path_km)
transects <- merge(transects, t.length[, c("site_id","tran_id","path_km")],
                   by = c("site_id","tran_id"), all.x = TRUE)

# ---- Transect-level summaries (replicates = transects) ----

# transect-level observation covariates from the 1–5 col. points
tran_covs <- transects %>%
  group_by(site_id, tran_id) %>%
  summarise(
    detection = sum(detection, na.rm = TRUE),
    time_sec = mean(time_sec, na.rm = TRUE),
    air_temp = mean(air_temp, na.rm = TRUE),
    gr_temp = mean(gr_temp,  na.rm = TRUE),
    wind = mean(avg_wind_spd, na.rm = TRUE),
    vpd = mean(vpd, na.rm = TRUE),
    rockiness = mean(rock_ind, na.rm = TRUE),
    shrub = mean(shr_pt_density, na.rm = TRUE),
    elev = mean(elev, na.rm = TRUE),
    prey = sum(pp_mam,pp_rep, pp_gbird, na.rm = TRUE),
    path_km = mean(path_km, na.rm = TRUE),
    .groups = "drop")


# ---- Transects to replicates within each site ----
tran_summary <- tran_covs %>%
  group_by(site_id) %>%
  arrange(tran_id, .by_group = TRUE) %>%
  mutate(rep_idx = row_number()) %>%
  ungroup()

all_sites <- sort(unique(tran_summary$site_id))

# max transects at a site
J_max <- tran_summary %>%
  count(site_id, name = "n_tran") %>%
  summarise(max(n_tran)) %>% pull()

# ---- Detection matrix ----
# (rows = sites, cols = transect replicates)

# long df
y_long <- expand.grid(site_id = all_sites, rep_idx = seq_len(J_max),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(tran_summary %>% select(site_id, rep_idx, detection),
            by = c("site_id","rep_idx")) %>%
  arrange(site_id, rep_idx)

# wide df
y_wide <- y_long %>%
  pivot_wider(names_from = rep_idx, values_from = detection, names_prefix = "rep_") %>%
  arrange(site_id)

y <- as.matrix(y_wide %>% select(starts_with("rep_")))
rownames(y) <- y_wide$site_id

# ---- Site Covariates (one row per site) ----
# transect-level summaries to site-level for occupancy covariates
site_covs <- tran_summary %>%
  group_by(site_id) %>%
  summarise(
    rockiness_mean = mean(rockiness, na.rm = TRUE),
    shrub_density_mean = mean(shrub, na.rm = TRUE),
    elevation_mean = mean(elev, na.rm = TRUE),
    prey_total = sum(prey, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(site_id)

# scale site covariates
site_covs$rockiness_z <- as.numeric(scale(site_covs$rockiness_mean))
site_covs$shrubs_z <- as.numeric(scale(site_covs$shrub_density_mean))
site_covs$elev_z <- as.numeric(scale(site_covs$elevation_mean))
site_covs$prey_z <- as.numeric(scale(site_covs$prey_total))

# align rows to y (sites as rownames)
site_covs <- site_covs %>% slice(match(rownames(y), site_id))

# ---- Transect level detection covariates ----
obs_base <- expand.grid(site_id = all_sites, rep_idx = seq_len(J_max),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(tran_summary %>%
              select(site_id, rep_idx, time_sec, air_temp, gr_temp, wind, vpd, rockiness) ,
            by = c("site_id","rep_idx")) %>%
  arrange(site_id, rep_idx)

obs_time <- obs_base %>% select(site_id, rep_idx, time_sec) %>%
  pivot_wider(names_from = rep_idx, values_from = time_sec, names_prefix = "rep_") %>% arrange(site_id)
obs_air  <- obs_base %>% select(site_id, rep_idx, air_temp) %>%
  pivot_wider(names_from = rep_idx, values_from = air_temp, names_prefix = "rep_") %>% arrange(site_id)
obs_gr   <- obs_base %>% select(site_id, rep_idx, gr_temp) %>%
  pivot_wider(names_from = rep_idx, values_from = gr_temp, names_prefix = "rep_") %>% arrange(site_id)
obs_wind <- obs_base %>% select(site_id, rep_idx, wind) %>%
  pivot_wider(names_from = rep_idx, values_from = wind, names_prefix = "rep_") %>% arrange(site_id)
obs_vpd  <- obs_base %>% select(site_id, rep_idx, vpd) %>%
  pivot_wider(names_from = rep_idx, values_from = vpd, names_prefix = "rep_") %>% arrange(site_id)
obs_rock <- obs_base %>% select(site_id, rep_idx, rockiness) %>%
  pivot_wider(names_from = rep_idx, values_from = rockiness, names_prefix = "rep_") %>% arrange(site_id)

# matrices aligned to y matrix
time_mat <- as.matrix(obs_time %>% select(starts_with("rep_")))[match(rownames(y), obs_time$site_id), ]
air_mat  <- as.matrix(obs_air  %>% select(starts_with("rep_")))[match(rownames(y), obs_air$site_id), ]
gr_mat   <- as.matrix(obs_gr   %>% select(starts_with("rep_")))[match(rownames(y), obs_gr$site_id), ]
wind_mat <- as.matrix(obs_wind %>% select(starts_with("rep_")))[match(rownames(y), obs_wind$site_id), ]
vpd_mat  <- as.matrix(obs_vpd  %>% select(starts_with("rep_")))[match(rownames(y), obs_vpd$site_id), ]
rock_mat <- as.matrix(obs_rock %>% select(starts_with("rep_")))[match(rownames(y), obs_rock$site_id), ]

# z-score scale each detection matrix
time_mat[] <- (time_mat - mean(time_mat, na.rm = TRUE)) / sd(time_mat, na.rm = TRUE)
air_mat[]  <- (air_mat  - mean(air_mat,  na.rm = TRUE)) / sd(air_mat,  na.rm = TRUE)
gr_mat[]   <- (gr_mat   - mean(gr_mat,   na.rm = TRUE)) / sd(gr_mat,   na.rm = TRUE)
wind_mat[] <- (wind_mat - mean(wind_mat, na.rm = TRUE)) / sd(wind_mat, na.rm = TRUE)
vpd_mat[]  <- (vpd_mat  - mean(vpd_mat,  na.rm = TRUE)) / sd(vpd_mat,  na.rm = TRUE)
rock_mat[] <- (rock_mat - mean(rock_mat, na.rm = TRUE)) / sd(rock_mat, na.rm = TRUE)

# ---- Prep Covariates for unmarked frame ----
# occupancy df
occ_df <- data.frame(
  rockiness = site_covs$rockiness_z,
  shrubs    = site_covs$shrubs_z,
  elev      = site_covs$elev_z,
  prey      = site_covs$prey_z)
rownames(occ_df) <- site_covs$site_id

# detection df
det_df <- list(
  time = time_mat,
  vpd  = vpd_mat,
  rock = rock_mat,
  t_air = air_mat,
  t_gr  = gr_mat,
  wind  = wind_mat)

# add to unmarked frame
umf <- unmarkedFrameOccu(y = y, siteCovs = occ_df, obsCovs = det_df)
summary(umf)

# calculate naive occupancy
naive_occ <- sum(ifelse(rowSums(y, na.rm = TRUE)>0,1,0)/nrow(y))
naive_occ
# naive occupancy is 0.83

# run the null model
null.model <- occu(formula = 
                     ~ 1 # detection
                   ~1, # occupancy
                   data = umf)

summary(null.model)
# null model occupancy = 1.81
# null model detection = -3.18

# transform probabilties from the log scale
p <-backTransform(null.model, type = "det")
p
# the probability of detecting a Gila is ~ 0.04 (4%) (0.0398) (SE = 0.0101) at a site

psi<- backTransform(null.model, type = "state")
psi
# there is roughly a 0.86 (86%) chance that a Gila occupies a site 

# write the global model
global.model.um <- occu(formula = 
                          ~vpd + rock + I(t_air^2) + t_air + I(t_gr^2) + t_gr + wind ~ 
                          ~ rockiness + shrubs + elev + prey,
                        data = umf)

summary(global.model.um)

# dredge & analyze models
mod.sel <- dredge(global.model.um)
mod.sel[1:10]    

# write model with occupancy constant (1)
occu.mb <- occu(formula = ~1 ~ rockiness + elev + prey, data = umf)
summary(occu.mb)

occu.psiPr <- predict(occu.mb, type = "state", appendData = TRUE)
head(occu.psiPr)
occu.psi <- occu(formula = ~rockiness + wind  ~rockiness + elev + prey, data = umf)
summary(occu.psi)

occu.psiPr <- predict(occu.psi,type = "state", appendData = TRUE)
head(occu.psiPr) 

# write model equations
#model_01 <- occu(~ rockiness ~ elev + shrubs, data = umf)
#model_02 <- occu(~ 1 ~ elev, data = umf)
#model_03 <- occu(~ rockiness + I(t_air^2) ~ prey, data = umf)


#print(summary(model_01))
#print(summary(model_02))
#print(summary(model_03))


m1 <- occu(~1 ~1, umf) # null model
summary(m1)

m2 <- occu(~rock ~1, umf) # effect of rockiness on detection when occupancy is constant
summary(m2)
confint(m2, type = "det")
# The confidence interval overlaps 0, meaning rockiness doesn’t strongly affect detection probability
intercept <- -3.367
slope <- 0.209
p0 <- plogis(intercept)           # rockiness = 0
p0
p1 <- plogis(intercept + slope)   # rockiness = 1
p1

  # Intercept (−3.391): baseline detection probability when rockiness = 0 (~3.3%)
  # Slope: rockiness (0.209): as rockiness increases (as it becomes sandier), detection probability slightly increases (up to ~4%) increase of ~0.8%
    # but this effect is not statistically significant (p = 0.21).
  # detection probability drops from ~3.3% to ~2.3% when rockiness increases by 1 unit.

m3 <- occu(~shrubs ~1, umf) # effect of shrub density on detection when occupancy is constant
summary(m3)
confint(m3, type = "det")
# The confidence interval overlaps 0, meaning shrub density doesn’t strongly affect detection probability
intercept <- -3.430
slope <- 0.496
p0 <- plogis(intercept)           # shrub density = 0
p0
p1 <- plogis(intercept + slope)   # shrub density = 1
p1
  # Intercept: -3.442: baseline dp when shrub density = 0 (~3.1%)
  # Slope: shrub density (0.504) dp increases as shrub density increases (~5%); (increase of ~1.9%) but the
    # effect is not statistically significant 

m4 <- occu(~prey ~1, umf) # effect of prey presence on detection when occupancy is constant
summary(m4)
confint(m4, type = "det")
# a statistically significant relationship: when more prey are present, detection probability increases.
intercept <- -3.219
slope <- 0.601
p0 <- plogis(intercept)           # prey presence = 0
p0
p1 <- plogis(intercept + slope)   # prey presence = 1
p1
  # Intercept: - 3.219: baseline dp at mean prey presence (~3.8%)
  # Slope: prey presence (0.601) dp increases as prey presence increases (~6.7%)

m5 <- occu(~elev ~1, umf) # effect of elevation on detection when occupancy is constant
summary(m5)
confint(m5, type = "det")
# The confidence interval overlaps 0, meaning elevation doesn’t strongly affect detection probability
intercept <- -3.20
slope <- 0.22
p0 <- plogis(intercept)           # elevation = 0
p0
p1 <- plogis(intercept + slope)   # elevation = 1
p1
  # Intercept: -3.20: baseline dp at mean elev = 0 (kind of) = (~3.9%)
  # Slope: elevation (0.22) dp increases with increasing elevation (?) slightly (~4.8%)

m6 <- occu(~vpd ~1, umf)
summary(m6)
# the confidence interval does not overlap 0
confint(m6, type = "det")
intercept <- -3.540 
slope <- 0.632
p0 <- plogis(intercept)           # elevation = 0
p0
p1 <- plogis(intercept + slope)   # elevation = 1
p1
  # Intercept: -3.54: baseline dp when elev = 0 (kind of) = (~3.9%)
  # Slope: elevation (0.22) dp increases with increasing elevation (?) slightly (~4.8%) (~0.9% increase)

m7 <- occu(~t_air + I(t_air^2) ~1, umf)
summary(m7)
confint(m7, type = "det")
  # baseline dp: ~4% at mean temperature
  # t_air (linear term) (positive slope = 0.0067); detection slightlyincreases with air temperature
  # t_air (quadratic term) (negative slope = -0.800); detection drops off again at high temperatures
  # (bell curve support?)
  # detection is likely highest at moderate temperatures; (low in cool mornings & hot afternoons)

m8 <- occu(~t_gr + I(t_gr^2) ~1, umf) 
summary(m8)
confint(m8, type = "det")
  # baseline dp: ~3% at mean temperature
  # t_gr (linear term) (positive slope = 0.248): weak positive slope; detection tends to increase slightly as ground temperature increases
  # t_gr (quadratice term) (negative slope = −0.256): slight negative; detection may decline again at higher ground temps slightly
  # Both effects confidence intervals include 0;
  # there’s no statistical evidence for a ground-temperature effect on detection

m9 <- occu(~wind ~1, umf)
summary(m9)
confint(m9, type = "det")
# does not overlap 0
intercept <- -3.68
slope <-  -1.47
p0 <- plogis(intercept)           # wind = 0
p0
p1 <- plogis(intercept + slope)   # wind = 1
p1
  # baseline dp: ~2.5% at mean wind
  # slope: wind(negative slope: -1.47) detection decreases to 0.05% as wind increases; windier days make detection harder



#### stop #####

m10 <- occu(~atemp.s + wind.s ~1, modocc)
summary(m9)
confint(m9, type = "det")

m11 <- occu(~atemp.s + I(atemp.s^2) + wind.s ~1, global.model.um)
summary(m10)
confint(m10, type = "det")

m12 <- occu(~interval.s ~1, global.model.um)
summary(m11)
confint(m11, type = "det")

m13 <- occu(~interval.s + I(interval.s^2) ~1, global.model.um)
summary(m12)
confint(m12, type = "det")

m14 <- occu(~interval.s + wind.s ~1, global.model.um)
summary(m13)
confint(m13, type = "det")

m15 <- occu(~interval.s + + I(interval.s^2) + wind.s ~1, global.model.um)
summary(m14)
confint(m14, type = "det")

##AIC selection for detection models
detmodels <- fitList('psi(.)p(.)' = m1,
                     'psi(.)p(tveg)' = m2,
                     'psi(.)p(trees)' = m3,
                     'psi(.)p(atemp)'= m4,
                     'psi(.)p(wind)'= m5,
                     'psi(.)p(stemp)' = m6,
                     'psi(.)p(stemp+wind)' = m6.a,
                     'psi(.)p(atemp+atemp^2)' = m7,
                     'psi(.)p(stemp+stemp^2)' = m8,
                     'psi(.)p(stemp+stemp^2+wind)' = m8.a,
                     'psi(.)p(atemp+wind)' = m9,
                     'psi(.)p(atemp+atemp^2+wind)' = m10,
                     'psi(.)p(time)' = m11,
                     'psi(.)p(time+time^2)' = m12,
                     'psi(.)p(time+wind)' = m13,
                     'psi(.)p(time+time^2+wind)' = m14
)




