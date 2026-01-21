library(dplyr)
library(tidyr)
library(spOccupancy)
library(coda)
library(stars)
library(MCMCvis)

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
# occupancy formula
global.occ.formula <- ~ rockiness + shrubs + elev + prey
# detection formula
global.det.formula <- ~ rock + I(t_air^2) + I(vpd^2) + rock + wind + I(t_gr^2)
# model type
cov.model <- "exponential"

# priors/initial values
# sigma.sq (spatial variance parameter)
# phi (spatial range parameter) the spatial range parameter phi is the most 
# sensitive to initial values. In general, the spatial range parameter will often have 
# poor mixing and take longer to converge than the rest of the parameters in the model
# as an initial value for the spatial range parameter phi, we compute the mean distance 
# between points in HBEF and then set it equal to 3 divided by this mean distance
# w (the latent spatial random effects at each site)
# and nu (spatial smoothness parameter): only specified if adopting a Matern covariance function

# pairwise distances between all transects (meters)
D <- as.matrix(dist(coords, method = "euclidean"))
off_diag <- D[upper.tri(D)]
min.dist <- min(off_diag, na.rm = TRUE)
max.dist <- max(off_diag, na.rm = TRUE)
mean.dist <- mean(off_diag, na.rm = TRUE)


# initial values
gm.inits <- list(
  alpha = 0,                           # occurrence coefficient (0 is the default; but setting explicitly for practice)
  beta = 0,                            # detection coefficient (0 is the default; but setting explicitly for practice)
  z = apply(y, 1, max, na.rm = TRUE),  # latent ("true") occurrence 
  sigma.sq  = 2,                       # spatial variance (relatively small bc all in PC)
  phi     = 3/mean.dist,         # spatial range (most sensitive) ;effective range is the average distance between transects
  w = rep(0, nrow(y)))                 # spatial random effects (0 is the default; but setting explicitly for practice)

# MCMC sampler
# we break up the total number of MCMC samples into a set of “batches”, where each batch has a specific number of MCMC samples
# the total number of MCMC samples is n.batch * batch.length. Typically, we set batch.length = 25 and then play around with n.batch until convergence of all model parameters is reached. 
# We recommend setting batch.length = 25 unless you have a specific reason to change it
batch.length <- 25
n.batch <- 400
n.burn <- 2000
n.thin <- 20
n.chains <- 3

# default accept.rate = 0.43
# The accept.rate argument specifies the ideal proportion of times we will accept the newly proposed values for these parameters. 
# Roberts and Rosenthal (2009) show that if we accept new values around 43% of the time, this will lead to optimal mixing and convergence of the MCMC chains.
# The values specified in the tuning argument helps control the initial values we will propose for both phi and nu
# The initial tuning value can be any value greater than 0, but we generally recommend starting the value out around 0.5. 
# This initial tuning value is closely related to the ideal (or target) acceptance rate we specified in accept.rate.
gm.tuning <- list(phi = 1) # copying example for now

# Priors
# We assume an inverse gamma prior for the spatial variance parameter sigma.sq (the tag of which is sigma.sq.ig), and uniform priors for the spatial decay parameter phi 
# and smoothness parameter nu (if using the Matern correlation function), with the associated tags phi.unif and nu.unif.
# Note that the priors for the spatial parameters in a spatially-explicit model must be at least weakly informative for the model to converge
gm.priors <- list(
  beta.normal  = list(mean = 0,  var = 2.72),   # default is 0 and 2.72, which corresponds to a relatively flat prior on the probability scale; 
  alpha.normal = list(mean = 0,  var = 2.72),   # such priors are an adequate choice when the goal is to specify relatively non-informative priors.
  sigma.sq.ig  = c(2, 1),
  # For the inverse-Gamma prior on the spatial variance, we typically set the shape parameter to 2 and the scale parameter equal to our best guess of the spatial variance.
  # The default prior hyperparameter values for the spatial variance σ2
  # are a shape parameter of 2 and a scale parameter of 1. This weakly informative prior suggests a prior mean of 1 for the spatial variance, which is a moderately 
  # small amount of spatial variation.
  phi.unif     = c(3/max.dist, 3/min.dist))
# For the spatial decay parameter, our default approach is to set the lower and upper bounds of the uniform prior based on the minimum and maximum distances between 
# sites in the data. More specifically, by default we set the lower bound to 3 / max and the upper bound to 3 / min, where min and max are the minimum and maximum 
# distances between sites in the data set, respectively. This equates to a vague prior that states the spatial autocorrelation in the data could only exist between sites
# that are very close together, or could span across the entire observed study area.

# The argument n.omp.threads specifies the number of threads to use for within-chain parallelization, while verbose specifies whether or not to print the progress of the 
# sampler. We highly recommend setting verbose = TRUE for all spatial models to ensure the adaptive MCMC is working as you want (and this is the reason for why this is the 
# default for this argument). The argument n.report specifies the interval to report the Metropolis-Hastings sampler acceptance rate. Note that n.report is specified in terms
# of batches, not the overall number of samples. Below we set n.report = 100, which will result in information on the acceptance rate and tuning parameters every 100th batch (not sample).

n.omp.threads <- 1
verbose <- TRUE
n.report <- 100 # Report progress at every 100th batch.

# The parameters NNGP, n.neighbors, and search.type relate to whether or not you want to fit the model with a Gaussian Process or with NNGP
out.sp <- spPGOcc(occ.formula = global.occ.formula, 
                  det.formula = global.det.formula, 
                  data = list(
                    y = y,
                    occ.covs = occ_df,   
                    det.covs = X.p,     
                    coords   = coords),
                  inits = gm.inits, 
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  priors = gm.priors, 
                  cov.model = cov.model, 
                  NNGP = TRUE, 
                  n.neighbors = 15,
                  tuning = gm.tuning, 
                  n.report = n.report, 
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = n.chains)

class(out.sp)
names(out.sp)
summary(out.sp)

# Rhat: R-hat compares the variance within each MCMC chain to the variance between different chains. The goal is for all chains to have converged to the same posterior distribution,
# meaning the within-chain and between-chain variances should be equal.
# Near 1: An R-hat value close to 1 (often with a threshold of 1.05 or 1.1) indicates that the chains have likely converged and 
# are reliably sampling from the target distribution.
# Greater than 1: Values significantly greater than 1 suggest that convergence has not been reached. The MCMC chains are likely
# not sampling from the same distribution, and you should run the simulation for more iterations.

ppc.sp.out <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2) # this time grouping by replicate, or survey occasion, (group = 2) instead of by site (group = 1).
summary(ppc.sp.out)

# binning the data across sites (group = 1) may help reveal whether the model fails to adequately represent variation in occurrence and detection probability across space
# binning the data across replicates (group = 2) may help reveal whether the model fails to adequately represent variation in detection probability across the different replicate surveys. 


waicOcc(out.sp)

# traceplot of occurrence coefficients
plot(out.sp, param = 'beta', density = FALSE)
# traceplot of detection coefficients
plot(out.sp, param = 'alpha', density = FALSE)

MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95),
         main = 'Occupancy Parameters')
MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95), 
         main = 'Detection Parameters')

