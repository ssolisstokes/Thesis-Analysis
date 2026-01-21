library(dplyr)
library(tidyr)
library(spOccupancy)
library(coda)
library(stars)
library(MCMCvis)
library(ggplot2)
library(tibble)

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
# aggregate transect-level summaries to site-level for occupancy covariates
site_covs <- tran_summary %>%
  group_by(site_id) %>%
  summarise(
    rockiness_mean      = mean(rockiness, na.rm = TRUE),
    shrub_density_mean  = mean(shrub, na.rm = TRUE),
    elevation_mean      = mean(elev, na.rm = TRUE),
    prey_total          = sum(prey, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(site_id)

# scale site covariates
site_covs$rockiness_z <- as.numeric(scale(site_covs$rockiness_mean))
site_covs$shrubs_z    <- as.numeric(scale(site_covs$shrub_density_mean))
site_covs$elev_z      <- as.numeric(scale(site_covs$elevation_mean))
site_covs$prey_z      <- as.numeric(scale(site_covs$prey_total))

# align rows to y (sites as rownames)
site_covs <- site_covs %>% slice(match(rownames(y), site_id))

# ---- Detection covariate matrices (N x J_max), aggregated at transect level ----
# start from tran_summary (transect-level covs), pad to full N x J_max, widen per covariate
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

# matrices aligned to y row order
time_mat <- as.matrix(obs_time %>% select(starts_with("rep_")))[match(rownames(y), obs_time$site_id), ]
air_mat  <- as.matrix(obs_air  %>% select(starts_with("rep_")))[match(rownames(y), obs_air$site_id), ]
gr_mat   <- as.matrix(obs_gr   %>% select(starts_with("rep_")))[match(rownames(y), obs_gr$site_id), ]
wind_mat <- as.matrix(obs_wind %>% select(starts_with("rep_")))[match(rownames(y), obs_wind$site_id), ]
vpd_mat  <- as.matrix(obs_vpd  %>% select(starts_with("rep_")))[match(rownames(y), obs_vpd$site_id), ]
rock_mat <- as.matrix(obs_rock %>% select(starts_with("rep_")))[match(rownames(y), obs_rock$site_id), ]

# z-score scale each detection matrix (across all cells)
time_mat[] <- (time_mat - mean(time_mat, na.rm = TRUE)) / sd(time_mat, na.rm = TRUE)
air_mat[]  <- (air_mat  - mean(air_mat,  na.rm = TRUE)) / sd(air_mat,  na.rm = TRUE)
gr_mat[]   <- (gr_mat   - mean(gr_mat,   na.rm = TRUE)) / sd(gr_mat,   na.rm = TRUE)
wind_mat[] <- (wind_mat - mean(wind_mat, na.rm = TRUE)) / sd(wind_mat, na.rm = TRUE)
vpd_mat[]  <- (vpd_mat  - mean(vpd_mat,  na.rm = TRUE)) / sd(vpd_mat,  na.rm = TRUE)
rock_mat[] <- (rock_mat - mean(rock_mat, na.rm = TRUE)) / sd(rock_mat, na.rm = TRUE)

# ---- Prep Covariates for spPGOcc ----
occ_df <- data.frame(
  rockiness = site_covs$rockiness_z,
  shrubs    = site_covs$shrubs_z,
  elev      = site_covs$elev_z,
  prey      = site_covs$prey_z)
rownames(occ_df) <- site_covs$site_id

cor(occ_df, use = "pairwise.complete.obs")
round(cor(occ_df), 2)
pairs(occ_df)
library(car)

vif(lm(rep(1, nrow(occ_df)) ~ ., data = occ_df))
# all VIFs <3, no collinearity issues with occupancy

X.p <- list(
  time = time_mat,
  vpd  = vpd_mat,
  rock = rock_mat,
  t_air = air_mat,
  t_gr  = gr_mat,
  wind  = wind_mat)


det_cov_df <- data.frame(
  rock = as.vector(X.p$rock),
  t_air = as.vector(X.p$t_air),
  t_gr  = as.vector(X.p$t_gr),
  vpd   = as.vector(X.p$vpd),
  wind  = as.vector(X.p$wind)
)


det_cov_df <- det_cov_df[complete.cases(det_cov_df), ]
round(cor(det_cov_df), 2)
pairs(det_cov_df)
vif(lm(rep(1, nrow(det_cov_df)) ~ ., data = det_cov_df))


# strong collinearity between temp/weather variables (expected)
# t_air, t_gr, and vpd are not statistically separable when together in a model
# WAIC may be high, but VIF shows caution in interp
# ---- Global model formulas, inits, priors, MCMC config ----
global.occ.formula <- ~ rockiness + shrubs + elev + prey
global.det.formula <- ~ rock + I(t_air^2) + t_air + vpd + wind + I(t_gr^2) + t_gr

gm.inits <- list(
  alpha = 0, 
  beta  = 0, 
  z = apply(y, 1, max, na.rm = TRUE))

gm.priors <- list(
  alpha.normal = list(mean = 0, var = 2.72), 
  beta.normal  = list(mean = 0, var = 2.72))

n.samples <- 5000
n.burn    <- 3000
n.thin    <- 2
n.chains  <- 3

out <- PGOcc(
  occ.formula = global.occ.formula,
  det.formula = global.det.formula, 
  data = list(
    y       = y,
    occ.covs = occ_df,
    det.covs = X.p),
  inits = gm.inits, 
  n.samples = n.samples, 
  priors = gm.priors, 
  n.omp.threads = 1, 
  verbose = TRUE, 
  n.report = 1000, 
  n.burn = n.burn, 
  n.thin = n.thin, 
  n.chains = n.chains)

names(out)
summary(out)
waicOcc(out)
# ---- Posterior Predictive Checks - Global model ----
plogis(-3.4642) # = 0.0303 

plot(out, 'beta', density = FALSE) # Occupancy parameters.
plot(out, 'alpha', density = FALSE) # Detection parameters.

ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
ppc.out2 <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.out)
summary(ppc.out2)

ppc.df <- data.frame(fit = ppc.out$fit.y, 
                     fit.rep = ppc.out$fit.y.rep, 
                     color = 'lightskyblue1')

ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'

plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

diff.fit <- ppc.out$fit.y.rep.group.quants[3, ] - ppc.out$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')






# model selection
n.samples <- 10000
n.burn    <- 5000
n.thin    <- 5
n.chains  <- 3

# occupancy formulas + global detection
om1 <- PGOcc(occ.formula = ~ 1,
                   det.formula = global.det.formula, 
                   data = list(
                     y = y,
                     occ.covs = occ_df,   
                     det.covs = X.p),    
                   inits = gm.inits, 
                   n.samples = n.samples, 
                   priors = gm.priors, 
                   n.omp.threads = 1, 
                   verbose = FALSE, 
                   n.report = 2000, 
                   n.burn = n.burn, 
                   n.thin = n.thin, 
                   n.chains = n.chains)

# habitat + elev
om2 <- PGOcc(occ.formula = ~ rockiness + shrubs + elev, 
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000,  
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# habitat + prey
om3 <- PGOcc(occ.formula = ~ rockiness + shrubs + prey, 
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# habitat
om4 <- PGOcc(occ.formula = ~ rockiness + shrubs, 
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rockiness + prey
om5 <- PGOcc(occ.formula = ~ rockiness + prey,
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# vegetation + prey
om6 <- PGOcc(occ.formula = ~ shrubs + prey,
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rockiness
om7 <- PGOcc(occ.formula = ~ rockiness,
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# shrub density
om8 <- PGOcc(occ.formula = ~ shrubs,
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# om_list <- list(om1, om2, om3, om4, om5, om6, om7,om8)

om1.waic <- waicOcc(om1)
om2.waic <- waicOcc(om2)
om3.waic <- waicOcc(om3)
om4.waic <- waicOcc(om4)
om5.waic <- waicOcc(om5)
om6.waic <- waicOcc(om6)
om7.waic <- waicOcc(om7)
om8.waic <- waicOcc(om8)

om1.df <- data.frame(model = "om1", t(om1.waic))
om2.df <- data.frame(model = "om2", t(om2.waic))
om3.df <- data.frame(model = "om3", t(om3.waic))
om4.df <- data.frame(model = "om4", t(om4.waic))
om5.df <- data.frame(model = "om5", t(om5.waic))
om6.df <- data.frame(model = "om6", t(om6.waic))
om7.df <- data.frame(model = "om7", t(om7.waic))
om8.df <- data.frame(model = "om8", t(om8.waic))

om.waic.df <- bind_rows(list(om1.df, om2.df, om3.df, om4.df, om5.df, om6.df, om7.df, om8.df))

om.waic_weights.tbl <- 
om.waic.df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta))) %>%# Burnham & Andserson 2002, 2004 converting model weights formula
    arrange(WAIC)


om.waic_weights.tbl
# MuMIn is built around likelihood-based models (lm/glm/glmer, unmarked::occu, etc.) and uses AIC/logLik/coef methods for those objects. Your PGOcc objects from spOccupancy are Bayesian MCMC fits with their own class ("PGOcc") and methods (waicOcc, plot, etc.), and there isn’t a ready-made MuMIn interface.

# model weights as objects
om1.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om1"]
om2.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om2"]
om3.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om3"]
om4.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om4"]
om5.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om5"]
om6.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om6"]
om7.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om7"]
om8.weight <- om.waic_weights.tbl$weight[om.waic_weights.tbl$model == "om8"]

# covariates in each model
  # rockiness: om2, om3, om4, om5, om7
  # shrub density: om2, om3, om4, om6, om8
  # elevation: om2
  # prey presence: om3, om5, om6

# calculating sum weight and variable importance
rockiness.sum_weight <- om2.weight + om3.weight + om4.weight + om5.weight + om7.weight
shrub.sum_weight <- om2.weight + om3.weight + om4.weight + om6.weight + om8.weight
elev.sum_weight <- om2.weight
prey.sum_weight <- om3.weight + om5.weight + om6.weight

var_importance <- data.frame(
  variable = c("rockiness","shrub","elev","prey"),
  sum_weight = c(rockiness.sum_weight, shrub.sum_weight, elev.sum_weight, prey.sum_weight))

var_importance
# pulling intercept/logit-scale coefficients
# om1 posterior means
om1.mat  <- as.matrix(om1$beta.samples)
om1.means <- colMeans(om1.mat)

om1_int <- om1.means["(Intercept)"]
om1_rock <- 0          
om1_shrubs <- 0         
om1_elev <- 0          
om1_prey <- 0          

# om2 posterior means
om2.mat  <- as.matrix(om2$beta.samples)
om2.means <- colMeans(om2.mat)

om2_int <- om2.means["(Intercept)"]
om2_rock <- om2.means["rockiness"]
om2_shrubs <- om2.means["shrubs"]
om2_elev <- om2.means["elev"]
om2_prey <- 0         

# om3 posterior means
om3.mat  <- as.matrix(om3$beta.samples)
om3.means <- colMeans(om3.mat)

om3_int <- om3.means["(Intercept)"]
om3_rock <- om3.means["rockiness"]
om3_shrubs <- om3.means["shrubs"]
om3_elev <- 0          
om3_prey <- om3.means["prey"]

# om4 posterior means
om4.mat  <- as.matrix(om4$beta.samples)
om4.means <- colMeans(om4.mat)

om4_int <- om4.means["(Intercept)"]
om4_rock <- om4.means["rockiness"]
om4_shrubs <- om4.means["shrubs"]
om4_elev <- 0          
om4_prey <- 0          

# om5 posterior means
om5.mat  <- as.matrix(om5$beta.samples)
om5.means <- colMeans(om5.mat)

om5_int <- om5.means["(Intercept)"]
om5_rock <- om5.means["rockiness"]
om5_shrubs <- 0         
om5_elev <- 0          
om5_prey <- om5.means["prey"]

# om6 posterior means
om6.mat  <- as.matrix(om6$beta.samples)
om6.means <- colMeans(om6.mat)

om6_int <- om6.means["(Intercept)"]
om6_rock <- 0       
om6_shrubs <- om6.means["shrubs"]
om6_elev <- 0      
om6_prey <- om6.means["prey"]

# om7 posterior means
om7.mat  <- as.matrix(om7$beta.samples)
om7.means <- colMeans(om7.mat)

om7_int <- om7.means["(Intercept)"]
om7_rock <- om7.means["rockiness"]
om7_shrubs <- 0          
om7_elev <- 0          
om7_prey <- 0          

# om8 posterior means
om8.mat  <- as.matrix(om8$beta.samples)
om8.means <- colMeans(om8.mat)

om8_int <- om8.means["(Intercept)"]
om8_rock <- 0       
om8_shrubs <- om8.means["shrubs"]
om8_elev <- 0       
om8_prey <- 0         


# multiply each coefficient by its model’s weight and sum across models
all.om_intercepts <- 
  (om1_int * om1.weight) + (om2_int * om2.weight) +
  (om3_int * om3.weight) + (om4_int * om4.weight) +
  (om5_int * om5.weight) + (om6_int * om6.weight) +
  (om7_int * om7.weight) + (om8_int * om8.weight)

all.om_rockiness <- 
  (om1_rock * om1.weight) + (om2_rock * om2.weight) +
  (om3_rock * om3.weight) + (om4_rock * om4.weight) +
  (om5_rock * om5.weight) + (om6_rock * om6.weight) +
  (om7_rock * om7.weight) + (om8_rock * om8.weight)

all.om_shrubs <- 
  (om1_shrubs * om1.weight) + (om2_shrubs * om2.weight) +
  (om3_shrubs * om3.weight) + (om4_shrubs * om4.weight) +
  (om5_shrubs * om5.weight) + (om6_shrubs * om6.weight) +
  (om7_shrubs * om7.weight) + (om8_shrubs * om8.weight)

all.om_elev <- 
  (om1_elev * om1.weight) + (om2_elev * om2.weight) +
  (om3_elev * om3.weight) + (om4_elev * om4.weight) +
  (om5_elev * om5.weight) + (om6_elev * om6.weight) +
  (om7_elev * om7.weight) + (om8_elev * om8.weight)

all.om_prey <- 
  (om1_prey * om1.weight) + (om2_prey * om2.weight) +
  (om3_prey * om3.weight) + (om4_prey * om4.weight) +
  (om5_prey * om5.weight) + (om6_prey * om6.weight) +
  (om7_prey * om7.weight) + (om8_prey * om8.weight)

om.modavg <- data.frame(
  avg_coefficients = c(all.om_intercepts, all.om_rockiness, all.om_shrubs, all.om_elev, all.om_prey))

om.modavg

# coefficient interpretation: (on the logit-scale, and each covariate is z-scored)
plogis(1.61844636)
# intercept: a site with average rockiness, shrub density, elevation and prey availability has an 
  # estimated occupancy of 0.834 (83%)
# rockiness:
# moderate effect
# sites that are mostly sandy are less likely to be occupied
# occupancy decreases as the substrate gets sandier; a 1 SD increase in rockiness reduces the log-odds of 
  # occupancy by 0.509 (half?)
  # this makes sense Gilas need rocky shelters for refugia as well as diverse habitat to
  # support prey, full sand dunes would not be appropriate GM habitat
# shrub density:
# strongest coefficient 
# sites with greater shrub density are more likely to be occupied
# a 1 SD increase in shrub density increases the log-odds of occupancy by 0.804
# elevation:
# very small/near zero
  # elevation has little to no effect on occupancy within our study sites
# prey:
# slight negative? due to correlation?
  # may be a better reflection of detection than occupancy
  # occupancy decreases as prey presence increases, slightly, but a weak effect

# k-fold cross validation of models
# k-folds = how many folds to split the data into
  # typically based on number of sites
  # for 6 sites, 3 or 2 (3 = 2 sites per fold, 2 = 3 sites per fold)
# k.fold.threads = how many folds to run in parallel/about speed not stats
  # based on how many cores your computer has
  # start with 1 and test from there
  #parallel::detectCores() to check cores; my current laptop has 14
  # can't run more threads than folds so 3 is the max for now

# best-fit occupancy formula
fit.occupancy.formula <- ~ rockiness + shrubs + prey

# detection formulas
#dm1 <- ~ 1
#dm2 <- ~ rock
#dm3 <- ~ t_air + I(t_air^2)
#dm4 <- ~ t_air + I(t_air^2) + vpd
#dm5 <- ~ t_air + I(t_air^2) + wind
#dm6 <- ~ t_air + I(t_air^2) + vpd + wind
#dm7 <- ~ t_gr + I(t_gr^2)
#dm8 <- ~ t_gr + I(t_gr^2) + vpd + wind
#dm9 <- ~ rock + t_air + I(t_air^2) + wind + vpd
#dm10 <- ~ rock  + vpd
#dm11 <- ~ rock  + t_air + I(t_air^2)
#dm12 <- ~ rock + I(t_air^2) + t_air + vpd + wind + I(t_gr^2) + t_gr, 

# null
dm1 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ 1, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rockiness
dm2 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ rock, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000,  
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# air temperature
dm3 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_air + I(t_air^2), 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# weather
dm4 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_air + I(t_air^2) + vpd, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# weather 2
dm5 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_air + I(t_air^2) + wind, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# weather 3
dm6 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_air + I(t_air^2) + vpd + wind, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# ground temp
dm7 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_gr + I(t_gr^2), 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# ground temp + weather
dm8 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ t_gr + I(t_gr^2) + vpd + wind, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rock and weather
dm9 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ rock + t_air + I(t_air^2) + wind + vpd, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rock and vpd
dm10 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ rock  + vpd, 
             data = list(
               y = y,
               occ.covs = occ_df,   
               det.covs = X.p),    
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

# rock and air
dm11 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock  + t_air + I(t_air^2), 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

# global 
dm12 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + I(t_air^2) + t_air + vpd + wind + I(t_gr^2) + t_gr, 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm1.waic <- waicOcc(dm1)
dm2.waic <- waicOcc(dm2)
dm3.waic <- waicOcc(dm3)
dm4.waic <- waicOcc(dm4)
dm5.waic <- waicOcc(dm5)
dm6.waic <- waicOcc(dm6)
dm7.waic <- waicOcc(dm7)
dm8.waic <- waicOcc(dm8)
dm9.waic <- waicOcc(dm9)
dm10.waic <- waicOcc(dm10)
dm11.waic <- waicOcc(dm11)
dm12.waic <- waicOcc(dm12)

dm1.df <- data.frame(model = "dm1", t(dm1.waic))
dm2.df <- data.frame(model = "dm2", t(dm2.waic))
dm3.df <- data.frame(model = "dm3", t(dm3.waic))
dm4.df <- data.frame(model = "dm4", t(dm4.waic))
dm5.df <- data.frame(model = "dm5", t(dm5.waic))
dm6.df <- data.frame(model = "dm6", t(dm6.waic))
dm7.df <- data.frame(model = "dm7", t(dm7.waic))
dm8.df <- data.frame(model = "dm8", t(dm8.waic))
dm9.df <- data.frame(model = "dm9", t(dm9.waic))
dm10.df <- data.frame(model = "dm10", t(dm10.waic))
dm11.df <- data.frame(model = "dm11", t(dm11.waic))
dm12.df <- data.frame(model = "dm12", t(dm12.waic))

dm.waic.df <- bind_rows(list(dm1.df, dm2.df, dm3.df, dm4.df, dm5.df, dm6.df, dm7.df, dm8.df, dm9.df, dm10.df, dm11.df, dm12.df))

dm.waic_weights.tbl <- 
  dm.waic.df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta))) %>%# Burnham & Anderson 2002, 2004 converting model weights formula
  arrange(WAIC)


 dm.waic_weights.tbl

# model averaging and variable importance weights
# model weights as objects
dm1.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm1"]
dm2.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm2"]
dm3.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm3"]
dm4.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm4"]
dm5.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm5"]
dm6.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm6"]
dm7.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm7"]
dm8.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm8"]
dm9.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm9"]
dm10.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm10"]
dm11.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm11"]
dm12.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm12"]

# detection formulas
#dm1 <- ~ 1
#dm2 <- ~ rock
#dm3 <- ~ t_air + I(t_air^2)
#dm4 <- ~ t_air + I(t_air^2) + vpd
#dm5 <- ~ t_air + I(t_air^2) + wind
#dm6 <- ~ t_air + I(t_air^2) + vpd + wind
#dm7 <- ~ t_gr + I(t_gr^2)
#dm8 <- ~ t_gr + I(t_gr^2) + vpd + wind
#dm9 <- ~ rock + t_air + I(t_air^2) + wind + vpd
#dm10 <- ~ rock  + vpd
#dm11 <- ~ rock  + t_air + I(t_air^2)

# covariates in each model
# rock: dm2, dm9, dm10, dm11, dm12
# air temp: dm3, dm4, dm5, dm6, dm9, dm11, dm12
# ground temp: dm7, dm8, dm12
# vpd: dm4, dm6, dm8, dm9, dm10, dm12
# wind: dm5, dm6, dm8, dm9, dm12

# calculating sum weight and variable importance
rock.sum_weight <- dm2.weight + dm9.weight + dm10.weight + dm11.weight + dm12.weight
t_air.sum_weight <- dm3.weight + dm4.weight + dm5.weight + dm6.weight + dm9.weight + dm11.weight + dm12.weight
t_gr.sum_weight <- dm7.weight + dm8.weight + dm12.weight
vpd.sum_weight <- dm4.weight + dm6.weight + dm8.weight + dm9.weight + dm10.weight + dm12.weight
wind.sum_weight <- dm5.weight + dm6.weight + dm8.weight + dm9.weight + dm12.weight
  
dm.var_importance <- data.frame(
  variable = c("rock","t_air","t_gr","vpd", "wind"),
  sum_weight = c(rock.sum_weight, t_air.sum_weight, t_gr.sum_weight, vpd.sum_weight,wind.sum_weight))

dm.var_importance


dm1.mat <- as.matrix(dm1$alpha.samples)
dm1.means <- colMeans(dm1.mat)

dm1_int   <- dm1.means["(Intercept)"]
dm1_rock  <- 0
dm1_tair  <- 0
dm1_tair2 <- 0
dm1_tgr   <- 0
dm1_tgr2  <- 0
dm1_vpd   <- 0
dm1_wind  <- 0


dm2.mat <- as.matrix(dm2$alpha.samples)
dm2.means <- colMeans(dm2.mat)

dm2_int   <- dm2.means["(Intercept)"]
dm2_rock  <- dm2.means["rock"]
dm2_tair  <- 0
dm2_tair2 <- 0
dm2_tgr   <- 0
dm2_tgr2  <- 0
dm2_vpd   <- 0
dm2_wind  <- 0


dm3.mat <- as.matrix(dm3$alpha.samples)
dm3.means <- colMeans(dm3.mat)

dm3_int   <- dm3.means["(Intercept)"]
dm3_tair  <- dm3.means["t_air"]
dm3_tair2 <- dm3.means["I(t_air^2)"]
dm3_rock  <- 0
dm3_tgr   <- 0
dm3_tgr2  <- 0
dm3_vpd   <- 0
dm3_wind  <- 0


dm4.mat <- as.matrix(dm4$alpha.samples)
dm4.means <- colMeans(dm4.mat)

dm4_int   <- dm4.means["(Intercept)"]
dm4_tair  <- dm4.means["t_air"]
dm4_tair2 <- dm4.means["I(t_air^2)"]
dm4_vpd   <- dm4.means["vpd"]
dm4_rock  <- 0
dm4_tgr   <- 0
dm4_tgr2  <- 0
dm4_wind  <- 0


dm5.mat <- as.matrix(dm5$alpha.samples)
dm5.means <- colMeans(dm5.mat)

dm5_int   <- dm5.means["(Intercept)"]
dm5_tair  <- dm5.means["t_air"]
dm5_tair2 <- dm5.means["I(t_air^2)"]
dm5_wind  <- dm5.means["wind"]
dm5_rock  <- 0
dm5_tgr   <- 0
dm5_tgr2  <- 0
dm5_vpd   <- 0


dm6.mat <- as.matrix(dm6$alpha.samples)
dm6.means <- colMeans(dm6.mat)

dm6_int   <- dm6.means["(Intercept)"]
dm6_tair  <- dm6.means["t_air"]
dm6_tair2 <- dm6.means["I(t_air^2)"]
dm6_vpd   <- dm6.means["vpd"]
dm6_wind  <- dm6.means["wind"]
dm6_rock  <- 0
dm6_tgr   <- 0
dm6_tgr2  <- 0


dm7.mat <- as.matrix(dm7$alpha.samples)
dm7.means <- colMeans(dm7.mat)

dm7_int   <- dm7.means["(Intercept)"]
dm7_tgr   <- dm7.means["t_gr"]
dm7_tgr2  <- dm7.means["I(t_gr^2)"]
dm7_rock  <- 0
dm7_tair  <- 0
dm7_tair2 <- 0
dm7_vpd   <- 0
dm7_wind  <- 0


dm8.mat <- as.matrix(dm8$alpha.samples)
dm8.means <- colMeans(dm8.mat)

dm8_int   <- dm8.means["(Intercept)"]
dm8_tgr   <- dm8.means["t_gr"]
dm8_tgr2  <- dm8.means["I(t_gr^2)"]
dm8_vpd   <- dm8.means["vpd"]
dm8_wind  <- dm8.means["wind"]
dm8_rock  <- 0
dm8_tair  <- 0
dm8_tair2 <- 0


dm9.mat <- as.matrix(dm9$alpha.samples)
dm9.means <- colMeans(dm9.mat)

dm9_int   <- dm9.means["(Intercept)"]
dm9_rock  <- dm9.means["rock"]
dm9_tair  <- dm9.means["t_air"]
dm9_tair2 <- dm9.means["I(t_air^2)"]
dm9_wind  <- dm9.means["wind"]
dm9_vpd   <- dm9.means["vpd"]
dm9_tgr   <- 0
dm9_tgr2  <- 0


dm10.mat <- as.matrix(dm10$alpha.samples)
dm10.means <- colMeans(dm10.mat)

dm10_int   <- dm10.means["(Intercept)"]
dm10_rock  <- dm10.means["rock"]
dm10_vpd   <- dm10.means["vpd"]
dm10_tair  <- 0
dm10_tair2 <- 0
dm10_tgr   <- 0
dm10_tgr2  <- 0
dm10_wind  <- 0


dm11.mat <- as.matrix(dm11$alpha.samples)
dm11.means <- colMeans(dm11.mat)

dm11_int   <- dm11.means["(Intercept)"]
dm11_rock  <- dm11.means["rock"]
dm11_tair  <- dm11.means["t_air"]
dm11_tair2 <- dm11.means["I(t_air^2)"]
dm11_tgr   <- 0
dm11_tgr2  <- 0
dm11_vpd   <- 0
dm11_wind  <- 0


dm12.mat <- as.matrix(dm12$alpha.samples)
dm12.means <- colMeans(dm12.mat)

dm12_int   <- dm12.means["(Intercept)"]
dm12_rock  <- dm12.means["rock"]
dm12_tair  <- dm12.means["t_air"]
dm12_tair2 <- dm12.means["I(t_air^2)"]
dm12_tgr   <- dm12.means["t_gr"]
dm12_tgr2  <- dm12.means["I(t_gr^2)"]
dm12_vpd   <- dm12.means["vpd"]
dm12_wind  <- dm12.means["wind"]


avg_int <-
  dm1_int*dm1.weight + dm2_int*dm2.weight + dm3_int*dm3.weight + dm4_int*dm4.weight +
  dm5_int*dm5.weight + dm6_int*dm6.weight + dm7_int*dm7.weight + dm8_int*dm8.weight +
  dm9_int*dm9.weight + dm10_int*dm10.weight + dm11_int*dm11.weight + dm12_int*dm12.weight

avg_rock <-
  dm1_rock*dm1.weight + dm2_rock*dm2.weight + dm3_rock*dm3.weight + dm4_rock*dm4.weight +
  dm5_rock*dm5.weight + dm6_rock*dm6.weight + dm7_rock*dm7.weight + dm8_rock*dm8.weight +
  dm9_rock*dm9.weight + dm10_rock*dm10.weight + dm11_rock*dm11.weight + dm12_rock*dm12.weight

avg_tair <-
  dm1_tair*dm1.weight + dm2_tair*dm2.weight + dm3_tair*dm3.weight + dm4_tair*dm4.weight +
  dm5_tair*dm5.weight + dm6_tair*dm6.weight + dm7_tair*dm7.weight + dm8_tair*dm8.weight +
  dm9_tair*dm9.weight + dm10_tair*dm10.weight + dm11_tair*dm11.weight + dm12_tair*dm12.weight

avg_tair2 <-
  dm1_tair2*dm1.weight + dm2_tair2*dm2.weight + dm3_tair2*dm3.weight + dm4_tair2*dm4.weight +
  dm5_tair2*dm5.weight + dm6_tair2*dm6.weight + dm7_tair2*dm7.weight + dm8_tair2*dm8.weight +
  dm9_tair2*dm9.weight + dm10_tair2*dm10.weight + dm11_tair2*dm11.weight + dm12_tair2*dm12.weight

avg_tgr <-
  dm1_tgr*dm1.weight + dm2_tgr*dm2.weight + dm3_tgr*dm3.weight + dm4_tgr*dm4.weight +
  dm5_tgr*dm5.weight + dm6_tgr*dm6.weight + dm7_tgr*dm7.weight + dm8_tgr*dm8.weight +
  dm9_tgr*dm9.weight + dm10_tgr*dm10.weight + dm11_tgr*dm11.weight + dm12_tgr*dm12.weight

avg_tgr2 <-
  dm1_tgr2*dm1.weight + dm2_tgr2*dm2.weight + dm3_tgr2*dm3.weight + dm4_tgr2*dm4.weight +
  dm5_tgr2*dm5.weight + dm6_tgr2*dm6.weight + dm7_tgr2*dm7.weight + dm8_tgr2*dm8.weight +
  dm9_tgr2*dm9.weight + dm10_tgr2*dm10.weight + dm11_tgr2*dm11.weight + dm12_tgr2*dm12.weight

avg_vpd <-
  dm1_vpd*dm1.weight + dm2_vpd*dm2.weight + dm3_vpd*dm3.weight + dm4_vpd*dm4.weight +
  dm5_vpd*dm5.weight + dm6_vpd*dm6.weight + dm7_vpd*dm7.weight + dm8_vpd*dm8.weight +
  dm9_vpd*dm9.weight + dm10_vpd*dm10.weight + dm11_vpd*dm11.weight + dm12_vpd*dm12.weight

avg_wind <-
  dm1_wind*dm1.weight + dm2_wind*dm2.weight + dm3_wind*dm3.weight + dm4_wind*dm4.weight +
  dm5_wind*dm5.weight + dm6_wind*dm6.weight + dm7_wind*dm7.weight + dm8_wind*dm8.weight +
  dm9_wind*dm9.weight + dm10_wind*dm10.weight + dm11_wind*dm11.weight + dm12_wind*dm12.weight


dm.modavg <- data.frame(
  variable = c("(Intercept)", "rock", "t_air", "t_air2", "t_gr", "t_gr2", "vpd", "wind"),
  avg_coeff = c(avg_int, avg_rock, avg_tair, avg_tair2, avg_tgr, avg_tgr2, avg_vpd, avg_wind)
)

dm.modavg

# new detection formulas after trimming variables due to collineratiry 
#dm13 <- ~ rock
#dm14 <- ~ t_gr + I(t_gr^2)
#dm15 <- ~ rock + wind
#dm16 <- ~ rock + t_gr + I(t_gr^2)
#dm17 <- ~ wind + t_gr + I(t_gr^2)
#dm18 <- ~ rock + wind + t_gr + I(t_gr^2)


dm13 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock, 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)


dm14 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ t_gr + I(t_gr^2), 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm15 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + wind, 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm16 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + t_gr + I(t_gr^2),
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm17 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ wind + t_gr + I(t_gr^2), 
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm18 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + wind + t_gr + I(t_gr^2),
              data = list(
                y = y,
                occ.covs = occ_df,   
                det.covs = X.p),    
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm13.waic <- waicOcc(dm13)
dm14.waic <- waicOcc(dm14)
dm15.waic <- waicOcc(dm15)
dm16.waic <- waicOcc(dm16)
dm17.waic <- waicOcc(dm17)
dm18.waic <- waicOcc(dm18)


dm13.df <- data.frame(model = "dm13", t(dm13.waic))
dm14.df <- data.frame(model = "dm14", t(dm14.waic))
dm15.df <- data.frame(model = "dm15", t(dm15.waic))
dm16.df <- data.frame(model = "dm16", t(dm16.waic))
dm17.df <- data.frame(model = "dm17", t(dm17.waic))
dm18.df <- data.frame(model = "dm18", t(dm18.waic))


dm.waic.df <- bind_rows(list(dm13.df,dm14.df, dm15.df, dm16.df, dm17.df, dm18.df ))

dm.waic_weights.tbl <- 
  dm.waic.df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta))) %>%# Burnham & Anderson 2002, 2004 converting model weights formula
  arrange(WAIC)


dm.waic_weights.tbl

dm13.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm13"]
dm14.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm14"]
dm15.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm15"]
dm16.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm16"]
dm17.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm17"]
dm18.weight <- dm.waic_weights.tbl$weight[dm.waic_weights.tbl$model == "dm18"]


# detection formulas
#dm13 <- ~ rock
#dm14 <- ~ t_gr + I(t_gr^2)
#dm15 <- ~ rock + wind
#dm16 <- ~ rock + t_gr + I(t_gr^2)
#dm17 <- ~ wind + t_gr + I(t_gr^2)
#dm18 <- ~ rock + wind + t_gr + I(t_gr^2)

# covariates in each model
# rock: dm13, dm15, dm16, dm18
# ground temp: dm14, dm16, dm17, dm18
# wind: dm15, dm17, dm18

# calculating sum weight and variable importance
rock.sum_weight <- dm13.weight + dm15.weight + dm16.weight + dm18.weight
t_gr.sum_weight <- dm14.weight + dm16.weight + dm17.weight + dm18.weight
wind.sum_weight <- dm15.weight + dm17.weight + dm18.weight

dm.var_importance <- data.frame(
  variable = c("rock","t_gr", "wind"),
  sum_weight = c(rock.sum_weight, t_gr.sum_weight,wind.sum_weight))

dm.var_importance


# model averaged coefficients
dm13.mat <- as.matrix(dm13$alpha.samples)
dm13.means <- colMeans(dm13.mat)

dm13_int  <- dm13.means["(Intercept)"]
dm13_rock <- dm13.means["rock"]
dm13_tgr  <- 0
dm13_tgr2 <- 0
dm13_wind <- 0

dm14.mat <- as.matrix(dm14$alpha.samples)
dm14.means <- colMeans(dm14.mat)

dm14_int  <- dm14.means["(Intercept)"]
dm14_rock <- 0
dm14_tgr  <- dm14.means["t_gr"]
dm14_tgr2 <- dm14.means["I(t_gr^2)"]
dm14_wind <- 0


dm15.mat <- as.matrix(dm15$alpha.samples)
dm15.means <- colMeans(dm15.mat)

dm15_int  <- dm15.means["(Intercept)"]
dm15_rock <- dm15.means["rock"]
dm15_tgr  <- 0
dm15_tgr2 <- 0
dm15_wind <- dm15.means["wind"]


dm16.mat <- as.matrix(dm16$alpha.samples)
dm16.means <- colMeans(dm16.mat)

dm16_int  <- dm16.means["(Intercept)"]
dm16_rock <- dm16.means["rock"]
dm16_tgr  <- dm16.means["t_gr"]
dm16_tgr2 <- dm16.means["I(t_gr^2)"]
dm16_wind <- 0


dm17.mat <- as.matrix(dm17$alpha.samples)
dm17.means <- colMeans(dm17.mat)

dm17_int  <- dm17.means["(Intercept)"]
dm17_rock <- 0
dm17_tgr  <- dm17.means["t_gr"]
dm17_tgr2 <- dm17.means["I(t_gr^2)"]
dm17_wind <- dm17.means["wind"]


dm18.mat <- as.matrix(dm18$alpha.samples)
dm18.means <- colMeans(dm18.mat)

dm18_int  <- dm18.means["(Intercept)"]
dm18_rock <- dm18.means["rock"]
dm18_tgr  <- dm18.means["t_gr"]
dm18_tgr2 <- dm18.means["I(t_gr^2)"]
dm18_wind <- dm18.means["wind"]

avg_int <-
  dm13_int*dm13.weight + dm14_int*dm14.weight + dm15_int*dm15.weight +
  dm16_int*dm16.weight + dm17_int*dm17.weight + dm18_int*dm18.weight

avg_rock <-
  dm13_rock*dm13.weight + dm14_rock*dm14.weight + dm15_rock*dm15.weight +
  dm16_rock*dm16.weight + dm17_rock*dm17.weight + dm18_rock*dm18.weight

avg_tgr <-
  dm13_tgr*dm13.weight + dm14_tgr*dm14.weight + dm15_tgr*dm15.weight +
  dm16_tgr*dm16.weight + dm17_tgr*dm17.weight + dm18_tgr*dm18.weight

avg_tgr2 <-
  dm13_tgr2*dm13.weight + dm14_tgr2*dm14.weight + dm15_tgr2*dm15.weight +
  dm16_tgr2*dm16.weight + dm17_tgr2*dm17.weight + dm18_tgr2*dm18.weight

avg_wind <-
  dm13_wind*dm13.weight + dm14_wind*dm14.weight + dm15_wind*dm15.weight +
  dm16_wind*dm16.weight + dm17_wind*dm17.weight + dm18_wind*dm18.weight



dm.modavg <- data.frame(
  variable = c("(Intercept)", "rock", "t_gr", "t_gr2", "wind"),
  avg_coeff = c(avg_int, avg_rock, avg_tgr, avg_tgr2, avg_wind)
)

dm.modavg
# interpretation
  # intercept: intercept indicates very low baseline detection probability when all covariates are at their (scaled) means
  # rockiness: Detection probability increases on sandier substrates
  # ground temp: bell-curve shape with positive (linear) and negative (quadratic)
    # Detection increases at moderate ground temperatures. Detection declines at high ground temperatures, likely due to: rapid track degradation
    #substrate crusting or slumping, increased wind activity later in the day
    # This pattern is biologically intuitive for track-based surveys and supports using ground temperature as a proxy for substrate condition rather than animal activity.
  # wind: strong negative effect

fit.detection.formula <- ~ rock + wind

global_model <- PGOcc(
  occ.formula = global.occ.formula,
  det.formula = global.det.formula, 
  data = list(
    y       = y,
    occ.covs = occ_df,
    det.covs = X.p),
  inits = gm.inits, 
  n.samples = n.samples, 
  priors = gm.priors, 
  n.omp.threads = 1, 
  verbose = FALSE, 
  n.report = 1000, 
  n.burn = n.burn, 
  n.thin = n.thin, 
  n.chains = n.chains)

best.fit_model <- PGOcc(occ.formula = fit.occupancy.formula,
      det.formula = fit.detection.formula, 
      data = list(
        y = y,
        occ.covs = occ_df,   
        det.covs = X.p),    
      inits = gm.inits, 
      n.samples = n.samples, 
      priors = gm.priors, 
      n.omp.threads = 1, 
      verbose = FALSE, 
      n.report = 2000, 
      n.burn = n.burn, 
      n.thin = n.thin, 
      n.chains = n.chains)

null_model <- PGOcc(occ.formula = ~ 1,
                        det.formula = fit.detection.formula, 
                        data = list(
                          y = y,
                          occ.covs = occ_df,   
                          det.covs = X.p),    
                        inits = gm.inits, 
                        n.samples = n.samples, 
                        priors = gm.priors, 
                        n.omp.threads = 1, 
                        verbose = FALSE, 
                        n.report = 2000, 
                        n.burn = n.burn, 
                        n.thin = n.thin, 
                        n.chains = n.chains)

summary(global_model)
summary(best.fit_model)
summary(null_model)
waicOcc(global_model)
waicOcc(best.fit_model)
waicOcc(null_model)


# k fold cross validation
global_cv  <- PGOcc(occ.formula = global.occ.formula,
                      det.formula = global.det.formula, 
                      data = list(
                        y = y,
                        occ.covs = occ_df,   
                        det.covs = X.p),    
                      inits = gm.inits, 
                      n.samples = n.samples, 
                      priors = gm.priors, 
                      n.omp.threads = 1, 
                      verbose = FALSE, 
                      n.report = 2000, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      k.fold = 3,
                      k.fold.threads = 3)

best.fit_cv  <- PGOcc(occ.formula = fit.occupancy.formula,
                        det.formula = fit.detection.formula, 
                        data = list(
                          y = y,
                          occ.covs = occ_df,   
                          det.covs = X.p),    
                        inits = gm.inits, 
                        n.samples = n.samples, 
                        priors = gm.priors, 
                        n.omp.threads = 1, 
                        verbose = FALSE, 
                        n.report = 2000, 
                        n.burn = n.burn, 
                        n.thin = n.thin, 
                        n.chains = n.chains,
                        k.fold = 3,
                        k.fold.threads = 3)

null_cv  <- PGOcc(occ.formula = ~1,
                      det.formula = fit.detection.formula, 
                      data = list(
                        y = y,
                        occ.covs = occ_df,   
                        det.covs = X.p),    
                      inits = gm.inits, 
                      n.samples = n.samples, 
                      priors = gm.priors, 
                      n.omp.threads = 1, 
                      verbose = FALSE, 
                      n.report = 2000, 
                      n.burn = n.burn, 
                      n.thin = n.thin, 
                      n.chains = n.chains,
                      k.fold = 3,
                      k.fold.threads = 3)

summary(global_cv)
summary(best.fit_cv)
summary(null_cv)

global_cv$k.fold.deviance
best.fit_cv$k.fold.deviance
null_cv$k.fold.deviance

# occupancy and detection estimates of best fit model
plogis(1.7390)
plogis(-3.4810)

plot(best.fit_model, 'beta', density = FALSE) # Occupancy parameters.
plot(best.fit_model, 'alpha', density = FALSE) # Detection parameters.

ppc.out <- ppcOcc(best.fit_model, fit.stat = 'freeman-tukey', group = 1)
ppc.out2 <- ppcOcc(best.fit_model, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.out)
summary(ppc.out2)
