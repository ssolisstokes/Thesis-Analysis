# ---- Libraries ----
library(dplyr)
library(tidyr)
library(spOccupancy)
library(car)

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


# calculate time since midnight
time <- strsplit(ifelse(is.na(gm$time), "NA:NA:NA", gm$time), ":", fixed = TRUE)
hh <- as.numeric(sapply(time, `[`, 1))
mm <- as.numeric(sapply(time, `[`, 2))
ss <- as.numeric(sapply(time, `[`, 3))
gm$time_sec <- hh*3600 + mm*60 + ss
gm$time_sec[!is.finite(gm$time_sec)] <- NA_real_

# replacement detection covariates
replacement.covs <- c("utm_x", "utm_y", "time", "detection","avg_wind_spd", "air_temp", "rel_hum", "gr_temp", "rock_ind", "shr_nw", "shr_ne", "shr_se", "shr_sw", "shr_pt_density", "vpd", "time_sec")

# function for replacing the detection point to appropriate time/collection point
gm <- gm %>%
  group_by(site_id, tran_id) %>%
  dplyr::group_modify(~ {
    df <- .x
    col.pts <- which(df$col_point %in% 1:5 & !is.na(df$time_sec))    # col_point 1-5
    det.pt   <- which(df$col_point == 6   & !is.na(df$time_sec))     # col_point 6 (detection point)
    if (length(col.pts) == 0 || length(det.pt) == 0) return(df)      # if no col_point 1-5 or col_point 6, do nothing ()
    
    for (i in det.pt) {                                               # for each col_point 6 row, find nearest-in-time col_point 1–5 and copy covariates
      dif   <- abs(df$time_sec[col.pts] - df$time_sec[i])
      nearest <- col.pts[which.min(dif)]
      df[nearest, replacement.covs] <- df[i, replacement.covs]
    }
    df
  }) %>%
  ungroup()


# keep only on-transect collection points 1–5 (removes off-transect detections)
transects <- subset(gm, col_point %in% 1:5)

# sanity check transects and detection points: 
# transects %>% filter(detection == 1)

# add the transect length (path_km)
transects <- merge(transects, t.length[, c("site_id","tran_id","path_km")],
                   by = c("site_id","tran_id"), all.x = TRUE)

# ---- Transect-level summaries (replicates = transects) ----

# transect-level observation covariates from the 1–5 col. points
tran.covs <- transects %>%
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

tran.summary <- tran.covs %>%
  group_by(site_id) %>%
  arrange(tran_id, .by_group = TRUE) %>%
  mutate(rep_idx = row_number()) %>%
  ungroup()

all.sites <- sort(unique(tran.summary$site_id))

# max transects at a site
tran.max <- tran.summary %>%
  count(site_id, name = "n_tran") %>%
  summarise(max(n_tran)) %>%
  pull()

# ---- Detection Matrix ----
# (rows = sites, cols = transect replicates)

# long df
y_long <- expand.grid(
  site_id = all.sites,
  rep_idx = seq_len(tran.max), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(tran.summary %>% 
            select(site_id, rep_idx, detection),
            by = c("site_id","rep_idx")) %>%
  arrange(site_id, rep_idx)

# wide df
y_wide <- y_long %>%
  pivot_wider(names_from = rep_idx, values_from = detection, names_prefix = "rep_") %>%
  arrange(site_id)

y <- as.matrix(y_wide %>% select(starts_with("rep_")))
rownames(y) <- y_wide$site_id

# ---- Site Covariates  ----
# transect-level summaries (one row per site)
site.covs <- tran.summary %>%
  group_by(site_id) %>%
  summarise(
    rockiness_mean = mean(rockiness, na.rm = TRUE),
    shrub_density_mean = mean(shrub, na.rm = TRUE),
    elevation_mean = mean(elev, na.rm = TRUE),
    prey_total = sum(prey, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(site_id)

# scale site covariates
site.covs$rockiness_z <- as.numeric(scale(site.covs$rockiness_mean))
site.covs$shrubs_z <- as.numeric(scale(site.covs$shrub_density_mean))
site.covs$elev_z <- as.numeric(scale(site.covs$elevation_mean))
site.covs$prey_z <- as.numeric(scale(site.covs$prey_total))

# align rows to y (sites as row names)
site.covs <- site.covs %>% slice(match(rownames(y), site_id))

# ---- Detection Covariates ----
# start from transect-level covs, pad to max transects per site
det.covs <- expand.grid(
  site_id = all.sites, 
  rep_idx = seq_len(tran.max),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  left_join(tran.summary %>%
              select(site_id, rep_idx, time_sec, air_temp, gr_temp, wind, vpd, rockiness),
            by = c("site_id","rep_idx")) %>%
  arrange(site_id, rep_idx)

# widen each covariate
time_wide <- det.covs %>% select(site_id, rep_idx, time_sec) %>%
  pivot_wider(names_from = rep_idx, values_from = time_sec, names_prefix = "rep_") %>% arrange(site_id)
a.temp_wide <- det.covs %>% select(site_id, rep_idx, air_temp) %>%
  pivot_wider(names_from = rep_idx, values_from = air_temp, names_prefix = "rep_") %>% arrange(site_id)
gr.temp_wide <- det.covs %>% select(site_id, rep_idx, gr_temp) %>%
  pivot_wider(names_from = rep_idx, values_from = gr_temp, names_prefix = "rep_") %>% arrange(site_id)
wind_wide <- det.covs %>% select(site_id, rep_idx, wind) %>%
  pivot_wider(names_from = rep_idx, values_from = wind, names_prefix = "rep_") %>% arrange(site_id)
vpd_wide <- det.covs %>% select(site_id, rep_idx, vpd) %>%
  pivot_wider(names_from = rep_idx, values_from = vpd, names_prefix = "rep_") %>% arrange(site_id)
rock_wide <- det.covs %>% select(site_id, rep_idx, rockiness) %>%
  pivot_wider(names_from = rep_idx, values_from = rockiness, names_prefix = "rep_") %>% arrange(site_id)

# covariate matricies; alligned to y matrix
time_mat <- as.matrix(time_wide %>% select(starts_with("rep_")))[match(rownames(y), time_wide$site_id), ]
a.temp_mat <- as.matrix(a.temp_wide  %>% select(starts_with("rep_")))[match(rownames(y), a.temp_wide$site_id), ]
gr.temp_mat <- as.matrix(gr.temp_wide   %>% select(starts_with("rep_")))[match(rownames(y), gr.temp_wide$site_id), ]
wind_mat <- as.matrix(wind_wide %>% select(starts_with("rep_")))[match(rownames(y), wind_wide$site_id), ]
vpd_mat <- as.matrix(vpd_wide  %>% select(starts_with("rep_")))[match(rownames(y), vpd_wide$site_id), ]
rock_mat <- as.matrix(rock_wide %>% select(starts_with("rep_")))[match(rownames(y), rock_wide$site_id), ]

# z-score scale each detection matrix; across all cells
time_mat[] <- (time_mat - mean(time_mat, na.rm = TRUE))/sd(time_mat, na.rm = TRUE)
a.temp_mat[] <- (a.temp_mat - mean(a.temp_mat, na.rm = TRUE))/sd(a.temp_mat,  na.rm = TRUE)
gr.temp_mat[] <- (gr.temp_mat - mean(gr.temp_mat, na.rm = TRUE))/sd(gr.temp_mat,   na.rm = TRUE)
wind_mat[] <- (wind_mat - mean(wind_mat, na.rm = TRUE))/sd(wind_mat, na.rm = TRUE)
vpd_mat[] <- (vpd_mat - mean(vpd_mat, na.rm = TRUE))/sd(vpd_mat,  na.rm = TRUE)
rock_mat[] <- (rock_mat - mean(rock_mat, na.rm = TRUE))/sd(rock_mat, na.rm = TRUE)

# ---- Prep Covariates for PGOcc ----
occ.cov_df <- data.frame(
  rockiness = site.covs$rockiness_z,
  shrubs    = site.covs$shrubs_z,
  elev      = site.covs$elev_z,
  prey      = site.covs$prey_z)
rownames(occ.cov_df) <- site.covs$site_id

det.cov_list <- list(
  time = time_mat,
  vpd  = vpd_mat,
  rock = rock_mat,
  t_air = a.temp_mat,
  t_gr  = gr.temp_mat,
  wind  = wind_mat)

## ---- Check Collinearity ----
# occupancy 
cor(occ.cov_df, use = "pairwise.complete.obs")
round(cor(occ.cov_df, use = "pairwise.complete.obs"), 2)
pairs(occ.cov_df )
vif(lm(rep(1, nrow(occ.cov_df)) ~ ., data = occ.cov_df))

# detection
# convert det_list to df 
det.cov_df <- data.frame(
  rock = as.vector(det.cov_list$rock),
  t_air = as.vector(det.cov_list$t_air),
  t_gr  = as.vector(det.cov_list$t_gr),
  vpd   = as.vector(det.cov_list$vpd),
  wind  = as.vector(det.cov_list$wind))

det.cov_df <- det.cov_df[complete.cases(det.cov_df), ]
round(cor(det.cov_df, use ="pairwise.complete.obs" ), 2)
pairs(det.cov_df)
vif(lm(rep(1, nrow(det.cov_df)) ~ ., data = det.cov_df))

# ---- Inits, Priors, MCMC Config ----
gm.inits <- list(
  alpha = 0, 
  beta  = 0, 
  z = apply(y, 1, max, na.rm = TRUE))

gm.priors <- list(
  alpha.normal = list(mean = 0, var = 2.72), 
  beta.normal  = list(mean = 0, var = 2.72))

n.samples <- 10000
n.burn    <- 5000
n.thin    <- 5
n.chains  <- 3

# ---- Global Model ----
global.occ.formula <- ~ rockiness + shrubs + elev + prey
global.det.formula <- ~ rock + I(t_air^2) + t_air + vpd + wind + I(t_gr^2) + t_gr

global_model <- PGOcc(
  occ.formula = global.occ.formula,
  det.formula = global.det.formula, 
  data = list(
    y = y,
    occ.covs = occ.cov_df,
    det.covs = det.cov_list),
  inits = gm.inits, 
  n.samples = n.samples, 
  priors = gm.priors, 
  n.omp.threads = 1, 
  verbose = TRUE,
  n.report = 2000, 
  n.burn = n.burn, 
  n.thin = n.thin, 
  n.chains = n.chains)

summary(global_model)
# check MCMC convergence; trace plots
plot(global_model, 'beta', density = FALSE) # occupancy parameters
plot(global_model, 'alpha', density = FALSE) # detection parameters

##  ---- Posterior Predictive Checks - Global model ----
# check occupancy and detection intercept
# plogis() 

# freeman-tukey GoF
global.site_ppc <- ppcOcc(global_model, fit.stat = 'freeman-tukey', group = 1)   # site
global.rep_ppc <- ppcOcc(global_model, fit.stat = 'freeman-tukey', group = 2)  # replicate
summary(global.site_ppc)
summary(global.rep_ppc)

# visualization 
# create df of group 1 (site); color for actual data set
global.site.ppc_df <- data.frame(fit = global.site_ppc$fit.y, 
                     fit.rep = global.site_ppc$fit.y.rep, 
                     color = 'lightskyblue1')
# set color for fit data
global.site.ppc_df$color[global.site.ppc_df$fit.rep > global.site.ppc_df$fit] <- 'lightsalmon'

# plot
plot(global.site.ppc_df$fit, global.site.ppc_df$fit.rep, bg = global.site.ppc_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True') 
  lines(global.site.ppc_df$fit, global.site.ppc_df$fit, col = 'black')

# create df of group 2 (replicates); color for actual data set
global.rep.ppc_df <- data.frame(fit = global.rep_ppc$fit.y, 
                                 fit.rep = global.rep_ppc$fit.y.rep, 
                                 color = 'lightskyblue1')
# set color for fit data
global.rep.ppc_df$color[global.rep.ppc_df$fit.rep > global.rep.ppc_df$fit] <- 'lightsalmon'

# plot
plot(global.rep.ppc_df$fit, global.rep.ppc_df$fit.rep, bg = global.rep.ppc_df$color, pch = 21, 
                              ylab = 'Fit', xlab = 'True') 
  lines(global.rep.ppc_df$fit, global.rep.ppc_df$fit, col = 'black')


# ---- Model Selection ----
## ---- Occupancy ----
### ---- Candidate Models ----
# occupancy formulas + global detection
om1 <- PGOcc(occ.formula = ~ 1,
             det.formula = global.det.formula, 
             data = list(
               y = y,
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
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

### ---- WAIC ----
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

om.waic_df <- bind_rows(list(om1.df, om2.df, om3.df, om4.df, om5.df, om6.df, om7.df, om8.df))

om.waic_tbl <- 
  om.waic_df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta)) # Burnham & Andserson 2002, 2004 model weights formula
    ) %>%  
  arrange(WAIC)


om.waic_tbl

### ---- Variable Importance & Model Averaging ----
#### ---- Variable Importance ----
# model weights as objects
om1.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om1"]
om2.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om2"]
om3.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om3"]
om4.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om4"]
om5.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om5"]
om6.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om6"]
om7.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om7"]
om8.weight <- om.waic_tbl$weight[om.waic_tbl$model == "om8"]

# covariates in each model
# rockiness: om2, om3, om4, om5, om7
# shrub density: om2, om3, om4, om6, om8
# elevation: om2
# prey presence: om3, om5, om6

# calculate sum weight and variable importance
rockiness.sum_weight <- om2.weight + om3.weight + om4.weight + om5.weight + om7.weight
shrub.sum_weight <- om2.weight + om3.weight + om4.weight + om6.weight + om8.weight
elev.sum_weight <- om2.weight
prey.sum_weight <- om3.weight + om5.weight + om6.weight

om.var.importance <- data.frame(
  variable = c("rockiness","shrub","elev","prey"),
  sum_weight = c(rockiness.sum_weight, shrub.sum_weight, elev.sum_weight, prey.sum_weight))

om.var.importance

#### ---- Model Averages -----
# om1 posterior means
om1.mat  <- as.matrix(om1$beta.samples)
om1.means <- colMeans(om1.mat)

om1_int <- om1.means["(Intercept)"]
om1_rockiness <- 0          
om1_shrubs <- 0         
om1_elev <- 0          
om1_prey <- 0          

# om2 posterior means
om2.mat  <- as.matrix(om2$beta.samples)
om2.means <- colMeans(om2.mat)

om2_int <- om2.means["(Intercept)"]
om2_rockiness <- om2.means["rockiness"]
om2_shrubs <- om2.means["shrubs"]
om2_elev <- om2.means["elev"]
om2_prey <- 0         

# om3 posterior means
om3.mat  <- as.matrix(om3$beta.samples)
om3.means <- colMeans(om3.mat)

om3_int <- om3.means["(Intercept)"]
om3_rockiness <- om3.means["rockiness"]
om3_shrubs <- om3.means["shrubs"]
om3_elev <- 0          
om3_prey <- om3.means["prey"]

# om4 posterior means
om4.mat  <- as.matrix(om4$beta.samples)
om4.means <- colMeans(om4.mat)

om4_int <- om4.means["(Intercept)"]
om4_rockiness <- om4.means["rockiness"]
om4_shrubs <- om4.means["shrubs"]
om4_elev <- 0          
om4_prey <- 0          

# om5 posterior means
om5.mat  <- as.matrix(om5$beta.samples)
om5.means <- colMeans(om5.mat)

om5_int <- om5.means["(Intercept)"]
om5_rockiness <- om5.means["rockiness"]
om5_shrubs <- 0         
om5_elev <- 0          
om5_prey <- om5.means["prey"]

# om6 posterior means
om6.mat  <- as.matrix(om6$beta.samples)
om6.means <- colMeans(om6.mat)

om6_int <- om6.means["(Intercept)"]
om6_rockiness <- 0       
om6_shrubs <- om6.means["shrubs"]
om6_elev <- 0      
om6_prey <- om6.means["prey"]

# om7 posterior means
om7.mat  <- as.matrix(om7$beta.samples)
om7.means <- colMeans(om7.mat)

om7_int <- om7.means["(Intercept)"]
om7_rockiness <- om7.means["rockiness"]
om7_shrubs <- 0          
om7_elev <- 0          
om7_prey <- 0          

# om8 posterior means
om8.mat  <- as.matrix(om8$beta.samples)
om8.means <- colMeans(om8.mat)

om8_int <- om8.means["(Intercept)"]
om8_rockiness <- 0       
om8_shrubs <- om8.means["shrubs"]
om8_elev <- 0       
om8_prey <- 0         


# multiply each coefficient by its model’s weight and sum across models
om.avg_int <- 
  (om1_int * om1.weight) + (om2_int * om2.weight) +
  (om3_int * om3.weight) + (om4_int * om4.weight) +
  (om5_int * om5.weight) + (om6_int * om6.weight) +
  (om7_int * om7.weight) + (om8_int * om8.weight)

om.avg_rockiness <- 
  (om1_rockiness * om1.weight) + (om2_rockiness * om2.weight) +
  (om3_rockiness * om3.weight) + (om4_rockiness * om4.weight) +
  (om5_rockiness * om5.weight) + (om6_rockiness * om6.weight) +
  (om7_rockiness * om7.weight) + (om8_rockiness * om8.weight)

om.avg_shrubs <- 
  (om1_shrubs * om1.weight) + (om2_shrubs * om2.weight) +
  (om3_shrubs * om3.weight) + (om4_shrubs * om4.weight) +
  (om5_shrubs * om5.weight) + (om6_shrubs * om6.weight) +
  (om7_shrubs * om7.weight) + (om8_shrubs * om8.weight)

om.avg_elev <- 
  (om1_elev * om1.weight) + (om2_elev * om2.weight) +
  (om3_elev * om3.weight) + (om4_elev * om4.weight) +
  (om5_elev * om5.weight) + (om6_elev * om6.weight) +
  (om7_elev * om7.weight) + (om8_elev * om8.weight)

om.avg_prey <- 
  (om1_prey * om1.weight) + (om2_prey * om2.weight) +
  (om3_prey * om3.weight) + (om4_prey * om4.weight) +
  (om5_prey * om5.weight) + (om6_prey * om6.weight) +
  (om7_prey * om7.weight) + (om8_prey * om8.weight)

# final data frame
om.modavg <- data.frame(
  avg_coefficients = c(om.avg_int, om.avg_rockiness, om.avg_shrubs, om.avg_elev, om.avg_prey))

om.modavg
### ----- Best-fit Occupancy Formula ----
# best-fit occupancy formula
fit.occupancy.formula <- ~ rockiness + shrubs + prey
## ---- Detection ----
### ---- Candidate Models ----
dm1 <- PGOcc(occ.formula = fit.occupancy.formula,
             det.formula = ~ 1, 
             data = list(
               y = y,
               occ.covs = occ.cov_df,
               det.covs = det.cov_list),
             inits = gm.inits, 
             n.samples = n.samples, 
             priors = gm.priors, 
             n.omp.threads = 1, 
             verbose = FALSE, 
             n.report = 2000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

dm2 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock, 
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)


dm3 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ t_gr + I(t_gr^2), 
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm4 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + wind, 
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm5 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + t_gr + I(t_gr^2),
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm6 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ wind + t_gr + I(t_gr^2), 
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)

dm7 <- PGOcc(occ.formula = fit.occupancy.formula,
              det.formula = ~ rock + wind + t_gr + I(t_gr^2),
              data = list(
                y = y,
                occ.covs = occ.cov_df,
                det.covs = det.cov_list),
              inits = gm.inits, 
              n.samples = n.samples, 
              priors = gm.priors, 
              n.omp.threads = 1, 
              verbose = FALSE, 
              n.report = 2000, 
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains)


### ---- WAIC ----
dm1.waic <- waicOcc(dm1)
dm2.waic <- waicOcc(dm2)
dm3.waic <- waicOcc(dm3)
dm4.waic <- waicOcc(dm4)
dm5.waic <- waicOcc(dm5)
dm6.waic <- waicOcc(dm6)
dm7.waic <- waicOcc(dm7)

dm1.df <- data.frame(model = "dm1", t(dm1.waic))
dm2.df <- data.frame(model = "dm2", t(dm2.waic))
dm3.df <- data.frame(model = "dm3", t(dm3.waic))
dm4.df <- data.frame(model = "dm4", t(dm4.waic))
dm5.df <- data.frame(model = "dm5", t(dm5.waic))
dm6.df <- data.frame(model = "dm6", t(dm6.waic))
dm7.df <- data.frame(model = "dm7", t(dm7.waic))

dm.waic_df <- bind_rows(list(dm1.df, dm2.df, dm3.df, dm4.df, dm5.df, dm6.df, dm7.df))

dm.waic_tbl <- 
  dm.waic_df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta))  # Burnham & Anderson 2002, 2004 model weights formula
    ) %>%
  arrange(WAIC)


dm.waic_tbl
### ---- Variable Importance & Model Averaging ----
#### ---- Variable Importance ----
# model weights as objects
dm1.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm1"]
dm2.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm2"]
dm3.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm3"]
dm4.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm4"]
dm5.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm5"]
dm6.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm6"]
dm7.weight <- dm.waic_tbl$weight[dm.waic_tbl$model == "dm7"]

# covariates in each model
# rock: dm2, dm4, dm5, dm7
# ground temp: dm3, dm5, dm6, dm7
# wind: dm4, dm6, dm7

rock.sum_weight <- dm2.weight + dm4.weight + dm5.weight + dm7.weight
t_gr.sum_weight <- dm3.weight + dm5.weight + dm6.weight + dm7.weight
wind.sum_weight <- dm4.weight + dm6.weight + dm7.weight

dm.var.importance <- data.frame(
  variable = c("rock","t_gr", "wind"),
  sum_weight = c(rock.sum_weight, t_gr.sum_weight,wind.sum_weight))

dm.var.importance

#### ---- Model Averages ----
# dm1 posterior means
dm1.mat <- as.matrix(dm1$alpha.samples)
dm1.means <- colMeans(dm1.mat)

dm1_int   <- dm1.means["(Intercept)"]
dm1_rock  <- 0
dm1_tgr   <- 0
dm1_tgr2  <- 0
dm1_wind  <- 0

# dm2 posterior means
dm2.mat <- as.matrix(dm2$alpha.samples)
dm2.means <- colMeans(dm2.mat)

dm2_int  <- dm2.means["(Intercept)"]
dm2_rock <- dm2.means["rock"]
dm2_tgr  <- 0
dm2_tgr2 <- 0
dm2_wind <- 0

# dm3 posterior means
dm3.mat <- as.matrix(dm3$alpha.samples)
dm3.means <- colMeans(dm3.mat)

dm3_int  <- dm3.means["(Intercept)"]
dm3_rock <- 0
dm3_tgr  <- dm3.means["t_gr"]
dm3_tgr2 <- dm3.means["I(t_gr^2)"]
dm3_wind <- 0

# dm4 posterior means
dm4.mat <- as.matrix(dm4$alpha.samples)
dm4.means <- colMeans(dm4.mat)

dm4_int  <- dm4.means["(Intercept)"]
dm4_rock <- dm4.means["rock"]
dm4_tgr  <- 0
dm4_tgr2 <- 0
dm4_wind <- dm4.means["wind"]

# dm5 posterior means
dm5.mat <- as.matrix(dm5$alpha.samples)
dm5.means <- colMeans(dm5.mat)

dm5_int  <- dm5.means["(Intercept)"]
dm5_rock <- dm5.means["rock"]
dm5_tgr  <- dm5.means["t_gr"]
dm5_tgr2 <- dm5.means["I(t_gr^2)"]
dm5_wind <- 0

# dm6 posterior means
dm6.mat <- as.matrix(dm6$alpha.samples)
dm6.means <- colMeans(dm6.mat)

dm6_int  <- dm6.means["(Intercept)"]
dm6_rock <- 0
dm6_tgr  <- dm6.means["t_gr"]
dm6_tgr2 <- dm6.means["I(t_gr^2)"]
dm6_wind <- dm6.means["wind"]

# dm7 posterior means
dm7.mat <- as.matrix(dm7$alpha.samples)
dm7.means <- colMeans(dm7.mat)

dm7_int  <- dm7.means["(Intercept)"]
dm7_rock <- dm7.means["rock"]
dm7_tgr  <- dm7.means["t_gr"]
dm7_tgr2 <- dm7.means["I(t_gr^2)"]
dm7_wind <- dm7.means["wind"]

# multiply each coefficient by its model’s weight and sum across models
dm.avg_int <-
  (dm1_int * dm1.weight) + (dm2_int * dm2.weight) + 
  (dm3_int * dm3.weight) + (dm4_int * dm4.weight) + 
  (dm5_int * dm5.weight) + (dm6_int * dm6.weight) + 
  (dm7_int * dm7.weight)

dm.avg_rock <-
  (dm1_rock * dm1.weight) + (dm2_rock * dm2.weight) + 
  (dm3_rock * dm3.weight) + (dm4_rock * dm4.weight) + 
  (dm5_rock * dm5.weight) + (dm6_rock * dm6.weight) + 
  (dm7_rock * dm7.weight)

dm.avg_tgr <-
  (dm1_tgr * dm1.weight) + (dm2_tgr * dm2.weight) + 
  (dm3_tgr * dm3.weight) + (dm4_tgr * dm4.weight) + 
  (dm5_tgr * dm5.weight) + (dm6_tgr * dm6.weight) + 
  (dm7_tgr * dm7.weight)

dm.avg_tgr2 <-
  (dm1_tgr2 * dm1.weight) + (dm2_tgr2 * dm2.weight) + 
  (dm3_tgr2 * dm3.weight) + (dm4_tgr2 * dm4.weight) + 
  (dm5_tgr2 * dm5.weight) + (dm6_tgr2 * dm6.weight) + 
  (dm7_tgr2 * dm7.weight)

dm.avg_wind <-
  (dm1_wind * dm1.weight) + (dm2_wind * dm2.weight) + 
  (dm3_wind * dm3.weight) + (dm4_wind * dm4.weight) + 
  (dm5_wind * dm5.weight) + (dm6_wind * dm6.weight) + 
  (dm7_wind * dm7.weight)

# final data frame
dm.modavg <- data.frame(
  avg_coefficients = c(dm.avg_int, dm.avg_rock, dm.avg_tgr, dm.avg_tgr2, dm.avg_wind))

dm.modavg
### ---- Best-fit Detection Formula ----
fit.detection.formula <- ~ rock + wind

# ---- Best-fit Model ----

best.fit_model <- PGOcc(occ.formula = fit.occupancy.formula,
                        det.formula = fit.detection.formula, 
                        data = list(
                          y = y,
                          occ.covs = occ.cov_df,
                          det.covs = det.cov_list),
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
                      occ.covs = occ.cov_df,
                      det.covs = det.cov_list),
                    inits = gm.inits, 
                    n.samples = n.samples, 
                    priors = gm.priors, 
                    n.omp.threads = 1, 
                    verbose = FALSE, 
                    n.report = 2000, 
                    n.burn = n.burn, 
                    n.thin = n.thin, 
                    n.chains = n.chains)


## ---- Best-fit vs Global vs Null ----
summary(global_model)
summary(best.fit_model)
summary(null_model)

### ---- WAIC ----
global.waic <- waicOcc(global_model)
best.fit.waic <- waicOcc(best.fit_model)
null.waic <- waicOcc(null_model)

global.waic_df <- data.frame(model = "global_model", t(global.waic))
best.fit.waic_df <- data.frame(model = "best.fit_model", t(best.fit.waic))
null.waic_df <- data.frame(model = "null_model", t(null.waic))

models.waic_df <- bind_rows(list(global.waic_df, best.fit.waic_df, null.waic_df))

models.waic_tbl <- 
  models.waic_df %>%
  mutate(
    elpd = as.numeric(elpd),
    pD = as.numeric(pD),
    WAIC = as.numeric(WAIC)) %>%
  mutate(
    delta = WAIC - min(WAIC),
    weight = exp(-0.5 * delta)/ sum (exp(-0.5 * delta))   # Burnham & Anderson 2002, 2004 model weights formula
    ) %>%
  arrange(WAIC)


models.waic_tbl
### ---- K-fold CV ----

global_cv  <- PGOcc(occ.formula = global.occ.formula,
                    det.formula = global.det.formula, 
                    data = list(
                      y = y,
                      occ.covs = occ.cov_df,
                      det.covs = det.cov_list),
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
                        occ.covs = occ.cov_df,
                        det.covs = det.cov_list),
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
                    occ.covs = occ.cov_df,
                    det.covs = det.cov_list),
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

global.kfold.dev <- global_cv$k.fold.deviance
best.fit.kfold.dev <- best.fit_cv$k.fold.deviance
null.kfold.dev <- null_cv$k.fold.deviance

kfold.dev_tbl <- data.frame(
  global = global.kfold.dev,
  best = best.fit.kfold.dev,
  null = null.kfold.dev)

kfold.dev_tbl

##  ---- Posterior Predictive Checks - Best-fit Model ----
# check occupancy and detection intercept
# plogis() 

# freeman-tukey GoF
fit.site_ppc <- ppcOcc(best.fit_model, fit.stat = 'freeman-tukey', group = 1)   # site
fit.rep_ppc <- ppcOcc(best.fit_model, fit.stat = 'freeman-tukey', group = 2)  # replicate
summary(fit.site_ppc)
summary(fit.rep_ppc)

# visualization 
# create df of group 1 (site); color for actual data set
fit.site.ppc_df <- data.frame(fit = fit.site_ppc$fit.y, 
                                 fit.rep = fit.site_ppc$fit.y.rep, 
                                 color = 'lightskyblue1')
# set color for fit data
fit.site.ppc_df$color[fit.site.ppc_df$fit.rep > fit.site.ppc_df$fit] <- 'lightsalmon'

# plot
plot(fit.site.ppc_df$fit, fit.site.ppc_df$fit.rep, bg = fit.site.ppc_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True') 
  lines(fit.site.ppc_df$fit, fit.site.ppc_df$fit, col = 'black')


# create df of group 2 (replicates); color for actual data set
fit.rep.ppc_df <- data.frame(fit = fit.rep_ppc$fit.y, 
                              fit.rep = fit.rep_ppc$fit.y.rep, 
                              color = 'lightskyblue1')
# set color for fit data
fit.rep.ppc_df$color[fit.rep.ppc_df$fit.rep > fit.rep.ppc_df$fit] <- 'lightsalmon'

# plot
plot(fit.rep.ppc_df$fit, fit.rep.ppc_df$fit.rep, bg = fit.rep.ppc_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True') 
  lines(fit.rep.ppc_df$fit, fit.rep.ppc_df$fit, col = 'black')

