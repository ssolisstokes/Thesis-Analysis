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

X.p <- list(
  time = time_mat,
  vpd  = vpd_mat,
  rock = rock_mat,
  t_air = air_mat,
  t_gr  = gr_mat,
  wind  = wind_mat)

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


out.small <- PGOcc(occ.formula = ~ 1,
                   det.formula = global.det.formula, 
                   data = list(
                     y = y,
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

waicOcc(out.small)
ppc.out.small <- ppcOcc(out.small, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.small)



# ---- Sites ----
# turn rownames into an explicit site_id column
occ_df$site_id <- rownames(occ_df)

# now make it a factor 
occ_df$site_factor <- factor(occ_df$site_id)


gm.inits <- list(
  alpha = 0, 
  beta  = 0, 
  z = apply(y, 1, max, na.rm = TRUE))

gm.priors <- list(
  alpha.normal = list(mean = 0, var = 2.72), 
  beta.normal  = list(mean = 0, var = 2.72))

n.samples <- 15000
n.burn    <- 3000
n.thin    <- 2
n.chains  <- 3

out.sites <- PGOcc(
  occ.formula = ~ site_factor,
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

summary(out.sites)


occ_formula <- eval(out.sites$call$occ.formula)


occ_df <- as.data.frame(occ_df)


X0 <- model.matrix(occ_formula, data = occ_df)

psi_pred <- predict(
  out.sites,
  X.0  = X0,
  type = "occupancy")

names(psi_pred)
psi_pred

# Posterior samples: rows = iterations, cols = sites
psi_samps <- psi_pred$psi.0.samples
dim(psi_samps)

# Attach site IDs to columns (make sure order matches occ_df)
colnames(psi_samps) <- occ_df$site_id

# Summarize mean and 95% credible interval per site
psi_df <- apply(psi_samps, 2, 
                function(x) {
                  c(psi_mean = mean(x),
                    psi_lcl  = quantile(x, 0.025),
                    psi_ucl  = quantile(x, 0.975))}) %>%
  t() %>% 
  as.data.frame()

psi_df$site_id <- rownames(psi_df)

psi_df <- psi_df %>%
  rename(psi_lcl = `psi_lcl.2.5%`, psi_ucl = `psi_ucl.97.5%`)


ggplot(psi_df, aes(x = reorder(site_id, psi_mean), y = psi_mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = psi_lcl, ymax = psi_ucl), width = 0.1) +
  coord_flip() +
  labs(
    x = "Site",
    y = "Occupancy probability (ψ)",
    title = "Posterior occupancy by site") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10))

## ---- Habitat Quality::Sites ----

high_qual <- c("paradise_canyon", "fort_pearce")
med_qual <- c("sun_river", "cove_wash")
low_qual <- c("white_reef", "turkey_farm")

# Posterior samples matrix: iterations × sites
psi_samps <- psi_pred$psi.0.samples
colnames(psi_samps) <- occ_df$site_id   # ensure names are attached

# Class-level posterior distributions
psi_high <- rowMeans(psi_samps[, high_qual])
psi_med  <- rowMeans(psi_samps[, med_qual])
psi_low  <- rowMeans(psi_samps[, low_qual])

psi_class_df <- data.frame(
  class = c("High", "Medium", "Low"),
  psi_mean = c(mean(psi_high), mean(psi_med), mean(psi_low)),
  psi_lcl  = c(quantile(psi_high, 0.025),
               quantile(psi_med, 0.025),
               quantile(psi_low, 0.025)),
  psi_ucl  = c(quantile(psi_high, 0.975),
               quantile(psi_med, 0.975),
               quantile(psi_low, 0.975)))

psi_class_df

prob_hi.med <- mean(psi_high > psi_med)
prob_hi.low <- mean(psi_high > psi_low)
prob_med.low <- mean(psi_med > psi_low)

prob_hi.med
prob_hi.low
prob_med.low


ggplot(psi_class_df, aes(x = class, y = psi_mean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = psi_lcl, ymax = psi_ucl), width = 0.1) +
  labs(
    x = "Habitat Quality",
    y = "Occupancy Probability (ψ)",
    title = " Predicted Occupancy by Habitat Quality") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text  = element_text(size = 11))


df_dens <- data.frame(
  high = psi_high,
  medium = psi_med,
  low = psi_low) %>%
  pivot_longer(everything(), names_to = "class", values_to = "psi")

ggplot(df_dens, aes(x = psi, fill = class)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Posterior distributions of occupancy by habitat class",
    x = "Occupancy probability (ψ)",
    y = "Density"
  ) +
  theme_bw()


# ---- Model Selection ----

# candidate occupancy formulas
occ_formulas <- list(
  om1 = ~ rockiness + shrubs + elev,      # habitat + elev
  om2 = ~ rockiness + shrubs + prey,     # habitat + prey
  om3 = ~ rockiness + shrubs,         # habitat
  om4 = ~ rockiness + prey,       # rockiness + prey
  om5 = ~ shrubs + prey,         # vegetation + prey
  om6 = ~ rockiness,                                 
  om7 = ~ shrubs)

# function to fit a PGOcc model with given occ formula
fit_pgocc <- function(occ_formula, det_formula = global.det.formula) {
  PGOcc(
    occ.formula   = occ_formula,
    det.formula   = det_formula, 
    data          = list(y = y,
                         occ.covs = occ_df,
                         det.covs = X.p),
    inits         = gm.inits, 
    n.samples     = n.samples, 
    priors        = gm.priors, 
    n.omp.threads = 1, 
    verbose       = FALSE,
    n.report      = 1000, 
    n.burn        = n.burn, 
    n.thin        = n.thin, 
    n.chains      = n.chains
  )
}


# fit all occupancy models
fits_occ <- lapply(occ_formulas, fit_pgocc)
names(fits_occ) <- names(occ_formulas)
print(fits_occ[1])
lapply(fits_occ, summary)
# check output of waic_list
# waic_list[[1]]
waic_list <- lapply(fits_occ, waicOcc)
waic_tbl <- tibble(
  model = names(waic_list),
  elpd  = sapply(waic_list, function(x) x["elpd"]),
  pD    = sapply(waic_list, function(x) x["pD"]),
  WAIC  = sapply(waic_list, function(x) x["WAIC"]))

# add delta & weights
waic_tbl <- waic_tbl %>%
  arrange(WAIC) %>%
  mutate(
    dWAIC  = WAIC - min(WAIC),
    weight = exp(-0.5 * dWAIC) / sum(exp(-0.5 * dWAIC)))

waic_tbl




































# ---- Detection Models -----
det_formulas <- list(
  d1 = global.det.formula,
  d2 = ~ t_air + I(t_air^2) + wind + vpd,
  d3 = ~ t_air + wind,
  d4 = ~ t_air,
  d5 = ~ 1
)

fits_det <- lapply(det_formulas, function(det_form) {
  PGOcc(
    occ.formula = best_occ_formula,      # from WAIC step
    det.formula = det_form,
    data        = list(y = y,
                       occ.covs = occ_df,
                       det.covs = X.p),
    inits       = gm.inits, 
    n.samples   = n.samples, 
    priors      = gm.priors, 
    n.omp.threads = 1,
    verbose     = FALSE,
    n.report    = 1000,
    n.burn      = n.burn, 
    n.thin      = n.thin, 
    n.chains    = n.chains
  )
})

waic_det <- lapply(fits_det, waicOcc)
# same WAIC table pattern as above

















