# running through occupancy model options
# option 1: the unmarked frame

# load the libraries

library(dplyr)
library(tidyr)
library(stringr)
library(unmarked)

# bring in my dataset
gm <- read.csv("data/data_gmocc_all-ut.csv")

# ---- Data Prep ----
# converting date to julian
gm$date <- mdy(gm$date)
gm$julian_day <- yday(gm$date)

# converting time to calc
hms_to_seconds <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  hhmmss <- strsplit(x, ":", fixed = TRUE)[[1]]
  as.numeric(hhmmss[1]) * 3600 + as.numeric(hhmmss[2]) * 60 + as.numeric(hhmmss[3])
}

gm <- gm %>% 
  mutate(time_sec = vapply(time, hms_to_seconds, numeric(1)))

planned  <- gm %>% 
  filter(col_point %in% 1:5)

enc_only <- gm %>% 
  filter(col_point == 6 & detection == 1) %>% 
  distinct(site_id, tran_id) %>% 
  mutate(encounter_flag = 1L)

tran_reps <- planned %>%
  group_by(site_id, tran_id) %>%
  summarise(
    det_planned       = as.integer(any(detection == 1, na.rm = TRUE)),
    air_temp_mean      = mean(air_temp,       na.rm = TRUE),
    gr_temp_mean      = mean(gr_temp, na.rm = TRUE),
    prey_presence = max(pp_total, na.rm = TRUE),
    wind_mean          = mean(avg_wind_spd,   na.rm = TRUE),
    time_mean_sec      = mean(time_sec,       na.rm = TRUE),
    vpd_mean           = mean(vpd,            na.rm = TRUE),
    rockiness_mean     = mean(rock_ind,       na.rm = TRUE),
    shrub_density_mean = mean(shr_pt_density, na.rm = TRUE),
    elev_mean          = mean(elev,           na.rm = TRUE),
    n_points_tran      = n_distinct(col_point),
    tran_utm_x         = mean(utm_x,          na.rm = TRUE),
    tran_utm_y         = mean(utm_y,          na.rm = TRUE),
    .groups = "drop") %>%
  left_join(enc_only, by = c("site_id","tran_id")) %>%
  mutate(encounter_flag = ifelse(is.na(encounter_flag), 0L, encounter_flag),
         det_transect   = pmin(1L, det_planned + encounter_flag)) %>%
  group_by(site_id) %>%
  arrange(tran_id, .by_group = TRUE) %>%
  mutate(rep_idx = row_number()) %>%
  ungroup()

max_reps <- tran_reps %>% count(site_id) %>% pull(n) %>% max()

y_wide <- tran_reps %>%
  select(site_id, rep_idx, det_transect) %>%
  complete(site_id, rep_idx = 1:max_reps) %>%
  arrange(site_id, rep_idx) %>%
  pivot_wider(names_from = rep_idx, values_from = det_transect, names_prefix = "rep_")

site_order <- unique(y_wide$site_id)

y <- y_wide %>% select(-site_id) %>% as.matrix()

rownames(y) <- site_order

site_covs <- tran_reps %>%
  group_by(site_id) %>%
  summarise(
    shrub_density_site_mean = mean(shrub_density_mean, na.rm = TRUE),
    elev_site_mean          = mean(elev_mean,          na.rm = TRUE),
    n_transects_site        = n(),
    .groups = "drop") %>%
  arrange(match(site_id, site_order))


stopifnot(identical(site_covs$site_id, rownames(y)))

siteCovs <- site_covs %>% 
  select(-site_id) %>% 
  mutate(across(where(is.numeric), scale))

to_obs_mat <- function(df, value_col) {
  w <- df %>%
    select(site_id, rep_idx, {{ value_col }}) %>%
    complete(site_id, rep_idx = 1:max_reps) %>%
    arrange(site_id, rep_idx) %>%
    pivot_wider(names_from = rep_idx, values_from = {{ value_col }})
  m <- as.matrix(w %>% select(-site_id))
  rownames(m) <- w$site_id
  scale(m)
}

obsCovs <- list(
  rockiness = to_obs_mat(tran_reps, rockiness_mean),
  time   = to_obs_mat(tran_reps, time_mean_sec),
  a.temp   = to_obs_mat(tran_reps, air_temp_mean),
  g.temp = to_obs_mat(tran_reps, gr_temp_mean),
  prey = to_obs_mat(tran_reps, prey_presence),
  wind   = to_obs_mat(tran_reps, wind_mean),
  vpd    = to_obs_mat(tran_reps, vpd_mean),
  effort = to_obs_mat(tran_reps, n_points_tran))

umf <- unmarkedFrameOccu(y = y, siteCovs = as.data.frame(siteCovs), obsCovs = obsCovs)

fit <- occu(
  ~ rockiness + time + a.temp + g.temp + wind + effort + vpd + prey
  ~ shrub_density_site_mean + elev_site_mean,
  data = umf)

fit1 <- occu(~ effort ~ elev_site_mean, data = umf)
fit2 <- occu(~ 1 ~ elev_site_mean, data = umf)
fit3 <- occu(~ effort + rockiness ~ elev_site_mean, data = umf)


print(summary(fit1))
print(summary(fit2))
print(summary(fit3))

print(summary(fit))
#saveRDS(list(umf=umf, fit=fit), file = "opt1_unmarked_transects_as_reps.rds")

summary(umf)         # see how many detections per site
plot(hist(rowSums(umf@y)))  # check detection frequencies
cor(as.data.frame(umf@obsCovs))  # check for collinearity
