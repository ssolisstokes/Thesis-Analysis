library(dplyr)
library(rsample)
library(purrr)
library(broom)
library(tidyverse)
library(ggplot2)
# bring in dataset
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

set.seed(456)

# Make a clean analysis df with explicit detected flag and factor for stratification
transects_bs <- transects %>%
  mutate(
    detected = detection == 1,                 # logical TRUE/FALSE
    detection_f = factor(detection, levels = c(0, 1)))

boot.resamples <- bootstraps(transects_bs, times = 2000, strata = detection_f)

# statistic: mean(rock_ind | detected) - mean(rock_ind | not detected)
#This function takes a 'split' object, extracts the analysis data,
# and calculates the mean of the 'rockiness' column
rockiness.means <- function(split) {
  df <- analysis(split)
  
  detections_m    <- mean(df$rock_ind[df$detected], na.rm = TRUE)
  nondetections_m <- mean(df$rock_ind[!df$detected], na.rm = TRUE)
  
  tibble(rock_diff = detections_m - nondetections_m)
}

# 3. Apply the function to all resamples using purrr::map
boot.results <- boot.resamples %>%
  mutate(estimates = map(splits, rockiness.means)) %>%
  unnest(estimates)

print(boot.results)

# Optional: summarize + CI
orig_diff <- with(transects_bs,
                  mean(rock_ind[detected], na.rm = TRUE) -
                    mean(rock_ind[!detected], na.rm = TRUE))

ci <- quantile(boot.results$rock_diff, c(0.025, 0.5, 0.975), na.rm = TRUE)

orig_diff
ci

ggplot(boot.results, aes(x = rock_diff)) +
  geom_histogram(bins = 30, fill = "tan", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Difference in mean sandiness (detected − not detected)",
    y = "Bootstrap frequency",
    title = "Bootstrap distribution of substrate effect on detection"
  ) +
  theme_classic()

  
ggplot(transects_bs, aes(x = factor(detected), y = rock_ind)) +
  geom_boxplot(fill = c("grey80", "tan")) +
  labs(
    x = "Track detected",
    y = "Sandiness (rock_ind)",
    title = "Substrate differences at detection vs non-detection points"
  ) +
  theme_classic()

ci_df <- tibble(
  ymin = quantile(boot.results$rock_diff, 0.025),
  ymed = quantile(boot.results$rock_diff, 0.5),
  ymax = quantile(boot.results$rock_diff, 0.975))

ci_df

#
ggplot(boot.results, aes(x = rock_diff)) +
  geom_histogram(
    bins = 30,
    fill = "#D2B48C",
    color = "black",
    alpha = 0.9 ) +
  geom_vline(
    xintercept = obs_diff,
    linetype = "dashed",
    linewidth = 1 ) +
  geom_vline(
    xintercept = c(ci_df$ymin, ci_df$ymax),
    linetype = "dotted",
    linewidth = 0.8) +
  annotate(
    "text",
    x = obs_diff,
    y = Inf,
    label = paste0(
      "Observed mean = ", round(obs_diff, 2),
      "\n95% CI: [", round(ci_df$ymin, 2), ", ", round(ci_df$ymax, 2), "]"),  
    vjust = 1.5,
    hjust = -0.05,
    size = 4) +
  labs(
    x = "Difference in mean rockiness (detected – not detected)",
    y = "Bootstrap frequency",
    title = "Bootstrap distribution of substrate effect on track detection",
    subtitle = "Positive values indicate sandier substrate at detection points"  ) +
  theme_minimal()


