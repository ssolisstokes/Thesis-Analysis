# ---- Libraries ----
library(dplyr)
library(rsample)
library(purrr)
library(broom)
library(tidyr)
library(ggplot2)

# ---- Data Wrangling ----
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

# sanity check
#transects %>% filter(detection == 1)

# add the transect length (path_km)
transects <- merge(transects, t.length[, c("site_id","tran_id","path_km")],
                   by = c("site_id","tran_id"), all.x = TRUE)

# mutate detection to be logical true/false
transects.bs <- transects %>%
  mutate(detection = as.logical(detection))

# sanity check
#transects.bs %>% filter(detection == TRUE)

# create df for detection covariates to run bootstrap
det.covs <- transects.bs %>%
  select(detection, rock_ind, air_temp, gr_temp, avg_wind_spd, vpd)

# ---- Bootstrap Analysis ----
# set seed
set.seed(456)

# use rsample::bootstraps
boot.resamples <- bootstraps(det.covs, times = 2000, strata = detection)

# function to run all bootstrap differences between detection/non-detection covariates
covariate.diffs <- function(split) {
  df <- analysis(split)
  
  tibble(
    rock.diff = mean(df$rock_ind[df$detection], na.rm = TRUE) - mean(df$rock_ind[!df$detection], na.rm = TRUE),
    air.diff = mean(df$air_temp[df$detection], na.rm = TRUE) - mean(df$air_temp[!df$detection], na.rm = TRUE),
    gr.diff = mean(df$gr_temp[df$detection], na.rm = TRUE) - mean(df$gr_temp[!df$detection], na.rm = TRUE),
    wind.diff = mean(df$avg_wind_spd[df$detection], na.rm = TRUE) - mean(df$avg_wind_spd[!df$detection], na.rm = TRUE),
    vpd.diff = mean(df$vpd[df$detection], na.rm = TRUE) - mean(df$vpd[!df$detection], na.rm = TRUE))
}

# apply function using purrr::map
boot.results <- boot.resamples %>%
  mutate(estimates = map(splits, covariate.diffs)) %>%
  unnest(estimates)

# check output
summary(boot.results)
print(boot.results)

# ---- Confidence Intervals ----
boot.confin_tbl <- boot.results %>%
  summarise(
    rock.min = quantile(rock.diff, 0.025, na.rm = TRUE),
    rock.med = quantile(rock.diff, 0.50, na.rm = TRUE),
    rock.max = quantile(rock.diff, 0.975, na.rm = TRUE),
    
    air.min = quantile(air.diff, 0.025, na.rm = TRUE),
    air.med = quantile(air.diff, 0.50, na.rm = TRUE),
    air.max = quantile(air.diff, 0.975, na.rm = TRUE),
    
    gr.min = quantile(gr.diff, 0.025, na.rm = TRUE),
    gr.med = quantile(gr.diff, 0.50, na.rm = TRUE),
    gr.max = quantile(gr.diff, 0.975, na.rm = TRUE),
    
    wind.min = quantile(wind.diff, 0.025, na.rm = TRUE),
    wind.med = quantile(wind.diff, 0.50, na.rm = TRUE),
    wind.max = quantile(wind.diff, 0.975, na.rm = TRUE),
    
    vpd.min = quantile(vpd.diff, 0.025, na.rm = TRUE),
    vpd.med = quantile(vpd.diff, 0.50,  na.rm = TRUE),
    vpd.max = quantile(vpd.diff, 0.975, na.rm = TRUE))

boot.confin_tbl

# reformat table 
boot.confin_tbl<- tibble(
  variable = c("rock_ind", "air_temp", "gr_temp", "avg_wind_spd", "vpd"),
  lcl = c(boot.confin_tbl$rock.min, boot.confin_tbl$air.min, boot.confin_tbl$gr.min, boot.confin_tbl$wind.min, boot.confin_tbl$vpd.min),
  med = c(boot.confin_tbl$rock.med, boot.confin_tbl$air.med, boot.confin_tbl$gr.med, boot.confin_tbl$wind.med, boot.confin_tbl$vpd.med),
  ucl = c(boot.confin_tbl$rock.max, boot.confin_tbl$air.max, boot.confin_tbl$gr.max, boot.confin_tbl$wind.max, boot.confin_tbl$vpd.max))

boot.confin_tbl

# ---- Plots ----
# function to plot bootstrap results
plot.hist <- function(df, colname, xlab) {
  
  ggplot(df, aes(x = .data[[colname]])) +
    geom_histogram(bins = 30, fill = "tan", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = paste("Bootstrap distribution:", xlab),
      x = paste0("Difference in mean ", xlab, " (detection vs non-detection)"),
      y = "Bootstrap Frequency") +
    theme_classic()
}

## ---- Histograms of Results ----
# apply function
plot.hist(boot.results, "rock.diff", "Substrate")
plot.hist(boot.results, "air.diff", "Air Temperature")
plot.hist(boot.results, "gr.diff", "Ground Temperature")
plot.hist(boot.results, "wind.diff", "Average Wind Speed")
plot.hist(boot.results, "vpd.diff", "Vapor Pressure Deficit")


# function to plot raw data
plot.box <- function(df, colname, xlab) {
  
  ggplot(df, aes(x = factor(detection), y = .data[[colname]])) +
    geom_boxplot(fill= c("grey80", "tan")) +
    labs(
      title = paste(xlab, "differences at detection vs non-detection points"),
      x = "Detection",
      y = xlab) +
    theme_classic()
}

## ---- Raw Data Plots ----
# apply function
plot.box(transects, "rock_ind", "Substrate")
plot.box(transects, "air_temp", "Air Temperature")
plot.box(transects, "gr_temp", "Ground Temperature")
plot.box(transects, "avg_wind_spd", "Average Wind Speed")
plot.box(transects, "vpd", "Vapor Pressure Deficit")

## ---- Two-sided bootstrap p-values ----
# calculate two-sided p-value for each variable
boot.pvals <- boot.results %>%
  summarise(
    rock.p = 2 * min(mean(rock.diff <= 0, na.rm = TRUE), mean(rock.diff >= 0, na.rm = TRUE)),
    air.p  = 2 * min(mean(air.diff <= 0, na.rm = TRUE), mean(air.diff >= 0, na.rm = TRUE)),
    gr.p   = 2 * min(mean(gr.diff <= 0, na.rm = TRUE), mean(gr.diff >= 0, na.rm = TRUE)),
    wind.p = 2 * min(mean(wind.diff <= 0, na.rm = TRUE), mean(wind.diff >= 0, na.rm = TRUE)),
    vpd.p  = 2 * min(mean(vpd.diff <= 0, na.rm = TRUE), mean(vpd.diff >= 0, na.rm = TRUE)))

# create summary table with CI tbl & p-values
boot.summary_tbl <- boot.confin_tbl %>%
  left_join(tibble(variable = c("rock_ind", "air_temp", "gr_temp", "avg_wind_spd", "vpd"),
      p_value = as.numeric(boot.pvals[1, ])),
    by = "variable")

boot.summary_tbl
