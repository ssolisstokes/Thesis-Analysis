## Calculating Total Transect Distances
# SSS - Fall 2025

# load the libraries
library(dplyr)
library(readr)
library(janitor)

# bring in the data 
dat <- read_csv("data/data_gmocc_all-ut.csv") %>%
  clean_names()

# ---- Calculating Distances ----
# exclude detection points from line path
dat <- dat %>%
  filter(!is.na(col_point), col_point != 6) 

# calculate transect length
transect_dist <- dat %>%
  arrange(site_id, tran_id, col_point) %>%                                        # order rows
  group_by(site_id, tran_id) %>%                                                                
  filter(n() >= 2) %>%                                                            # exclude off-transect detections 
  mutate(
    dx = utm_x - lag(utm_x),               # calculating the difference in X from the previous point
    dy = utm_y - lag(utm_y),               # calculating the difference in Y from the previous point
    seg_m = sqrt(dx^2 + dy^2)) %>%         # calculating the pythagorean distance for the line segment between points
  summarise(
    n_points   = n(),
    n_segments = sum(!is.na(seg_m)),
    path_m     = sum(seg_m, na.rm = TRUE),                                       # sum of all segment lengths (actual path)
    straight_m = sqrt((last(utm_x) - first(utm_x))^2 +
                        (last(utm_y) - first(utm_y))^2),                         # start to end straight line
    curve  = ifelse(straight_m > 0, path_m / straight_m, NA_real_),              # how curvy the path is; curviness ratio; path/straight
    .groups = "drop") %>%          
  mutate(
    path_km = path_m / 1000,                                                     # m to km convert
    straight_km = straight_m / 1000)

write_csv(transect_dist, "output/gmocc_tran-distance.csv")

# ---- Summary Tables ----
## ---- overall distance summary ----
overall_summary <- transect_dist %>%
  summarise(
    n_transects = n(),
    avg_path_km = mean(path_km, na.rm = TRUE),
    total_path_km = round(sum(path_km, na.rm = TRUE), 0),
    median_path_km = median(path_km, na.rm = TRUE),
    .groups = "drop")

overall_summary
write_csv(overall_summary, "output/gmocc_tran-dist_sum-table.csv")

## ---- by-site distance summary ----
by_site_summary <- transect_dist %>%
  group_by(site_id) %>%
  summarise(
    n_transects = n(),
    avg_path_km = mean(path_km, na.rm = TRUE),
    total_path_km = sum(path_km, na.rm = TRUE),
    median_path_km = median(path_km, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(desc(total_path_km))                            

by_site_summary
write_csv(by_site_summary, "output/gmocc_tran-dist_sum-table-sites.csv")