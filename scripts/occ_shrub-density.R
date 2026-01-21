## Calculating Shrub Density
# SSS - Fall 2025

# load libraries
library(dplyr)

# bring in dataset
data <- read.csv("data/data_gmocc_all-ut.csv")
names(data) <- trimws(names(data))

# ---- Writing Shrub Subset ----
# select just shrub data
shrub <- data %>%
  select("site_id","tran_id","date","col_point","shr_nw","shr_ne","shr_sw","shr_se") %>%
  rename(
    `1_dist` = shr_nw,
    `2_dist` = shr_ne,
    `3_dist` = shr_sw,
    `4_dist` = shr_se)

write.csv(shrub, "data/data_gmocc_shrubbery-ut.csv")

# ---- Calculating Density ----
shrub_density <- shrub %>%
  mutate(across(ends_with("dist"), ~ .x / 100)) %>%                          # convert cm to                               
  rowwise() %>%                                                              # one mean per row across all distance columns
  mutate(point_mean_dist = mean(c_across(ends_with("dist")), na.rm = TRUE),  # calculate mean distance at each point
         point_mean_area = point_mean_dist^2,                                # calculate mean area at each point                    
         point_density = ifelse(point_mean_area > 0,1 /                      # calculate density at each point
                                point_mean_area, NA_real_)) %>% 
  ungroup()
                    
# summarize by transect
shrubden_by_transect <- shrub_density %>%
  group_by(tran_id) %>%
  mutate(
    shr_m_area   = round(mean(point_mean_area, na.rm = TRUE), 3),
    shr_med_area = round(median(point_mean_area, na.rm = TRUE), 3),
    shr_m_density   = round(mean(point_density, na.rm = TRUE), 3),
    shr_med_density = round(median(point_density, na.rm = TRUE), 3))

# ---- Exporting -----
# export to new csv
write.csv(shrubden_by_transect, "data/data_gmocc_shrub-density-ut.csv")

