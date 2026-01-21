## Data Cleaning & Subsetting
# SSS - Fall 2025
library(dplyr)
library(janitor)

# load in data
data_gmocc <- read.csv("data/data_gmocc_all-ut.csv")
names(data_gmocc) <- trimws(names(data_gmocc))

# look at the data
glimpse(data_gmocc)

# ---- Checking for transect errors ----
# looking for transects that have extra/not enough collection points
errors <- 
data_gmocc %>%
  group_by(tran_id) %>%
  summarise(n_points = n_distinct(col_point)) %>%
  filter(n_points != 5)                                       # n = 6 for on transect detections / n = 1 for off transect detections

# summary of total transects
sum_table <-
  data_gmocc %>%
  count(tran_id)
sum_table

# ---- Summarizing PP & adding to master data ----
# adding the total pp columns
data_gmocc <- data_gmocc %>%
  group_by(tran_id) %>%
  mutate(pp_total = sum(pp_mam, pp_rep, pp_gbird))

# updating all_data csv
write.csv(data_gmocc, "data/data_gmocc_all-ut.csv", row.names = FALSE)

# ---- Exporting Subset CSVs ----
# create off transect detection csv
data_gmocc_off_tran <-
  data_gmocc %>%
  filter(grepl("OTD", tran_id))

write.csv(data_gmocc_off_tran, "data/data_gmocc_off-sign-ut.csv", row.names = FALSE)

# create on transect detection csv
data_gmocc_on_tran <-
  data_gmocc %>%
  filter(!grepl("OTD", tran_id) & detection == 1)

write.csv(data_gmocc_on_tran, "data/data_gmocc_on-sign-ut.csv", row.names = FALSE)
