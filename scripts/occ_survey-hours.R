# Calculating Survey Effort/Hours
# SSS - Fall 2025

# load in libraries
library(dplyr)
library(lubridate)
library(hms)
library(tidyr)

# bring in dataset
data <- read.csv("data/data_gmocc_all-ut.csv")

names(data) <- trimws(names(data))

# ---- Data Wrangling ----
# select just survey hour data
hours <- data %>%
  select(site_id, tran_id, date, obs_1, obs_2, obs_3, obs_4, time, col_point) %>%
  filter(!grepl("OTD", tran_id, ignore.case = TRUE)) %>%                                    # exclude off-transect times
  filter(col_point %in% c(1, 5))                                                            # only start/end

# make the data wider to prep for the calculations
hours <- hours %>%
  pivot_wider(
    names_from = col_point,
    values_from = time,
    names_prefix = "col_point_") %>%
  rename(start_time = col_point_1,
         end_time   = col_point_5)
  

# ---- Calculating Total Survey Hours ----
survey_hours <- hours %>%
  mutate(
    start_time = as_hms(parse_date_time(start_time, orders = c("H:M", "H:M:S"))),           # Convert to hms times
    end_time   = as_hms(parse_date_time(end_time,   orders = c("H:M", "H:M:S"))),
    duration_hr = as.numeric(end_time - start_time, units = "hours"),                       # duration in hours
    n_observers = rowSums(!is.na(select(., obs_1:obs_4)) & select(., obs_1:obs_4) != ""),   # count observers
    total_hours = duration_hr * n_observers)                                                # person-hours

# save this df as a csv
write.csv(survey_hours, "data/data_gmocc_survey-hours-ut.csv")
# ---- Summaries ----
overall_total <- survey_hours %>%
  summarise(
    total_person_hours = round(sum(total_hours, na.rm = TRUE),0),        # adjusted by num observers ; essentially total human labor hours
    total_effort_hours = round(sum(duration_hr, na.rm = TRUE),0))        # actual survey duration ; time in field regardless of how many ppl

hours_by_site <- survey_hours %>%
  group_by(site_id) %>%
  summarise(
    total_person_hours = round(sum(total_hours, na.rm = TRUE), 0),
    total_effort_hours = round(sum(duration_hr, na.rm = TRUE), 0),
    avg_observers_per_survey = mean(n_observers, na.rm = TRUE),   
    n_surveys = n(),
    .groups = "drop")

# ---- Exporting ----
write.csv(overall_total, "output/gmocc_hours_sum-table.csv")
write.csv(hours_by_site, "output/gmocc_hours_sum-table-sites.csv")

