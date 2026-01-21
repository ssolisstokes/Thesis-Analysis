## Variable Density Plots
# SSS - Fall 2025

# load in libraries
library(tidyverse)
library(janitor)
library(ggplot2)
library(patchwork)
library(fs)

# bring in data set
data_gmocc <- read.csv("data/data_gmocc_all-ut.csv")
names(data_gmocc) <- trimws(names(data_gmocc))

# create output folder
dir_create("images/density-plots/all-data")

# columns NOT to include in all-vars density facet
no_den <- c("site_id","tran_id","date","col_point","utm_x","utm_y",
            "detection","shr_ne","shr_nw","shr_se","shr_sw",
           "shr_m_area", "pp_mam", "pp_gbird", "pp_rep", "med_vpd", "shr_m_density")

# variable list (numeric only, excluding above and prey presence)
density_vars <- data_gmocc %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()
# --------- All Sites -------
## ---- All Variables ----
# create a pivot longer table to plot all variables at once
data_gmocc.long <- data_gmocc |>
  pivot_longer(
    cols = all_of(density_vars),
    names_to = "variable",
    values_to = "value")

# create the plot
plot_all.vars <-
  ggplot(data_gmocc.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +       # facet wrap is how we see all at once
  labs(title = "Density Plots of All Variables",
       x = "Value", y = "Density") +
  theme_minimal()

plot_all.vars
ggsave("images/density-plots/all-data/density_all-vars.png", plot = plot_all.vars, width = 8, height = 10, dpi = 300)

## ---- Individual Variables ----

### ---- Rockiness ----
plot_rock <-
  ggplot(data_gmocc, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Rockiness",
    x = "Rockiness Index",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density-rockiness.png", width = 8, height = 10, dpi = 300)

# show densities at detections/non-detections
plot_rock.det <-
ggplot(data_gmocc, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Rockiness by Detection Status",
    x = "Rockiness Index Values",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_rockiness-det.png", width = 8, height = 10, dpi = 300)

### ---- Air Temperature ----
plot_air.temp <-
  ggplot(data_gmocc, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Air Temp (C)",
    x = "Air Temp (C)",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_air-temp.png", width = 8, height = 10, dpi = 300)

# show densities at detections/non-detections
plot_air.temp.det <-
  ggplot(data_gmocc, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Air Temp by Detection Status",
    x = "Air Temp (C)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_air-temp-det.png", width = 8, height = 10, dpi = 300)

### ---- Ground Temperature ----
plot_gr.temp <- 
  ggplot(data_gmocc, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Ground Temp (C)",
    x = "Ground Temp (C)",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_gr-temp.png", width = 8, height = 10, dpi = 300)

# show densities at detections/non-detections
plot_gr.temp.det <-
  ggplot(data_gmocc, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Ground Temp by Detection Status",
    x = "Ground Temp (C)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_gr-temp-det.png", width = 8, height = 10, dpi = 300)

### ---- Average Wind Speed ----
plot_aws <- 
ggplot(data_gmocc, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Avg. Wind Speed (kmh)",
    x = "Average Wind Speed (kmh)",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_aws.png", width = 8, height = 10, dpi = 300)

# show densities at detections/non-detections
plot_aws.det <-
  ggplot(data_gmocc, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Avg. Wind Speed by Detection Status",
    x = "Average Wind Speed (kmh)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_aws-det.png", width = 8, height = 10, dpi = 300)

### ---- Relative Humidity ----
plot_rel.hum <- 
ggplot(data_gmocc, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Relative Humidity (%)",
    x = "Relative Humidity (%)",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_rel-hum.png", width = 8, height = 10, dpi = 300)

# show densities at detections/non-detections
plot_rel.hum.det <-
  ggplot(data_gmocc, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Relative Humidity by Detection Status",
    x = "Relative Humidity (%)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_rel-hum-det.png", width = 8, height = 10, dpi = 300)

### ---- Elevation ----
plot_elev <- 
ggplot(data_gmocc, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Elevation",
    x = "Elevation (m)",
    y= "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_elev.png")

# show densities at detections/non-detections
plot_elev.det <-
  ggplot(data_gmocc, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Elevation by Detection Status",
    x = "Elevation (m)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_elev-det.png", width = 8, height = 10, dpi = 300)

### ---- Shrub density ----
plot_shrub <-
  ggplot(data_gmocc, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_shrub.png",
        plot = plot_shrub, width = 6, height = 4, dpi = 300)

plot_shrub.det <-
  ggplot(data_gmocc, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_shrub-det.png",
      plot = plot_shrub.det, width = 6, height = 4, dpi = 300)

### ---- Prey Presence ----
plot_pp <-
  ggplot(data_gmocc, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_pp.png",
       plot = plot_pp, width = 6, height = 4, dpi = 300)

plot_pp.det <-
  ggplot(data_gmocc, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_pp-det.png",
       plot = plot_pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
plot_vpd <-
  ggplot(data_gmocc, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_vpd.png",
       plot = plot_vpd, width = 6, height = 4, dpi = 300)

plot_vpd.det <-
  ggplot(data_gmocc, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/all-data/density_vpd-det.png",
       plot = plot_vpd.det, width = 6, height = 4, dpi = 300)

# --------- Site Specific -------
## ----- Cove Wash ------
# select site
cowa <- data_gmocc %>%
  filter(site_id == "cove_wash")

# create output folder
dir_create("images/density-plots/cove-wash")

# variable list (numeric only, excluding above and prey presence)
cowa.density_vars <- cowa %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()

### ---- All Variables ----
# long table for all-vars plot
cowa.long <- cowa %>%
  pivot_longer(
    cols = all_of(cowa.density_vars),
    names_to = "variable",
    values_to = "value")

# ALL variables density (facetted)
cowa.plot_all.vars <-
  ggplot(cowa.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Density Plots of All Variables — Cove Wash",
       x = "Value", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_all-vars.png",
       plot = cowa.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness ----
cowa.plot_rock <-
  ggplot(cowa, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Rockiness — Cove Wash",
       x = "Rockiness Index", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_rockiness.png",
       plot = cowa.plot_rock, width = 6, height = 4, dpi = 300)

cowa.plot_rock.det <-
  ggplot(cowa, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Rockiness by Detection Status — Cove Wash",
       x = "Rockiness Index", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_rockiness-det.png",
       plot = cowa.plot_rock.det, width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
cowa.plot_air.temp <-
  ggplot(cowa, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Air Temp (°C) — Cove Wash",
       x = "Air Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_air-temp.png",
       plot = cowa.plot_air.temp, width = 6, height = 4, dpi = 300)

cowa.plot_air.temp.det <-
  ggplot(cowa, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Air Temp by Detection Status — Cove Wash",
       x = "Air Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_air-temp-det.png",
       plot = cowa.plot_air.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Ground Temperature ----
cowa.plot_gr.temp <-
  ggplot(cowa, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Ground Temp (°C) — Cove Wash",
       x = "Ground Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_gr-temp.png",
       plot = cowa.plot_gr.temp, width = 6, height = 4, dpi = 300)

cowa.plot_gr.temp.det <-
  ggplot(cowa, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Ground Temp by Detection Status — Cove Wash",
       x = "Ground Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_gr-temp-det.png",
       plot = cowa.plot_gr.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Average Wind Speed ----
cowa.plot_aws <-
  ggplot(cowa, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Avg. Wind Speed (km/h) — Cove Wash",
       x = "Average Wind Speed (km/h)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_aws.png",
       plot = cowa.plot_aws, width = 6, height = 4, dpi = 300)

cowa.plot_aws.det <-
  ggplot(cowa, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Avg. Wind Speed by Detection Status — Cove Wash",
       x = "Average Wind Speed (km/h)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_aws-det.png",
       plot = cowa.plot_aws.det, width = 6, height = 4, dpi = 300)

#### ---- Relative Humidity ----
cowa.plot_rel.hum <-
  ggplot(cowa, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Relative Humidity (%) — Cove Wash",
       x = "Relative Humidity (%)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_rel-hum.png",
       plot = cowa.plot_rel.hum, width = 6, height = 4, dpi = 300)

cowa.plot_rel.hum.det <-
  ggplot(cowa, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Relative Humidity by Detection Status — Cove Wash",
       x = "Relative Humidity (%)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_rel-hum-det.png",
       plot = cowa.plot_rel.hum.det, width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
cowa.plot_elev <-
  ggplot(cowa, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Elevation — Cove Wash",
       x = "Elevation (m)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_elev.png",
       plot = cowa.plot_elev, width = 6, height = 4, dpi = 300)

cowa.plot_elev.det <-
  ggplot(cowa, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Elevation by Detection Status — Cove Wash",
       x = "Elevation (m)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_elev-det.png",
       plot = cowa.plot_elev.det, width = 6, height = 4, dpi = 300)

#### ---- Shrub Density ----
cowa.plot_shrub <-
  ggplot(cowa, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density — Cove Wash",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_shrub.png",
       plot = cowa.plot_shrub, width = 6, height = 4, dpi = 300)

cowa.plot_shrub.det <-
  ggplot(cowa, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status — Cove Wash",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_shrub-det.png",
       plot = cowa.plot_shrub.det, width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
cw.plot_pp <-
  ggplot(cowa, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - Cove Wash",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_pp.png",
       plot = cw.plot_pp, width = 6, height = 4, dpi = 300)

cw.pp.det <-
  ggplot(cowa, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - Cove Wash",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_pp-det.png",
      plot = cw.pp.det, width = 6, height = 4, dpi = 300)

#### ---- Vapor Pressure Deficit ----
cw.plot_vpd <-
  ggplot(cowa, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - Cove Wash",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_vpd.png",
       plot = cw.plot_vpd, width = 6, height = 4, dpi = 300)

cw.plot_vpd.det <-
  ggplot(cowa, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - Cove Wash",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/cove-wash/cw-density_vpd-det.png",
       plot = cw.plot_vpd.det, width = 6, height = 4, dpi = 300)

## ---- White Reef ----
# select site
whre <- data_gmocc %>%
  filter(site_id == "white_reef")

# create output folder
dir_create("images/density-plots/white-reef")

# variable list (numeric only, excluding above and prey presence)
whre.density_vars <- whre %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()

### ---- All Variables ----
# long table for all-vars plot
whre.long <- whre %>%
  pivot_longer(
    cols = all_of(whre.density_vars),
    names_to = "variable",
    values_to = "value")

whre.plot_all.vars <-
  ggplot(whre.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Density Plots of All Variables — White Reef",
       x = "Value", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_all-vars.png",
       plot = whre.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness ----
whre.plot_rock <-
  ggplot(whre, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Rockiness — White Reef",
       x = "Rockiness Index", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_rockiness.png",
       plot = whre.plot_rock, width = 6, height = 4, dpi = 300)

whre.plot_rock.det <-
  ggplot(whre, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Rockiness by Detection Status — White Reef",
       x = "Rockiness Index", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_rockiness-det.png",
       plot = whre.plot_rock.det, width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
whre.plot_air.temp <-
  ggplot(whre, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Air Temp (°C) — White Reef",
       x = "Air Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_air-temp.png",
       plot = whre.plot_air.temp, width = 6, height = 4, dpi = 300)

whre.plot_air.temp.det <-
  ggplot(whre, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Air Temp by Detection Status — White Reef",
       x = "Air Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_air-temp-det.png",
       plot = whre.plot_air.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Ground temperature ----
whre.plot_gr.temp <-
  ggplot(whre, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Ground Temp (°C) — White Reef",
       x = "Ground Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_gr-temp.png",
       plot = whre.plot_gr.temp, width = 6, height = 4, dpi = 300)

whre.plot_gr.temp.det <-
  ggplot(whre, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Ground Temp by Detection Status — White Reef",
       x = "Ground Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_gr-temp-det.png",
       plot = whre.plot_gr.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Average wind speed ----
whre.plot_aws <-
  ggplot(whre, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Avg. Wind Speed (km/h) — White Reef",
       x = "Average Wind Speed (km/h)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_aws.png",
       plot = whre.plot_aws, width = 6, height = 4, dpi = 300)

whre.plot_aws.det <-
  ggplot(whre, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Avg. Wind Speed by Detection Status — White Reef",
       x = "Average Wind Speed (km/h)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_aws-det.png",
       plot = whre.plot_aws.det, width = 6, height = 4, dpi = 300)

#### ---- Relative humidity ----
whre.plot_rel.hum <-
  ggplot(whre, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Relative Humidity (%) — White Reef",
       x = "Relative Humidity (%)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_rel-hum.png",
       plot = whre.plot_rel.hum, width = 6, height = 4, dpi = 300)

whre.plot_rel.hum.det <-
  ggplot(whre, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Relative Humidity by Detection Status — White Reef",
       x = "Relative Humidity (%)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_rel-hum-det.png",
       plot = whre.plot_rel.hum.det, width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
whre.plot_elev <-
  ggplot(whre, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Elevation — White Reef",
       x = "Elevation (m)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_elev.png",
       plot = whre.plot_elev, width = 6, height = 4, dpi = 300)

whre.plot_elev.det <-
  ggplot(whre, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Elevation by Detection Status — White Reef",
       x = "Elevation (m)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_elev-det.png",
       plot = whre.plot_elev.det, width = 6, height = 4, dpi = 300)

#### ---- Shrub density ----
whre.plot_shrub <-
  ggplot(whre, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density — White Reef",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_shrub.png",
       plot = whre.plot_shrub, width = 6, height = 4, dpi = 300)

whre.plot_shrub.det <-
  ggplot(whre, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status — White Reef",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_shrub-det.png",
       plot = whre.plot_shrub.det, width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
wr.plot_pp <-
  ggplot(whre, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - White Reef",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_pp.png",
       plot = wr.plot_pp, width = 6, height = 4, dpi = 300)

wr.pp.det <-
  ggplot(whre, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - White Reef",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_pp-det.png",
       plot = wr.pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
wr.plot_vpd <-
  ggplot(whre, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - White Reef",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_vpd.png",
       plot = wr.plot_vpd, width = 6, height = 4, dpi = 300)

wr.plot_vpd.det <-
  ggplot(whre, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - White Reef",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/white-reef/wr-density_vpd-det.png",
       plot = wr.plot_vpd.det, width = 6, height = 4, dpi = 300)

## ----- Turkey Farm  ------
# selecting site
tufa <- data_gmocc %>%
  filter(site_id == "turkey_farm")

# create output folder
dir_create("images/density-plots/turkey-farm")

# create the df for the variables
tufa.density_vars <- tufa %>%
  select(where(is.numeric) & 
           !all_of(no_den)) %>% 
  names()

### ---- All Variables ----
# create a pivot longer table to plot all variables at once
tufa.long <- tufa %>%
  pivot_longer(
    cols = all_of(tufa.density_vars),
    names_to = "variable",
    values_to = "value")

# create the plot
tufa.plot_all.vars <-
  ggplot(tufa.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +      
  labs(title = "Density Plots of All Variables - Turkey Farm",
       x = "Value", y = "Density") +
  theme_minimal()

tufa.plot_all.vars
ggsave("images/density-plots/turkey-farm/tf-density_all-vars.png", plot = tufa.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness----
tufa.plot_rock <-
  ggplot(tufa, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Rockiness - Turkey Farm",
    x = "Rockiness Index",
    y= "Density") +
  theme_minimal()

tufa.plot_rock
ggsave("images/density-plots/turkey-farm/tf-density_rockiness.png", width = 6, height = 4, dpi = 300)

# show densities at detections/non-detections
tufa.plot_rock.det <-
  ggplot(tufa, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Rockiness by Detection Status - Turkey Farm",
    x = "Rockiness Index Values",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_rock.det
ggsave("images/density-plots/turkey-farm/tf-density_rockiness-det.png", width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
tufa.plot_air.temp <-
  ggplot(tufa, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Air Temp (C) - Turkey Farm",
    x = "Air Temp (C)",
    y= "Density") +
  theme_minimal()

tufa.plot_air.temp
ggsave("images/density-plots/turkey-farm/tf-density_air-temp.png", width = 6, height = 4, dpi = 300)

# show densities at detections/non-detections
tufa.plot_air.temp.det <-
  ggplot(tufa, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Air Temp by Detection Status - Turkey Farm",
    x = "Air Temp (C)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_air.temp.det
ggsave("images/density-plots/turkey-farm/tf-density_air-temp-det.png", width = 6, height = 4, dpi = 300)

#### ---- Ground temperature ----
tufa.plot_gr.temp <- 
  ggplot(tufa, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Ground Temp (C) - Turkey Farm",
    x = "Ground Temp (C)",
    y= "Density") +
  theme_minimal()

tufa.plot_gr.temp
ggsave("images/density-plots/turkey-farm/tf-density_gr-temp.png", width = 6, height = 4, dpi = 300)

# show densities at detections/non-detections
tufa.plot_gr.temp.det <-
  ggplot(tufa, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Ground Temp by Detection Status - Turkey Farm",
    x = "Ground Temp (C)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_gr.temp.det
ggsave("images/density-plots/turkey-farm/tf-density_gr-temp-det.png", width = 6, height = 4, dpi = 300)

#### ---- Average wind speed ----
tufa.plot_aws <- 
  ggplot(tufa, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Avg. Wind Speed (kmh) - Turkey Farm",
    x = "Average Wind Speed (kmh)",
    y= "Density") +
  theme_minimal()

tufa.plot_aws
ggsave("images/density-plots/turkey-farm/tf-density_aws.png", width = 6, height = 4, dpi = 300) 

# show densities at detections/non-detections
tufa.plot_aws.det <-
  ggplot(tufa, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Avg. Wind Speed by Detection Status - Turkey Farm",
    x = "Average Wind Speed (kmh)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_aws.det
ggsave("images/density-plots/turkey-farm/tf-density_aws-det.png", width = 6, height = 4, dpi = 300)

#### ---- Relative humidity ----
tufa.plot_rel.hum <- 
  ggplot(tufa, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Relative Humidity (%) - Turkey Farm",
    x = "Relative Humidity (%)",
    y= "Density") +
  theme_minimal()

tufa.plot_rel.hum
ggsave("images/density-plots/turkey-farm/tf-density_rel-hum.png", width = 6, height = 4, dpi = 300) 

# show densities at detections/non-detections
tufa.plot_rel.hum.det <-
  ggplot(tufa, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Relative Humidity by Detection Status - Turkey Farm",
    x = "Relative Humidity (%)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_rel.hum.det
ggsave("images/density-plots/turkey-farm/tf-density_rel-hum-det.png", width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
tufa.plot_elev <- 
  ggplot(tufa, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Elevation - Turkey Farm",
    x = "Elevation (m)",
    y= "Density") +
  theme_minimal()

tufa.plot_elev
ggsave("images/density-plots/turkey-farm/tf-density_elev.png", width = 6, height = 4, dpi = 300)

# show densities at detections/non-detections
tufa.plot_elev.det <-
  ggplot(tufa, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Elevation by Detection Status - Turkey Farm",
    x = "Elevation (m)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_elev.det
ggsave("images/density-plots/turkey-farm/tf-density_elev-det.png", width = 6, height = 4, dpi = 300)

#### ---- Shrub density ----
tufa.plot_shrub <- 
  ggplot(tufa, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(
    title = "Density Plot of Shrub Density - Turkey Farm",
    x = "Shrub Density (m^2)",
    y= "Density") +
  theme_minimal()

tufa.plot_shrub
ggsave("images/density-plots/turkey-farm/tf-density_shrub.png", width = 6, height = 4, dpi = 300)

# show densities at detections/non-detections
tufa.plot_shrub.det <-
  ggplot(tufa, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection", "Detection")) +
  labs(
    title = "Shrub Density by Detection Status - Turkey Farm",
    x = "Shrub Density (m^2)",
    y = "Density",
    fill = "Detection") +
  theme_minimal()

tufa.plot_shrub.det
ggsave("images/density-plots/turkey-farm/tf-density_shrub-det.png", width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
tf.plot_pp <-
  ggplot(tufa, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - Turkey Farm",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/turkey-farm/tf-density_pp.png",
       plot = tf.plot_pp, width = 6, height = 4, dpi = 300)

tf.pp.det <-
  ggplot(tufa, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - Turkey Farm",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/turkey-farm/tf-density_pp-det.png",
       plot = tf.pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
tf.plot_vpd <-
  ggplot(tufa, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - Turkey Farm",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/turkey-farm/tf-density_vpd.png",
       plot = tf.plot_vpd, width = 6, height = 4, dpi = 300)

tf.plot_vpd.det <-
  ggplot(tufa, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - Turkey Farm",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/turkey-farm/tf-density_vpd-det.png",
       plot = tf.plot_vpd.det, width = 6, height = 4, dpi = 300)

## ----- Sun River  ------
# select site
suri <- data_gmocc %>%
  filter(site_id == "sun_river")

# create output folder
dir_create("images/density-plots/sun-river")

# variable list (numeric only, excluding above and prey presence)
suri.density_vars <- suri %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()
### ---- All Variables ----
# long table for all-vars plot
suri.long <- suri %>%
  pivot_longer(
    cols = all_of(suri.density_vars),
    names_to = "variable",
    values_to = "value")

# ALL variables density (facetted)
suri.plot_all.vars <-
  ggplot(suri.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Density Plots of All Variables — Sun River",
       x = "Value", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_all-vars.png",
       plot = suri.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness ----
suri.plot_rock <-
  ggplot(suri, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Rockiness — Sun River",
       x = "Rockiness Index", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_rockiness.png",
       plot = suri.plot_rock, width = 6, height = 4, dpi = 300)

suri.plot_rock.det <-
  ggplot(suri, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Rockiness by Detection Status — Sun River",
       x = "Rockiness Index", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_rockiness-det.png",
       plot = suri.plot_rock.det, width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
suri.plot_air.temp <-
  ggplot(suri, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Air Temp (°C) — Sun River",
       x = "Air Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_air-temp.png",
       plot = suri.plot_air.temp, width = 6, height = 4, dpi = 300)

suri.plot_air.temp.det <-
  ggplot(suri, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Air Temp by Detection Status — Sun River",
       x = "Air Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_air-temp-det.png",
       plot = suri.plot_air.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Ground temperature ----
suri.plot_gr.temp <-
  ggplot(suri, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Ground Temp (°C) — Sun River",
       x = "Ground Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_gr-temp.png",
       plot = suri.plot_gr.temp, width = 6, height = 4, dpi = 300)

suri.plot_gr.temp.det <-
  ggplot(suri, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Ground Temp by Detection Status — Sun River",
       x = "Ground Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_gr-temp-det.png",
       plot = suri.plot_gr.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Average wind speed ----
suri.plot_aws <-
  ggplot(suri, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Avg. Wind Speed (km/h) — Sun River",
       x = "Average Wind Speed (km/h)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_aws.png",
       plot = suri.plot_aws, width = 6, height = 4, dpi = 300)

suri.plot_aws.det <-
  ggplot(suri, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Avg. Wind Speed by Detection Status — Sun River",
       x = "Average Wind Speed (km/h)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_aws-det.png",
       plot = suri.plot_aws.det, width = 6, height = 4, dpi = 300)

#### ---- Relative humidity ----
suri.plot_rel.hum <-
  ggplot(suri, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Relative Humidity (%) — Sun River",
       x = "Relative Humidity (%)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_rel-hum.png",
       plot = suri.plot_rel.hum, width = 6, height = 4, dpi = 300)

suri.plot_rel.hum.det <-
  ggplot(suri, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Relative Humidity by Detection Status — Sun River",
       x = "Relative Humidity (%)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_rel-hum-det.png",
       plot = suri.plot_rel.hum.det, width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
suri.plot_elev <-
  ggplot(suri, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Elevation — Sun River",
       x = "Elevation (m)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_elev.png",
       plot = suri.plot_elev, width = 6, height = 4, dpi = 300)

suri.plot_elev.det <-
  ggplot(suri, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Elevation by Detection Status — Sun River",
       x = "Elevation (m)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_elev-det.png",
       plot = suri.plot_elev.det, width = 6, height = 4, dpi = 300)

#### ---- Shrub density ----
suri.plot_shrub <-
  ggplot(suri, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density — Sun River",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_shrub.png",
       plot = suri.plot_shrub, width = 6, height = 4, dpi = 300)

suri.plot_shrub.det <-
  ggplot(suri, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status — Sun River",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_shrub-det.png",
       plot = suri.plot_shrub.det, width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
sr.plot_pp <-
  ggplot(suri, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - Sun River",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_pp.png",
       plot = sr.plot_pp, width = 6, height = 4, dpi = 300)

sr.pp.det <-
  ggplot(suri, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - Sun River",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_pp-det.png",
       plot = sr.pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
sr.plot_vpd <-
  ggplot(suri, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - Sun River",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_vpd.png",
       plot = sr.plot_vpd, width = 6, height = 4, dpi = 300)

sr.plot_vpd.det <-
  ggplot(suri, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - Sun River",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/sun-river/sr-density_vpd-det.png",
       plot = sr.plot_vpd.det, width = 6, height = 4, dpi = 300)

## ----- Paradise Canyon  ------
# select site
paca <- data_gmocc %>%
  filter(site_id == "paradise_canyon")

# create output folder
dir_create("images/density-plots/paradise-canyon")

# variable list (numeric only, excluding above and prey presence)
paca.density_vars <- paca %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()
### ---- All Variables ----
# long table for all-vars plot
paca.long <- paca %>%
  pivot_longer(
    cols = all_of(paca.density_vars),
    names_to = "variable",
    values_to = "value")

# ALL variables density (facetted)
paca.plot_all.vars <-
  ggplot(paca.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Density Plots of All Variables — Paradise Canyon",
       x = "Value", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_all-vars.png",
       plot = paca.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness ----
paca.plot_rock <-
  ggplot(paca, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Rockiness — Paradise Canyon",
       x = "Rockiness Index", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_rockiness.png",
       plot = paca.plot_rock, width = 6, height = 4, dpi = 300)

paca.plot_rock.det <-
  ggplot(paca, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Rockiness by Detection Status — Paradise Canyon",
       x = "Rockiness Index", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_rockiness-det.png",
       plot = paca.plot_rock.det, width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
paca.plot_air.temp <-
  ggplot(paca, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Air Temp (°C) — Paradise Canyon",
       x = "Air Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_air-temp.png",
       plot = paca.plot_air.temp, width = 6, height = 4, dpi = 300)

paca.plot_air.temp.det <-
  ggplot(paca, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Air Temp by Detection Status — Paradise Canyon",
       x = "Air Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_air-temp-det.png",
       plot = paca.plot_air.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Ground temperature ----
paca.plot_gr.temp <-
  ggplot(paca, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Ground Temp (°C) — Paradise Canyon",
       x = "Ground Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_gr-temp.png",
       plot = paca.plot_gr.temp, width = 6, height = 4, dpi = 300)

paca.plot_gr.temp.det <-
  ggplot(paca, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Ground Temp by Detection Status — Paradise Canyon",
       x = "Ground Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_gr-temp-det.png",
       plot = paca.plot_gr.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Average wind speed ----
paca.plot_aws <-
  ggplot(paca, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Avg. Wind Speed (km/h) — Paradise Canyon",
       x = "Average Wind Speed (km/h)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_aws.png",
       plot = paca.plot_aws, width = 6, height = 4, dpi = 300)

paca.plot_aws.det <-
  ggplot(paca, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Avg. Wind Speed by Detection Status — Paradise Canyon",
       x = "Average Wind Speed (km/h)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_aws-det.png",
       plot = paca.plot_aws.det, width = 6, height = 4, dpi = 300)

#### ---- Relative humidity ----
paca.plot_rel.hum <-
  ggplot(paca, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Relative Humidity (%) — Paradise Canyon",
       x = "Relative Humidity (%)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_rel-hum.png",
       plot = paca.plot_rel.hum, width = 6, height = 4, dpi = 300)

paca.plot_rel.hum.det <-
  ggplot(paca, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Relative Humidity by Detection Status — Paradise Canyon",
       x = "Relative Humidity (%)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_rel-hum-det.png",
       plot = paca.plot_rel.hum.det, width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
paca.plot_elev <-
  ggplot(paca, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Elevation — Paradise Canyon",
       x = "Elevation (m)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_elev.png",
       plot = paca.plot_elev, width = 6, height = 4, dpi = 300)

paca.plot_elev.det <-
  ggplot(paca, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Elevation by Detection Status — Paradise Canyon",
       x = "Elevation (m)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_elev-det.png",
       plot = paca.plot_elev.det, width = 6, height = 4, dpi = 300)

#### ---- Shrub density ----
paca.plot_shrub <-
  ggplot(paca, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density — Paradise Canyon",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_shrub.png",
       plot = paca.plot_shrub, width = 6, height = 4, dpi = 300)

paca.plot_shrub.det <-
  ggplot(paca, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status — Paradise Canyon",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_shrub-det.png",
       plot = paca.plot_shrub.det, width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
pc.plot_pp <-
  ggplot(paca, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - Paradise Canyon",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_pp.png",
       plot = pc.plot_pp, width = 6, height = 4, dpi = 300)

pc.pp.det <-
  ggplot(paca, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - Paradise Canyon",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_pp-det.png",
       plot = pc.pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
pc.plot_vpd <-
  ggplot(paca, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - Paradise Canyon",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_vpd.png",
       plot = pc.plot_vpd, width = 6, height = 4, dpi = 300)

pc.plot_vpd.det <-
  ggplot(paca, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - Paradise Canyon",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/paradise-canyon/pc-density_vpd-det.png",
       plot = pc.plot_vpd.det, width = 6, height = 4, dpi = 300)

## ----- Fort Pearce ------
# select site
fope <- data_gmocc %>%
  filter(site_id == "fort_pearce")

# create output folder
dir_create("images/density-plots/fort-pearce")

# variable list (numeric only, excluding above and prey presence)
fope.density_vars <- fope %>%
  select(where(is.numeric) &
           !all_of(no_den)) %>%
  names()

### ---- All Variables ----
# long table for all-vars plot
fope.long <- fope %>%
  pivot_longer(
    cols = all_of(fope.density_vars),
    names_to = "variable",
    values_to = "value")

# ALL variables density (facetted)
fope.plot_all.vars <-
  ggplot(fope.long, aes(x = value)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(title = "Density Plots of All Variables — Fort Pearce",
       x = "Value", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_all-vars.png",
       plot = fope.plot_all.vars, width = 8, height = 10, dpi = 300)

### ---- Individual Variables ----
#### ---- Rockiness ----
fope.plot_rock <-
  ggplot(fope, aes(x = rock_ind)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Rockiness — Fort Pearce",
       x = "Rockiness Index", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_rockiness.png",
       plot = fope.plot_rock, width = 6, height = 4, dpi = 300)

fope.plot_rock.det <-
  ggplot(fope, aes(x = rock_ind, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Rockiness by Detection Status — Fort Pearce",
       x = "Rockiness Index", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_rockiness-det.png",
       plot = fope.plot_rock.det, width = 6, height = 4, dpi = 300)

#### ---- Air temperature ----
fope.plot_air.temp <-
  ggplot(fope, aes(x = air_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Air Temp (°C) — Fort Pearce",
       x = "Air Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_air-temp.png",
       plot = fope.plot_air.temp, width = 6, height = 4, dpi = 300)

fope.plot_air.temp.det <-
  ggplot(fope, aes(x = air_temp, fill = factor(detection))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Air Temp by Detection Status — Fort Pearce",
       x = "Air Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_air-temp-det.png",
       plot = fope.plot_air.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Ground temperature ----
fope.plot_gr.temp <-
  ggplot(fope, aes(x = gr_temp)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Ground Temp (°C) — Fort Pearce",
       x = "Ground Temp (°C)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_gr-temp.png",
       plot = fope.plot_gr.temp, width = 6, height = 4, dpi = 300)

fope.plot_gr.temp.det <-
  ggplot(fope, aes(x = gr_temp, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Ground Temp by Detection Status — Fort Pearce",
       x = "Ground Temp (°C)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_gr-temp-det.png",
       plot = fope.plot_gr.temp.det, width = 6, height = 4, dpi = 300)

#### ---- Average wind speed ----
fope.plot_aws <-
  ggplot(fope, aes(x = avg_wind_spd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Avg. Wind Speed (km/h) — Fort Pearce",
       x = "Average Wind Speed (km/h)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_aws.png",
       plot = fope.plot_aws, width = 6, height = 4, dpi = 300)

fope.plot_aws.det <-
  ggplot(fope, aes(x = avg_wind_spd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Avg. Wind Speed by Detection Status — Fort Pearce",
       x = "Average Wind Speed (km/h)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_aws-det.png",
       plot = fope.plot_aws.det, width = 6, height = 4, dpi = 300)

#### ---- Relative humidity ----
fope.plot_rel.hum <-
  ggplot(fope, aes(x = rel_hum)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Relative Humidity (%) — Fort Pearce",
       x = "Relative Humidity (%)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_rel-hum.png",
       plot = fope.plot_rel.hum, width = 6, height = 4, dpi = 300)

fope.plot_rel.hum.det <-
  ggplot(fope, aes(x = rel_hum, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Relative Humidity by Detection Status — Fort Pearce",
       x = "Relative Humidity (%)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_rel-hum-det.png",
       plot = fope.plot_rel.hum.det, width = 6, height = 4, dpi = 300)

#### ---- Elevation ----
fope.plot_elev <-
  ggplot(fope, aes(x = elev)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Elevation — Fort Pearce",
       x = "Elevation (m)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_elev.png",
       plot = fope.plot_elev, width = 6, height = 4, dpi = 300)

fope.plot_elev.det <-
  ggplot(fope, aes(x = elev, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Elevation by Detection Status — Fort Pearce",
       x = "Elevation (m)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_elev-det.png",
       plot = fope.plot_elev.det, width = 6, height = 4, dpi = 300)

#### ---- Shrub density ----
fope.plot_shrub <-
  ggplot(fope, aes(x = shr_pt_density)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Shrub Density — Fort Pearce",
       x = "Shrub Density (m^2)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_shrub.png",
       plot = fope.plot_shrub, width = 6, height = 4, dpi = 300)

fope.plot_shrub.det <-
  ggplot(fope, aes(x = shr_pt_density, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Shrub Density by Detection Status — Fort Pearce",
       x = "Shrub Density (m^2)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_shrub-det.png",
       plot = fope.plot_shrub.det, width = 6, height = 4, dpi = 300)

#### ---- Prey Presence ----
fp.plot_pp <-
  ggplot(fope, aes(x = pp_total)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Prey Presence - Fort Pearce",
       x = "Prey Presence (count)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_pp.png",
       plot = fp.plot_pp, width = 6, height = 4, dpi = 300)

fp.pp.det <-
  ggplot(fope, aes(x = pp_total, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Prey Presence by Detection Status - Fort Pearce",
       x = "Prey Presence (count)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_pp-det.png",
       plot = fp.pp.det, width = 6, height = 4, dpi = 300)

### ---- Vapor Pressure Deficit ----
fp.plot_vpd <-
  ggplot(fope, aes(x = vpd)) +
  geom_density(fill = "darkorange2", color = "darkorange2", alpha = 0.7) +
  labs(title = "Density Plot of Vapor Pressure Deficit - Fort Pearce",
       x = "Vapor Pressure Deficit (kPa)", y = "Density") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_vpd.png",
       plot = fp.plot_vpd, width = 6, height = 4, dpi = 300)

fp.plot_vpd.det <-
  ggplot(fope, aes(x = vpd, fill = factor(detection))) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("0" = "grey", "1" = "darkorange2"),
                    labels = c("No Detection","Detection")) +
  labs(title = "Vapor Pressure Deficit by Detection Status - Fort Pearce",
       x = "Vapor Pressure Deficit (kPa)", y = "Density", fill = "Detection") +
  theme_minimal()

ggsave("images/density-plots/fort-pearce/fp-density_vpd-det.png",
       plot = fp.plot_vpd.det, width = 6, height = 4, dpi = 300)