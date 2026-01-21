## Checking UTM Accuracy
# SSS - Fall 2025

# Load required libraries
library(sf)
library(tidyverse)
library(leaflet)
library(janitor)

#import data
transect_points <- read_csv("data/data_gmocc_all-ut.csv") %>%
  janitor::clean_names()

unique(transect_points$site_id) #pick your site

# --- Fort Pearce ----
## --- shapefiles/geometry ----
fp.transect_points <- transect_points %>%
  filter(site_id == "fort_pearce") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
fp.transect_points_sf <- fp.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
fp.transect_lines_sf <- fp.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = fp.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = fp.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# --- White Reef ----
## --- shapefiles/geometry ----
wr.transect_points <- transect_points %>%
  filter(site_id == "white_reef") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
wr.transect_points_sf <- wr.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
wr.transect_lines_sf <- wr.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = wr.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = wr.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# --- Cove Wash ----
## --- shapefiles/geometry ----
cw.transect_points <- transect_points %>%
  filter(site_id == "cove_wash") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
cw.transect_points_sf <- cw.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
cw.transect_lines_sf <- cw.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = cw.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = cw.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# --- Paradise Canyon----
## --- shapefiles/geometry ----
pc.transect_points <- transect_points %>%
  filter(site_id == "paradise_canyon") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
pc.transect_points_sf <- pc.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
pc.transect_lines_sf <- pc.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = pc.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = pc.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# --- Sun River ----
## --- shapefiles/geometry ----
sr.transect_points <- transect_points %>%
  filter(site_id == "sun_river") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
sr.transect_points_sf <- sr.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
sr.transect_lines_sf <- sr.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = sr.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = sr.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# --- Turkey Farm ----
## --- shapefiles/geometry ----
tf.transect_points <- transect_points %>%
  filter(site_id == "turkey_farm") %>%
  filter(!is.na(utm_x) | !is.na(utm_y))

# Convert to sf points (Utah Zone 12N - EPSG:32612, WGS84 - EPSG:4326)
tf.transect_points_sf <- tf.transect_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

# Create transect line geometries by grouping and connecting points
tf.transect_lines_sf <- tf.transect_points_sf %>%
  group_by(tran_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addPolylines(data = tf.transect_lines_sf,                                         # Add transect lines 
               color = "blue",
               weight = 3,
               label = ~tran_id,
               group = "Transects") %>%
  addCircleMarkers(data = tf.transect_points_sf,                                    # Add points as markers
                   radius = 3,
                   color = "red",
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Collection Points") %>%
  addLayersControl(                                                               # Layer controls
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Collection Points"),
    options = layersControlOptions(collapsed = FALSE))

# ---- Detection Points ----
data_occ <- read_csv("data/data_gmocc_all-ut.csv")

# create detection dataset
det_points <- 
  data_occ %>%
    filter(grepl("1", detection))

## --- shapefiles/geometry ----
# get the points into sf
det_points_sf <- det_points %>%
  st_as_sf(coords = c("utm_x", "utm_y"), crs = 32612) %>%
  st_transform(crs = 4326)

## ---- Leaflet Map ----
  leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  addCircleMarkers(data = det_points_sf,                                            # Add points as markers
                   radius = 3,
                   color = "darkorange",
                   label = ~paste("Detection Point", tran_id),
                   group = "detection points") %>%                                  # Layer controls
  addLayersControl(
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Detection Points"),
    options = layersControlOptions(collapsed = FALSE))

# ---- Combining Transects ----
# bind all transect lines together
transect_lines_all <-
  bind_rows(fp.transect_lines_sf,
            wr.transect_lines_sf,
            tf.transect_lines_sf,
            cw.transect_lines_sf,
            pc.transect_lines_sf,
            sr.transect_lines_sf)

# bind all points together
transect_points_all <-
  bind_rows(fp.transect_points_sf,
            wr.transect_points_sf,
            tf.transect_points_sf,
            cw.transect_points_sf,
            pc.transect_points_sf,
            sr.transect_points_sf)

## ---- Leaflet Map ----
leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri Satellite") %>%
  # Add transect lines
  addPolylines(data = transect_lines_all,
               color = "black",
               weight = 3,
               opacity = 1,
               label = ~tran_id,
               group = "Transect Lines") %>%
  # Add points as markers
  addCircleMarkers(data = transect_points_all,
                   radius = 3,
                   color = "white",
                   opacity = 1,
                   label = ~paste(tran_id, "Collection Point", col_point),
                   group = "Transect Points") %>%
  addCircleMarkers(data = det_points_sf,
                   radius = 5,
                   color = "darkorange",
            
                   opacity = 1,
                   label = ~paste("Detection Point", tran_id),
                   group = "Detection Points") %>%
  # Layer controls 
  addLayersControl(
    baseGroups = c("OpenStreetMap", "Esri Satellite"),
    overlayGroups = c("Transects", "Transect points", "Detection points"),
    options = layersControlOptions(collapsed = FALSE))


  