# code to calculate TWI from elevation data

# load libraries------------------

library(sf)
library(terra)
library(tidyverse)
library(raster)
library(tmap)
library(spatialEco)

library(whitebox)
# whitebox::install_whitebox()
whitebox::wbt_init()


# load data-----------------------
# elevation
# although .rdata, this needs to be read using read.table - thanks Bier!

hkk_elev <- read.table("data/HKK-other/CTFSElev_hkk.rdata")

# convert to topographic wetness index
# following code from https://vt-hydroinformatics.github.io/rgeoraster.html


theme_set(theme_classic())

# convert to raster

# first convert hkk_elev$col.elev to matrix
elev_mat <- matrix(hkk_elev$col.elev, nrow = 201, ncol = 101, byrow = T)

# then convert to raster
elev_rast <- raster(elev_mat,
  xmn = 0, xmx = 50, ymn = 0, ymx = 100,
  crs = "+proj=utm +zone=47 +datum=WGS84"
)
plot(elev_rast)

writeRaster(elev_rast, "data/HKK-other/HKKDEMsm_5m.tif", overwrite = TRUE)

tm_shape(elev_rast) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE) +
  tm_scale_bar()

# hillshade

wbt_hillshade(
  dem = "data/HKK-other/HKKDEMsm_5m.tif",
  output = "data/HKK-other/HKK_hillshade.tif",
  azimuth = 115
)

hillshade <- raster("data/HKK-other/HKK_hillshade.tif")

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar()

wbt_breach_depressions_least_cost(
  dem = "data/HKK-other/HKKDEMsm_5m.tif",
  output = "data/HKK-other/HKKDEMsm_breached.tif",
  dist = 5,
  fill = TRUE
)

wbt_fill_depressions_wang_and_liu(
  dem = "data/HKK-other/HKKDEMsm_breached.tif",
  output = "data/HKK-other/HKKDEMsm_filled_breached.tif"
)

filled_breached <- raster("data/HKK-other/HKKDEMsm_filled_breached.tif")

## What did this do?
difference <- elev_rast - filled_breached

difference[difference == 0] <- NA

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar() +
  tm_shape(difference) +
  tm_raster(style = "cont", legend.show = TRUE) +
  tm_scale_bar()


# TWI

wbt_d_inf_flow_accumulation(
  input = "data/HKK-other/HKKDEMsm_filled_breached.tif",
  output = "data/HKK-other/HKK_DinfFAsca.tif",
  out_type = "Specific Contributing Area"
)

wbt_slope(
  dem = "data/HKK-other/HKKDEMsm_filled_breached.tif",
  output = "data/HKK-other/demslope.tif",
  units = "degrees"
)

wbt_wetness_index(
  sca = "data/HKK-other/HKK_DinfFAsca.tif",
  slope = "data/HKK-other/demslope.tif",
  output = "data/HKK-other/TWI.tif"
)

library(raster)
twi <- raster("data/HKK-other/TWI.tif")

str(twi)

# find what percent cells have negative values
sum(twi[] < 0) / length(twi[]) * 100

min(twi[], na.rm = TRUE)
max(twi[], na.rm = TRUE)


# twi[twi > 0] <- NA

# tm_shape(hillshade)+
#   tm_raster(style = "cont",palette = "-Greys", legend.show = FALSE)+
tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()

png("doc/display/explore/twi.png", width = 4, height = 8, units = "in", res = 300)
tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 1) +
  tm_scale_bar()
dev.off()


# make TWI with SRTM--------------------

# load data-----------------------
# elevation


# convert to topographic wetness index
# following code from https://vt-hydroinformatics.github.io/rgeoraster.html

hkk_elev <- raster("data/HKK-other/N15E099.SRTMGL1.hgt/N15E099.hgt")

theme_set(theme_classic())

# read hkk shapefile
hkkwls <- read_sf("data/HKK-other/HKK_boundary/shapefile/WDPA_WDOECM_Feb2025_Public_1407_shp-polygons.shp")

# HKK plot as point
# from the HKK site data repo (15.6324 N, 99.217 E)
hkk_point <- st_point(c(99.217, 15.6324))

# Convert hkk_point to an sf object
hkk_point_sf <- st_sfc(hkk_point, crs = st_crs(hkkwls))
hkk_point_sf <- st_sf(geometry = hkk_point_sf)

writeRaster(hkk_elev, "data/HKK-other/HKKDEM_30m.tif", overwrite = TRUE)


png("doc/display/explore/hkk_elev_30m.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(hkk_elev) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE) +
  tm_shape(hkkwls) +
  tm_borders(lwd = 4) +
  tm_shape(hkk_point_sf) +
  tm_symbols(size = 2, shape = 15, border.lwd = 2, col = NA) +
  tm_scale_bar()
dev.off()

# hillshade

wbt_hillshade(
  dem = "data/HKK-other/HKKDEM_30m.tif",
  output = "data/HKK-other/HKK_hillshade_30m.tif",
  azimuth = 115
)

hillshade <- raster("data/HKK-other/HKK_hillshade_30m.tif")

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar()

wbt_breach_depressions_least_cost(
  dem = "data/HKK-other/HKKDEM_30m.tif",
  output = "data/HKK-other/HKKDEM_30m_breached.tif",
  dist = 5,
  fill = TRUE
)

wbt_fill_depressions_wang_and_liu(
  dem = "data/HKK-other/HKKDEM_30m_breached.tif",
  output = "data/HKK-other/HKKDEM_30m_filled_breached.tif"
)

filled_breached <- raster("data/HKK-other/HKKDEM_30m_filled_breached.tif")

## What did this do?
difference <- hkk_elev - filled_breached

difference[difference == 0] <- NA

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar() +
  tm_shape(difference) +
  tm_raster(style = "cont", legend.show = TRUE) +
  tm_scale_bar()


# TWI

wbt_d_inf_flow_accumulation(
  input = "data/HKK-other/HKKDEM_30m_filled_breached.tif",
  output = "data/HKK-other/HKK_DinfFAsca_30m.tif",
  out_type = "Specific Contributing Area"
)

wbt_slope(
  dem = "data/HKK-other/HKKDEM_30m_filled_breached.tif",
  output = "data/HKK-other/demslope_30m.tif",
  units = "degrees"
)

wbt_wetness_index(
  sca = "data/HKK-other/HKK_DinfFAsca_30m.tif",
  slope = "data/HKK-other/demslope_30m.tif",
  output = "data/HKK-other/TWI_30m.tif"
)

library(raster)
twi <- raster("data/HKK-other/TWI_30m.tif")

str(twi)

# find what percent cells have negative values
sum(twi[] < 0) / length(twi[]) * 100

min(twi[], na.rm = TRUE)
max(twi[], na.rm = TRUE)


# twi[twi > 0] <- NA

# tm_shape(hillshade)+
#   tm_raster(style = "cont",palette = "-Greys", legend.show = FALSE)+
tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()

png("doc/display/explore/twi_30m.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 1) +
  tm_scale_bar()
dev.off()


# repeat with clipped raster -----------

# Define the size of the square buffer (e.g., 1000 meters)
buffer_size <- 0.05

# Create a bounding box around the point
bbox <- st_bbox(hkk_point_sf)
bbox <- bbox + c(-buffer_size, -buffer_size, buffer_size, buffer_size)

# Convert the bounding box to an sf object
bbox_sf <- st_as_sfc(bbox)

# Convert the bounding box to a Spatial object
bbox_sp <- as(bbox_sf, "Spatial")

# Clip the raster using the bounding box
clipped_raster <- crop(mask(hkk_elev, bbox_sp), bbox_sp)


# Plot the clipped raster
png("doc/display/explore/hkk_elev_30m_clipped.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(clipped_raster) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE) +
  # tm_shape(hkkwls) +
  # tm_borders(lwd = 4) +
  tm_shape(hkk_point_sf) +
  tm_symbols(size = 2, shape = 15, border.lwd = 2, col = NA) +
  tm_scale_bar()
dev.off()

writeRaster(clipped_raster, "data/HKK-other/HKKDEM_30m_clipped.tif", overwrite = TRUE)

# hillshade

wbt_hillshade(
  dem = "data/HKK-other/HKKDEM_30m_clipped.tif",
  output = "data/HKK-other/HKK_hillshade_30m_clipped.tif",
  azimuth = 115
)

hillshade <- raster("data/HKK-other/HKK_hillshade_30m_clipped.tif")

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar()

wbt_breach_depressions_least_cost(
  dem = "data/HKK-other/HKKDEM_30m_clipped.tif",
  output = "data/HKK-other/HKKDEM_30m_clipped_breached.tif",
  dist = 5,
  fill = TRUE
)

wbt_fill_depressions_wang_and_liu(
  dem = "data/HKK-other/HKKDEM_30m_clipped_breached.tif",
  output = "data/HKK-other/HKKDEM_30m_clipped_filled_breached.tif"
)

filled_breached <- raster("data/HKK-other/HKKDEM_30m_clipped_filled_breached.tif")

## What did this do?
difference <- clipped_raster - filled_breached

difference[difference == 0] <- NA

tm_shape(hillshade) +
  tm_raster(style = "cont", palette = "-Greys", legend.show = FALSE) +
  tm_scale_bar() +
  tm_shape(difference) +
  tm_raster(style = "cont", legend.show = TRUE) +
  tm_scale_bar()

# TWI

wbt_d_inf_flow_accumulation(
  input = "data/HKK-other/HKKDEM_30m_clipped_filled_breached.tif",
  output = "data/HKK-other/HKK_DinfFAsca_30m_clipped.tif",
  out_type = "Specific Contributing Area"
)

wbt_slope(
  dem = "data/HKK-other/HKKDEM_30m_clipped_filled_breached.tif",
  output = "data/HKK-other/demslope_30m_clipped.tif",
  units = "degrees"
)

wbt_wetness_index(
  sca = "data/HKK-other/HKK_DinfFAsca_30m_clipped.tif",
  slope = "data/HKK-other/demslope_30m_clipped.tif",
  output = "data/HKK-other/TWI_30m_clipped.tif"
)

library(raster)
twi <- raster("data/HKK-other/TWI_30m_clipped.tif")

str(twi)

# find what percent cells have negative values
sum(twi[] < 0) / length(twi[]) * 100

min(twi[], na.rm = TRUE)
max(twi[], na.rm = TRUE)


# twi[twi > 0] <- NA

# tm_shape(hillshade)+
#   tm_raster(style = "cont",palette = "-Greys", legend.show = FALSE)+
tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()


# clip further to hkk_point with x distance of 0.005 degrees and y-dist of 0.01 degrees
# Define the size of the square buffer (e.g., 1000 meters)
buffer_size_x <- 0.005
buffer_size_y <- 0.01

# Create a bounding box around the point
bbox <- st_bbox(hkk_point_sf)
bbox <- bbox + c(-buffer_size_x, -buffer_size_y, buffer_size_x, 2 * buffer_size_y)

# Convert the bounding box to an sf object
bbox_sf <- st_as_sfc(bbox)

# Convert the bounding box to a Spatial object
bbox_sp <- as(bbox_sf, "Spatial")

# Clip the raster using the bounding box
clipped_raster <- crop(mask(twi, bbox_sp), bbox_sp)

# plot clipped raster
png("doc/display/explore/twi_30m_clipped_to_plot.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(clipped_raster) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE) +
  # tm_shape(hkk_point_sf) +
  # tm_symbols(size = 2, shape = 15, border.lwd = 2, col = NA) +
  tm_scale_bar()
dev.off()


# TPI ---------------------------
# library(spatialEco)

# read elev_rast as terra SpatRaster

elev_rast_terra <- terra::rast("data/HKK-other/HKKDEMsm_5m.tif")
tpi_5 <- spatialEco::tpi(elev_rast_terra, scale = 5, win = "circle")

png("doc/display/explore/tpi_5.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(tpi_5) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()
dev.off()

tpi_3 <- spatialEco::tpi(elev_rast_terra, scale = 3, win = "circle")
png("doc/display/explore/tpi_3.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(tpi_3) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()
dev.off()

tpi_7 <- spatialEco::tpi(elev_rast_terra, scale = 7, win = "circle")
png("doc/display/explore/tpi_7.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(tpi_7) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()
dev.off()

tpi_1 <- spatialEco::tpi(elev_rast_terra, scale = 1, win = "circle")
png("doc/display/explore/tpi_1.png", width = 8, height = 8, units = "in", res = 300)
tm_shape(tpi_1) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE, alpha = 0.5) +
  tm_scale_bar()
dev.off()

# write these as rasters
terra::writeRaster(tpi_5, "data/HKK-other/TPI_5.tif", overwrite = TRUE)
terra::writeRaster(tpi_3, "data/HKK-other/TPI_3.tif", overwrite = TRUE)
terra::writeRaster(tpi_7, "data/HKK-other/TPI_7.tif", overwrite = TRUE)
terra::writeRaster(tpi_1, "data/HKK-other/TPI_1.tif", overwrite = TRUE)
