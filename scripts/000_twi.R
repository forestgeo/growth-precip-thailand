# code to calculate TWI from elevation data

# load libraries------------------

library(sf)
library(terra)
library(tidyverse)
library(raster)
library(tmap)

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

tm_shape(twi) +
  tm_raster(style = "cont", palette = "PuOr", legend.show = TRUE) +
  tm_scale_bar()
