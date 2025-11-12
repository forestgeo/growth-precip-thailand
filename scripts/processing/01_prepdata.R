# Code to prep data for analysing the sensitivity of growth to drought events in HKK

rm(list = ls())
# load libraries---------------------------
library(tidyverse)
library(lubridate)
library(httr)
library(janitor)
library(raster)
library(ggplot2)

# reading the GitHub PAT (stored privately and locally)
# pat<-as.character(read.table("../data/HKK-dendro/git_pat.txt")[1,1])

# remotes:::download("DendroByBand.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-dendrobands/main/data/modified_data/DendroByBand.RData",
#                    auth_token = pat)


# tree x time variables---------------------------
load("data/HKK-dendro/DendroByBand.RData")

dendrobyband <- merge(dendrobyband, trees, by = "Tag")

hkk_OBSID$year <- year(hkk_OBSID$Date2)
hkk_OBSID$month <- month(hkk_OBSID$Date2)
hkk_OBSID$day <- day(hkk_OBSID$Date2)

# correcting two sets of wrong dates

hkk_OBSID$year[hkk_OBSID$year == 2209] <- 2009
hkk_OBSID$year[hkk_OBSID$month == 1 & hkk_OBSID$Cno == 5] <- 2011

dendrobyband <- merge(dendrobyband, hkk_OBSID, by = c("OBSID", "Cno", "Date2", "Tag"))

dendrobyband$Date2 <- as.Date(paste0(dendrobyband$year, "-", dendrobyband$month, "-", dendrobyband$day))


# calculate dbh by stem-------------------------------

# remove all data with decimal errors, misplaced decimal

dendro.data <- dendrobyband %>%
  dplyr::rename(sp = SPCODE.UPDATE) %>%
  mutate(year = year(Date2)) %>%
  filter(
    is.na(DBH1_flag), DBH1_error == "none",
    win.decimal.error == "none", is.na(win.flag)
  ) %>%
  group_by(Cno, sp, Tag) %>%
  dplyr::summarise(
    calcDBH = calcDBH[which(n.dendro == min(n.dendro))],
    Date2 = Date2[which(n.dendro == min(n.dendro))],
    n.dendro = min(n.dendro), # which dendroband is it from?
    Cii = max(Cii, na.rm = T)
  )

dbh.data <- dendrobyband %>%
  dplyr::rename(sp = SPCODE.UPDATE) %>%
  mutate(year = year(Date2)) %>%
  filter(is.na(dbh.flag), dbh.decimal.error == "none") %>%
  group_by(Cno, sp, Tag) %>%
  dplyr::summarise(
    dbh = DbhDendro[which(n.dendro == min(n.dendro))],
    Date2 = Date2[which(n.dendro == min(n.dendro))],
    n.dendro = min(n.dendro),
    Cii = max(Cii, na.rm = T)
  )

# calculate long-term growth rates-----------------------------

# remotes:::download("../data/HKK-dendro/stem1.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem1.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem2.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem2.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem3.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem3.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem4.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem4.rdata",
#                    auth_token = pat)

load("data/HKK-dendro/stem1.RData")
load("data/HKK-dendro/stem2.RData")
load("data/HKK-dendro/stem3.RData")
load("data/HKK-dendro/stem4.RData")

hkk_census <- bind_rows(hkk.stem1, hkk.stem2, hkk.stem3, hkk.stem4)

# hkk_census <- hkk_census %>%
#   filter(status=="A")%>%
#   filter(tag %in% dendro.data$Tag)%>%
#   group_by(tag, StemTag)%>%
#   dplyr::mutate(dbhmin1 = dplyr::lag(dbh, n=1, default=NA, order_by=StemTag),
#                 inc_annual = (dbh-dbhmin1)/5)%>%
#   filter(!is.na(dbhmin1))%>%
#   #average increment in cm
#   dplyr::summarise(avg_inc = mean(inc_annual/10, na.rm=T),
#                    sd_inc = sd(inc_annual/10, na.rm=T))


# #since dendrobands are installed on trees with tag==StemTag
# hkk_census<-hkk_census %>%
# filter(tag==StemTag)

# #make the colnames compatible
# hkk_census<-hkk_census %>%
# rename(Tag=tag)

# #calculate increments-----------------------------

# #add long-term increment to dendro data
# dendro.data<-merge(dendro.data, hkk_census, by="Tag", all.x=T)
# dbh.data<-merge(dbh.data, hkk_census, by="Tag", all.x=T)

# making increments for each year

alt_census <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)

dendro_inc <- dendro.data %>%
  filter(Cno %in% alt_census) %>%
  arrange(Cno) %>%
  group_by(Tag) %>% # grouping by treeID to calculate increment
  dplyr::mutate(
    calcDBH_min1 = dplyr::lag(calcDBH, n = 1, default = NA, order_by = Tag),
    n.dendro_diff = n.dendro - dplyr::lag(n.dendro, n = 1, default = NA, order_by = Tag),
    cii_min1 = dplyr::lag(Cii, n = 1, default = NA, order_by = Tag),
    tdif = difftime(Date2, dplyr::lag(Date2, n = 1, default = NA, order_by = Tag), units = "days"),
    inc = calcDBH - calcDBH_min1,
    inc_annual = (inc / as.numeric(tdif)) * 365,
    inc.time = Cno - dplyr::lag(Cno, n = 1, default = NA, order_by = Tag),
    # removing increment value if the previous measure was before the last year
    inc_annual = ifelse(inc.time == 2, inc_annual, NA),
    treeID = paste0("HKK.", Tag)
  ) %>%
  filter(!is.na(inc)) %>%
  filter(!is.na(inc_annual)) %>%
  ungroup() %>%
  dplyr::mutate(
    diainc.scaled = scale(inc_annual),
    out.inc = ifelse(abs(diainc.scaled) > 3, 1, 0)
  ) %>%
  filter(out.inc == 0)


dbh_inc <- dbh.data %>%
  filter(Cno %in% alt_census) %>%
  arrange(Cno) %>%
  group_by(Tag) %>% # grouping by treeID to calculate increment
  dplyr::mutate(
    dbh_min1 = dplyr::lag(dbh, n = 1, default = NA, order_by = Tag),
    n.dendro_diff = n.dendro - dplyr::lag(n.dendro, n = 1, default = NA, order_by = Tag),
    cii_min1 = dplyr::lag(Cii, n = 1, default = NA, order_by = Tag),
    tdif = difftime(Date2, dplyr::lag(Date2, n = 1, default = NA, order_by = Tag), units = "days"),
    inc = dbh - dbh_min1,
    inc_annual = (inc / as.numeric(tdif)) * 365,
    inc.time = Cno - dplyr::lag(Cno, n = 1, default = NA, order_by = Tag),
    # removing increment value if the previous measure was before the last year
    inc_annual = ifelse(inc.time == 2, inc_annual, NA),
    treeID = paste0("HKK.", Tag)
  ) %>%
  filter(!is.na(inc)) %>%
  filter(!is.na(inc_annual)) %>%
  ungroup() %>%
  dplyr::mutate(
    diainc.scaled = scale(inc_annual),
    out.inc = ifelse(abs(diainc.scaled) > 3, 1, 0)
  ) %>%
  filter(out.inc == 0)

# identiy outliers------------------------------

# merge dbh and dendro data by Tag and Cno
inc_all <- merge(dendro_inc, dbh_inc, by = c("sp", "Tag", "Cno"), all = T)

# Perpendicular distance from point 'a' to a line with 'slope' and 'intercept'
dist_point_line <- function(a, slope, intercept) {
  b <- c(1, intercept + slope)
  c <- c(-intercept / slope, 0)
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1, v2)
  return(abs(det(m)) / sqrt(sum(v1 * v1)))
}

inc_all <- inc_all %>%
  dplyr::rename(
    inc_annual.dendro = inc_annual.x,
    inc_annual.dbh = inc_annual.y,
    band.diff.dendro = n.dendro_diff.x,
    band.diff.dbh = n.dendro_diff.y
  ) %>%
  # split into three size classes
  dplyr::mutate(
    size_class = ifelse(dbh < 35, "small",
      ifelse(dbh < 70, "medium", "large")
    ),
    size_class = factor(size_class, levels = c("small", "medium", "large"))
    # find the deviation of the increment from the 1:1 line
  ) %>%
  dplyr::select(
    sp, Tag, Cno, inc_annual.dendro, inc_annual.dbh, dbh,
    size_class, band.diff.dbh, band.diff.dendro
  )

# using the distance to the line
dev.1 <- mapply(function(x) dist_point_line(c(inc_all$inc_annual.dendro[x], inc_all$inc_annual.dbh[x]), 1, 0), x = 1:nrow(inc_all))
inc_all$dev.1 <- dev.1

# remove outliers-----------------------------------

# outlier dataframe
outliers <- inc_all %>% filter(dev.1 >= 0.5)
# write.csv(outliers, "growth-precip-thailand/data/HKK-dendro/inc_outliers.csv", row.names=F)

# clean dataframe

dendro_inc <- merge(dendro_inc, inc_all %>% dplyr::select(Tag, Cno, sp, dev.1),
  by = c("Tag", "Cno", "sp"), all.x = T
)

dendro_inc_clean <- dendro_inc %>% filter(dev.1 < 0.5)

# calculate average increment for each tree from dendro measures

dendro_inc_clean <- dendro_inc_clean %>%
  group_by(Tag) %>%
  dplyr::mutate(avg_inc = mean(inc_annual, na.rm = T))

# how many avg_inc is -ve?
dendro_inc_clean %>%
  filter(avg_inc < 0) %>%
  nrow()
nrow(dendro_inc_clean)

# which species are these?
neg_species <- dendro_inc_clean %>%
  filter(avg_inc < 0) %>%
  group_by(sp) %>%
  dplyr::summarise(n.trees.neg = length(unique(Tag))) %>%
  ungroup() %>%
  arrange(n.trees.neg)

# total 72 trees with a maximum 5 trees from a species is removed with this approach
# LAGETO and SACCLI have 5 trees each with negative growth rates

# drop these trees
dendro_inc_clean <- dendro_inc_clean %>% filter(avg_inc >= 0)

# what is the distribution of avg_inc?
tree_incs <- dendro_inc_clean %>%
  group_by(Tag) %>%
  dplyr::summarise(avg_inc = mean(avg_inc, na.rm = T))

ggplot(tree_incs, aes(x = avg_inc)) +
  geom_density() +
  theme_minimal()

# how many trees with avg_inc < 0.001?
# sum(tree_incs$avg_inc<0.001) #58 trees
# nrow(tree_incs)

# remove these trees from the dendro_inc_clean dataset
dendro_inc_clean <- dendro_inc_clean %>% filter(avg_inc >= 0.01)

# calculate drought sensitivity-----------------------------

# here we define drought years as 2010 and 2015

sensitivity <- dendro_inc_clean %>%
  # filter(Cno %in% c(5, 15, 25))%>%
  group_by(Tag) %>% # grouping by treeID to calculate increment
  dplyr::mutate(
    # yr=ifelse(Cno==5, 2010, ifelse(Cno==15, 2015, 2020)),
    yr = 2008 + ((Cno - 1) / 2),
    # sens.div = inc_annual/avg_inc,
    # sens.dif = inc_annual-avg_inc,
    inc_annual_zero = ifelse(inc_annual < 0, 0, inc_annual),
    sens.prop = (inc_annual - avg_inc) / avg_inc,
    sens.prop.zero = (inc_annual_zero - avg_inc) / avg_inc
  ) %>%
  ungroup() %>%
  filter(!is.infinite(sens.prop)) %>%
  group_by(Cno) %>%
  # making quantiles for each year
  dplyr::mutate(
    # quantiles_div=ntile(resistance.div, 4),
    # quantiles_dif=ntile(resistance.dif, 4),
    quantiles_prop = ntile(sens.prop, 4)
  )

# merge with tree attributes

sensitivity <- merge(sensitivity, dplyr::select(trees, -"SPCODE.UPDATE"), by = "Tag", all.x = T)

tree.time <- sensitivity %>%
  dplyr::select(Tag, treeID, Cno, yr, inc_annual, inc_annual_zero, sens.prop, sens.prop.zero, avg_inc, cii_min1, calcDBH_min1) %>%
  rename(avg_inc_tree = avg_inc)

# tree table---------------------------

# pick trees in the tree.time dataset

trees <- trees %>% filter(Tag %in% tree.time$Tag)

# twi values
twi <- raster::raster("data/HKK-other/TWI.tif")

# pick twi values for each tree
trees$twi <- raster::extract(twi, 0.1 * trees[, c("Y", "X")]) # multiply by 0.1 to get the correct resolution

# tpi
tpi <- raster::raster("data/HKK-other/TPI_3.tif")
trees$tpi <- raster::extract(tpi, 0.1 * trees[, c("Y", "X")])

# plot twi and tpi for all numbers

# tpi_7 <- ggplot(trees, aes(x = twi, y = -tpi)) +
#   geom_point() +
#   geom_abline() +
#   theme_bw()

# library(patchwork)
# layout <- "
# AB
# CD
# "
# png("doc/display/explore/twi_tpi.png", width = 5, height = 5, units = "in", res = 300)
# tpi_1 + tpi_3 + tpi_5 + tpi_7 + plot_layout(design = layout) +
#   plot_annotation(title = "TWI and TPI for HKK trees")
# dev.off()

# number of trees left out with each window size of tpi
# 1- 64, 3 - 259, 5 - 483, 7 - 721

# habitat---------------------------

# read habitat data
habitat <- read.csv("data/HKK-other/HKK_habtype(in).csv")
# convert habitat dataframe to raster
habitat_raster <- raster::rasterFromXYZ(habitat[, c("x", "y", "hab")])

# pick habitat values for each tree
trees$habitat <- raster::extract(habitat_raster, trees[, c("X", "Y")])

# Neighborhood Crowding Index--------------------------

# read NCI values calculated using 00_nci.R

nci <- read.csv("data/HKK-dendro/nci.csv")

nci <- nci %>%
  rename(Tag = tag)

trees <- merge(trees, dplyr::select(nci, -"X"), by = "Tag", all.x = T)

tree_vars <- trees %>%
  rename(Species = SPCODE.UPDATE) %>%
  dplyr::select(Tag, Species, X, Y, twi, tpi, habitat, nci_15)


# species table---------------------------

# read in the deciduousness values
dec_williams <- read.csv("data/HKK-dendro/species.csv")

dec_williams <- dec_williams %>%
  filter(n.tree >= 20)


# filter trees based on these species

tree_vars <- tree_vars %>%
  filter(Species %in% dec_williams$sp)

# filter observations based on these trees

tree.time <- tree.time %>%
  filter(Tag %in% tree_vars$Tag)

# species wise growth rates---------------------------

# use census data for these species
hkk_census_growth <- hkk_census %>%
  filter(sp %in% dec_williams$sp) %>%
  filter(status == "A") %>%
  group_by(sp, tag) %>%
  dplyr::mutate(
    dbhmin1 = dplyr::lag(dbh, n = 1, default = NA, order_by = tag),
    inc_annual = (dbh - dbhmin1) / 5
  ) %>%
  filter(!is.na(dbhmin1)) %>%
  ungroup() %>%
  group_by(sp) %>%
  # average increment in cm
  dplyr::summarise(
    avg_inc = mean(inc_annual, na.rm = T),
    median_inc = median(inc_annual, na.rm = T)
  )

# max DBH and occupancy for each species from the census

hkk.stem4$habitat <- raster::extract(habitat_raster, hkk.stem4[c("gx", "gy")])
hkk.stem4$twi <- raster::extract(twi, 0.1 * hkk.stem4[, c("gy", "gx")]) # multiply by 0.1 to get the correct resolution
hkk.stem4$tpi <- raster::extract(tpi, 0.1 * hkk.stem4[, c("gy", "gx")]) # multiply by 0.1 to get the correct resolution
head(hkk.stem4)

# plot SACCLI trees on twi to make sure that the scale is right

twi_df <- as.data.frame(twi, xy = TRUE)

# png("doc/display/explore/twi_SACCLI.png", width = 5, height = 8, units = "in", res = 300)
# ggplot() +
#   geom_raster(data = twi_df, aes(x = x, y = y, fill = TWI)) +
#   scale_fill_distiller(palette = "Spectral", direction = 1) +
#   geom_point(data = hkk.stem4[hkk.stem4$sp == "SACCLI", ], aes(x = 0.1 * gy, y = 0.1 * gx), pch = 16, alpha = 0.5) +
#   coord_fixed()
# dev.off()

# # plot distribution of twi across top 4 species
# png("doc/display/explore/twi_dist_top4.png", width = 5, height = 5, units = "in", res = 300)
# ggplot(hkk.stem4 %>% filter(sp %in% c("SACCLI", "HOPEOD", "POLYVI", "TETRNU")), aes(x = twi, fill = sp)) +
#   geom_density(alpha = 0.5) +
#   scale_fill_viridis_d() +
#   theme_minimal() +
#   labs(title = "Distribution of top 4 species in TWI space (HKK)")
# dev.off()

# plot distribution as violin plot for top 10 species
dec_william_sort <- dec_williams %>%
  arrange(williams_dec)

# png("doc/display/explore/twi_violin_top10.png", width = 5, height = 5, units = "in", res = 300)
# ggplot(hkk.stem4 %>% filter(sp %in% dec_william_sort$sp), aes(x = factor(sp, levels = dec_william_sort$sp), y = twi, fill = sp)) +
#   geom_violin() +
#   geom_boxplot(width = 0.1) +
#   scale_fill_viridis_d() +
#   theme_minimal() +
#   labs(title = "Distribution of top 10 species in TWI space (HKK)")
# dev.off()

sp.traits <- hkk.stem4 %>%
  filter(status == "A") %>%
  filter(sp %in% dec_williams$sp) %>%
  pivot_wider(
    names_from = habitat,
    values_from = habitat,
    names_prefix = "hab"
  ) %>%
  group_by(sp) %>%
  dplyr::summarise(
    maxDBH = max(dbh, na.rm = T),
    # habitat occupancy
    hab1 = ifelse(is.na(sum(hab1, na.rm = T)) == T, 0, 1),
    hab2 = ifelse(is.na(sum(hab2, na.rm = T)) == T, 0, 1),
    hab3 = ifelse(is.na(sum(hab3, na.rm = T)) == T, 0, 1),
    hab4 = ifelse(is.na(sum(hab4, na.rm = T)) == T, 0, 1),
    hab5 = ifelse(is.na(sum(hab5, na.rm = T)) == T, 0, 1),
    hab6 = ifelse(is.na(sum(hab6, na.rm = T)) == T, 0, 1),
    # twi mean
    twi_mean = mean(twi, na.rm = T),
    twi_median = median(twi, na.rm = T),
    # twi sd
    twi_sd = sd(twi, na.rm = T),
    # twi weighted mean with dbh
    twi_dbh_mean = mean(twi * dbh, na.rm = T),
    twi_dbh_median = median(twi * dbh, na.rm = T),
    twi_dbh_sd = sd(twi * dbh, na.rm = T),
    # twi weighted mean with agb
    twi_agb_mean = mean(twi * agb, na.rm = T),
    twi_agb_median = median(twi * agb, na.rm = T),
    twi_agb_sd = sd(twi * agb, na.rm = T),
    # tpi mean
    tpi_mean = mean(tpi, na.rm = T),
    tpi_median = median(tpi, na.rm = T),
    # twi sd
    tpi_sd = sd(tpi, na.rm = T),
    # twi weighted mean with dbh
    tpi_dbh_mean = mean(tpi * dbh, na.rm = T),
    tpi_dbh_median = median(tpi * dbh, na.rm = T),
    tpi_dbh_sd = sd(tpi * dbh, na.rm = T),
    # twi weighted mean with agb
    tpi_agb_mean = mean(tpi * agb, na.rm = T),
    tpi_agb_median = median(tpi * agb, na.rm = T),
    tpi_agb_sd = sd(tpi * agb, na.rm = T)
  )

sp.traits <- merge(sp.traits, hkk_census_growth, by = "sp")


sp_vars <- merge(sp.traits, dplyr::select(dec_williams, c("sp", "spfull", "williams_dec")), by = "sp")

# make an RData object with these dataframes
saveRDS(list(tree_vars = tree_vars, tree.time = tree.time, sp_vars = sp_vars), "data/HKK-dendro/sensitivity_data.RData")


rm(list = ls())
datasets <- readRDS("data/HKK-dendro/sensitivity_data.RData")

tree.time <- datasets$tree.time
tree_vars <- datasets$tree_vars
sp_vars <- datasets$sp_vars

colnames(sp_vars)[1] <- "Species"

# merge tree_vars and sp_vars
tree_vars <- merge(tree_vars, sp_vars, by = "Species", all.x = TRUE)
# merge tree_vars and tree.time
tree.time <- merge(tree.time, tree_vars, by = "Tag", all.x = TRUE)

colnames(tree.time)

# standardise the variables within species and across all--------------------
yrs <- c(2010, 2015, 2020)

tree.time <- tree.time %>%
  dplyr::mutate(
    calcDBH_min1_scaled = scale(calcDBH_min1, center = TRUE, scale = TRUE),
    cii_min1_scaled = scale(cii_min1, center = TRUE, scale = TRUE),
    cii_min1_scaled = ifelse(cii_min1_scaled == "NaN", 0, cii_min1_scaled),
    twi_scaled = scale(twi, center = TRUE, scale = TRUE),
    tpi_scaled = scale(tpi, center = TRUE, scale = TRUE),
    median_inc_scaled = scale(median_inc, center = TRUE, scale = TRUE),
    avg_inc_tree_scaled = scale(avg_inc_tree, center = TRUE, scale = TRUE)
  ) %>%
  group_by(Species) %>%
  # scale while retaining original values
  dplyr::mutate(
    calcDBH_min1_scaled_sp = scale(calcDBH_min1, center = TRUE, scale = TRUE),
    cii_min1_scaled_sp = scale(cii_min1, center = TRUE, scale = TRUE),
    cii_min1_scaled_sp = ifelse(cii_min1_scaled_sp == "NaN", 0, cii_min1_scaled_sp),
    twi_scaled_sp = scale(twi, center = TRUE, scale = TRUE),
    tpi_scaled_sp = scale(tpi, center = TRUE, scale = TRUE)
  ) %>%
  ungroup() %>%
  # remove large outliers for each year
  group_by(yr) %>%
  # find sens.prop values that are 3 sds from the mean
  #    dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 3 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 3 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
  dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 4 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 4 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
  filter(!is.na(sens.prop) & !is.na(cii_min1) & !is.na(calcDBH_min1) & !is.na(twi)) %>%
  ungroup()

# make a plot with anomalies instead of growth increments

# first summarise anomalies by species
tree.time.anom <- tree.time %>%
  group_by(spname, yr) %>%
  dplyr::summarise(
    median_inc = median(inc_annual, na.rm = T)
  ) %>%
  ungroup() %>%
  group_by(spname) %>%
  dplyr::mutate(
    mean = mean(median_inc, na.rm = T),
    sd = sd(median_inc, na.rm = T),
    anomaly = (median_inc - mean) / sd
  ) %>%
  ungroup()

tree.time.anom.all <- tree.time %>%
  group_by(yr) %>%
  dplyr::summarise(
    median_inc = median(inc_annual, na.rm = T)
  ) %>%
  dplyr::mutate(
    mean = mean(median_inc, na.rm = T),
    sd = sd(median_inc, na.rm = T),
    anomaly = (median_inc - mean) / sd
  ) %>%
  ungroup() %>%
  dplyr::mutate(spname = "all")

tree.time.anom <- bind_rows(tree.time.anom, tree.time.anom.all)

# bar plot of anomalies

spagplot_top10_anom <- ggplot() +
  # species plots
  geom_bar(
    data = tree.time.anom %>% filter(spname %in% top_10_sp$spfull),
    aes(x = yr, y = anomaly, fill = spname), stat = "identity",
    position = "dodge"
  ) +
  # mean of all trees
  geom_bar(
    data = tree.time.anom %>% filter(spname == "all"),
    aes(x = yr, y = anomaly), fill = NA, col = "black", linewidth = 1, stat = "identity",
    position = "dodge"
  ) +
  # make all years show on x-axis
  scale_x_continuous(breaks = unique(tree.time$yr)) +
  scale_fill_viridis_d() +
  geom_hline(yintercept = c(-2, -1, 0, 1, 2), linetype = "dashed") +
  # add text on these lines
  geom_text(aes(
    x = c(2010, 2015, 2020), y = 1,
    label = c("ENSO drought", "ENSO drought", "drought")
  ), hjust = 0.8, vjust = -0.2, angle = 90) +
  guides(fill = guide_legend("species"), nrow = 3) +
  xlab("year") +
  ylab("anomaly") +
  # ggtitle("growth increments for top 10 species") +
  theme_bw() +
  # theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

spagplot_top10_anom

png("doc/display/growth_anom_fullseries.png", width = 8, height = 8, units = "in", res = 300)
spagplot_top10_anom
dev.off()


# SI figure - map of dendroband trees ---------------

dendro_tree_map_sp <- ggplot() +
  geom_point(data = tree.time %>% filter(yr %in% c(2010, 2015, 2020)), aes(x = Y, y = X, size = calcDBH_min1), alpha = 0.3) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(spfull ~ yr, nrow = 4) +
  # theme(
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.ticks.y = element_blank()
  # ) +
  coord_fixed()

png("doc/display/dendro_tree_map_sp.png", width = 36, height = 16, units = "in", res = 300)
dendro_tree_map_sp
dev.off()


dendro_tree_map <- ggplot() +
  geom_point(data = tree.time %>% filter(yr %in% c(2010, 2015, 2020)), aes(x = Y, y = X, size = calcDBH_min1, col = spfull), alpha = 0.5) +
  scale_color_viridis_d() +
  theme_bw() +
  facet_wrap(~yr) +
  xlab("metres east") +
  ylab("metres north") +
  guides(
    size = guide_legend(title = "DBH"),
    color = guide_legend(title = "Species")
  ) +
  # theme(
  #     legend.position = "none",
  #     axis.text.x = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.ticks.y = element_blank()
  # ) +
  coord_fixed()

png("doc/display/dendro_tree_map.png", width = 24, height = 14, units = "in", res = 300)
dendro_tree_map
dev.off()


colnames(tree.time)

tree.time.public <- tree.time %>%
  dplyr::select(Tag, treeID, Species, spfull, williams_dec, yr, sens.prop, cii_min1, calcDBH_min1_scaled, calcDBH_min1_scaled_sp, twi, twi_scaled, twi_scaled_sp, tpi, tpi_scaled, tpi_scaled_sp) %>%
  filter(yr %in% yrs) %>%
  arrange(yr)

colnames(tree.time.public)
unique(tree.time.public$yr)

write.csv(tree.time.public, "data/dendro/sensitivity_dataset.csv", row.names = F)

median_incs <- tree.time %>%
  group_by(yr) %>%
  dplyr::summarise(median_inc_mm = median(inc_annual) * 10)

median_incs

median_inc_sp <- tree.time %>%
  group_by(yr, Species) %>%
  dplyr::summarise(median_inc = median(inc_annual) * 10) %>%
  ungroup() %>%
  group_by(yr) %>%
  dplyr::summarise(median_inc = median(median_inc))

median_inc_sp

median_incs$median_inc_sp_mm <- median_inc_sp$median_inc

median_incs

sp_median <- tree.time %>%
  group_by(Species) %>%
  dplyr::summarise(median_inc = median(inc_annual) * 10)

median_incs_all <- data.frame(yr = "all", median_inc_mm = median(tree.time$inc_annual, na.rm = T) * 10, median_inc_sp_mm = median(sp_median$median_inc))

median_incs <- rbind(median_incs, median_incs_all)
median_incs

write.csv(median_incs, "data/dendro/summaries_dataset.csv", row.names = F)
