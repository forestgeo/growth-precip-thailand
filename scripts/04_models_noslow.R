# run the same models removing extremely slow growing trees from the sample


# Load required libraries---------------------
library(tidyverse)
library(brms)

# load data--------------------------
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

tree.time <- tree.time %>%
    dplyr::mutate(
        calcDBH_min1_scaled = scale(calcDBH_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled = scale(cii_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled = ifelse(cii_min1_scaled == "NaN", 0, cii_min1_scaled),
        twi_scaled = scale(twi, center = TRUE, scale = TRUE),
        median_inc_scaled = scale(median_inc, center = TRUE, scale = TRUE)
    ) %>%
    group_by(Species) %>%
    # scale while retaining original values
    dplyr::mutate(
        calcDBH_min1_scaled_sp = scale(calcDBH_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled_sp = scale(cii_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled_sp = ifelse(cii_min1_scaled_sp == "NaN", 0, cii_min1_scaled_sp),
        twi_scaled_sp = scale(twi, center = TRUE, scale = TRUE)
    ) %>%
    ungroup() %>%
    # remove large outliers for each year
    group_by(yr) %>%
    # find sens.prop values that are 3 sds from the mean
    dplyr::mutate(sens.prop = ifelse(sens.prop > mean(sens.prop, na.rm = TRUE) + 3 * sd(sens.prop, na.rm = TRUE), NA, sens.prop)) %>%
    filter(!is.na(sens.prop) & !is.na(cii_min1) & !is.na(calcDBH_min1) & !is.na(twi)) %>%
    ungroup()
