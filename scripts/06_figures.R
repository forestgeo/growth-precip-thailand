# script to code for figures

# load libraries
library(tidyverse)
library(lubridate)

# Fig 1
# climate data

spei <- read.csv("data/climate/SPEI_HKK_from_GEE.csv")
ffstation <- read.csv("data/climate/WeatherData_ALL_forest_fire_research_station.csv")

enso <- read.table("data/climate/ENSO_index_meiv2.data", header = F, skip = 1)

# calculate vpd using ffstation data
ffstation <- ffstation %>%
    mutate(
        vpdmax = 0.6108 * exp(17.27 * TempMax / (TempMax + 237.3)) * (1 - (Relative.Humidity / 100)),
        vpdmin = 0.6108 * exp(17.27 * TempMin / (TempMin + 237.3)) * (1 - (Relative.Humidity / 100))
    )

ffstation_monthly <- ffstation %>%
    group_by(YearII, Month) %>%
    dplyr::summarise(
        precipitation = sum(Precipitation),
        dry_days = sum(Precipitation == 0),
        tmax = max(TempMax),
        vpdmax = mean(vpdmax),
        vpdmin = mean(vpdmin)
    ) %>%
    rename(year = YearII, month = Month) %>%
    dplyr::mutate(year = as.character(year)) %>%
    select(year, month, precipitation, dry_days, tmax, vpdmax, vpdmin)

ffstation_lt <- ffstation_monthly %>%
    filter(year <= 2021, year >= 2009) %>%
    ungroup() %>%
    group_by(month) %>%
    dplyr::mutate(
        precipitation = mean(precipitation, na.rm = T),
        dry_days = mean(dry_days, na.rm = T),
        tmax = mean(tmax, na.rm = T),
        vpdmax = mean(vpdmax, na.rm = T),
        vpdmin = mean(vpdmin, na.rm = T),
        year = "long-term"
    ) %>%
    select(year, month, precipitation, dry_days, tmax, vpdmax, vpdmin)

# bind these two dataframes
ffstation_full <- bind_rows(ffstation_monthly, ffstation_lt) %>%
    filter(year %in% c("long-term", "2010", "2015"))


# spei

spei_month <- spei %>%
    dplyr::mutate(
        Date = as.Date(system.time_start, format = "%B %d, %Y"),
        year = year(Date),
        month = month(Date)
    ) %>%
    rename(spei_val = SPEI_01_month) %>%
    select(year, month, spei_val)

spei_lt <- spei_month %>%
    filter(year <= 2021, year >= 2009) %>%
    ungroup() %>%
    group_by(month) %>%
    dplyr::mutate(
        spei_val = mean(spei_val, na.rm = T),
        year = "long-term"
    ) %>%
    select(year, month, spei_val)

spei_month <- spei_month %>%
    dplyr::mutate(year = as.character(year)) %>%
    filter(year %in% c("long-term", "2010", "2015"))

spei_full <- bind_rows(spei_month, spei_lt)

# merge data
clim <- merge(ffstation_full, spei_full, by = c("year", "month"))

clim <- clim %>%
    pivot_longer(
        cols = c(precipitation, dry_days, tmax, vpdmax, vpdmin, spei_val),
        names_to = "climvar",
        values_to = "value"
    ) %>%
    mutate(variable = factor(climvar, levels = c("precipitation", "dry_days", "tmax", "vpdmax", "vpdmin", "spei_val")))


# plot

varnames <- as_labeller(c(
    precipitation = "Precipitation (mm)",
    dry_days = "Dry days",
    tmax = "Max temperature (°C)",
    vpdmax = "VPD max (kPa)",
    vpdmin = "VPD min (kPa)",
    spei_val = "SPEI"
))

climplot <- ggplot() +
    geom_line(data = clim, aes(x = month, y = value, color = year, linetype = year), linewidth = 1.5) +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames) +
    theme_bw() +
    scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # make every month show up on the x-axis
    scale_x_continuous(breaks = 1:12) +
    # add rectangles for dry season
    geom_rect(
        data = data.frame(xmin = c(11, 1), xmax = c(12, 4), ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
    ) +
    guides(linetype = "none") +
    labs(x = "Month", y = "Value", color = "Year") +
    scale_color_manual(values = c("long-term" = "grey40", "2010" = "indianred2", "2015" = "indianred4"))

climplot

png("doc/display/Fig1.png", width = 8, height = 4, units = "in", res = 300)
climplot
dev.off()

# figure 2 - growth increments ENSO plot + sensitivity raw distributions

# Load required libraries---------------------
library(tidyverse)

# load data--------------------------
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
        median_inc_scaled = scale(median_inc, center = TRUE, scale = TRUE),
        avg_inc_tree_scaled = scale(avg_inc_tree, center = TRUE, scale = TRUE)
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
    #    dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 3 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 3 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
    dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 4 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 4 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
    filter(!is.na(sens.prop) & !is.na(cii_min1) & !is.na(calcDBH_min1) & !is.na(twi)) %>%
    ungroup()

# plot the sensitivity values for all trees in 2010 and 2015
yrs <- c(2010, 2015)

sens.all <- ggplot(data = tree.time %>% filter(yr %in% yrs), aes(x = sens.prop)) +
    geom_density(alpha = 0.5) +
    scale_fill_viridis_d() +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
    facet_wrap(~yr) +
    # xlim(-5, 5)+
    labs(title = "Distribution of drought sensitivities \nfor all trees", x = "Sensitivity", y = "Density") +
    theme_bw()

# png("doc/display/sens_all.png", width = 6, height = 4, units = "in", res = 300)
# sens.all
# dev.off()

# run the model for the top 10 species for 2015

top_10_sp <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)

spagplot_top10 <- ggplot() +
    # species plots
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(
            x = yr, y = median_inc,
            group = Species, col = Species
        )
    ) +
    # mean of all trees
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 2
    ) +
    # add points of these
    geom_point(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 3
    ) +
    # mean of species
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
    # make all years show on x-axis
    scale_x_continuous(breaks = unique(tree.time$yr)) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(x = c(2010, 2015), y = 0.55, label = "ENSO drought"), hjust = 0.8, vjust = -0.2, angle = 90) +
    guides(col = guide_legend("species")) +
    xlab("year") +
    ylab("annualised diameter increment (cm)") +
    # ggtitle("growth increments for top 10 species") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# png("doc/display/spaghetti_top10_new.png", width = 4, height = 4, units = "in", res = 300)
# spagplot_top10
# dev.off()

library(gridExtra)
png("doc/display/Fig2.png", width = 4, height = 6, units = "in", res = 300)
grid.arrange(spagplot_top10, sens.all, nrow = 2, heights = c(1.2, 1))
dev.off()
