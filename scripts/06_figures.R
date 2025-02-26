# script to code for figures

# load libraries
library(tidyverse)
library(lubridate)
library(ggplot2)


# figure 1 ------------------------------------------------

# Fig 1
# climate data

# read data

spei <- read.csv("data/climate/SPEI_HKK_from_GEE.csv")
ffstation <- read.csv("data/climate/WeatherData_ALL_forest_fire_research_station.csv")

enso <- read.table("data/climate/ENSO_index_meiv2.data", header = F, skip = 1)
head(ffstation)
# calculate vpd using ffstation data
ffstation <- ffstation %>%
    mutate(
        vpdmax = 0.6108 * exp(17.27 * TempMax / (TempMax + 237.3)) * (1 - (Relative.Humidity / 100)),
        vpdmin = 0.6108 * exp(17.27 * TempMin / (TempMin + 237.3)) * (1 - (Relative.Humidity / 100)),
        vpdmean = 0.6108 * exp(17.27 * TempMean / (TempMean + 237.3)) * (1 - (Relative.Humidity / 100)),
        dry_days = ifelse(Precipitation == 0, 1, 0)
    )

# make a figure with rolling monthly means
ffstation_monthly_rl <- ffstation %>%
    pivot_longer(cols = c(Precipitation, dry_days, TempMax, vpdmax, vpdmin, vpdmean), names_to = "climvar", values_to = "value") %>%
    select(YearII, Month, Date, climvar, value) %>%
    arrange(YearII, Month, Date, climvar) %>%
    group_by(climvar) %>%
    dplyr::mutate(
        rlsum = zoo::rollsum(value, k = 30, fill = NA, align = "center"),
        rlmean = zoo::rollmean(value, k = 30, fill = NA, align = "center"),
        rlval = ifelse(climvar == c("Precipitation"), rlsum, rlmean)
    ) %>%
    rename(year = YearII, month = Month) %>%
    dplyr::mutate(year = as.character(year)) %>%
    dplyr::mutate(
        Date2 = as.Date(paste0(year, "-", month, "-", Date)),
        day_of_year = as.numeric(format(Date2, "%j"))
    ) %>%
    select(-c(rlsum, rlmean))

ffstation_lt_rl <- ffstation_monthly_rl %>%
    filter(year <= 2021, year >= 2009) %>%
    ungroup() %>%
    group_by(day_of_year, climvar) %>%
    dplyr::summarise(
        rlse = sd(rlval, na.rm = T) / sqrt(sum(!is.na(rlval))),
        rlval = mean(rlval, na.rm = T)
    ) %>%
    ungroup() %>%
    dplyr::mutate(year = "long-term")
select(year, month, Date, climvar, rlval, rlse, day_of_year)

# bind these two dataframes
ffstation_full_rl <- bind_rows(ffstation_monthly_rl, select(ffstation_lt_rl, -("rlse"))) %>%
    filter(year %in% c("long-term", "2010", "2015")) %>%
    filter(climvar %in% c("Precipitation", "dry_days", "TempMax", "vpdmax"))

varnames <- as_labeller(c(
    Precipitation = "Precipitation (mm)",
    dry_days = "Dry days",
    TempMax = "Max temperature (°C)",
    vpdmax = "VPD max (kPa)",
    vpdmin = "VPD min (kPa)",
    vpdmean = "VPD mean (kPa)",
    spei_val = "SPEI"
))


climplot <- ggplot() +
    geom_rect(
        data = data.frame(xmin = c(305, 1), xmax = c(366, 120), ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
    ) +
    geom_line(data = ffstation_full_rl, aes(x = day_of_year, y = rlval, color = year, linetype = year), linewidth = 0.8) +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames, strip.position = "left") +
    theme_bw() +
    scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
    geom_ribbon(data = ffstation_lt_rl %>% filter(climvar %in% c("Precipitation", "dry_days", "TempMax", "vpdmax")), aes(x = day_of_year, ymin = rlval - rlse, ymax = rlval + rlse), alpha = 0.3) +
    # text for wet and dry season
    annotate("text", x = -Inf, y = -Inf, label = "dry season", col = "black", hjust = -0.5, vjust = -1, fontface = "italic", size = 1.15) +
    annotate("text", x = -Inf, y = -Inf, label = "wet season", col = "black", hjust = -2.5, vjust = -1, fontface = "italic", size = 1.15) +
    guides(linetype = "none") +
    labs(x = "Day of year", y = "", color = "Year") +
    scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4")) +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )

# how to change the order of layers to plot geom_rect first
# https://stackoverflow.com/questions/75007050/manually-adjusting-the-z-order-of-geom-layers

climplot

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


# spei plot with bars
speiplot <- ggplot() +
    geom_rect(data = data.frame(xmin = c(10.5, 0.5), xmax = c(12.5, 4.5), ymin = -Inf, ymax = Inf), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3) +
    geom_bar(data = spei_full, aes(x = month, y = spei_val, fill = year), position = "dodge", stat = "identity") +
    scale_fill_manual(values = c("long-term" = "grey40", "2010" = "indianred2", "2015" = "indianred4")) +
    theme_bw() +
    labs(x = "Month", y = "SPEI", fill = "Year") +
    scale_x_continuous(breaks = 1:12) +
    # text for wet and dry season
    annotate("text", x = -Inf, y = -Inf, label = "dry season", col = "black", hjust = -0.5, vjust = -1, fontface = "italic", size = 2) +
    annotate("text", x = -Inf, y = -Inf, label = "wet season", col = "black", hjust = -2.5, vjust = -1, fontface = "italic", size = 2) +
    guides(linetype = "none") +
    geom_hline(yintercept = c(0, -1, -2), linetype = "dashed") +
    guides(fill = "none")

library(patchwork)
png("doc/display/Fig1.png", width = 8.5, height = 4, units = "in", res = 300)
climplot + speiplot + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "a")
dev.off()


# monthly means as supplementary plot-------------
ffstation_monthly <- ffstation %>%
    group_by(YearII, Month) %>%
    dplyr::summarise(
        precipitation = sum(Precipitation),
        dry_days = sum(Precipitation == 0),
        tmax = max(TempMax),
        vpdmax = mean(vpdmax),
        vpdmin = mean(vpdmin),
        vpsmean = mean(vpdmean)
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
    geom_rect(
        data = data.frame(xmin = c(11, 1), xmax = c(12, 4), ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
    ) +
    geom_line(data = clim, aes(x = month, y = value, color = year, linetype = year), linewidth = 1.5) +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames) +
    theme_bw() +
    scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # make every month show up on the x-axis
    scale_x_continuous(breaks = 1:12) +
    # add rectangles for dry season
    guides(linetype = "none") +
    labs(x = "Month", y = "Value", color = "Year") +
    scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4"))

climplot

png("doc/display/env_vars_SI.png", width = 8, height = 4, units = "in", res = 300)
climplot
dev.off()


# figure 2 - growth increments ENSO plot + sensitivity raw distributions-------------------------

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
    # facet_wrap(~yr) +
    facet_wrap(~yr, ncol = 1) +
    # xlim(-5, 5)+
    labs(x = "Drought sensitivity", y = "Density") +
    theme_bw() +
    # make the strip title bold
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

sens.all


# get the top 10 species for 2015

top_10_sp <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species, spfull) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)

# top_10_sp$Species
# top_10_sp$spname <- c("Miliusa horsfieldii", "Hopea odorata", "Polyalthia viridis", "Tetrameles nudiflora", "Vatica harmandiana", "Alphonsea ventriculosa", "Garcinia speciosa", "Dipterocarpus alatus", "Baccaurea ramiflora", "Garuga pinnata")

tree.time$spname <- top_10_sp$spfull[match(tree.time$Species, top_10_sp$Species)]

spagplot_top10 <- ggplot() +
    # species plots
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, spname) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(
            x = yr, y = median_inc,
            group = spname, col = spname
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
            group_by(yr, spname) %>%
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
    guides(col = guide_legend("species"), nrow = 3) +
    xlab("year") +
    ylab("annualised diameter increment (cm)") +
    # ggtitle("growth increments for top 10 species") +
    theme_bw() +
    # theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

spagplot_top10

library(patchwork)
png("doc/display/Fig2.png", width = 10, height = 6, units = "in", res = 300)
spagplot_top10 + sens.all + plot_annotation(tag_levels = "a") +
    # plot_layout(guides = "collect", widths = c(1, 2.2)) & theme(
    plot_layout(guides = "collect") & theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.text = element_text(face = "italic")
)
dev.off()

# figure 3 --------------------------------------------

# related issue - https://github.com/forestgeo/growth-precip-thailand/issues/12


# coefs plot

coefs_df <- readRDS("results/models/non_negative/sensitivity_model_intercept.RData")
pred_df <- readRDS("results/models/non_negative/predictions_intercept.RData")

# plot random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)

# add intercept to ranef_df
# make intercept vector
intercepts <- coefs_df %>%
    filter(param %in% "Intercept") %>%
    dplyr::select(median) %>%
    pull(median)
intercepts <- rep(intercepts, each = nrow(ranef_df) / 2)
ranef_df <- ranef_df %>%
    dplyr::mutate(
        intercept = intercepts + median,
        lwr = intercepts + lwr,
        upr = intercepts + upr
    )

# plot intercepts against species characteristics

# join ranef_df with sp_vars
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

# lms
coefs_dec_lms <- ranef_df %>%
    filter(yr %in% yrs) %>%
    filter(Species != "ALPHVE") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(lm(median ~ williams_dec, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))

coefs_dec_lms

# make plot with just deciduousness
dec_intercept_plot <- ggplot(ranef_df, aes(x = williams_dec, y = intercept, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    geom_smooth(method = "lm", col = "grey40") +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(name ~ factor(yr, levels = c(2010, 2015)), scales = "free_x") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015)),
        # scales = "free",
        ncol = 2
    ) +
    # add text with p value
    geom_text(data = coefs_dec_lms %>% filter(term == "williams_dec"), inherit.aes = F, aes(x = c(3, 3), y = c(1.5, 1), label = paste("p = ", round(p.value, 2)), hjust = 0, vjust = 0)) +
    labs(x = "Deciduousness", y = "Predicted sensitivity") +
    guides(color = "none") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )


# panel with TWI/deciduousness difference

new_preds <- readRDS("results/models/non_negative/new_preds_isoclines_nore.RDS")

new_preds_df <- do.call(rbind, new_preds)
coefs_df <- do.call(rbind, coefs)

iso_plot <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = twi, fill = Estimate)
) +
    geom_tile() +
    labs(x = "Deciduousness", y = "Topographic Wetness Index", fill = "Sensitivity") +
    # geom_contour(aes(z = Estimate), colour = "black") +
    # facet_grid(yr ~ Species) +
    facet_wrap(~yr) +
    theme_bw() +
    scale_fill_gradient2() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )
iso_plot


library(patchwork)
png("doc/display/Fig3.png", width = 8, height = 8, units = "in", res = 300)
dec_intercept_plot / iso_plot + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.8, 1))
dev.off()



# supplementary plot of other values ------------------------------------


# make long df for plotting
ranef_df_long <- ranef_df %>%
    pivot_longer(c("twi_sd", "maxDBH", "williams_dec")) %>%
    # rename the variables
    dplyr::mutate(name = case_when(
        # name == "twi_median" ~ "median TWI",
        name == "twi_sd" ~ "sd(TWI)",
        name == "maxDBH" ~ "maximum DBH",
        name == "williams_dec" ~ "deciduousness"
    ))


# lms
coefs_dec_lms <- ranef_df_long %>%
    filter(yr %in% yrs) %>%
    # filter(Species != "ALPHVE") %>%
    nest_by(yr, name) %>%
    dplyr::mutate(mod = list(lm(median ~ value, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))

coefs_dec_lms

# make plot with just deciduousness
dec_intercept_plot <- ggplot(ranef_df_long, aes(x = value, y = intercept)) +
    geom_point() +
    geom_smooth(method = "lm", col = "grey40") +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(name ~ factor(yr, levels = c(2010, 2015)), scales = "free_x") +
    facet_wrap(name ~ factor(yr, levels = c(2010, 2015)),
        scales = "free", ncol = 2, strip.position = "left"
    ) +
    # add text with p value
    geom_text(data = coefs_dec_lms %>% filter(term == "value"), aes(x = 3, y = 1, label = paste("p = ", round(p.value, 2)), hjust = 0, vjust = 0)) +
    labs(x = "species trait value", y = "intercept") +
    guides(color = "none") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )


library(patchwork)
png("doc/display/Fig_trait_intercept_plot.png", width = 6, height = 6, units = "in", res = 300)
# sensplot + sp_intercept_plot + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1.8, 1))
dec_intercept_plot
dev.off()

# SI figure - species all species sensitivities---------------

# order of species abundance
tree.time$Species <- factor(tree.time$Species, levels = table(tree.time$Species) %>% sort(decreasing = T) %>% names())

# plot of median sensitivities across species for 2010 and 2015
sensplot <- ggplot(data = tree.time %>% filter(yr %in% c(2010, 2015)), aes(x = Species, y = sens.prop, color = factor(yr))) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    geom_boxplot(fill = NA) +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("2010" = "indianred2", "2015" = "indianred4")) +
    labs(
        # title = "Distribution of drought sensitivities \nfor all species",
        x = "Species", y = "Sensitivity"
    ) +
    theme_bw() +
    guides(color = guide_legend(title = "Year")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("doc/display/Fig_SI_species_sensitivities.png", width = 8, height = 4, units = "in", res = 300)
sensplot
dev.off()


# # fig 4 - DAG + coefs --------------------------------------------

# dag_img <- magick::image_read("doc/display/hypotheses.png")

# print(dag_img)
# dag_img

# dag_gg <- magick::image_ggplot(dag_img, interpolate = F)

# # read coefs
# coefs_df <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")

# # plot coefs
# colours <- c("#e15f41", "#546de5", "#f7b731")

# par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

# coefs_all_sp <- ggplot(
#     data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
#     aes(
#         x = param, y = median,
#         # col = factor(signif, levels = c("neg", "pos", "no"))
#         col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
#     )
# ) +
#     geom_point() +
#     scale_x_discrete(labels = par_names) +
#     # make error bars with narrow heads
#     geom_errorbar(aes(
#         ymin = lwr, ymax = upr,
#         # col = factor(signif, levels = c("neg", "pos", "no"))
#         col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
#     ), width = 0.1) +
#     # scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
#     scale_color_manual(values = rep(colours, 2), drop = FALSE) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     # facet_grid(~param, scales = "free", labeller = par_names) +
#     facet_grid(~ factor(yr, levels = c(2010, 2015)), scales = "free") +
#     labs(title = "", x = "", y = "coefficient") +
#     guides(color = "none") +
#     theme_bw() +
#     coord_flip()

# coefs_all_sp
# library(patchwork)



# png("doc/display/Fig4.png", width = 6, height = 4, units = "in", res = 300)
# (dag_gg / coefs_all_sp) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(2, 1))
# dev.off()


# figure 4 ----------------------------------
# read coefs
# coefs_df <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")
coefs <- readRDS("results/models/orderedcii/coefs_rel_spre.rds")
coefs_df <- do.call(rbind, coefs)

# par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

par_names <- as_labeller(c("b_sensprop_calcDBH_min1_scaled" = "DBH effect", "bsp_sensprop_mocii_min1" = "exposure effect", "b_sensprop_twi_scaled" = "wetness effect"))

coefs_all_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_twi_scaled")),
    aes(
        x = factor(param, levels = c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_twi_scaled")),
        y = median,
        ymin = lwr, ymax = upr
        # col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
    )
) +
    geom_pointrange() +
    scale_x_discrete(labels = par_names) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015)), scales = "free") +
    labs(title = "", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12)) +
    coord_flip()


coefs_all_sp

# plot of TWI slopes
coefs_twi <- coefs_df %>%
    # filter(grepl("b_twi_scaled_sp", param))
    filter(grepl("b_sensprop_twi_scaled", param))

head(coefs_twi)

coefs_sp_twi <- coefs_df %>%
    filter(grepl("r_Species", param)) %>%
    filter(grepl("twi", param)) %>%
    filter(grepl("cii|DBH|Intercept", param) == FALSE) %>%
    dplyr::mutate(
        Species = gsub("r_Species__sensprop\\[|\\,twi_scaled\\]", "", param)
    ) %>%
    group_by(yr) %>%
    dplyr::mutate(
        total_effect = median + coefs_twi$median[match(yr, coefs_twi$yr)],
        total_lwr = lwr + coefs_twi$median[match(yr, coefs_twi$yr)],
        total_upr = upr + coefs_twi$median[match(yr, coefs_twi$yr)]
    )


# merge species vars
coefs_sp_twi <- merge(coefs_sp_twi, sp_vars, by = "Species", all.x = TRUE)

# cors
twi_cor <- coefs_sp_twi %>%
    group_by(yr) %>%
    dplyr::summarise(
        cor = cor.test(williams_dec, total_effect)[4]$estimate,
        cor_p = cor.test(williams_dec, total_effect)[3]$p.value
    )

twi_cor
# plot coefs

library(ggpubr)
twi_slopes_plot <- ggplot(
    data = coefs_sp_twi,
    # aes(x = reorder(Species, median), y = median),
    aes(
        x = williams_dec, y = total_effect,
        ymin = total_lwr, ymax = total_upr
    ),
    order = median
) +
    # geom_point(size = 3, alpha = 0.7) +
    # # geom_smooth(method = "lm", col = "grey40") +
    # geom_errorbar(aes(ymin = total_lwr, ymax = total_upr),
    #     width = 0.005, linewidth = 1, alpha = 0.7
    # ) +
    geom_pointrange(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015)),
        # scales = "free",
        ncol = 2
    ) +
    geom_smooth(
        data = coefs_sp_twi %>% filter(yr == 2015),
        method = "lm", col = "grey40"
    ) +
    stat_cor(method = "pearson", label.x = 0.5, label.y = 0.5) +
    labs(x = "Deciduousness", y = "wetness effect") +
    # coord_flip()+
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

twi_slopes_plot

# preds_df <- readRDS("results/models/non_negative/predictions_spre.RData")

preds <- readRDS("results/models/orderedcii/pred_sens_rel_spre.rds")
preds_df <- do.call(rbind, preds)

brbg.5 <- RColorBrewer::brewer.pal(5, "BrBG")
brbg.3 <- RColorBrewer::brewer.pal(3, "BrBG")

# plot predictions
# pred_plot <- ggplot(data = preds_df, aes(x = twi_scaled_sp, y = median)) +
pred_plot <- ggplot(data = preds_df, aes(x = twi_scaled, y = median)) +
    geom_smooth(aes(group = Species, col = williams_dec), method = "lm", alpha = 0.2, se = F) +
    geom_smooth(method = "lm", col = "black", linewidth = 2) +
    # geom_point(aes(x = twi_scaled, y = sens.prop), alpha = 0.1) +
    scale_color_gradient(low = brbg.5[5], high = brbg.5[1]) +
    # scale_color_gradient(low = brbg.3[3], high = brbg.7[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015)),
        # scales = "free",
        ncol = 2
    ) +
    labs(x = "Topographic Wetness Index", y = "Predicted sensitivity") +
    guides(color = guide_legend(title = "Deciduousness")) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

pred_plot

str(preds)

library(patchwork)
png("doc/display/Fig4.png", width = 16, height = 8, units = "in", res = 300)
(coefs_all_sp / twi_slopes_plot) | pred_plot # & plot_annotation(tag_levels = "a") & plot_layout(widths = c(1, 1.5))
dev.off()

# fig 5-----------------------------------

fullmed_img <- magick::image_read("doc/display/full_mediation.png")
fullmed_gg <- magick::image_ggplot(fullmed_img, interpolate = F)

parmed_img <- magick::image_read("doc/display/hypotheses_grey.png")
parmed_gg <- magick::image_ggplot(parmed_img, interpolate = F)

# coefficient plots

coefs_fullmed <- readRDS(paste0(models_dir, "/coefs_fullmediation.rds"))
coefs_parmed <- readRDS(paste0(models_dir, "/coefs_partialmed.rds"))

coefs_fullmed <- do.call(rbind, coefs_fullmed)
coefs_parmed <- do.call(rbind, coefs_parmed)

pars.keep <- c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_twi_scaled")

par_names <- as_labeller(c(
    "b_sensprop_calcDBH_min1_scaled" = "DBH effect",
    "bsp_sensprop_mocii_min1" = "CII effect",
    "b_sensprop_twi_scaled" = "TWI effect"
))

# plot the coefficients
coefs_plot_fullmed <- ggplot(coefs_fullmed %>% filter(param %in% pars.keep), aes(x = param, y = median, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    facet_wrap(~yr) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    scale_x_discrete(labels = par_names) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Coefficients for the full mediation model",
        x = "Parameter",
        y = "Coefficient"
    ) +
    coord_flip()

coefs_plot_fullmed

coefs_plot_parmed <- ggplot(coefs_parmed %>% filter(param %in% pars.keep), aes(x = param, y = median, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    facet_wrap(~yr) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    scale_x_discrete(labels = par_names) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Coefficients for the partial mediation model",
        x = "Parameter",
        y = "Coefficient"
    ) +
    coord_flip()

png("doc/display/Fig5.png", width = 12, height = 8, units = "in", res = 300)
(coefs_plot_parmed + parmed_gg) / (coefs_plot_fullmed + fullmed_gg)
dev.off()

# plot conditional effects of cii

# https://discourse.mc-stan.org/t/conditional-effects-plot-for-monotonic-predictors-with-logit-link/35850/2

# read fits
# fits <- readRDS("results/models/orderedcii/fits_partialmed_rel.rds")
fits <- readRDS("results/models/orderedcii/fits_rel_spre.rds")

library(tidybayes)
get_variables(fits[[1]])

cii_fit1 <- fits[[1]] %>%
    spread_draws(
        b_sensprop_Intercept, bsp_sensprop_mocii_min1,
        simo_sensprop_mocii_min11[i]
    ) %>%
    dplyr::mutate(yr = 2010)

cii_fit2 <- fits[[2]] %>%
    spread_draws(
        b_sensprop_Intercept, bsp_sensprop_mocii_min1,
        simo_sensprop_mocii_min11[i]
    ) %>%
    dplyr::mutate(yr = 2015)

cii_fit <- bind_rows(cii_fit1, cii_fit2)

cii_fit <- cii_fit %>%
    group_by(yr) %>%
    dplyr::mutate(
        # D is equal to number of categories minus 1
        D = length(unique(i) - 1)
    ) %>%
    group_by(.chain, .iteration) %>%
    # add a row within each of these groups

    dplyr::mutate(
        cumsumi = cumsum(simo_sensprop_mocii_min11),
        post_mu = b_sensprop_Intercept + (bsp_sensprop_mocii_min1 * D * cumsumi)
    )

# add a row per chain and iteration with i = 0
cii_fit_add <- cii_fit %>%
    group_by(yr, .chain, .iteration, .draw) %>%
    dplyr::summarise(
        b_sensprop_Intercept = b_sensprop_Intercept[1],
        bsp_sensprop_mocii_min1 = bsp_sensprop_mocii_min1[1],
        i = 0,
        simo_sensprop_mocii_min11 = 0,
        yr = yr[1],
        cumsumi = 0,
        post_mu = b_sensprop_Intercept
    )

cii_fit <- bind_rows(cii_fit, cii_fit_add)


# for 2015, the categories are 2, 3, 4, 5
cii_fit <- cii_fit %>%
    mutate(i = ifelse(yr == 2015, i + 1, i))


p_manual_ce <- ggplot(data = cii_fit %>% filter(i < 5), aes(y = post_mu, x = factor(i))) +
    scale_x_discrete(labels = c("1", "2", "3", "4", "5")) +
    stat_halfeye() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~yr) +
    labs(x = "CII", y = "Sensitivity") +
    theme_bw()

png("doc/display/Fig5_alternate.png", width = 12, height = 8, units = "in", res = 300)
(coefs_plot_parmed + parmed_gg) / p_manual_ce
dev.off()


# DAGS with values

dag_2010 <- magick::image_read("doc/display/dag_2010_vals.png")
dag_2015 <- magick::image_read("doc/display/dag_2015_vals.png")

dag_2010_gg <- magick::image_ggplot(dag_2010, interpolate = F)
dag_2015_gg <- magick::image_ggplot(dag_2015, interpolate = F)

png("doc/display/Fig5_alternate2.png", width = 8, height = 8, units = "in", res = 300)
(dag_2010_gg + dag_2015_gg) / p_manual_ce
dev.off()


# supplementary figure of CII effects by species----------------

cii_fit_sp1 <- fits[[1]] %>%
    spread_draws(
        b_sensprop_Intercept, bsp_sensprop_mocii_min1,
        simo_sensprop_mocii_min11[i],
        r_Species__sensprop[condition, mocii_min1]
    ) %>%
    dplyr::mutate(yr = 2010)

cii_fit_sp2 <- fits[[2]] %>%
    spread_draws(
        b_sensprop_Intercept, bsp_sensprop_mocii_min1,
        simo_sensprop_mocii_min11[i],
        r_Species__sensprop[condition, mocii_min1]
    ) %>%
    dplyr::mutate(yr = 2015)

cii_fit_sp <- bind_rows(cii_fit_sp1, cii_fit_sp2)

cii_fit_sp <- cii_fit_sp %>%
    group_by(yr) %>%
    dplyr::mutate(
        # D is equal to number of categories minus 1
        D = length(unique(i) - 1)
    ) %>%
    group_by(.chain, .iteration, condition) %>%
    # add a row within each of these groups

    dplyr::mutate(
        cumsumi = cumsum(simo_sensprop_mocii_min11),
        post_mu = b_sensprop_Intercept + (bsp_sensprop_mocii_min1 * D * cumsumi)
    )

# TODO - where to include r_Species__sensprop in this calculation?
