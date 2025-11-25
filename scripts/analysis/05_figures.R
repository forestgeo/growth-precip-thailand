# script to code for figures

# load libraries ---------------------------
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidybayes)

# figure 1 ------------------------------------------------

# Fig 1
# climate data

# read data

spei <- read.csv("data/climate/SPEI_HKK_from_GEE.csv")

# how to change the order of layers to plot geom_rect first
# https://stackoverflow.com/questions/75007050/manually-adjusting-the-z-order-of-geom-layers

# spei

spei_month <- spei %>%
    dplyr::mutate(
        Date = as.Date(system.time_start, format = "%B %d, %Y"),
        year = year(Date),
        year = as.character(year),
        month = month(Date)
    ) %>%
    filter(year %in% c("2010", "2015", "2020")) %>%
    pivot_longer(cols = c(SPEI_01_month, SPEI_03_month, SPEI_06_month, SPEI_12_month), names_to = "spei_var", values_to = "spei_val") %>%
    select(year, month, spei_var, spei_val)


# figure 1 alternate with 3 years chirps and era5land-----------------------

varnames <- as_labeller(c(
    SPEI_01_month = "1 month",
    SPEI_03_month = "3 month",
    SPEI_06_month = "6 month",
    SPEI_12_month = "12 month"
))


bounds <- spei_month %>%
    # pivot_wider(names_from = year, values_from = spei_val) %>%
    mutate(
        # ymax = pmax(0, spei_val),
        ymax = -1,
        ymin = pmin(-1, spei_val),
        fill = ymin < -1
    )

bounds

bounds2 <- spei_month %>%
    dplyr::mutate(
        x = month,
        y = pmin(-1, spei_val)
    )

bounds2

df <- data.frame(
    x = rep(9.5, 4), y = c(0.75, -0.7, -1.5, -3),
    label = c("wetter", "drier", "drought", "severe\ndrought")
)

cols <- viridis::viridis(4, option = "F")


# spei plot with lines
speiplot_line <- ggplot() +
    # geom_rect(data = data.frame(xmin = c(10.5, 0.5), xmax = c(12.5, 4.5), ymin = -Inf, ymax = Inf), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3) +
    geom_rect(
        data = data.frame(xmin = 4.5, xmax = 10.5, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    # geom_bar(data = spei_month %>% filter(year != "long-term"), aes(x = month, y = spei_val, fill = year), position = "dodge", stat = "identity") +
    geom_line(data = spei_month, aes(
        x = month, y = spei_val,
        #  linetype = spei_var,
        color = year
    ), linewidth = 1.5) +
    # add text for wet, dry drought regions
    # annotate("text",
    #     x = c(12, 12, 12, 12), y = c(0.75, -0.75, -1.5, -2.2),
    #     label = c("wetter", "drier", "drought", "severe drought"),
    #     col = "black", hjust = -0.5, vjust = -1, fontface = "italic", size = 4
    # ) +
    # geom_text(inherit.aes = F, aes(
    #     x = c(10, 10, 10, 10), y = c(0.75, -0.75, -1.5, -2.2),
    #     label = c("wetter", "drier", "drought", "severe drought"), fontface = "italic"
    # )) +
    geom_text(
        data = df, aes(x = x, y = y, label = label),
        col = "black", hjust = 0, vjust = 0, fontface = "italic", size = 3
    ) +
    ggtitle("SPEI") +
    # scale_color_manual(values = c("2010" = "indianred2", "2015" = "indianred4")) +
    # viridis::scale_color_viridis(discrete = T, option="H") +
    scale_color_manual(values = cols[1:3]) +
    facet_wrap(~spei_var,
        # scales = "free_y",
        strip.position = "left", nrow = 4, labeller = varnames
    ) +
    ylim(-3, 1.5) +
    theme_bw() +
    labs(x = "Month", y = "SPEI", fill = "Year") +
    scale_x_continuous(breaks = 1:12) +
    guides(linetype = "none") +
    guides(color = "none") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    ) +
    geom_hline(yintercept = c(0, -1, -2), linetype = "dashed") +
    guides(fill = "none")

speiplot_line


# read chirps data
chirps <- read.csv("data/climate/CHIRPS_DAILY_HKK.csv")
era5land <- read.csv("data/climate/ERA5Land_Daily_HKK.csv")


# bind chirps and era5land data
clim_sat <- merge(chirps, era5land, by = "date")
# make rolling means for dry days and precipitation
climsat_rlmean <- clim_sat %>%
    dplyr::rename(
        tmax = temperature_2m_max,
        tmean = temperature_2m
    ) %>%
    dplyr::mutate(
        dry_day = ifelse(precipitation == 0, 1, 0),
        date = as.Date(date, format = "%Y-%m-%d"),
        year = year(date),
        month = month(date),
        day = day(date),
        vpdmax = 0.6108 * exp((17.27 * (tmax - 273.3)) / (tmax)) * (1 - Rh_100),
    ) %>%
    select(date, year, month, day, precipitation, dry_day, tmax, vpdmax) %>%
    pivot_longer(cols = c(precipitation, dry_day, tmax, vpdmax), names_to = "climvar", values_to = "value") %>%
    select(year, month, day, date, climvar, value) %>%
    arrange(date, climvar) %>%
    group_by(climvar) %>%
    dplyr::mutate(
        rlsum = zoo::rollsum(value, k = 30, fill = NA, align = "center"),
        rlmean = zoo::rollmean(value, k = 30, fill = NA, align = "center"),
        rlval = ifelse(climvar == c("precipitation"), rlsum, rlmean)
    ) %>%
    # compute anomalies
    group_by(climvar, month, day) %>%
    dplyr::mutate(
        long.term = mean(rlval, na.rm = TRUE),
        long.term.sd = sd(rlval, na.rm = TRUE),
        anomaly = (rlval - long.term) / long.term.sd
    ) %>%
    dplyr::mutate( # day of year
        doy = as.numeric(format(date, "%j"))
    )

head(climsat_rlmean)

varnames <- as_labeller(c(
    precipitation = "Precipitation (mm)",
    dry_day = "Dry days",
    tmax = "Max temperature (K)",
    vpdmax = "VPD max (kPa)"
))

# plot values

climsat_plot <- ggplot(climsat_rlmean %>% filter(year %in% c("2010", "2015", "2020")), col = year) +
    geom_rect(
        data = data.frame(xmin = 121, xmax = 304, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    # geom_bar(aes(x = month, y = value, fill = factor(year)), position = "dodge", stat = "identity") +
    geom_line(aes(x = doy, y = rlval, color = factor(year)), linewidth = 0.8) +
    geom_line(aes(x = doy, y = long.term), col = "grey40", lty = "dashed") +
    geom_ribbon(aes(x = doy, ymin = long.term - long.term.sd, ymax = long.term + long.term.sd), col = "grey60", alpha = 0.2) +
    facet_wrap(~climvar,
        scales = "free_y", ncol = 1,
        strip.position = "left",
        labeller = varnames
    ) +
    theme_bw() +
    labs(x = "Day of year", y = "", color = "Year") +
    # scale_x_continuous(breaks = 1:12) +
    guides(linetype = "none", color = "none") +
    ggtitle("Climate variables") +
    # scale_color_manual(values = c("2010" = "indianred2", "2015" = "indianred4")) +
    # viridis::scale_color_viridis(discrete = T) +
    scale_color_manual(values = cols[1:3]) +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )

climsat_plot

climsat_anomaly_plot <- ggplot(climsat_rlmean %>% filter(year %in% c("2010", "2015", "2020")), col = year) +
    geom_rect(
        data = data.frame(xmin = 121, xmax = 304, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    # geom_bar(aes(x = month, y = anomaly, fill = factor(year)), position = "dodge", stat = "identity") +
    geom_line(aes(x = doy, y = anomaly, color = factor(year)), linewidth = 0.8) +
    facet_wrap(~climvar,
        scales = "free_y",
        strip.position = "left", ncol = 1,
        labeller = varnames
    ) +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Day of year", y = "", color = "Year") +
    # scale_x_continuous(breaks = 1:12) +
    guides(linetype = "none") +
    ggtitle("Anomalies") +
    # viridis::scale_color_viridis(discrete = T) +
    # scale_color_manual(values = c("2010" = "indianred2", "2015" = "indianred4")) +
    scale_color_manual(values = cols[1:3]) +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )
climsat_anomaly_plot
speiplot_line


png("doc/display/climsat_plot.png", width = 6, height = 8, units = "in", res = 300)
climsat_plot + climsat_anomaly_plot
dev.off()

layout <- "
ABC
ABC
"

png("doc/display/Fig1.png", width = 8, height = 8, units = "in", res = 300)
climsat_plot + climsat_anomaly_plot + speiplot_line + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

# SI figures with remote data----------------------

varnames <- as_labeller(c(
    precipitation = "Precipitation (mm)",
    dry_day = "Dry days",
    tmax = "Max temperature (°C)",
    vpdmax = "VPD max (kPa)",
    vpdmin = "VPD min (kPa)",
    vpdmean = "VPD mean (kPa)",
    spei_val = "SPEI"
))


climplot_SI <- ggplot(
    climsat_rlmean %>%
        filter(climvar %in% c("precipitation", "tmax", "dry_day", "vpdmax")) %>%
        filter(year >= 2009)
) +
    annotate(
        geom = "rect",
        xmin = as.Date("2015-01-01"), xmax = as.Date("2016-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[2], alpha = 0.5, col = cols[2]
    ) +
    annotate(
        geom = "rect",
        xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[1], alpha = 0.5, col = cols[1]
    ) +
    annotate(
        geom = "rect",
        xmin = as.Date("2020-01-01"), xmax = as.Date("2021-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[3], alpha = 0.5, col = cols[3]
    ) +
    geom_line(aes(x = date, y = rlval), size = 0.15) +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames, strip.position = "left", nrow = 4) +
    # annotate("rect",
    #     fill = "indianred2", alpha = 0.5,
    #     xmin = 2010, xmax = 2011, ymin = -Inf, ymax = Inf
    # ) +
    theme_bw() +
    ylab("") +
    xlab("Date") +
    # show all years on x-axis
    scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )

climplot_SI

png("doc/display/climvars_fullseries_remote.png", width = 8, height = 6, units = "in", res = 300)
climplot_SI
dev.off()

# plot anomalies

climplot_anomalies_SI <- ggplot(
    climsat_rlmean %>%
        filter(climvar %in% c("precipitation", "tmax", "dry_day", "vpdmax"), year >= 2009) %>%
        filter(year != "long-term.sd" & year != "long-term")
) +
    # add rectangles for 2010 and 2015
    # geom_rect(
    #     aes(xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf),
    #     fill = cols[1], alpha = 0.5, col = cols[1]
    # ) +
    annotate(
        geom = "rect",
        xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[1], alpha = 0.5, col = cols[1]
    ) +
    annotate(
        geom = "rect",
        xmin = as.Date("2015-01-01"), xmax = as.Date("2016-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[2], alpha = 0.5, col = cols[2]
    ) +
    annotate(
        geom = "rect",
        xmin = as.Date("2020-01-01"), xmax = as.Date("2021-01-01"), ymin = -Inf, ymax = Inf,
        fill = cols[3], alpha = 0.5, col = cols[3]
    ) +
    geom_line(aes(x = date, y = anomaly), size = 0.15) +
    facet_wrap(~climvar,
        scales = "free_y",
        labeller = varnames, strip.position = "left", nrow = 4
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    ylab("") +
    xlab("Date") +
    # show all years on x-axis
    scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )

climplot_anomalies_SI

png("doc/display/climvars_anomalies_fullseries_remote.png", width = 8, height = 6, units = "in", res = 300)
climplot_anomalies_SI
dev.off()

# figure 2 - growth increments ENSO plot + sensitivity raw distributions-------------------------

# load data--------------------------
tree.time <- read.csv("data/dendro/sensitivity_dataset.csv")
median_incs <- read.csv("data/dendro/summaries_dataset.csv")

# plot the sensitivity values for all trees in 2010 and 2015
yrs <- c(2010, 2015, 2020)

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

cols <- viridis::viridis(4, option = "F")

sens.all.2 <- ggplot(data = tree.time %>% filter(yr %in% yrs), aes(x = sens.prop)) +
    geom_density(alpha = 0.5, aes(col = factor(yr), fill = factor(yr))) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
    # facet_wrap(~yr) +
    facet_wrap(~ factor(yr), ncol = 3) +
    # xlim(-5, 5)+
    labs(x = "Drought sensitivity", y = "Density") +
    theme_bw() +
    guides(color = "none", fill = "none") +
    # make the strip title bold
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )


sens.all.2

sens.all.3 <- ggplot(data = tree.time %>% filter(yr %in% yrs), aes(x = sens.prop)) +
    geom_density(alpha = 0.5, aes(col = factor(yr), fill = factor(yr))) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
    # facet_wrap(~factor(yr), ncol = 1) +
    # xlim(-5, 5)+
    labs(x = "Drought sensitivity", y = "Density") +
    theme_bw() +
    # guides(color="none", fill="none")+
    # make the strip title bold
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
    )

sens.all.3

png("doc/display/sens_all_1panel.png", width = 4, height = 4, units = "in", res = 300)
sens.all.3
dev.off()

# get the top 10 species for 2015

top_10_sp <- tree.time %>%
    filter(yr == 2015) %>%
    group_by(Species, spfull) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)

# top_10_sp$Species
# top_10_sp$spname <- c("Miliusa horsfieldii", "Hopea odorata", "Polyalthia viridis", "Tetrameles nudiflora", "Vatica harmandiana", "Alphonsea ventriculosa", "Garcinia speciosa", "Dipterocarpus alatus", "Baccaurea ramiflora", "Garuga pinnata")

tree.time$spname <- top_10_sp$spfull[match(tree.time$Species, top_10_sp$Species)]

spagplot_nosp <- ggplot() +
    # mean of all trees
    # geom_line(
    #     data = tree.time %>%
    #         # filter(Species %in% top_10_sp$Species) %>%
    #         group_by(yr) %>%
    #         dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
    #     aes(x = yr, y = median_inc), col = "black", size = 2
    # ) +
    geom_line(
        data = median_incs, aes(x = yr, y = median_inc_mm), col = "black", size = 2
    ) +
    # add points of these
    # geom_point(
    #     data = tree.time %>%
    #         group_by(yr) %>%
    #         dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
    #     aes(x = yr, y = median_inc), col = "black", size = 3
    # )
    geom_point(data = median_incs, aes(x = yr, y = median_inc_mm), col = "black", size = 3) +
    # mean of species
    # geom_line(
    #     data = tree.time %>%
    #         group_by(yr, spname) %>%
    #         dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)) %>%
    #         ungroup() %>%
    #         group_by(yr) %>% dplyr::summarise(median_inc = median(median_inc, na.rm = T)),
    #     aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    # ) +
    geom_line(data = median_incs, aes(x = yr, y = median_inc_sp_mm), col = "grey40", size = 0.8) +
    # make all years show on x-axis
    scale_x_continuous(limits = c(2009, 2022), breaks = 2009:2022) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015, 2020), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(
        x = c(2010, 2015, 2020), y = 3.5,
        label = "drought",
    ), color = cols[1:3], hjust = 0.8, vjust = -0.3, angle = 90) +
    # guides(col = guide_legend("species"), nrow = 3) +
    guides(col = "none") +
    xlab("year") +
    ylab("annualised diameter increment (mm)") +
    # ggtitle("growth increments for top 10 species") +
    # xlim(2009, 2024) +
    theme_bw() +
    # theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        # legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

spagplot_nosp


layout <- "
AAABB
"

layout <- "
AAA
AAA
BBB
"

png("doc/display/Fig2_alt2.png", width = 8, height = 8, units = "in", res = 300)
spagplot_nosp + sens.all.2 + plot_layout(design = layout) + plot_annotation(tag_levels = "a")
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
intercepts <- rep(intercepts, each = nrow(ranef_df) / 3)
ranef_df <- ranef_df %>%
    dplyr::mutate(
        intercept = intercepts + median,
        lwr = intercepts + lwr,
        upr = intercepts + upr
    )

# plot intercepts against species characteristics

sp_vars <- tree.time %>%
    group_by(Species, spfull) %>%
    dplyr::summarise(williams_dec = first(williams_dec))

# join ranef_df with sp_vars
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

yrs <- c(2010, 2015, 2020)

# lms
coefs_dec_lms <- ranef_df %>%
    filter(yr %in% yrs) %>%
    filter(Species != "ALPHVE") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(cor.test(data$williams_dec, data$intercept))) %>%
    # dplyr::mutate(mod = list(cor.test(data$maxDBH, data$intercept))) %>%
    dplyr::reframe(broom::tidy(mod))

coefs_dec_lms
colnames(ranef_df)

# make plot with just deciduousness
dec_intercept_plot <- ggplot(ranef_df, aes(x = williams_dec, y = intercept, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    geom_smooth(
        data = ranef_df %>% filter(yr %in% c(2015, 2020)),
        method = "lm", col = "grey40"
    ) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(name ~ factor(yr, levels = c(2010, 2015)), scales = "free_x") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015, 2020)),
        # scales = "free",
        ncol = 3
    ) +
    # add text with p value
    geom_text(
        # data = coefs_dec_lms %>% filter(term == "williams_dec"),
        data = coefs_dec_lms,
        inherit.aes = F,
        aes(
            x = c(2, 2, 2), y = c(0.75, 0.75, 0.75),
            label = paste("r = ", round(estimate, 2), "\np = ", round(p.value, 2)), hjust = 0, vjust = 0
        )
    ) +
    labs(x = "Deciduousness", y = "Species sensitivity") +
    guides(color = "none") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )
dec_intercept_plot


# Fig 3 panel with TWI/deciduousness difference--------------

new_preds <- readRDS("results/models/non_negative/new_preds_isoclines_twi.RDS")

new_preds_df <- do.call(rbind, new_preds)
# coefs_df <- do.call(rbind, coefs)

iso_plot_twi <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = twi, fill = Estimate)
) +
    geom_tile() +
    labs(
        x = "Deciduousness",
        fill = "Sensitivity"
    ) +
    ylab("Topographic Wetness Index\n drier \u2194 wetter") +
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

iso_plot_twi

png("doc/display/Fig3.png", width = 8, height = 8, units = "in", res = 300)
dec_intercept_plot / iso_plot_twi + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.8, 1))
dev.off()


# SI figure with points and deciduousness--------------------

pred_df <- readRDS("results/models/non_negative/predictions_intercept.RData")
str(pred_df)

fits <- readRDS("results/models/non_negative/fits_interceptonly.RDS")
res1 <- residuals(fits[[1]])
res2 <- residuals(fits[[2]])
res3 <- residuals(fits[[3]])

residuals <- rbind(res1, res2, res3)

# bind residuals to preds
pred_res <- cbind(pred_df, residuals)

# plot residuals with variables

# first make long df
pred_res_long <- pred_res %>%
    pivot_longer(cols = c(cii_min1, calcDBH_min1_scaled, twi_scaled, williams_dec)) %>%
    dplyr::select(Tag, treeID, yr, sens.prop, median, lwr, upr, Estimate, Est.Error, Q2.5, Q97.5, name, value) %>%
    dplyr::mutate(name = ifelse(name == "calcDBH_min1_scaled", "scaled DBH",
        ifelse(name == "cii_min1", "CII",
            ifelse(name == "twi_scaled", "scaled TWI", "deciduousness")
        )
    ))

res_cors <- pred_res_long %>%
    nest_by(yr, name) %>%
    dplyr::mutate(mod = list(lm(Estimate ~ value, data = data))) %>%
    dplyr::reframe(broom::tidy(mod)) %>%
    filter(term == "value", p.value < 0.05)


# plot these

plot_res <- ggplot(pred_res_long, aes(x = value, y = Estimate)) +
    geom_point(alpha = 0.3) +
    geom_smooth(
        data = pred_res_long %>% filter(paste0(yr, name) %in% paste0(res_cors$yr, res_cors$name)),
        col = "grey40", method = "lm"
    ) +
    ylab("Residual") +
    facet_grid(name ~ yr) +
    theme_bw()

png("doc/display/residuals.png", width = 6, height = 8, units = "in", res = 300)
plot_res
dev.off()

# Figure 3 panel with TPI -----------------------------------------

# coefs plot
new_preds <- readRDS("results/models/non_negative/new_preds_isoclines_tpi.RDS")

new_preds_df <- do.call(rbind, new_preds)
# coefs_df <- do.call(rbind, coefs)

iso_plot_tpi <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = -tpi, fill = Estimate)
) +
    geom_tile() +
    labs(x = "Deciduousness", y = "-Topographic Position Index\n wetter-->", fill = "Sensitivity") +
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
# iso_plot

layout <- "
A
B
C
"

png("doc/display/FigSI_TPI_preds.png", width = 6, height = 2.5, units = "in", res = 300)
iso_plot_tpi
dev.off()

# check models
coefs_df_twi <- readRDS("results/models/non_negative/coefs_isoclines_twi.RDS")
coefs_df_tpi <- readRDS("results/models/non_negative/coefs_isoclines_tpi.RDS")

coefs_df_twi <- do.call(rbind, coefs_df_twi)
coefs_df_tpi <- do.call(rbind, coefs_df_tpi)

coefs_df_twi$sig <- ifelse(coefs_df_twi$lwr > 0 | coefs_df_twi$upr < 0, "sig", "not sig")
coefs_df_tpi$sig <- ifelse(coefs_df_tpi$lwr > 0 | coefs_df_tpi$upr < 0, "sig", "not sig")

coefs_df_twi_tpi <- readRDS("results/models/non_negative/coefs_isoclines_twi_tpi.RDS")
coefs_df_twi_tpi <- do.call(rbind, coefs_df_twi_tpi)
coefs_df_twi_tpi$sig <- ifelse(coefs_df_twi_tpi$lwr > 0 | coefs_df_twi_tpi$upr < 0, "sig", "not sig")




# supplementary plot of other values ------------------------------------

# lms
coefs_dec_lms <- ranef_df %>%
    filter(yr %in% yrs) %>%
    # filter(Species != "ALPHVE", Species != "MACASI") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(lm(median ~ williams_dec, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))

coefs_dec_lms

# SI figure - species all species sensitivities---------------

# order of species abundance
tree.time$Species <- factor(tree.time$Species, levels = table(tree.time$Species) %>% sort(decreasing = T) %>% names())

# plot of median sensitivities across species for 2010 and 2015
sensplot <- ggplot(data = tree.time %>% filter(yr %in% c(2010, 2015, 2020)), aes(x = Species, y = sens.prop, color = factor(yr))) +
    geom_point(position = position_jitterdodge(), alpha = 0.5) +
    geom_boxplot(fill = NA) +
    # scale_color_viridis_d() +
    # scale_color_manual(values = c("2010" = "indianred2", "2015" = "indianred4")) +
    scale_color_manual(values = cols[1:3]) +
    ylim(-5, 5) +
    labs(
        # title = "Distribution of drought sensitivities \nfor all species",
        x = "Species", y = "Sensitivity"
    ) +
    theme_bw() +
    guides(color = guide_legend(title = "Year")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("doc/display/Fig_SI_species_sensitivities.png", width = 10, height = 4, units = "in", res = 300)
sensplot
dev.off()


# figure 4 ----------------------------------
# read coefs
# coefs_df <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")
coefs <- readRDS("results/models/orderedcii/coefs_rel_spre.rds")
coefs_df <- do.call(rbind, coefs)

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
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
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
    facet_wrap(~ factor(yr, levels = c(2010, 2015, 2020)),
        # scales = "free",
        ncol = 3
    ) +
    labs(x = "Scaled Topographic Wetness Index\n drier \u2194 wetter", y = "Predicted sensitivity") +
    guides(color = guide_legend(title = "Deciduousness")) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

pred_plot

str(preds)

layout <- "
A
B
"

png("doc/display/Fig4.png", width = 8, height = 8, units = "in", res = 300)
coefs_all_sp + pred_plot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

# alternate figure 4 with TPI--------------
coefs <- readRDS("results/models/orderedcii_tpi/coefs_rel_spre.rds")
coefs_df <- do.call(rbind, coefs)

# par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

par_names <- as_labeller(c("b_sensprop_calcDBH_min1_scaled" = "DBH effect", "bsp_sensprop_mocii_min1" = "exposure effect", "b_sensprop_tpi_scaled" = "wetness effect"))

coefs_all_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_tpi_scaled")),
    aes(
        x = factor(param, levels = c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_tpi_scaled")),
        y = median,
        ymin = lwr, ymax = upr
        # col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
    )
) +
    geom_pointrange() +
    scale_x_discrete(labels = par_names) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12)) +
    coord_flip()


coefs_all_sp

preds <- readRDS("results/models/orderedcii_tpi/pred_sens_rel_spre.rds")
preds_df <- do.call(rbind, preds)

brbg.5 <- RColorBrewer::brewer.pal(5, "BrBG")
brbg.3 <- RColorBrewer::brewer.pal(3, "BrBG")

# plot predictions
# pred_plot <- ggplot(data = preds_df, aes(x = twi_scaled_sp, y = median)) +
pred_plot <- ggplot(data = preds_df, aes(x = tpi_scaled, y = median)) +
    geom_smooth(aes(group = Species, col = williams_dec), method = "lm", alpha = 0.2, se = F) +
    geom_smooth(method = "lm", col = "black", linewidth = 2) +
    # geom_point(aes(x = twi_scaled, y = sens.prop), alpha = 0.1) +
    scale_color_gradient(low = brbg.5[5], high = brbg.5[1]) +
    # scale_color_gradient(low = brbg.3[3], high = brbg.7[1]) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015, 2020)),
        # scales = "free",
        ncol = 3
    ) +
    labs(x = "Topographic Position Index", y = "Predicted sensitivity") +
    guides(color = guide_legend(title = "Deciduousness")) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

pred_plot

layout <- "
A
B
"

png("doc/display/Fig4_tpi_3yrs.png", width = 8, height = 8, units = "in", res = 300)
coefs_all_sp + pred_plot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

# SI figure with random intercept only----------------

coefs <- readRDS("results/models/orderedcii/coefs_tree.rds")
coefs_df <- do.call(rbind, coefs)

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
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    theme(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 12)) +
    coord_flip()

coefs_all_sp

preds <- readRDS("results/models/orderedcii/pred_sens_tree.rds")
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
    facet_wrap(~ factor(yr, levels = c(2010, 2015, 2020)),
        # scales = "free",
        ncol = 3
    ) +
    labs(x = "Scaled Topographic Wetness Index", y = "Predicted sensitivity") +
    guides(color = guide_legend(title = "Deciduousness")) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

pred_plot

layout <- "
AA
BB
"


png("doc/display/FigSI_coefs_nospre.png", width = 8, height = 8, units = "in", res = 300)
coefs_all_sp + pred_plot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()


# fig 5-----------------------------------

# plot conditional effects of cii

# https://discourse.mc-stan.org/t/conditional-effects-plot-for-monotonic-predictors-with-logit-link/35850/2

# read fits
# fits <- readRDS("results/models/orderedcii/fits_partialmed_rel.rds")
fits <- readRDS("results/models/orderedcii/fits_rel_spre.rds")


# get_variables(fits[[1]])

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

cii_fit3 <- fits[[3]] %>%
    spread_draws(
        b_sensprop_Intercept, bsp_sensprop_mocii_min1,
        simo_sensprop_mocii_min11[i]
    ) %>%
    dplyr::mutate(yr = 2020)

cii_fit <- bind_rows(cii_fit1, cii_fit2, cii_fit3)

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
    theme_bw() +
    theme(strip.background = element_blank())

p_manual_ce

# DAGS with values

dag_2010 <- magick::image_read("doc/display/dag_2010_vals.png")
dag_2015 <- magick::image_read("doc/display/dag_2015_vals.png")
dag_2020 <- magick::image_read("doc/display/dag_2020_vals.png")

dag_2010_gg <- magick::image_ggplot(dag_2010, interpolate = F)
dag_2015_gg <- magick::image_ggplot(dag_2015, interpolate = F)
dag_2020_gg <- magick::image_ggplot(dag_2020, interpolate = F)

layout <- "
AABBCC
AABBCC
AABBCC
DDDDDD
DDDDDD
"

png("doc/display/Fig5.png", width = 8, height = 6, units = "in", res = 300)
dag_2010_gg + dag_2015_gg + dag_2020_gg + p_manual_ce + plot_layout(design = layout) + plot_annotation(tag_levels = "a")
dev.off()

# cors between species median sensitivities ----------

# tree.time.sp <- tree.time %>%
#     group_by(williams_dec, spfull) %>%
#     dplyr::summarise(
#         inc_annual_med = round(median(inc_annual, na.rm = T), 2),
#         inc_annual_sd = round(sd(inc_annual, na.rm = T), 2)
#     ) %>%
#     arrange(inc_annual_med)

sens.sp <- tree.time %>%
    filter(yr %in% c(2010, 2015, 2020)) %>%
    group_by(williams_dec, spfull, yr) %>%
    dplyr::summarise(
        sens_med = round(median(sens.prop, na.rm = T), 2),
        sens_mean = round(mean(sens.prop, na.rm = T), 2),
        sens_sd = round(sd(sens.prop, na.rm = T), 2)
    ) %>%
    arrange(sens_med)

# spread this for each year

sens.sp.wide <- sens.sp %>%
    ungroup() %>%
    dplyr::select(williams_dec, spfull, yr, sens_med) %>%
    pivot_wider(names_from = yr, values_from = sens_med, names_prefix = "sens.") %>%
    dplyr::mutate(dec = ifelse(williams_dec >= 2, "dec", "evg")) %>%
    dplyr::mutate(dec = ifelse(is.na(dec), "evg", dec))


# paired correlation of 2010 and 2015 sensitivities across species

sens.sp.cor_10_15 <- cor.test(sens.sp.wide$sens.2010, sens.sp.wide$sens.2015)
sens.sp.cor_15_20 <- cor.test(sens.sp.wide$sens.2015, sens.sp.wide$sens.2020)
sens.sp.cor_10_20 <- cor.test(sens.sp.wide$sens.2010, sens.sp.wide$sens.2020)

# plot sensivity cor by species group

brbg.5 <- RColorBrewer::brewer.pal(5, "BrBG")
brbg.3 <- RColorBrewer::brewer.pal(3, "BrBG")

sens.sp_plot_1 <- ggplot(sens.sp.wide, aes(x = sens.2010, y = sens.2015, group = dec)) +
    geom_point(aes(col = williams_dec)) +
    # geom_smooth(aes(col = dec), method = "lm") +
    scale_color_gradient(low = brbg.5[5], high = brbg.5[1]) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    # add species names to points with the name starting at the point
    geom_text(aes(label = spfull, col = williams_dec), hjust = 0, vjust = 0) +
    # add correlation values
    geom_text(aes(
        label = paste0("r = ", round(sens.sp.cor_10_15$estimate, 2), "\np = ", round(sens.sp.cor_10_15$p.value, 2)),
        x = -0.7, y = 0
    )) +
    xlim(c(-1, 1)) +
    guides(color = guide_legend(title = "Deciduousness")) +
    labs(x = "species drought sensitivity in 2010", y = "species drought sensitivity in 2015") +
    theme_bw()

sens_sp_plot_2 <- ggplot(sens.sp.wide, aes(x = sens.2015, y = sens.2020, group = dec)) +
    geom_point(aes(col = williams_dec)) +
    # geom_smooth(aes(col = dec), method = "lm") +
    scale_color_gradient(low = brbg.5[5], high = brbg.5[1]) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    # add species names to points with the name starting at the point
    geom_text(aes(label = spfull, col = williams_dec), hjust = 0, vjust = 0) +
    # add correlation values
    geom_text(aes(
        label = paste0("r = ", round(sens.sp.cor_15_20$estimate, 2), "\np = ", round(sens.sp.cor_15_20$p.value, 2)),
        x = -0.7, y = 0.5
    )) +
    xlim(c(-1, 0.7)) +
    guides(color = guide_legend(title = "Deciduousness")) +
    labs(x = "species drought sensitivity in 2015", y = "species drought sensitivity in 2020") +
    theme_bw()

sens_sp_plot_3 <- ggplot(sens.sp.wide, aes(x = sens.2010, y = sens.2020, group = dec)) +
    geom_point(aes(col = williams_dec)) +
    # geom_smooth(aes(col = dec), method = "lm") +
    scale_color_gradient(low = brbg.5[5], high = brbg.5[1]) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    # add species names to points with the name starting at the point
    geom_text(aes(label = spfull, col = williams_dec), hjust = 0, vjust = 0) +
    # add correlation values
    geom_text(aes(
        label = paste0("r = ", round(sens.sp.cor_10_20$estimate, 2), "\np = ", round(sens.sp.cor_10_20$p.value, 2)),
        x = -0.7, y = 0.5
    )) +
    xlim(c(-1, 1.2)) +
    guides(color = guide_legend(title = "Deciduousness")) +
    labs(x = "species drought sensitivity in 2010", y = "species drought sensitivity in 2020") +
    theme_bw()

sens.sp.wide$williams_dec[sens.sp.wide$spfull == "Alphonsea ventricosa"] <- 1

# plot difference between 2010 and 2015 as a range plot

# reorder species by deciduousness
sens.sp.wide$spfull <- factor(sens.sp.wide$spfull, levels = sens.sp.wide$spfull[order(sens.sp.wide$williams_dec)])

sens_rangeplot <- ggplot(sens.sp.wide, aes(x = spfull)) +
    # geom_linerange(aes(ymin = sens.2010, ymax = sens.2015), alpha = 0.5) +
    geom_point(aes(y = sens.2010), col = cols[1]) +
    geom_point(aes(y = sens.2015), col = cols[2]) +
    geom_point(aes(y = sens.2020), col = cols[3]) +
    coord_flip() +
    xlab("species in order of deciduousness") +
    ylab("species drought sensitivity across droughts") +
    theme_bw()

sens_rangeplot

layout <- "
AB
CD
"

png("doc/display/Fig_SI_sensitivity_cor_dec.png", width = 16, height = 16, units = "in", res = 300)
sens.sp_plot_1 + sens_sp_plot_2 + sens_sp_plot_3 + sens_rangeplot + plot_layout(design = layout)
dev.off()

## SI plot - correlation between 2010 and 2015

forcor <- tree.time %>%
    # filter(yr %in% c(2010, 2015, 2020)) %>%
    dplyr::select(Tag, Species, yr, sens.prop) %>%
    group_by(Tag) %>%
    pivot_wider(
        names_from = yr,
        names_prefix = "sens_",
        values_from = sens.prop
    )

# how many complete observations
nrow(forcor[complete.cases(forcor), ])

pval_1_2 <- cor.test(forcor$sens_2010, forcor$sens_2015, use = "pairwise.complete.obs")[3]$p.value
rval_1_2 <- cor.test(forcor$sens_2010, forcor$sens_2015, use = "pairwise.complete.obs")[4]$estimate
pval_1_3 <- cor.test(forcor$sens_2010, forcor$sens_2020, use = "pairwise.complete.obs")[3]$p.value
rval_1_3 <- cor.test(forcor$sens_2010, forcor$sens_2020, use = "pairwise.complete.obs")[4]$estimate
pval_2_3 <- cor.test(forcor$sens_2015, forcor$sens_2020, use = "pairwise.complete.obs")[3]$p.value
rval_2_3 <- cor.test(forcor$sens_2015, forcor$sens_2020, use = "pairwise.complete.obs")[4]$estimate


# plot the correlation
corplot_1 <- ggplot(forcor, aes(x = sens_2010, y = sens_2015)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey40") +
    # geom_smooth(method="lm", col="grey20", linetype)+
    xlab("Sensitivity in 2010") +
    ylab("Sensitivity in 2015") +
    # add text with r value
    geom_text(
        data = data.frame(x = -3, y = 3, label = paste0("R = ", round(rval_1_2, 2), "\np = ", round(pval_1_2, 3))),
        aes(x = x, y = y, label = label), hjust = 0.5, vjust = 0.5
    ) +
    theme_bw()

corplot_2 <- ggplot(forcor, aes(x = sens_2010, y = sens_2020)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey20") +
    geom_smooth(method = "lm", col = "grey20") +
    xlab("Sensitivity in 2010") +
    ylab("Sensitivity in 2020") +
    # add text with r value
    geom_text(
        data = data.frame(x = -3, y = 3, label = paste0("R = ", round(rval_1_3, 2), "\np = ", round(pval_1_3, 3))),
        aes(x = x, y = y, label = label), hjust = 0.5, vjust = 0.5
    ) +
    theme_bw()

corplot_3 <- ggplot(forcor, aes(x = sens_2015, y = sens_2020)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey20") +
    geom_smooth(method = "lm", col = "grey20") +
    xlab("Sensitivity in 2015") +
    ylab("Sensitivity in 2020") +
    # add text with r value
    geom_text(
        data = data.frame(x = -3, y = 3, label = paste0("R = ", round(rval_2_3, 2), "\np = ", round(pval_2_3, 3))),
        aes(x = x, y = y, label = label), hjust = 0.5, vjust = 0.5
    ) +
    theme_bw()

png("doc/display/sens_corr.png", width = 12, height = 4, units = "in", res = 300)
corplot_1 + corplot_2 + corplot_3 + plot_layout(ncol = 3)
dev.off()

# SI figure with distribution of deciduousness------------------
sp_vars <- tree.time %>%
    group_by(Species) %>%
    dplyr::summarise(williams_dec = first(williams_dec)) %>%
    ungroup()

# make a histogram of williams_dec

dec_hist <- ggplot(sp_vars, aes(x = williams_dec)) +
    geom_histogram(bins = 10, alpha = 0.7) +
    ylab("Number of species") +
    xlab("Deciduousness") +
    scale_y_continuous(breaks = 1:10) +
    theme_bw()

png("doc/display/Fig_dec_dist.png", width = 4, height = 4, units = "in", res = 300)
dec_hist
dev.off()

# conditional dependencies-------------------------
cond_dep_all <- tree.time %>%
    # filter(yr %in% yrs) %>%
    group_by(yr) %>%
    pivot_longer(cols = c("calcDBH_min1_scaled", "cii_min1"), names_to = "var", values_to = "value") %>%
    dplyr::mutate(varnames = ifelse(var == "calcDBH_min1_scaled", "DBH", "CII"))

# these variables are conditionally independent for the most part

# plot the conditional independencies

cond_dep_all_plot <- ggscatter(
    data = cond_dep_all,
    x = "value", y = "twi_scaled",
    add = "reg.line", conf.int = TRUE, alpha = 0.3,
    cor.coef = TRUE, cor.method = "spearman",
    xlab = "variable", ylab = "scaled TWI"
) +
    facet_wrap(varnames ~ yr, scales = "free")

cond_dep_all_plot

png("doc/display/cond_dep_alltrees.png", width = 8, height = 8, units = "in", res = 300)
cond_dep_all_plot
dev.off()

# SI conditional dependencies -----------------------

## test conditional independences for each species-----------------------------

# DBH | TWI
# CII | TWI

# 2015 data only

cond_dep <- tree.time %>%
    filter(yr == 2015) %>%
    group_by(Species) %>%
    dplyr::summarise(
        DBH_TWI = cor.test(calcDBH_min1_scaled, twi_scaled)[4]$estimate,
        DBH_TWI_p = cor.test(calcDBH_min1_scaled, twi_scaled)[3]$p.value,
        CII_TWI = cor.test(cii_min1, twi_scaled)[4]$estimate,
        CII_TWI_p = cor.test(cii_min1, twi_scaled)[3]$p.value
    )

# these variables are conditionally independent for the most part

# plot the conditional independencies

cond_dep_dbh_twi <- ggscatter(
    data = tree.time %>% filter(yr == 2015),
    x = "calcDBH_min1_scaled", y = "twi_scaled",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "spearman",
    xlab = "DBH scaled", ylab = "TWI scaled"
) +
    facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

png("doc/display/cond_dep_dbh_twi.png", width = 12, height = 12, units = "in", res = 300)
cond_dep_dbh_twi
dev.off()

# CII | TWI
cond_dep_cii_twi <- ggscatter(
    data = tree.time %>% filter(yr == 2015),
    x = "cii_min1", y = "twi_scaled",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "pearson",
    xlab = "CII", ylab = "TWI scaled"
) +
    facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

png("doc/display/cond_dep_cii_twi.png", width = 12, height = 12, units = "in", res = 300)
cond_dep_cii_twi
dev.off()

# raw distribution plots--------------------------

# plot of sensitivity against CII
sens_cii <- ggplot(
    tree.time,
    aes(x = factor(cii_min1), y = sens.prop)
) +
    geom_jitter(alpha = 0.2) +
    geom_hline(yintercept = 0, col = "grey20") +
    geom_boxplot(alpha = 0.5) +
    # geom_smooth(method = "lm", col = "grey40") +
    facet_wrap(~yr) +
    xlab("CII") +
    ylab("Sensitivity") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

sens_cii

sens_dbh_cor <- tree.time %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(cor.test(data$calcDBH_min1_scaled, data$sens.prop))) %>%
    dplyr::reframe(broom::tidy(mod))

# make a supsmu df

smu_dbh <- tree.time %>%
    ungroup() %>%
    group_by(yr) %>%
    dplyr::mutate(
        x = stats::supsmu(x = calcDBH_min1_scaled, y = sens.prop)$x,
        y = stats::supsmu(x = calcDBH_min1_scaled, y = sens.prop)$y
    )

sens_dbh <- ggplot(
    tree.time,
    aes(x = calcDBH_min1_scaled, y = sens.prop)
) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        data = tree.time %>% filter(yr == 2010),
        method = "lm", col = "darkblue"
    ) +
    geom_line(data = smu_dbh, aes(x = x, y = y), col = "red") +
    facet_wrap(~yr) +
    xlab("scaled DBH") +
    ylab("Sensitivity") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

sens_dbh

sens_twi_cor <- tree.time %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(cor.test(data$twi_scaled, data$sens.prop))) %>%
    dplyr::reframe(broom::tidy(mod))

smu_twi <- tree.time %>%
    ungroup() %>%
    group_by(yr) %>%
    dplyr::reframe(
        x = stats::supsmu(x = twi_scaled, y = sens.prop)$x,
        y = stats::supsmu(x = twi_scaled, y = sens.prop)$y
    )

smu_twi

sens_twi <- ggplot(
    tree.time,
    aes(x = twi_scaled, y = sens.prop)
) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        data = tree.time %>%
            filter(yr == 2020),
        method = "lm", col = "darkblue"
    ) +
    geom_line(data = smu_twi, aes(x = x, y = y), col = "red") +
    facet_wrap(~yr) +
    xlab("scaled TWI") +
    ylab("Sensitivity") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

sens_twi

png("doc/display/sens_var_rels.png", width = 10, height = 12, units = "in", res = 300)
sens_cii + sens_dbh + sens_twi + plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a")
dev.off()
