# script to code for figures

# load libraries ---------------------------
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)

# figure 1 ------------------------------------------------

# Fig 1
# climate data

# read data

spei <- read.csv("data/climate/SPEI_HKK_from_GEE.csv")
ffstation <- read.csv("data/climate/WeatherData_ALL_forest_fire_research_station.csv")

# enso <- read.table("data/climate/ENSO_index_meiv2.data", header = F, skip = 1)
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

# ffstation_lt_rl <- ffstation_monthly_rl %>%
#     filter(year <= 2021, year >= 2009) %>%
#     ungroup() %>%
#     group_by(day_of_year, climvar) %>%
#     dplyr::summarise(
#         # rlse = sd(rlval, na.rm = T) / sqrt(sum(!is.na(rlval))),
#         rlval = mean(rlval, na.rm = T)
#     ) %>%
#     ungroup() %>%
#     dplyr::mutate(year = "long-term")
# # select(year, month, Date, climvar, rlval, rlse, day_of_year)

ffstation_lt_rl <- ffstation_monthly_rl %>%
    filter(year <= 2021, year >= 2009) %>%
    ungroup() %>%
    group_by(day_of_year, climvar) %>%
    dplyr::summarise(
        rlsd = sd(rlval, na.rm = T),
        rlval = mean(rlval, na.rm = T)
    ) %>%
    ungroup() %>%
    dplyr::mutate(year = "long-term")

ffstation_lt_rl_long <- ffstation_lt_rl %>%
    pivot_longer(cols = c(rlval, rlsd), names_to = "stat", values_to = "rlval") %>%
    dplyr::mutate(year = ifelse(stat == "rlsd", "long-term.sd", year)) %>%
    select(-stat)

# bind these two dataframes
ffstation_full_rl <- bind_rows(ffstation_monthly_rl, ffstation_lt_rl) %>%
    filter(year %in% c("long-term", "2010", "2015")) %>%
    filter(climvar %in% c("Precipitation", "dry_days", "TempMax", "vpdmax"))

ffstation_full_rl_long <- bind_rows(ffstation_monthly_rl, ffstation_lt_rl_long) %>%
    # filter(year %in% c("long-term", "long-term.sd", "2010", "2015")) %>%
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
    # geom_rect(
    #     data = data.frame(xmin = c(305, 1), xmax = c(366, 120), ymin = -Inf, ymax = Inf),
    #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
    # ) +
    geom_rect(
        data = data.frame(xmin = 121, xmax = 304, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    geom_line(data = ffstation_full_rl, aes(x = day_of_year, y = rlval, color = year, linetype = year), linewidth = 0.8) +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames, strip.position = "left", nrow = 4) +
    theme_bw() +
    scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
    geom_ribbon(data = ffstation_lt_rl %>% filter(climvar %in% c("Precipitation", "dry_days", "TempMax", "vpdmax")), aes(x = day_of_year, ymin = rlval - rlsd, ymax = rlval + rlsd), alpha = 0.3) +
    # # text for wet and dry season
    # annotate("text", x = -Inf, y = -Inf, label = "dry season", col = "black", hjust = -0.5, vjust = -1, fontface = "italic", size = 1.15) +
    # annotate("text", x = -Inf, y = -Inf, label = "wet season", col = "black", hjust = -2.5, vjust = -1, fontface = "italic", size = 1.15) +
    guides(linetype = "none") +
    ggtitle("Climate variables") +
    labs(x = "Day of year", y = "", color = "Year") +
    scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4")) +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside"
    )

# how to change the order of layers to plot geom_rect first
# https://stackoverflow.com/questions/75007050/manually-adjusting-the-z-order-of-geom-layers

climplot

# calculate anomalies
ffstation_anomalies_full <- ffstation_full_rl_long %>%
    group_by(climvar, day_of_year) %>%
    dplyr::mutate(
        long.term = rlval[year == "long-term"],
        long.term.sd = rlval[year == "long-term.sd"],
        anomaly = (rlval - long.term) / long.term.sd
    ) %>%
    ungroup()

ffstation_anomalies <- ffstation_anomalies_full %>%
    filter(year %in% c("2010", "2015"))

# plot anomalies

climplot_anom <- ggplot() +
    # geom_rect(
    #     data = data.frame(xmin = c(305, 1), xmax = c(366, 120), ymin = -Inf, ymax = Inf),
    #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
    # ) +
    geom_rect(
        data = data.frame(xmin = 121, xmax = 304, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    geom_line(data = ffstation_anomalies %>% filter(year %in% c("2010", "2015", "long-term")), aes(x = day_of_year, y = anomaly, color = year), linewidth = 0.8) +
    # geom_bar(data = ffstation_anomalies, aes(x = day_of_year, y = anomaly, fill = year), stat = "identity", position = "dodge") +
    facet_wrap(~climvar, scales = "free_y", labeller = varnames, strip.position = "left", nrow = 4) +
    theme_bw() +
    ggtitle("Anomalies") +
    # scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
    # # text for wet and dry season
    # annotate("text", x = -Inf, y = -Inf, label = "dry season", col = "black", hjust = -0.5, vjust = -1, fontface = "italic", size = 1.15) +
    # annotate("text", x = -Inf, y = -Inf, label = "wet season", col = "black", hjust = -2.5, vjust = -1, fontface = "italic", size = 1.15) +
    guides(linetype = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Day of year", y = "", color = "Year") +
    # scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4")) +
    scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4")) +
    guides(color = "none") +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_blank()
    )

climplot_anom


# spei

spei_month <- spei %>%
    dplyr::mutate(
        Date = as.Date(system.time_start, format = "%B %d, %Y"),
        year = year(Date),
        year = as.character(year),
        month = month(Date)
    ) %>%
    # filter(year %in% c("2010", "2015")) %>%
    filter(year %in% c("2010", "2015", "2020")) %>%
    pivot_longer(cols = c(SPEI_01_month, SPEI_03_month, SPEI_06_month, SPEI_12_month), names_to = "spei_var", values_to = "spei_val") %>%
    select(year, month, spei_var, spei_val)

# spei_lt <- spei_month %>%
#     filter(year <= 2021, year >= 2009) %>%
#     ungroup() %>%
#     group_by(month) %>%
#     dplyr::mutate(
#         spei_val = mean(spei_val, na.rm = T),
#         year = "long-term"
#     ) %>%
#     select(year, month, spei_val)

# spei_month <- spei_month %>%
#     dplyr::mutate(year = as.character(year)) %>%
#     filter(year %in% c("long-term", "2010", "2015"))

# spei_full <- bind_rows(spei_month, spei_lt)


# spei plot with bars
speiplot <- ggplot() +
    # geom_rect(data = data.frame(xmin = c(10.5, 0.5), xmax = c(12.5, 4.5), ymin = -Inf, ymax = Inf), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3) +
    geom_rect(
        data = data.frame(xmin = 4.5, xmax = 10.5, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    geom_bar(data = spei_month %>% filter(year %in% c("2010", "2015")), aes(x = month, y = spei_val, fill = year), position = "dodge", stat = "identity") +
    scale_fill_manual(values = c("long-term" = "grey40", "2010" = "indianred2", "2015" = "indianred4")) +
    theme_bw() +
    labs(x = "Month", y = "SPEI", fill = "Year") +
    scale_x_continuous(breaks = 1:12) +
    guides(linetype = "none") +
    geom_hline(yintercept = c(0, -1, -2), linetype = "dashed") +
    guides(fill = "none")



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
    x = rep(9.5, 4), y = c(0.75, -0.7, -1.5, -2.5),
    label = c("wetter", "drier", "drought", "severe drought")
)

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
    scale_color_manual(values = c("long-term" = "grey40", "2010" = "indianred2", "2015" = "indianred4")) +
    # ggh4x::stat_difference(data = spei_month, aes(x = month, ymin = -1, ymax = spei_val), levels = c("above", "below")) +
    # scale_fill_manual(limits = c("above", "below"), values = c("grey40", "grey60")) +
    facet_wrap(~spei_var,
        # scales = "free_y",
        strip.position = "left", nrow = 4, labeller = varnames
    ) +
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


library(patchwork)
# png("doc/display/Fig1.png", width = 8.5, height = 4, units = "in", res = 300)
# climplot + speiplot + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "a")
# dev.off()

# layout <- "
# AAAACC
# BBBBCC
# "

# png("doc/display/Fig1_anom.png", width = 12, height = 4, units = "in", res = 300)
# # wrap_elements(climplot / climplot_anom) + speiplot + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "a")
# climplot + climplot_anom + speiplot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a")
# dev.off()

layout <- "
ABCCC
ABCCC
"

png("doc/display/Fig1.png", width = 10, height = 8, units = "in", res = 300)
climplot + climplot_anom + speiplot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

layout <- "
ABC
ABC
"

png("doc/display/Fig1_new.png", width = 8, height = 8, units = "in", res = 300)
climplot + climplot_anom + speiplot_line + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

## SI figures of climvars and anomalies------

varnames <- as_labeller(c(
    Precipitation = "Precipitation (mm)",
    dry_days = "Dry days",
    TempMax = "Max temperature (°C)",
    vpdmax = "VPD max (kPa)",
    vpdmin = "VPD min (kPa)",
    vpdmean = "VPD mean (kPa)",
    spei_val = "SPEI"
))


climplot_SI <- ggplot(
    ffstation_monthly_rl %>%
        filter(climvar %in% c("Precipitation", "TempMax", "dry_days", "vpdmax")) %>%
        filter(year >= 2009)
) +
    # add rectangles for 2010 and 2015
    # geom_rect(
    #     data = data.frame(xmin = decimal_date(as.Date("2010-01-01")), xmax = decimal_date(as.Date("2011-01-01")), ymin = -Inf, ymax = Inf),
    #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #     # fill = c("indianred2", "indianred4"),
    #     fill = "indianred2",
    #     alpha = 0.3
    # ) +
    geom_rect(
        aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2016-01-01"), ymin = -Inf, ymax = Inf),
        fill = "indianred4", alpha = 0.3, col = "indianred4"
    ) +
    geom_rect(
        aes(xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf),
        fill = "indianred2", alpha = 0.3, col = "indianred2"
    ) +
    geom_line(aes(x = Date2, y = rlval), size = 0.15) +
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

png("doc/display/climvars_fullseries.png", width = 8, height = 6, units = "in", res = 300)
climplot_SI
dev.off()

# plot anomalies

climplot_anomalies_SI <- ggplot(
    ffstation_anomalies_full %>%
        filter(climvar %in% c("Precipitation", "TempMax", "dry_days", "vpdmax"), year >= 2009) %>%
        filter(year != "long-term.sd" & year != "long-term")
) +
    # add rectangles for 2010 and 2015
    geom_rect(
        aes(xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf),
        fill = "indianred2", alpha = 0.3, col = "indianred2"
    ) +
    geom_rect(
        aes(xmin = as.Date("2015-01-01"), xmax = as.Date("2016-01-01"), ymin = -Inf, ymax = Inf),
        fill = "indianred4", alpha = 0.3, col = "indianred4"
    ) +
    geom_line(aes(x = Date2, y = anomaly), size = 0.15) +
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

png("doc/display/climvars_anomalies_fullseries.png", width = 8, height = 6, units = "in", res = 300)
climplot_anomalies_SI
dev.off()

# # monthly means as supplementary plot-------------
# ffstation_monthly <- ffstation %>%
#     group_by(YearII, Month) %>%
#     dplyr::summarise(
#         precipitation = sum(Precipitation),
#         dry_days = sum(Precipitation == 0),
#         tmax = max(TempMax),
#         vpdmax = mean(vpdmax),
#         vpdmin = mean(vpdmin),
#         vpsmean = mean(vpdmean)
#     ) %>%
#     rename(year = YearII, month = Month) %>%
#     dplyr::mutate(year = as.character(year)) %>%
#     select(year, month, precipitation, dry_days, tmax, vpdmax, vpdmin)

# ffstation_lt <- ffstation_monthly %>%
#     filter(year <= 2021, year >= 2009) %>%
#     ungroup() %>%
#     group_by(month) %>%
#     dplyr::mutate(
#         precipitation = mean(precipitation, na.rm = T),
#         dry_days = mean(dry_days, na.rm = T),
#         tmax = mean(tmax, na.rm = T),
#         vpdmax = mean(vpdmax, na.rm = T),
#         vpdmin = mean(vpdmin, na.rm = T),
#         year = "long-term"
#     ) %>%
#     select(year, month, precipitation, dry_days, tmax, vpdmax, vpdmin)

# # bind these two dataframes
# ffstation_full <- bind_rows(ffstation_monthly, ffstation_lt) %>%
#     filter(year %in% c("long-term", "2010", "2015"))



# # merge data
# clim <- merge(ffstation_full, spei_full, by = c("year", "month"))

# clim <- clim %>%
#     pivot_longer(
#         cols = c(precipitation, dry_days, tmax, vpdmax, vpdmin, spei_val),
#         names_to = "climvar",
#         values_to = "value"
#     ) %>%
#     mutate(variable = factor(climvar, levels = c("precipitation", "dry_days", "tmax", "vpdmax", "vpdmin", "spei_val")))


# # plot

# varnames <- as_labeller(c(
#     precipitation = "Precipitation (mm)",
#     dry_days = "Dry days",
#     tmax = "Max temperature (°C)",
#     vpdmax = "VPD max (kPa)",
#     vpdmin = "VPD min (kPa)",
#     spei_val = "SPEI"
# ))

# climplot <- ggplot() +
#     geom_rect(
#         data = data.frame(xmin = c(11, 1), xmax = c(12, 4), ymin = -Inf, ymax = Inf),
#         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3
#     ) +
#     geom_line(data = clim, aes(x = month, y = value, color = year, linetype = year), linewidth = 1.5) +
#     facet_wrap(~climvar, scales = "free_y", labeller = varnames) +
#     theme_bw() +
#     scale_linetype_manual(values = c("long-term" = "longdash", "2010" = "solid", "2015" = "solid")) +
#     # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     # make every month show up on the x-axis
#     scale_x_continuous(breaks = 1:12) +
#     # add rectangles for dry season
#     guides(linetype = "none") +
#     labs(x = "Month", y = "Value", color = "Year") +
#     scale_color_manual(values = c("long-term" = "grey60", "2010" = "indianred2", "2015" = "indianred4"))

# climplot

# png("doc/display/env_vars_SI.png", width = 8, height = 4, units = "in", res = 300)
# climplot
# dev.off()

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

climsat_plot <- ggplot(climsat_rlmean %>% filter(year %in% c("2010", "2015", "2020", "long-term")), col = year) +
    geom_rect(
        data = data.frame(xmin = 121, xmax = 304, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "lightblue", alpha = 0.3
    ) +
    # geom_bar(aes(x = month, y = value, fill = factor(year)), position = "dodge", stat = "identity") +
    geom_line(aes(x = doy, y = rlval, color = factor(year)), linewidth = 0.8) +
    facet_wrap(~climvar,
        scales = "free_y", ncol = 1,
        strip.position = "left",
        labeller = varnames
    ) +
    theme_bw() +
    labs(x = "Dy of year", y = "", color = "Year") +
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

library(patchwork)
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


# correlation between remote and ground data and anomalies----------------------
head(ffstation_anomalies_full)
head(climsat_rlmean)

ff_tomerge <- ffstation_anomalies_full %>%
    dplyr::rename(date = Date2) %>%
    dplyr::mutate(
        climvar = ifelse(climvar == "Precipitation", "precipitation", climvar),
        climvar = ifelse(climvar == "dry_day", "dry_days", climvar),
        climvar = ifelse(climvar == "TempMax", "tmax", climvar)
    ) %>%
    filter(year %in% c("2010", "2015")) %>%
    select(date, climvar, rlval, anomaly) %>%
    filter(climvar %in% c("precipitation", "dry_days", "tmax", "vpdmax")) %>%
    dplyr::mutate(source = "station")
head(ff_tomerge)
sum(ffstation_anomalies_full$year == "2015")

trial <- ffstation_anomalies_full %>%
    filter(year == "2015")

climsat_tomerge <- climsat_rlmean %>%
    ungroup() %>%
    filter(year %in% c("2010", "2015")) %>%
    select(date, climvar, rlval, anomaly) %>%
    dplyr::mutate(
        climvar = ifelse(climvar == "dry_day", "dry_days", climvar),
        rlval = ifelse(climvar == "tmax", rlval - 273.15, rlval) # convert tmax from K to C
    ) %>%
    dplyr::mutate(source = "remote")

head(climsat_tomerge)
max(ff_tomerge$rlval[ff_tomerge$climvar == "precipitation"], na.rm = TRUE)

nrow(ff_tomerge)
nrow(climsat_tomerge)
clim_merge <- bind_rows(ff_tomerge, climsat_tomerge)

# for each climvar calculate the correlation between station and remote data
clim_corr <- clim_merge %>%
    filter(!is.na(rlval)) %>%
    filter(!is.na(anomaly)) %>%
    group_by(climvar, year(date)) %>%
    dplyr::summarise(
        cor_rlval = cor.test(rlval[source == "station"], rlval[source == "remote"])$estimate,
        cor_rlval_p = cor.test(rlval[source == "station"], rlval[source == "remote"])$p.value,
        cor_anomaly = cor.test(anomaly[source == "station"], anomaly[source == "remote"])$estimate,
        cor_anomaly_p = cor.test(anomaly[source == "station"], anomaly[source == "remote"])$p.value
    ) %>%
    ungroup()

clim_corr

write.csv(clim_corr, "doc/display/clim_corr_remote_station.csv", row.names = FALSE)

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
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(
            x = yr, y = median_inc,
            group = spname, col = spname
        ), alpha = 0.7
    ) +
    # mean of all trees
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 2
    ) +
    # add points of these
    geom_point(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 3
    ) +
    # mean of species
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr, spname) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
    # make all years show on x-axis
    scale_x_continuous(breaks = unique(tree.time$yr)) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015, 2020), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(
        x = c(2010, 2015, 2020), y = 5.5,
        label = c("ENSO drought", "ENSO drought", "drought")
    ), hjust = 0.8, vjust = -0.2, angle = 90) +
    guides(col = guide_legend("species"), nrow = 3) +
    xlab("year") +
    ylab("annualised diameter increment (mm)") +
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

spagplot_top10_names <- ggplot() +
    # species plots
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, spname) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(
            x = yr, y = median_inc,
            group = spname, col = spname
        ), alpha = 0.7
    ) +
    # mean of all trees
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 2
    ) +
    # add points of these
    geom_point(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 3
    ) +
    # mean of species
    geom_line(
        data = tree.time %>%
            # filter(Species %in% top_10_sp$Species) %>%
            group_by(yr, spname) %>%
            dplyr::summarise(median_inc = median(inc_annual * 10, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
    # make all years show on x-axis
    scale_x_continuous(limits = c(2009, 2024), breaks = 2009:2024) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015, 2020), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(
        x = c(2010, 2015, 2020), y = 5.5,
        label = c("ENSO drought", "ENSO drought", "drought")
    ), hjust = 0.8, vjust = -0.2, angle = 90) +
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
    ) +
    # add species names at x = 2022, y = last value, in corresponding colour
    geom_text(
        data = tree.time %>%
            filter(Species %in% top_10_sp$Species) %>%
            group_by(spname) %>%
            arrange(yr) %>%
            summarise(
                y = median(inc_annual[yr == max(yr)] * 10, na.rm = TRUE),
                col = unique(spname)
            ),
        aes(x = 2022, y = y, label = spname, color = col),
        hjust = 0, size = 3, show.legend = FALSE
    )

spagplot_top10_names


layout <- "
AAABB
"

library(patchwork)
# png("doc/display/Fig2.png", width = 10, height = 8, units = "in", res = 300)
png("doc/display/Fig2.png", width = 10, height = 8, units = "in", res = 300)
spagplot_top10 + sens.all + plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect", design = layout) & theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.text = element_text(face = "italic")
)
dev.off()

layout <- "
AAA
AAA
BBB
"

png("doc/display/Fig2_alt.png", width = 8, height = 8, units = "in", res = 300)
spagplot_top10_names + sens.all.2 + plot_layout(design = layout) + plot_annotation(tag_levels = "a")
dev.off()



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

## SI plots - negative growth increments
neg.incs.sp <- tree.time %>%
    # filter(inc_annual < 0) %>%
    group_by(yr, Species) %>%
    dplyr::summarise(
        negs = sum(inc_annual < 0, na.rm = T),
        n = n()
    ) %>%
    ungroup() %>%
    dplyr::mutate(
        Species = factor(Species, levels = unique(Species))
    )

neg.incs <- neg.incs.sp %>%
    group_by(yr) %>%
    dplyr::summarise(
        negs = sum(negs, na.rm = T),
        n = sum(n)
    ) %>%
    ungroup() %>%
    dplyr::mutate(
        Species = "all",
        forcol = ifelse(yr == 2010, "2010",
            ifelse(yr == 2015, "2015",
                ifelse(yr == "2020", "2020", "other")
            )
        ),
        forcol = factor(forcol, levels = c("2010", "2015", "2020", "other"))
    )

head(neg.incs)

neg_incs_plot <- ggplot() +
    geom_bar(
        data = neg.incs,
        aes(x = yr, y = negs / n, fill = forcol),
        stat = "identity"
    ) +
    scale_x_continuous(breaks = unique(tree.time$yr)) +
    # change the colour of 2010 and 2015 bars
    # scale_fill_manual(values = c("indianred2", "indianred4", "grey40")) +
    ylab("Proportion of negative growth increments") +
    xlab("Year") +
    guides(fill = guide_legend("year"), nrow = 1) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
png("doc/display/neg_incs.png", width = 6, height = 4, units = "in", res = 300)
neg_incs_plot
dev.off()

## SI plot - correlation between 2010 and 2015

forcor <- tree.time %>%
    filter(yr %in% c(2010, 2015, 2020)) %>%
    select(Tag, Species, yr, sens.prop) %>%
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


# vals<- cor.test(forcor$sens_2010, forcor$sens_2015, use = "pairwise.complete.obs")
# vals[4]

# plot the correlation
corplot_1 <- ggplot(forcor, aes(x = sens_2010, y = sens_2015)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
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
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
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
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
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

# conditional dependencies-------------------------

cond_dep_all <- tree.time %>%
    # filter(Cno == 15) %>%
    filter(Cno %in% c(5, 15, 25)) %>%
    group_by(Cno) %>%
    pivot_longer(cols = c("calcDBH_min1", "cii_min1"), names_to = "var", values_to = "value") %>%
    dplyr::mutate(varnames = ifelse(var == "calcDBH_min1", "DBH", "CII"))

# these variables are conditionally independent for the most part

# plot the conditional independencies

library(ggpubr)
cond_dep_all_plot <- ggscatter(
    data = cond_dep_all,
    x = "value", y = "twi",
    add = "reg.line", conf.int = TRUE, alpha = 0.3,
    cor.coef = TRUE, cor.method = "spearman",
    xlab = "variable", ylab = "TWI"
) +
    facet_wrap(varnames ~ yr, scales = "free")

cond_dep_all_plot

png("doc/display/cond_dep_alltrees.png", width = 8, height = 8, units = "in", res = 300)
cond_dep_all_plot
dev.off()

# raw distribution plots--------------------------

# plot of sensitivity against CII
sens_cii <- ggplot(
    tree.time %>% filter(yr %in% c(2010, 2015, 2020)),
    aes(x = factor(cii_min1), y = sens.prop)
) +
    geom_boxplot() +
    geom_jitter(alpha = 0.3) +
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

sens_dbh <- ggplot(
    tree.time %>% filter(yr %in% c(2010, 2015, 2020)),
    aes(x = calcDBH_min1, y = sens.prop)
) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", col = "grey40") +
    facet_wrap(~yr) +
    xlab("DBH (cm)") +
    ylab("Sensitivity") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

sens_dbh

sens_twi <- ggplot(
    tree.time %>% filter(yr %in% c(2010, 2015, 2020)),
    aes(x = twi, y = sens.prop)
) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", col = "grey40") +
    facet_wrap(~yr) +
    xlab("TWI") +
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


library(patchwork)
# png("doc/display/Fig3.png", width = 8, height = 8, units = "in", res = 300)
png("doc/display/Fig3.png", width = 8, height = 8, units = "in", res = 300)
dec_intercept_plot / iso_plot_twi + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.8, 1))
dev.off()

# SI figure with points and deciduousness--------------------

# preds <- readRDS("results/models/non_negative/pred_isoclines_twi.RDS")
# pred_df <- do.call(rbind, preds)

# fits <- readRDS("results/models/non_negative/fits_isoclines_twi.RDS")
# library(brms)

# res1 <- residuals(fits[[1]])
# res2 <- residuals(fits[[2]])
# res3 <- residuals(fits[[3]])

# head(res1)
# head(pred_df)

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
    select(Tag, treeID, yr, sens.prop, median, lwr, upr, Estimate, Est.Error, Q2.5, Q97.5, name, value) %>%
    dplyr::mutate(name = ifelse(name == "calcDBH_min1_scaled", "scaled DBH",
        ifelse(name == "cii_min1", "CII",
            ifelse(name == "twi_scaled", "scaled TWI", "deciduousness")
        )
    ))

# plot these


plot_res <- ggplot(pred_res_long, aes(x = value, y = Estimate)) +
    geom_point(alpha = 0.3) +
    geom_smooth() +
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
    facet_wrap(name ~ factor(yr, levels = c(2010, 2015, 2020)),
        scales = "free", ncol = 3, strip.position = "left"
    ) +
    # add text with p value
    geom_text(data = coefs_dec_lms %>% filter(term == "value"), aes(x = 1, y = 0.5, label = paste("p = ", round(p.value, 2)), hjust = 0, vjust = 0)) +
    labs(x = "species trait value", y = "intercept") +
    guides(color = "none") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )


library(patchwork)
png("doc/display/Fig_trait_intercept_plot.png", width = 8, height = 6, units = "in", res = 300)
# sensplot + sp_intercept_plot + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1.8, 1))
dec_intercept_plot
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
    facet_wrap(~ factor(yr, levels = c(2010, 2015, 2020)),
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

library(patchwork)
# png("doc/display/Fig4.png", width = 16, height = 8, units = "in", res = 300)
# (coefs_all_sp / twi_slopes_plot) | pred_plot & plot_annotation(tag_levels = "a")
# dev.off()


layout <- "
AACCC
BBCCC
"

png("doc/display/Fig4.png", width = 10, height = 6, units = "in", res = 300)
coefs_all_sp + twi_slopes_plot + pred_plot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()

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

library(patchwork)


layout <- "
A
B
"

png("doc/display/Fig4_tpi_3yrs.png", width = 8, height = 8, units = "in", res = 300)
coefs_all_sp + pred_plot + plot_layout(design = layout, guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))
dev.off()



# fig 5-----------------------------------

# fullmed_img <- magick::image_read("doc/display/full_mediation.png")
# fullmed_gg <- magick::image_ggplot(fullmed_img, interpolate = F)

# parmed_img <- magick::image_read("doc/display/hypotheses_grey.png")
# parmed_gg <- magick::image_ggplot(parmed_img, interpolate = F)

# # coefficient plots

# models_dir <- "results/models/orderedcii"

# coefs_fullmed <- readRDS(paste0(models_dir, "/coefs_fullmediation.rds"))
# coefs_parmed <- readRDS(paste0(models_dir, "/coefs_partialmed.rds"))

# coefs_fullmed <- do.call(rbind, coefs_fullmed)
# coefs_parmed <- do.call(rbind, coefs_parmed)

# pars.keep <- c("b_sensprop_calcDBH_min1_scaled", "bsp_sensprop_mocii_min1", "b_sensprop_twi_scaled")

# par_names <- as_labeller(c(
#     "b_sensprop_calcDBH_min1_scaled" = "DBH effect",
#     "bsp_sensprop_mocii_min1" = "CII effect",
#     "b_sensprop_twi_scaled" = "TWI effect"
# ))

# # plot the coefficients
# coefs_plot_fullmed <- ggplot(coefs_fullmed %>% filter(param %in% pars.keep), aes(x = param, y = median, ymin = lwr, ymax = upr)) +
#     geom_pointrange() +
#     facet_wrap(~yr) +
#     theme_bw() +
#     theme(strip.background = element_blank()) +
#     scale_x_discrete(labels = par_names) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     labs(
#         title = "Coefficients for the full mediation model",
#         x = "Parameter",
#         y = "Coefficient"
#     ) +
#     coord_flip()

# coefs_plot_fullmed

# coefs_plot_parmed <- ggplot(coefs_parmed %>% filter(param %in% pars.keep), aes(x = param, y = median, ymin = lwr, ymax = upr)) +
#     geom_pointrange() +
#     facet_wrap(~yr) +
#     theme_bw() +
#     theme(strip.background = element_blank()) +
#     scale_x_discrete(labels = par_names) +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     labs(
#         title = "Coefficients for the partial mediation model",
#         x = "Parameter",
#         y = "Coefficient"
#     ) +
#     coord_flip()

# # png("doc/display/Fig5.png", width = 12, height = 8, units = "in", res = 300)
# # (coefs_plot_parmed + parmed_gg) / (coefs_plot_fullmed + fullmed_gg)
# # dev.off()

# plot conditional effects of cii

# https://discourse.mc-stan.org/t/conditional-effects-plot-for-monotonic-predictors-with-logit-link/35850/2

# read fits
# fits <- readRDS("results/models/orderedcii/fits_partialmed_rel.rds")
fits <- readRDS("results/models/orderedcii/fits_rel_spre.rds")

library(tidybayes)
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

# png("doc/display/Fig5_alternate.png", width = 12, height = 8, units = "in", res = 300)
# (coefs_plot_parmed + parmed_gg) / p_manual_ce
# dev.off()


# DAGS with values

dag_2010 <- magick::image_read("doc/display/dag_2010_vals.png")
dag_2015 <- magick::image_read("doc/display/dag_2015_vals.png")
dag_2020 <- magick::image_read("doc/display/dag_2020_vals.png")

dag_2010_gg <- magick::image_ggplot(dag_2010, interpolate = F)
dag_2015_gg <- magick::image_ggplot(dag_2015, interpolate = F)
dag_2020_gg <- magick::image_ggplot(dag_2020, interpolate = F)

# layout <- "
# AAABBB
# AAABBB
# AAABBB
# #CCCC#
# #CCCC#
# "
# png("doc/display/Fig5.png", width = 8, height = 6, units = "in", res = 300)
# dag_2010_gg + dag_2015_gg + p_manual_ce + plot_layout(design = layout) + plot_annotation(tag_levels = "a")
# dev.off()

layout <- "
AABBCC
AABBCC
#DDDD#
#DDDD#
"
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


# # plot of change in effect sizes

# # read coefs
# coefs_spre <- readRDS("results/models/orderedcii/coefs_rel_spre.rds")
# coefs_sp <- readRDS("results/models/non_negative/sensitivity_model_intercept.RData")
# coefs_isocline <- readRDS("results/models/non_negative/coefs_isoclines_nore.rds")

# coefs_isocline <- do.call(rbind, coefs_isocline)
# coefs_isocline

# cors between species median sensitivities ----------

library(tidyverse)
tree.time <- readRDS("data/HKK-dendro/sensitivity_data_formodels.RData")


tree.time.sp <- tree.time %>%
    group_by(williams_dec, spfull) %>%
    dplyr::summarise(
        inc_annual_med = round(median(inc_annual, na.rm = T), 2),
        inc_annual_sd = round(sd(inc_annual, na.rm = T), 2)
    ) %>%
    arrange(inc_annual_med)

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

# do cor.test by deciduousness - this is not significant

# sens.sp.cor_dec <- sens.sp.wide %>%
#     group_by(dec) %>%
#     dplyr::summarise(
#         cor = cor.test(sens.2010, sens.2015)[4]$estimate,
#         cor_p = cor.test(sens.2010, sens.2015)[3]$p.value
#     )

# sens.sp.cor_dec

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
    labs(x = "2010 sensitivity", y = "2015 sensitivity") +
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
    labs(x = "2015 sensitivity", y = "2020 sensitivity") +
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
    labs(x = "2010 sensitivity", y = "2020 sensitivity") +
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
    ylab("sensitivity range between droughts") +
    theme_bw()

sens_rangeplot

library(patchwork)

layout <- "
AB
CD
"

png("doc/display/Fig_SI_sensitivity_cor_dec.png", width = 16, height = 16, units = "in", res = 300)
sens.sp_plot_1 + sens_sp_plot_2 + sens_sp_plot_3 + sens_rangeplot + plot_layout(design = layout)
dev.off()


# SI conditional dependencies -----------------------

## test conditional independences for each species-----------------------------

# DBH | TWI
# CII | TWI

# 2015 data only

cond_dep <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        DBH_TWI = cor.test(calcDBH_min1, twi)[4]$estimate,
        DBH_TWI_p = cor.test(calcDBH_min1, twi)[3]$p.value,
        CII_TWI = cor.test(cii_min1, twi)[4]$estimate,
        CII_TWI_p = cor.test(cii_min1, twi)[3]$p.value
    )

# these variables are conditionally independent for the most part

# plot the conditional independencies

library(ggpubr)
cond_dep_dbh_twi <- ggscatter(
    data = tree.time %>% filter(Cno == 15),
    x = "calcDBH_min1", y = "twi",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "spearman",
    xlab = "DBH", ylab = "TWI"
) +
    facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

png("doc/display/explore/cond_dep_dbh_twi.png", width = 12, height = 12, units = "in", res = 300)
cond_dep_dbh_twi
dev.off()

# CII | TWI
cond_dep_cii_twi <- ggscatter(
    data = tree.time %>% filter(Cno == 15),
    x = "cii_min1", y = "twi",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "pearson",
    xlab = "CII", ylab = "TWI"
) +
    facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

png("doc/display/explore/cond_dep_cii_twi.png", width = 12, height = 12, units = "in", res = 300)
cond_dep_cii_twi
dev.off()

# plot conditional dependencies across all trees

cond_dep_cii_twi_all <- ggscatter(
    data = tree.time %>% filter(yr %in% c(2010, 2015)),
    x = "cii_min1", y = "twi",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "pearson",
    xlab = "CII", ylab = "TWI",
    alpha = 0.5
) +
    facet_wrap(~yr)

cond_dep_cii_twi_all

cond_dep_dbh_twi_all <- ggscatter(
    data = tree.time %>% filter(yr %in% c(2010, 2015)),
    x = "calcDBH_min1", y = "twi",
    add = "reg.line", conf.int = TRUE,
    cor.coef = TRUE, cor.method = "pearson",
    xlab = "DBH", ylab = "TWI",
    alpha = 0.5
) +
    facet_wrap(~yr)

cond_dep_dbh_twi_all
library(patchwork)
png("doc/display/explore/cond_dep_alltrees.png", width = 8, height = 8, units = "in", res = 300)
cond_dep_dbh_twi_all + cond_dep_cii_twi_all + plot_layout(ncol = 1)
dev.off()
