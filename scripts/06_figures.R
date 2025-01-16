# script to code for figures

# load libraries
library(tidyverse)
library(lubridate)

# Fig 1
# climate data

spei <- read.csv("data/climate/SPEI_HKK_from_GEE.csv")
ffstation <- read.csv("data/climate/WeatherData_ALL_forest_fire_research_station.csv")

ffstation_monthly <- ffstation %>%
    group_by(YearII, Month) %>%
    dplyr::summarise(
        precipitation = sum(Precipitation),
        dry_days = sum(Precipitation == 0),
        tmax = max(TempMax),
    ) %>%
    rename(year = YearII, month = Month) %>%
    dplyr::mutate(year = as.character(year)) %>%
    select(year, month, precipitation, dry_days, tmax)

ffstation_lt <- ffstation_monthly %>%
    filter(year <= 2021, year >= 2009) %>%
    ungroup() %>%
    group_by(month) %>%
    dplyr::mutate(
        precipitation = mean(precipitation, na.rm = T),
        dry_days = mean(dry_days, na.rm = T),
        tmax = mean(tmax, na.rm = T),
        year = "long-term"
    ) %>%
    select(year, month, precipitation, dry_days, tmax)

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
        cols = c(precipitation, dry_days, tmax, spei_val),
        names_to = "climvar",
        values_to = "value"
    ) %>%
    mutate(variable = factor(climvar, levels = c("precipitation", "dry_days", "tmax", "spei_val")))


# plot

climplot <- ggplot() +
    geom_line(clim, aes(x = month, y = value, color = year), linewidth = 2) +
    facet_wrap(~variable, scales = "free_y") +
    theme_bw() +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # make every month show up on the x-axis
    scale_x_continuous(breaks = 1:12) +

    # add rectangles for dry season
    geom_rect(
        data = data.frame(xmin = c(11, 1), xmax = c(12, 4), ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.2
    ) +
    labs(x = "Month", y = "Value", color = "Year") +
    scale_color_manual(values = c("long-term" = "black", "2010" = "gold", "2015" = "magenta"))

climplot

png("doc/display/Fig1.png", width = 8, height = 8, units = "in", res = 300)
climplot
dev.off()
