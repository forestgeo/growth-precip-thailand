# load libraries
library(tidyverse)
library(brms)
library(ggplot2)

# read data
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

colours <- c("#e15f41", "#546de5", "#f7b731")

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


## SI plot - correlation between 2010 and 2015

forcor <- tree.time %>%
    filter(yr %in% c(2010, 2015, 2020)) %>%
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


# vals<- cor.test(forcor$sens_2010, forcor$sens_2015, use = "pairwise.complete.obs")
# vals[4]

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


# conditional dependencies-------------------------

cond_dep_all <- tree.time %>%
    # filter(Cno == 15) %>%
    filter(yr %in% yrs) %>%
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
    filter(yr %in% c(2010, 2015, 2020)) %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(cor.test(data$calcDBH_min1, data$sens.prop))) %>%
    dplyr::reframe(broom::tidy(mod))

# make a supsmu df

smu_dbh <- tree.time %>%
    ungroup() %>%
    filter(yr %in% yrs) %>%
    group_by(yr) %>%
    dplyr::mutate(
        x = stats::supsmu(x = calcDBH_min1, y = sens.prop)$x,
        y = stats::supsmu(x = calcDBH_min1, y = sens.prop)$y
    )

sens_dbh <- ggplot(
    tree.time %>% filter(yr %in% c(2010, 2015, 2020)),
    aes(x = calcDBH_min1, y = sens.prop)
) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        data = tree.time %>% filter(yr == 2010),
        method = "lm", col = "darkblue"
    ) +
    geom_line(data = smu_dbh, aes(x = x, y = y), col = "red") +
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


sens_twi_cor <- tree.time %>%
    filter(yr %in% c(2010, 2015, 2020)) %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(cor.test(data$twi, data$sens.prop))) %>%
    dplyr::reframe(broom::tidy(mod))

smu_twi <- tree.time %>%
    ungroup() %>%
    filter(yr %in% yrs) %>%
    group_by(yr) %>%
    dplyr::reframe(
        x = stats::supsmu(x = twi, y = sens.prop)$x,
        y = stats::supsmu(x = twi, y = sens.prop)$y
    )

smu_twi

str(tree.time$sens.prop[tree.time$yr == 2010])
trial <- stats::supsmu(x = tree.time$twi[tree.time$yr == 2010], y = tree.time$sens.prop[tree.time$yr == 2010])
str(trial)

sens_twi <- ggplot(
    tree.time %>% filter(yr %in% c(2010, 2015, 2020)),
    aes(x = twi, y = sens.prop)
) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        data = tree.time %>%
            filter(yr == 2020),
        method = "lm", col = "darkblue"
    ) +
    geom_line(data = smu_twi, aes(x = x, y = y), col = "red") +
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

# SI figure with distribution of deciduousness------------------

str(sp_vars)

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

str(tree_vars)
tree_n <- tree_vars %>%
    group_by(Species) %>%
    dplyr::summarise(n = n())

str(tree_n)
sp_n <- merge(sp_vars, tree_n, by = "Species")

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

png("doc/display/cond_dep_dbh_twi.png", width = 12, height = 12, units = "in", res = 300)
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

png("doc/display/cond_dep_cii_twi.png", width = 12, height = 12, units = "in", res = 300)
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
png("doc/display/cond_dep_alltrees.png", width = 8, height = 8, units = "in", res = 300)
cond_dep_dbh_twi_all + cond_dep_cii_twi_all + plot_layout(ncol = 1)
dev.off()
