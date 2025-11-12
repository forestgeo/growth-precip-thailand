# load libraries
library(tidyverse)
library(brms)
library(ggplot2)
library(ggpubr)
library(patchwork)

# read data
rm(list = ls())
tree.time <- read.csv("data/dendro/sensitivity_dataset.csv")
median_incs <- read.csv("data/dendro/summaries_dataset.csv")

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

# str(tree.time$sens.prop[tree.time$yr == 2010])
# trial <- stats::supsmu(x = tree.time$twi[tree.time$yr == 2010], y = tree.time$sens.prop[tree.time$yr == 2010])
# str(trial)

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
