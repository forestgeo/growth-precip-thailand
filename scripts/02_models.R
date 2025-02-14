# models to assess sensitivity of tree growth to precipitation

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

colours <- c("#e15f41", "#546de5", "#f7b731")

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

saveRDS(tree.time, "data/HKK-dendro/sensitivity_data_formodels.RData")
# plot the sensitivity values for all trees in 2010 and 2015
yrs <- c(2010, 2015)

sens.all <- ggplot(data = tree.time %>% filter(yr %in% yrs), aes(x = sens.prop)) +
    geom_density(alpha = 0.5) +
    scale_fill_viridis_d() +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
    facet_wrap(~yr) +
    # xlim(-5, 5)+
    labs(title = "Distribution of drought sensitivities for all trees", x = "Sensitivity", y = "Density") +
    theme_bw()

png("doc/display/sens_all.png", width = 6, height = 4, units = "in", res = 300)
sens.all
dev.off()

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
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 2
    ) +
    # add points of these
    geom_point(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 3
    ) +
    # mean of species
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
    # add points of these
    # geom_point(
    #     data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
    #         dplyr::summarise(median_inc = median(inc_annual, na.rm = T)) %>%
    #         ungroup() %>%
    #         group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
    #     aes(x = yr, y = median_inc), col = "grey40", size = 3
    # ) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(x = c(2010, 2015), y = 0.55, label = "ENSO drought"), hjust = 0.8, vjust = -0.2, angle = 90) +
    guides(col = guide_legend("species")) +
    xlab("year") +
    ylab("diameter increment (cm)") +
    ggtitle("growth increments for top 10 species") +
    theme_bw()

png("doc/display/spaghetti_top10_new.png", width = 4, height = 4, units = "in", res = 300)
spagplot_top10
dev.off()




# plot the distribution of sensitivities for the top 10 species
yrs <- c(2010, 2015)
sens_top10 <- ggplot(data = tree.time %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = sens.prop, fill = Species)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~yr) +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(title = "Distribution of drought sensitivities for top 10 species", x = "Species", y = "Sensitivity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

png("doc/display/sens_top10.png", width = 6, height = 4, units = "in", res = 300)
sens_top10
dev.off()


# pull out the species with rings
ring_sp <- tree.time %>%
    filter(Species %in% c("AFZEXY", "NEOLOB", "TOONCI", "CHUKTA", "MELIAZ"))

remove_sp <- tree.time %>%
    dplyr::mutate(Species = "ALL")

ring_sp <- rbind(ring_sp, remove_sp)

library(viridis)

ring_sp$Species <- factor(ring_sp$Species, levels = c("AFZEXY", "NEOLOB", "TOONCI", "CHUKTA", "MELIAZ", "ALL"))

# plot the distribution of sensitivities for the species with rings
sens_ring <- ggplot(data = ring_sp %>% filter(yr %in% yrs), aes(
    x = Species, y = sens.prop,
    fill = Species
)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~yr) +
    # make color scale viridis and one grey
    scale_fill_manual(values = c(viridis(3), "grey80")) +
    # make y axis percent
    scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    # add n for the number of trees
    geom_text(
        # make n for each species
        data = ring_sp %>% filter(yr %in% yrs) %>% group_by(Species, yr) %>% dplyr::summarise(n = n()),
        aes(label = paste0("n = ", n), x = Species, y = 5), size = 5.5, hjust = 0.5, vjust = 0
    ) +
    labs(title = "Distribution of drought sensitivities for species with rings", x = "Sensitivity", y = "Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

png("doc/display/sens_ring.png", width = 12, height = 8, units = "in", res = 300)
sens_ring
dev.off()

# Model 0: intercept only model--------------------------------------------

intercept_model <- bf(sens.prop ~ 1 + (1 | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(intercept_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    # fit <- brm(intercept_model, data = tree.time%>%filter(yr==yrs[i]), family = skew_normal(), iter=3000, warmup=1000, chains = 4, cores = 4)
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# save the fits
saveRDS(fits, "results/models/non_negative/fits_interceptonly.RDS")


# save the coefs and predictions
coefs_df <- do.call(rbind, coefs)
pred_df <- do.call(rbind, pred)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_intercept.RData")

saveRDS(pred_df, "results/models/non_negative/predictions_intercept.RData")

coefs_df <- readRDS("results/models/non_negative/sensitivity_model_intercept.RData")
pred_df <- readRDS("results/models/non_negative/predictions_intercept.RData")

# plot random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)


# plot the random effects
ranef_plot <- ggplot(data = ranef_df, aes(x = Species, y = median, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/ranefs_intercept.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

library(patchwork)
png("results/plots/non_negative/pp_intercept.png", width = 8, height = 4, res = 300, units = "in")
pp_check(fits[[1]]) + pp_check(fits[[2]])
dev.off()

# plot predicted vs observed distributions within species

# make long df for plotting
pred_obs_sp <- pred_df %>%
    pivot_longer(c("sens.prop", "median")) %>%
    dplyr::mutate(name = ifelse(name == "sens.prop", "obs", "pred"))

top10_predobs <- ggplot(data = pred_obs_sp %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = value, fill = Species, color = name)) +
    geom_violin() + # geom_boxplot(width=0.1, fill="white") +
    facet_wrap(~yr) +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # scale_color_manual(values=c())+
    # make y axis percent
    scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(title = "Distribution of sensitivities for top 10 species", x = "Sensitivity", y = "Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

# plot violin plots of observed sensitivity with a line for predictions

# add intercept to ranef_df
# make intercept vector
intercepts <- coefs_df %>%
    filter(param %in% "Intercept") %>%
    select(median) %>%
    pull(median)
intercepts <- rep(intercepts, each = nrow(ranef_df) / 2)
ranef_df <- ranef_df %>%
    mutate(
        intercept = intercepts + median,
        lwr = intercepts + lwr,
        upr = intercepts + upr
    )

top10_predobs <- ggplot(data = pred_df %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = sens.prop, fill = Species)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~yr) +
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    # add lines for each species predicted median
    # geom_hline(data=ranef_df %>% filter(Species %in% top_10_sp$Species), aes(x=Species, yintercept=median), col="grey20", linewidth=1.2) +
    geom_segment(
        data = ranef_df %>% filter(Species %in% top_10_sp$Species), aes(x = rep(1:10, 2) - 0.5, xend = rep(1:10, 2) + 0.5, y = intercepts, yend = intercepts),
        # data = data.frame(x = 1:3, y = c(5, 6, 7)),
        colour = "red", linetype = "dashed"
    ) +
    labs(title = "Distribution of sensitivities for top 10 species", x = "Sensitivity", y = "Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

top10_predobs


png("results/plots/non_negative/ranefs_intercept_bysp_new.png", width = 12, height = 8, units = "in", res = 300)
top10_predobs
dev.off()


# plot intercepts against species characteristics

# join ranef_df with sp_vars
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

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

# sp_intercept_plot <- ggplot(ranef_df_long, aes(x = value, y = intercept, col = factor(signif, levels = c("neg", "pos", "no")))) +
sp_intercept_plot <- ggplot(ranef_df_long, aes(x = value, y = intercept)) +
    geom_point() +
    geom_smooth(aes(x = value, y = intercept), inherit.aes = F, method = "lm") +
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(name ~ factor(yr, levels = c(2010, 2015)), scales = "free", ncol = 2) +
    labs(title = "Species intercepts", x = "species trait value", y = "intercept") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/sp_intercept_plot.png", width = 4, height = 6, units = "in", res = 300)
sp_intercept_plot
dev.off()

# plot the range of sensitivities for each species
sp_intercept_range_plot <- ggplot(ranef_df_long, aes(x = value, y = upr - lwr)) +
    geom_point() +
    geom_smooth(aes(x = value, y = upr - lwr), inherit.aes = F, method = "lm") +
    facet_wrap(name ~ factor(yr, levels = c(2010, 2015)), scales = "free", ncol = 2) +
    labs(title = "Species intercept range", x = "species trait value", y = "intercept range") +
    theme_bw()


# run intercept only model for all years

yrs <- sort(unique(tree.time$yr))

for (i in 1:length(yrs)) {
    fit <- brm(intercept_model, data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# save the coefs and predictions
coefs_df <- do.call(rbind, coefs)
pred_df <- do.call(rbind, pred)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_intercept_all13yrs.RData")

saveRDS(pred_df, "data/HKK-dendro/predictions_intercept_all13yrs.RData")

# plot random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)


# plot the random effects
ranef_plot <- ggplot(data = ranef_df, aes(x = Species, y = median, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = yrs), scales = "free") +
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("doc/display/ranefs_intercept_all13yrs.png", width = 12, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

# Model 1 - all species separately-----------------------------------------
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

## run the models-----------------------------


# model from DAG

sp_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp)

# make a function to run the model and return coefs

run_model <- function(data, model) {
    fit <- brm(model, data = data, family = gaussian(), chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("mean", "lwr", "upr")
    post_sum$param <- rownames(post_sum)

    # predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and n.trees columns
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("mean", "lwr", "upr")
    # add this to the observation data
    pred_sum <- cbind(data, pred_sum)

    return(list(post_sum, pred_sum))
}

# TODO: update with predictions and save predictions

# run this function for all 30 species in 2010 and 2015
# year<-2010
for (year in yrs) {
    tree.time.yr <- tree.time %>% filter(yr == year)

    all_sp <- tree.time %>%
        filter(Cno == 15) %>%
        group_by(Species) %>%
        dplyr::summarise(
            n = n()
        ) %>%
        arrange(desc(n))

    coefs <- list()
    preds <- list()

    for (i in 1:nrow(all_sp)) {
        print(paste0("Running model for", year, " for species ", i, " : ", all_sp$Species[i]))
        sp <- all_sp$Species[i]
        data <- tree.time.yr %>% filter(Species == sp)
        results <- run_model(data, sp_model)
        coefs[[i]] <- results[[1]]
        preds[[i]] <- results[[2]]
    }

    # results


    # unlist coefs, make a df and add species names
    coefs_df <- do.call(rbind, coefs)
    head(coefs_df)


    # add species names by repeating each element of top_10_sp$Species 8 times
    coefs_df$Species <- rep(all_sp$Species, each = 8)

    # add a column for significance
    coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
        ifelse(coefs_df$lwr > 0, "pos", "neg")
    )

    saveRDS(list(coefs_df, all_sp), paste0("results/models/non_negative/sensitivity_model_", year, ".RData"))

    preds_df <- do.call(rbind, preds)
    saveRDS(preds_df, paste0("results/models/non_negative/predictions_", year, ".RData"))
}
# plot the mean and 95% credible interval for each parameter

# read the RDS file
# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[1]]
# all_sp<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[2]]

# read both rds files and make a df for coefs and another for species
coefs_df <- list()
all_sp <- list()
preds_df <- list()
for (yr in yrs) {
    i <- which(yrs == yr)
    coefs_df[[i]] <- readRDS(paste0("results/models/non_negative/sensitivity_model_", yr, ".RData"))[[1]]
    coefs_df[[i]]$yr <- rep(yr, nrow(coefs_df[[i]]))
    all_sp[[i]] <- readRDS(paste0("results/models/non_negative/sensitivity_model_", yr, ".RData"))[[2]]
    all_sp[[i]]$yr <- rep(yr, nrow(all_sp[[i]]))

    preds_df[[i]] <- readRDS(paste0("results/models/non_negative/predictions_", yr, ".RData"))
}

str(coefs_df)

coefs_df <- do.call(rbind, coefs_df)
# all_sp<-do.call(rbind, all_sp)

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect", "2010" = "2010", "2015" = "2015"))

`%nin%` <- Negate(`%in%`)
coefs_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp"), Species %nin% c("NEOLOB", "IRVIMA")),
    aes(x = factor(Species, levels = rev(all_sp[[1]]$Species)), y = mean, col = factor(signif, levels = c("neg", "pos", "no")))
) +
    geom_point() +
    # make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(yr ~ param, scales = "free", labeller = par_names) +
    labs(title = "Effect of parameters on growth sensitivity", x = "Species", y = "Mean") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_sp

# write these as pngs
png(paste0("results/plots/non_negative/coefs_sp.png"), width = 8, height = 8, units = "in", res = 300)
coefs_sp
dev.off()

# do this by species
coefs_sp_param <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(x = param, y = mean, col = factor(signif, levels = c("neg", "pos", "no")))
) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(yr ~ factor(Species, levels = rev(all_sp[[1]]$Species)), scales = "free") + # , labeller = par_names) +
    labs(title = "Effect of parameters on growth sensitivity", x = "Species", y = "Mean") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/coefs_sp_separate.png", width = 12, height = 4, units = "in", res = 300)
coefs_sp_param
dev.off()

# plot these coefs against species variables

# merge the sp_vars with the coefs_df
coefs_df <- merge(coefs_df, sp_vars, by = "Species", all.x = TRUE)

# make a long df for plotting
coefs_df_long <- coefs_df %>%
    filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")) %>%
    dplyr::mutate(param = ifelse(param == "b_calcDBH_min1_scaled_sp", "DBH effect",
        ifelse(param == "b_cii_min1_scaled_sp", "CII effect", "TWI effect")
    )) %>%
    pivot_longer(c("maxDBH", "williams_dec", "median_inc"), names_to = "sp_vars")

# plot the effects of the species variables on the coefs
coefs_spvars_plot <- ggplot(coefs_df_long %>% filter(yr == 2010, Species %nin% c("NEOLOB", "IRVIMA")), aes(x = value, y = mean, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    geom_smooth(method = "lm", col = "grey") +
    scale_color_manual(values = c("red", "blue", "black")) +
    facet_grid(param ~ sp_vars, scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Effect of species variables on growth sensitivity", x = "Species variable", y = "Mean") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/coefs_spvars.png", width = 8, height = 8, units = "in", res = 300)
coefs_spvars_plot
dev.off()

# predicted vs observed plots

preds_df <- do.call(rbind, preds_df)
head(preds_df)

preds_plot <- ggplot(data = preds_df %>% filter(Species %nin% c("NEOLOB", "IRVIMA")), aes(x = sens.prop, y = mean, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_viridis_d() +
    facet_grid(. ~ yr, scales = "free") +
    labs(title = "Predicted vs observed growth sensitivity", x = "Observed", y = "Predicted") +
    theme_bw()

preds_plot

png("results/plots/non_negative/preds_plot_spsep.png", width = 8, height = 4, units = "in", res = 300)
preds_plot
dev.off()



# Model 2 : tree is a tree model --------------------------------------------------

# first test conditional dependencies DBH | TWI and CII | TWI
# make a long dataframe for plotting conditional independencies

cond_dep_all <- tree.time %>%
    # filter(Cno == 15) %>%
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

png("doc/display/cond_dep_alltrees.png", width = 8, height = 8, units = "in", res = 300)
cond_dep_all_plot
dev.off()



# first find sensitivity relative to the mean growth across all trees of that species

# #make a new column for the mean growth rate of each species
# tree.time<-tree.time %>% group_by(Species, yr) %>%
# dplyr::mutate(sp_mean_inc = mean(inc_annual, na.rm = TRUE))%>%
# ungroup()%>%
# dplyr::mutate(sens_prop_sp = (inc_annual - sp_mean_inc)/sp_mean_inc)

## define and run model -------------------------

tree_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + cii_min1 + twi_scaled + (1 | Species))

# run the model

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(tree_model,
        # data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(),
        data = tree.time %>% filter(yr == yrs[i]), family = gaussian(),
        chains = 4, iter = 4000, warmup = 2000, cores = 4
    )
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

saveRDS(fits, "results/models/non_negative/fits_treeisatree.RDS")

summary(fits[[1]])

png("results/plots/non_negative/fits_treeisatree.png")
plot(fits[[1]], variable = "^b", regex = T)
dev.off()

library(patchwork)
png("results/plots/non_negative/pp_treeisatree.png", width = 8, height = 4, units = "in", res = 300)
pp_check(fits[[1]]) + pp_check(fits[[2]])
dev.off()

# save the coefs and predictions

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_treeisatree.RData")

coefs_df <- readRDS("results/models/non_negative/sensitivity_model_treeisatree.RData")

saveRDS(pred_df, "results/models/non_negative/predictions_treeisatree.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_tree <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1", "b_twi_scaled")),
    aes(
        x = param, y = median,
        # col = factor(signif, levels = c("neg", "pos", "no"))
        col = factor(param, levels = c("b_calcDBH_min1_scaled", "b_twi_scaled", "b_cii_min1"))
    )
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(param, levels = c("b_calcDBH_min1_scaled", "b_twi_scaled", "b_cii_min1"))), width = 0.1) +
    # scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    scale_color_manual(values = rep(colours, 2), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_tree

# write these as pngs
png("results/plots/non_negative/coefs_tree.png", width = 6, height = 4, units = "in", res = 300)
coefs_tree
dev.off()

# make plots for the species random effects

# first subset the coefs_df for the random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
# add a column with the species names by removing the "r_Species[" and ",Intercept]" from the param column
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)

# plot the random effects
ranef_plot <- ggplot(data = ranef_df, aes(x = Species, y = median, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/ranefs_tree.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

colnames(pred_df)

# plot predictions against observations
pred_plot <- ggplot(data = pred_df, aes(x = sens.prop, y = median, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    scale_color_viridis_d() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    labs(title = "Predictions vs observations", x = "Observations", y = "Predictions") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/pred_vs_obs_treeisatree.png", width = 8, height = 4, units = "in", res = 300)
pred_plot
dev.off()

# plotting the random effects against species characteristics----------

# merge the sp_vars with the coefs_df
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

head(ranef_df)

# plot ranefs against maxDBH, deciduousness and growth rate

# first make long dataframe for plotting
ranef_long_2015 <- ranef_df %>%
    filter(yr == 2015) %>%
    pivot_longer(cols = c("maxDBH", "williams_dec", "median_inc"), names_to = "sp_var", values_to = "value")

# make labels
sp_var_names <- as_labeller(c("maxDBH" = "maxDBH", "williams_dec" = "deciduousness", "median_inc" = "growth rate"))

ranef_sp_plot <- ggplot(data = ranef_long_2015, aes(
    x = value, y = mean,
    color = factor(signif, levels = c("neg", "no", "pos"))
)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color = factor(signif, levels = c("neg", "no", "pos"))), width = 0.1) +
    scale_color_manual(values = c("red", "grey40", "blue")) +
    ggtitle("Species variables and random effect") +
    ylab("species random effect") +
    facet_wrap(~sp_var, scales = "free", labeller = sp_var_names) +
    guides(color = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw()

png("results/models/non_negative/ranef_sp_plot.png", width = 10, height = 4, units = "in", res = 300)
ranef_sp_plot
dev.off()

# plot model predictions


# model that controls for individual growth rate----------------------------------------

tree_model_gr <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + cii_min1_scaled + twi_scaled + avg_inc_tree_scaled + (1 | Species))

coefs <- list()
pred <- list()

for (i in 1:2) {
    # yrs <- c(2010, 2015, 2020)
    fit <- brm(tree_model_gr, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# save the coefs and predictions

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_treeisatreegrowthrate.RData")

# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model_treeisatree.RData")

saveRDS(pred_df, "results/models/non_negative/predictions_treeisatree_growthrate.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect", "b_avg_inc_tree_scaled" = "growth rate effect"))

coefs_tree <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled", "b_avg_inc_tree_scaled")),
    aes(x = param, y = median, col = factor(signif, levels = c("neg", "pos", "no")))
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_tree

# write these as pngs
png("results/plots/non_negative/coefs_treegrowthrate.png", width = 6, height = 4, units = "in", res = 300)
coefs_tree
dev.off()

# make plots for the species random effects

# first subset the coefs_df for the random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
# add a column with the species names by removing the "r_Species[" and ",Intercept]" from the param column
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)

# plot the random effects
ranef_plot <- ggplot(data = ranef_df, aes(x = Species, y = median, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/ranefs_treegrowthrate_allyrs.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()


# plot predictions against observations
pred_plot <- ggplot(data = pred_df, aes(x = sens.prop, y = median, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    scale_color_viridis_d() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    labs(title = "Predictions vs observations", x = "Observations", y = "Predictions") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/pred_vs_obs_treeisatreegrowthrate.png", width = 8, height = 4, units = "in", res = 300)
pred_plot
dev.off()


# Model3: model that uses predictors scaled at the species level ----------------------------

tree_model_rel <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    #    fit <- brm(tree_model_rel, data = tree.time%>%filter(yr==yrs[i]), family = skew_normal(), iter=3000, warmup = 1000, chains = 4, cores = 4)
    fit <- brm(tree_model_rel, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# save models
# summary(fits[[1]]) #Rhats are fine
saveRDS(fits, "results/models/non_negative/fits_treeisatreerel.RDS")

library(patchwork)
# plot posterior predictive checks
png("results/plots/non_negative/pp_treeisatreerel.png", width = 8, height = 4, res = 300, units = "in")
pp_check(fits[[1]]) + pp_check(fits[[2]])
dev.off()

# plot chains

p1 <- plot(fits[[1]], variable = "^b", regex = T)
p2 <- plot(fits[[2]], variable = "^b", regex = T)
p <- cbind(p1, p2)
plot(p)

library(gtable)
library(gridExtra)

png("doc/display/diagnostics_treeisatreerel.png", width = 8, height = 4, units = "in", res = 300)

dev.off()



# save the coefs and predictions

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_treeisatreerel.RData")

# coefs_df <- readRDS("results/models/non_negative/sensitivity_model_treeisatreerel.RData")

saveRDS(pred_df, "results/models/non_negative/predictions_treeisatree_rel.RData")

# pred_df <- readRDS("results/models/non_negative/predictions_treeisatree_rel.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_tree <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(
        x = param, y = median,
        # col = factor(signif, levels = c("neg", "pos", "no"))
        col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
    )
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(
        aes(
            ymin = lwr, ymax = upr,
            # col = factor(signif, levels = c("neg", "pos", "no"))
            col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
        ),
        width = 0.1
    ) +
    # scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    scale_color_manual(values = rep(colours, 2), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_tree

# write these as pngs
png("results/plots/non_negative/coefs_treerel.png", width = 6, height = 4, units = "in", res = 300)
coefs_tree
dev.off()

# plot conditional effects of cii against sensitivity

cond_plot <- conditional_effects(fits[[1]], "cii_min1_scaled_sp", categorical = FALSE, points = TRUE)
cond_plot <- cond_plot + theme_bw()

cond_plot

# predicted plots
predcii_plot <- ggplot(pred_df, aes(x = cii_min1_scaled_sp, y = median, col = Species)) +
    geom_smooth(method = "lm", se = FALSE) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr, col=Species), alpha=0.2)+
    geom_point(aes.inherit = F, aes(x = cii_min1_scaled_sp, y = sens.prop), alpha = 0.1) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col=Species), width=0.1)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    scale_color_viridis_d() +
    labs(title = "Predicted sensitivity vs CII", x = "CII", y = "Predicted sensitivity") +
    # guides(color="none")+
    theme_bw()

png("results/plots/non_negative/pred_cii_treerel_withpts.png", width = 8, height = 4, units = "in", res = 300)
predcii_plot
dev.off()

# make plots for the species random effects

# first subset the coefs_df for the random effects
ranef_df <- coefs_df %>% filter(grepl("r_Species", param))
# add a column with the species names by removing the "r_Species[" and ",Intercept]" from the param column
ranef_df$Species <- gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)

# plot the random effects
ranef_plot <- ggplot(data = ranef_df, aes(x = Species, y = median, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/ranefs_treerel.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

# plot predictions against observations
pred_plot <- ggplot(data = pred_df, aes(x = sens.prop, y = median, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    scale_color_viridis_d() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    labs(title = "Predictions vs observations", x = "Observations", y = "Predictions") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/pred_vs_obs_treeisatreerel.png", width = 12, height = 4, units = "in", res = 300)
pred_plot
dev.off()

# plot species variables against random effects----------------

# merge the sp_vars with the coefs_df
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

head(ranef_df)

# plot ranefs against maxDBH, deciduousness and growth rate

# first make long dataframe for plotting
ranef_long <- ranef_df %>% # filter(yr==2015) %>%
    pivot_longer(cols = c("maxDBH", "williams_dec", "median_inc"), names_to = "sp_var", values_to = "value")

# make labels
sp_var_names <- as_labeller(c(
    "maxDBH" = "maxDBH", "williams_dec" = "deciduousness", "median_inc" = "growth rate", "2010" =
        "2010", "2015" = "2015", "2020" = "2020"
))

ranef_sp_plot <- ggplot(data = ranef_long, aes(
    x = value, y = median,
    color = factor(signif, levels = c("neg", "no", "pos"))
)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color = factor(signif, levels = c("neg", "no", "pos"))), width = 0.1) +
    scale_color_manual(values = c("red", "grey40", "blue"), drop = FALSE) +
    ggtitle("Species variables and random effect") +
    ylab("species random effect") +
    facet_grid(yr ~ sp_var, scales = "free", labeller = sp_var_names) +
    guides(color = "none") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw()

png("doc/display/ranef_sp_plot_tree_rel.png", width = 10, height = 8, units = "in", res = 300)
ranef_sp_plot
dev.off()


# Model 5: species random effects on slope ------------------------------------------------------

twi_slope_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 + twi_scaled_sp | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(twi_slope_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_twi_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(x = param, y = median, col = factor(signif, levels = c("neg", "pos", "no")))
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_twi_sp

png("results/plots/non_negative/coefs_twi_sp.png", width = 8, height = 4, units = "in", res = 300)
coefs_twi_sp
dev.off()

pp_check(fits[[1]]) + pp_check(fits[[2]])

# Model 6: Model with random slope on CII

cii_slope_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 + cii_min1_scaled_sp | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(cii_slope_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

# plot the coefs

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_cii_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(x = param, y = median, col = factor(signif, levels = c("neg", "pos", "no")))
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_cii_sp

png("results/plots/non_negative/coefs_cii_sp.png", width = 6, height = 4, units = "in", res = 300)
coefs_cii_sp
dev.off()

# predicted plots
predcii_plot <- ggplot(pred_df, aes(x = cii_min1_scaled_sp, y = median, col = Species)) +
    geom_smooth(method = "lm", se = FALSE) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr, col=Species), alpha=0.2)+
    geom_point(aes(x = cii_min1_scaled_sp, y = sens.prop, col = Species), alpha = 0.3) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col=Species), width=0.1)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    scale_color_viridis_d() +
    labs(title = "Predicted sensitivity vs CII", x = "CII", y = "Predicted sensitivity") +
    # guides(color="none")+
    theme_bw()

png("results/plots/non_negative/pred_cii_tree_ciislope2.png", width = 8, height = 4, units = "in", res = 300)
predcii_plot
dev.off()


# pred vs obs

# plot predictions against observations
pred_plot <- ggplot(data = pred_df, aes(x = sens.prop, y = median, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    scale_color_viridis_d() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    labs(title = "Predictions vs observations", x = "Observations", y = "Predictions") +
    guides(color = "none") +
    theme_bw()

png("results/plots/non_negative/pred_vs_obs_cii_sp.png", width = 8, height = 4, units = "in", res = 300)
pred_plot
dev.off()

# Model 6: species random effects on all slopes

all_slope_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 + cii_min1_scaled_sp + calcDBH_min1_scaled_sp + twi_scaled_sp | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(all_slope_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    fits[[i]] <- fit
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_spre.RData")
saveRDS(pred_df, "results/models/non_negative/predictions_spre.RData")

# read the saved files
# coefs_df <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")
# pred_df <- readRDS("results/models/non_negative/predictions_spre.RData")

# plot the coefs

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_all_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(
        x = param, y = median,
        # col = factor(signif, levels = c("neg", "pos", "no"))
        col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
    )
) +
    geom_point() +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(aes(
        ymin = lwr, ymax = upr,
        # col = factor(signif, levels = c("neg", "pos", "no"))
        col = factor(param, levels = c("b_calcDBH_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1_scaled_sp"))
    ), width = 0.1) +
    # scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    scale_color_manual(values = rep(colours, 2), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_all_sp

png("results/plots/non_negative/coefs_sp_slopes.png", width = 6, height = 4, units = "in", res = 300)
coefs_all_sp
dev.off()

# make long dataframe to plot predicted effects of cii, dbh and twi

pred_coef_long <- pred_df %>%
    pivot_longer(c("cii_min1_scaled_sp", "calcDBH_min1_scaled_sp", "twi_scaled_sp"), names_to = "predictor")

# make labels
par_names <- as_labeller(c("cii_min1_scaled_sp" = "CII effect", "calcDBH_min1_scaled_sp" = "DBH effect", "twi_scaled_sp" = "TWI effect", "2010" = "2010", "2015" = "2015", "2020" = "2020"))

# predicted plots
pred_coefs_plot <- ggplot(pred_coef_long, aes(x = value, y = median, col = Species)) +
    geom_point(aes(x = value, y = sens.prop, col = Species), alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr, col=Species), alpha=0.2)+
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col=Species), width=0.1)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(predictor ~ yr, scales = "free", labeller = par_names) +
    scale_color_viridis_d() +
    labs(title = "Predicted sensitivity vs predictor", x = "value", y = "Predicted sensitivity") +
    # guides(color="none")+
    theme_bw()

pred_coefs_plot

png("results/plots/non_negative/pred_coefs_spre_slope.png", width = 6, height = 6, units = "in", res = 300)
pred_coefs_plot
dev.off()

# pred plot with cii

pred_cii_plot <- ggplot(pred_df, aes(x = cii_min1_scaled_sp, y = median, col = Species)) +
    geom_smooth(method = "lm", se = FALSE) +
    # geom_ribbon(aes(ymin = lwr, ymax = upr, col=Species), alpha=0.2)+
    geom_point(aes(x = cii_min1_scaled_sp, y = sens.prop, col = Species), alpha = 0.1) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col=Species), width=0.1)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    scale_color_viridis_d() +
    labs(title = "Predicted sensitivity vs CII", x = "CII", y = "Predicted sensitivity") +
    # guides(color="none")+
    theme_bw()

png("results/plots/non_negative/pred_cii_spre.png", width = 8, height = 4, units = "in", res = 300)
pred_cii_plot
dev.off()



# Combined plots---------------------------------------
# make a plot with three different coefs together

coefs_df_tree <- readRDS("results/models/non_negative/sensitivity_model_treeisatree.RData")
coefs_df_tree_rel <- readRDS("results/models/non_negative/sensitivity_model_treeisatreerel.RData")
coefs_df_spre <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")

# make a long dataframe for plotting
coefs_long <- rbind(
    coefs_df_tree %>% mutate(model = "varying intercepts"),
    coefs_df_tree_rel %>% mutate(model = "varying intercepts + \nspecies scaling"),
    coefs_df_spre %>% mutate(model = "varying slopes + \nspecies scaling")
) %>%
    # filter only the parameters of interest
    filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp", "b_cii_min1", "b_calcDBH_min1_scaled", "b_twi_scaled")) %>%
    mutate(param_name = case_when(
        grepl("DBH", param) ~ "DBH effect",
        grepl("cii", param) ~ "CII effect",
        grepl("twi", param) ~ "TWI effect"
    ))


# make plot

coefs_all_plot <- ggplot(coefs_long, aes(
    x = factor(param_name, levels = c("DBH effect", "CII effect", "TWI effect")), y = median,
    col = factor(param_name, levels = c("DBH effect", "TWI effect", "CII effect"))
)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, linewidth = 1.3) +
    facet_grid(model ~ yr, scales = "free") +
    scale_color_manual(values = colours, drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/non_negative/coefs_all_plot.png", width = 8, height = 6, units = "in", res = 300)
coefs_all_plot
dev.off()


# mediation effect models------------------------------------------------

# first test for the effect of DBH and TWI on sensitivity

no_cii_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + twi_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()

yrs <- c(2010, 2015)

for (i in 1:length(yrs)) {
    fit <- brm(no_cii_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_nocii.RData")
coefs_df <- readRDS("results/models/non_negative/sensitivity_model_nocii.RData")

saveRDS(pred_df, "results/models/non_negative/predictions_nocii.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_tree_no_cii <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    # aes(x = param, y = median, col = factor(signif, levels = c("neg", "pos", "no")))
    aes(x = param, y = median, col = param)
) +
    geom_point(size = 4) +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = param), width = 0.1, linewidth = 1.3) +
    scale_color_manual(values = rep(colours, 2), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015)), scales = "free") +
    labs(title = "Effect of parameters on growth sensitivity with CII removed", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_tree_no_cii

# then test the effect of DBH on CII

dbh_cii_model <- bf(cii_min1_scaled_sp ~ 1 + calcDBH_min1_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()



for (i in 1:length(yrs)) {
    # yrs <- c(2010, 2015, 2020)
    fit <- brm(dbh_cii_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("median", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    post_sum$yr <- rep(yrs[i], nrow(post_sum))
    coefs[[i]] <- post_sum

    # make predictions
    preds <- posterior_predict(fit)
    # this makes a dataframe with 4000 rows (chains* sampling iterations) and 1449 columns (number of trees)
    pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sum) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
}

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

# save these
saveRDS(coefs_df, "results/models/non_negative/sensitivity_model_dbh_cii.RData")
saveRDS(pred_df, "results/models/non_negative/predictions_dbh_cii.RData")

coefs_df <- readRDS("results/models/non_negative/sensitivity_model_dbh_cii.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_dbh_cii <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
    aes(x = param, y = median)
) +
    geom_point(col = colours[1], size = 4) +
    scale_x_discrete(labels = par_names) +
    # make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr), col = colours[1], width = 0.1, linewidth = 1.3) +
    # scale_color_manual(values = c("red", "blue", "grey40"), drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~ factor(yr, levels = c(2010, 2015, 2020)), scales = "free") +
    labs(title = "Effect of DBH on CII", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_dbh_cii

# save these together as a png
library(gridExtra)
png("results/plots/non_negative/coefs_mediation.png", width = 8, height = 4, units = "in", res = 300)
grid.arrange(coefs_tree_no_cii, coefs_dbh_cii, ncol = 1, heights = c(2, 1))
dev.off()


# plot species random effect slope

# read coefs
coefs_df <- readRDS("results/models/non_negative/sensitivity_model_spre.RData")

head(coefs_df)
# species random effects on TWI slope

# mean slope
coefs_twi <- coefs_df %>%
    filter(grepl("b_twi_scaled_sp", param))

# random effect on slope
coefs_sp_twi <- coefs_df %>%
    filter(grepl("r_Species", param)) %>%
    filter(grepl("twi", param)) %>%
    filter(grepl("cii|DBH|Intercept", param) == FALSE) %>%
    dplyr::mutate(
        Species = gsub("r_Species\\[|\\,twi_scaled_sp\\]", "", param)
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
    aes(x = williams_dec, y = total_effect),
    order = median
) +
    geom_point() +
    # geom_smooth(method = "lm", col = "grey40") +
    geom_errorbar(aes(ymin = total_lwr, ymax = total_upr), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015)),
        scales = "free", ncol = 2
    ) +
    stat_cor(method = "pearson", label.x = 0.5, label.y = 0.5) +
    labs(x = "Deciduousness", y = "TWI slope") +
    # coord_flip()+
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

twi_slopes_plot

png("results/plots/non_negative/twislope_dec.png", width = 8, height = 4, units = "in", res = 300)
twi_slopes_plot
dev.off()


# repeat for cii slopes

# mean slope
coefs_cii <- coefs_df %>%
    filter(grepl("b_cii_min1_scaled_sp", param))

# random effect on slope
coefs_sp_cii <- coefs_df %>%
    filter(grepl("r_Species", param)) %>%
    filter(grepl("cii", param)) %>%
    filter(grepl("twi|DBH|Intercept", param) == FALSE) %>%
    dplyr::mutate(
        Species = gsub("r_Species\\[|\\,cii_min1_scaled_sp\\]", "", param)
    ) %>%
    group_by(yr) %>%
    dplyr::mutate(
        total_effect = median + coefs_cii$median[match(yr, coefs_cii$yr)],
        total_lwr = lwr + coefs_cii$median[match(yr, coefs_cii$yr)],
        total_upr = upr + coefs_cii$median[match(yr, coefs_cii$yr)]
    )

# merge species vars
coefs_sp_cii <- merge(coefs_sp_cii, sp_vars, by = "Species", all.x = TRUE)

# cors
cii_cor <- coefs_sp_cii %>%
    group_by(yr) %>%
    dplyr::summarise(
        cor_dec = cor.test(williams_dec, total_effect)[4]$estimate,
        cor_dec_p = cor.test(williams_dec, total_effect)[3]$p.value,
        cor_dbh = cor.test(maxDBH, total_effect)[4]$estimate,
        cor_dbh_p = cor.test(maxDBH, total_effect)[3]$p.value
    )

cii_cor
# plot coefs

library(ggpubr)
cii_slopes_plot <- ggplot(
    data = coefs_sp_cii,
    # aes(x = reorder(Species, median), y = median),
    aes(x = williams_dec, y = total_effect),
    order = median
) +
    geom_point() +
    # geom_smooth(method = "lm", col = "grey40") +
    geom_errorbar(aes(ymin = total_lwr, ymax = total_upr), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ factor(yr, levels = c(2010, 2015)),
        scales = "free", ncol = 2
    ) +
    stat_cor(method = "pearson", label.x = 0.5, label.y = 0.5) +
    labs(x = "Deciduousness", y = "CII slope") +
    # coord_flip()+
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 12)
    )

cii_slopes_plot

png("results/plots/non_negative/ciislope_dec.png", width = 8, height = 4, units = "in", res = 300)
cii_slopes_plot
dev.off()
