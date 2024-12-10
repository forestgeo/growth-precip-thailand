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
    # mean of species
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
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


# spagplot using zero growth model
spagplot_top10_zero <- ggplot() +
    # species plots
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual_zero, na.rm = T)),
        aes(
            x = yr, y = median_inc,
            group = Species, col = Species
        )
    ) +
    # mean of all trees
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr) %>%
            dplyr::summarise(median_inc = median(inc_annual_zero, na.rm = T)),
        aes(x = yr, y = median_inc), col = "black", size = 2
    ) +
    # mean of species
    geom_line(
        data = tree.time %>% filter(Species %in% top_10_sp$Species) %>% group_by(yr, Species) %>%
            dplyr::summarise(median_inc = median(inc_annual_zero, na.rm = T)) %>%
            ungroup() %>%
            group_by(yr) %>% dplyr::summarise(median_inc = mean(median_inc, na.rm = T)),
        aes(x = yr, y = median_inc), col = "grey40", size = 0.8
    ) +
    scale_color_viridis_d() +
    geom_vline(xintercept = c(2010, 2015), linetype = "dashed") +
    # add text on these lines
    geom_text(aes(x = c(2010, 2015), y = 0.55, label = "ENSO drought"), hjust = 0.8, vjust = -0.2, angle = 90) +
    guides(col = guide_legend("species")) +
    xlab("year") +
    ylab("diameter increment (cm)") +
    ggtitle("growth increments for top 10 species with zero growth assumption") +
    theme_bw()

png("doc/display/spaghetti_top10_zero_growth.png", width = 4, height = 4, units = "in", res = 300)
spagplot_top10_zero
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
    labs(title = "Distribution of sensitivities for top 10 species", x = "Sensitivity", y = "Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

png("doc/display/sens_top10.png", width = 6, height = 4, units = "in", res = 300)
sens_top10
dev.off()

# sensitivities for top 10 species using zero growth model
sens_top10_zero <- ggplot(data = tree.time %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = sens.prop.zero, fill = Species)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~yr) +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(title = "Distribution of sensitivities for top 10 species with zero growth assumption", x = "Sensitivity", y = "Density") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

png("doc/display/sens_top10_zero.png", width = 6, height = 4, units = "in", res = 300)
sens_top10_zero
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
sens_ring <- ggplot(data = ring_sp, aes(
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
        data = ring_sp %>% group_by(Species, yr) %>% dplyr::summarise(n = n()),
        aes(label = paste0("n = ", n), x = Species, y = 5), size = 5.5, hjust = 0.5, vjust = 0
    ) +
    labs(title = "Distribution of sensitivities for species with rings", x = "Sensitivity", y = "Density") +
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
saveRDS(fits, "data/HKK-dendro/fits_interceptonly.RDS")


# save the coefs and predictions
coefs_df <- do.call(rbind, coefs)
pred_df <- do.call(rbind, pred)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_intercept.RData")

saveRDS(pred_df, "data/HKK-dendro/predictions_intercept.RData")

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

png("doc/display/ranefs_intercept_allyrs.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

library(patchwork)
png("doc/display/pp_intercept.png", width = 8, height = 4, res = 300, units = "in")
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

png("doc/display/ranefs_intercept_bysp_new.png", width = 12, height = 8, units = "in", res = 300)
top10_predobs
dev.off()

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
    return(post_sum)
}

# run this function for all 30 species in 2010 and 2015

for (yr in yrs) {
    tree.time.yr <- tree.time %>% filter(yr == yr)
    # tree.time.2015<-tree.time %>% filter(Cno == 15)

    all_sp <- tree.time %>%
        filter(Cno == 15) %>%
        group_by(Species) %>%
        dplyr::summarise(
            n = n()
        ) %>%
        arrange(desc(n))

    coefs <- list()


    for (i in 1:nrow(all_sp)) {
        print(paste0("Running model for species ", i, " : ", all_sp$Species[i]))
        sp <- all_sp$Species[i]
        data <- tree.time.yr %>% filter(Species == sp)
        post_sum <- run_model(data, sp_model)
        coefs[[i]] <- post_sum
    }


    # unlist coefs, make a df and add species names
    coefs_df <- do.call(rbind, coefs)
    head(coefs_df)

    # add species names by repeating each element of top_10_sp$Species 8 times
    coefs_df$Species <- rep(all_sp$Species, each = 8)

    # add a column for significance
    coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
        ifelse(coefs_df$lwr > 0, "pos", "neg")
    )

    saveRDS(list(coefs_df, all_sp), paste0("data/HKK-dendro/sensitivity_model_", yr, "_.RData"))
}
# plot the mean and 95% credible interval for each parameter

# read the RDS file
# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[1]]
# all_sp<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[2]]

# read both rds files and make a df for coefs and another for species
coefs_df <- list()
all_sp <- list()
for (yr in yrs) {
    i <- which(yrs == yr)
    coefs_df[[i]] <- readRDS(paste0("data/HKK-dendro/sensitivity_model_", yr, "_.RData"))[[1]]
    coefs_df[[i]]$yr <- rep(yr, nrow(coefs_df[[i]]))
    all_sp[[i]] <- readRDS(paste0("data/HKK-dendro/sensitivity_model_", yr, "_.RData"))[[2]]
    all_sp[[i]]$yr <- rep(yr, nrow(all_sp[[i]]))
}

coefs_df <- do.call(rbind, coefs_df)
# all_sp<-do.call(rbind, all_sp)

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect", "2010" = "2010", "2015" = "2015"))

coefs_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled_sp", "b_cii_min1_scaled_sp", "b_twi_scaled_sp")),
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
png(paste0("doc/display/coefs_sp.png"), width = 8, height = 8, units = "in", res = 300)
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

png("doc/display/coefs_sp_separate_nolowgrowth.png", width = 12, height = 4, units = "in", res = 300)
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
coefs_spvars_plot <- ggplot(coefs_df_long %>% filter(yr == 2010), aes(x = value, y = mean, col = factor(signif, levels = c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = factor(signif, levels = c("neg", "pos", "no"))), width = 0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    facet_grid(param ~ sp_vars, scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Effect of species variables on growth sensitivity", x = "Species variable", y = "Mean") +
    guides(color = "none") +
    theme_bw()

coefs_spvars_plot

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
        data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(),
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

saveRDS(fits, "data/HKK-dendro/fits_treeisatree.RDS")

summary(fits[[1]])

png("doc/display/fits_treeisatree.png")
plot(fits[[1]], variable = "^b", regex = T)
dev.off()

library(patchwork)
png("doc/display/pp_treeisatree.png", width = 8, height = 4, units = "in", res = 300)
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

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_treeisatree.RData")

# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model_treeisatree.RData")

saveRDS(pred_df, "data/HKK-dendro/predictions_treeisatree.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_tree <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1", "b_twi_scaled")),
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
png("doc/display/coefs_tree_allyrs.png", width = 6, height = 4, units = "in", res = 300)
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

png("doc/display/ranefs_tree_allyrs.png", width = 8, height = 8, units = "in", res = 300)
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

png("doc/display/pred_vs_obs_treeisatree.png", width = 12, height = 4, units = "in", res = 300)
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

png("doc/display/ranef_sp_plot.png", width = 10, height = 4, units = "in", res = 300)
ranef_sp_plot
dev.off()

# plot model predictions


# model that controls for individual growth rate----------------------------------------

tree_model_gr <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + cii_min1_scaled + twi_scaled + avg_inc_tree_scaled + (1 | Species))

coefs <- list()
pred <- list()

for (i in 1:3) {
    yrs <- c(2010, 2015, 2020)
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

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_treeisatreegrowthrate.RData")

# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model_treeisatree.RData")

saveRDS(pred_df, "data/HKK-dendro/predictions_treeisatree_growthrate.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect", "b_median_inc_scaled" = "growth rate effect"))

coefs_tree <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled", "b_median_inc_scaled")),
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
png("doc/display/coefs_treegrowthrate_allyrs.png", width = 6, height = 4, units = "in", res = 300)
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

png("doc/display/ranefs_treegrowthrate_allyrs.png", width = 8, height = 8, units = "in", res = 300)
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

png("doc/display/pred_vs_obs_treeisatreegrowthrate.png", width = 12, height = 4, units = "in", res = 300)
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
saveRDS(fits, "data/HKK-dendro/fits_treeisatreerel.RDS")

library(patchwork)
# plot posterior predictive checks
png("doc/display/pp_treeisatreerel.png", width = 8, height = 4, res = 300, units = "in")
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

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_treeisatreerel.RData")

coefs_df <- readRDS("data/HKK-dendro/sensitivity_model_treeisatreerel.RData")

saveRDS(pred_df, "data/HKK-dendro/predictions_treeisatree_rel.RData")

pred_df <- readRDS("data/HKK-dendro/predictions_treeisatree_rel.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_tree <- ggplot(
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

coefs_tree

# write these as pngs
png("doc/display/coefs_treerel_allyrs.png", width = 6, height = 4, units = "in", res = 300)
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
    # geom_point()+
    # geom_errorbar(aes(ymin = lwr, ymax = upr, col=Species), width=0.1)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    scale_color_viridis_d() +
    labs(title = "Predicted sensitivity vs CII", x = "CII", y = "Predicted sensitivity") +
    # guides(color="none")+
    theme_bw()

png("doc/display/pred_cii_treerel.png", width = 8, height = 4, units = "in", res = 300)
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

png("doc/display/ranefs_treerel_allyrs.png", width = 8, height = 8, units = "in", res = 300)
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

png("doc/display/pred_vs_obs_treeisatreerel.png", width = 12, height = 4, units = "in", res = 300)
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

# Models 4: mediation effect models------------------------------------------------

# first test for the effect of DBH and TWI on sensitivity

no_cii_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + twi_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()


for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
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

# saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_treeisatreerel.RData")

# coefs_df<-readRDS("data/HKK-dendro/sensitivity_model_treeisatree.RData")

# saveRDS(pred_df, "data/HKK-dendro/predictions_treeisatree_rel.RData")

# first make labels
par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

coefs_tree_no_cii <- ggplot(
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

coefs_tree_no_cii

# then test the effect of DBH on CII

dbh_cii_model <- bf(cii_min1_scaled_sp ~ 1 + calcDBH_min1_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()



for (i in 1:3) {
    yrs <- c(2010, 2015, 2020)
    fit <- brm(dbh_cii_model, data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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

coefs_dbh_cii <- ggplot(
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
    labs(title = "Effect of DBH on CII", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_dbh_cii

# save these together as a png
library(gridExtra)
png("doc/display/coefs_mediation.png", width = 8, height = 4, units = "in", res = 300)
grid.arrange(coefs_tree_no_cii, coefs_dbh_cii, ncol = 1, heights = c(2, 1))
dev.off()

# Model 5: species random effects on slope ------------------------------------------------------

twi_slope_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 + twi_scaled_sp | Species))

coefs <- list()
pred <- list()


for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(twi_slope_model, data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
