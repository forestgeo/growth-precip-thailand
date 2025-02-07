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

# colours <- c("#e15f41", "#546de5", "#f7b731")

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

# species median and lms------------------------



sens_sp <- tree.time %>%
    filter(yr %in% yrs) %>%
    group_by(Species, yr) %>%
    dplyr::summarise(
        sens_median = median(sens.prop, na.rm = TRUE),
        sens_sp_mean = mean(sens.prop, na.rm = TRUE),
        sens_lwr = quantile(sens.prop, probs = 0.05, na.rm = TRUE),
        sens_upr = quantile(sens.prop, probs = 0.95, na.rm = TRUE),
        sens_boot = boot::boot(sens.prop, statistic = function(x, i) median(x[i]), R = 1000)$t0 # ,
        # sens_boot_median = sens_boot$t0
        # sens_boot_lwr = sens_boot %>% map_dbl(~quantile(., probs = 0.05)),
        # sens_boot_upr = sens_boot %>% map_dbl(~quantile(., probs = 0.95))
        # sens_boot_ci = boot::boot.ci(sens_boot, 0.05, 0.95)
    ) %>%
    ungroup()

# b <- boot::boot(tree.time$sens.prop, function(x, i) median(x[i]), R = 1000)
# b$t0
# boot::boot.ci(b)

head(sens_sp)

sens.sp <- merge(sens_sp, sp_vars, by = "Species", all.x = TRUE)

# plot median sensitivity against deciduousness
sens_sp_median <- ggplot(sens.sp, aes(x = williams_dec, y = sens_median)) +
    geom_point() +
    geom_errorbar(aes(ymin = sens_lwr, ymax = sens_upr), width = 0.1) +
    theme_bw() +
    facet_wrap(~yr) +
    # scale_color_manual(values = colours) +
    labs(x = "Deciduousness", y = "Sensitivity", title = "Sensitivity vs Deciduousness") +
    theme(legend.position = "none")

sens_sp_median

# make lms for each year
sens_sp_lms <- sens.sp %>%
    filter(Species != "ALPHVE") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(lm(sens_median ~ williams_dec + twi_sd + maxDBH, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))

sens_sp_lms

sens_sp_lms <- sens.sp %>%
    filter(Species != "ALPHVE") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(lm(sens_sp_mean ~ williams_dec, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))


sens_sp_lms <- tree.time %>%
    filter(yr %in% yrs) %>%
    filter(Species != "ALPHVE") %>%
    nest_by(yr) %>%
    dplyr::mutate(mod = list(lme4::lmer(sens.prop ~ williams_dec + twi_sd + maxDBH + (1 | Species), data = data))) %>%
    dplyr::reframe(broom.mixed::tidy(mod))
sens_sp_lms

dec_model <- bf(sens_median ~ 1 + williams_dec)

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(dec_model, data = sens.sp %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
    # pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
    pred_sum$yr <- rep(yrs[i], nrow(pred_sum))
    pred[[i]] <- pred_sum
}

coefs_df <- do.call(rbind, coefs)
coefs_df

# plot coefs
coefs_plot <- ggplot(coefs_df %>% filter(param == "b_williams_dec"), aes(x = yr, y = median)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1) +
    # facet_wrap(~yr) +
    theme_bw() +
    labs(x = "Deciduousness", y = "Sensitivity", title = "Sensitivity vs Deciduousness") +
    theme(legend.position = "none")

coefs_plot

# model--------------------------------------------
library(brms)
isocline_model <- bf(sens.prop ~ 1 + twi + williams_dec + (1 + twi | Species))

isocline_model_twi <- bf(sens.prop ~ 1 + twi + williams_dec + twi:williams_dec)


coefs <- list()
pred <- list()
fits <- list()
new_preds <- list()

# make new data
sps <- unique(tree.time$Species)
sps <- setdiff(sps, c("ALPHVE", "MACASI"))
# sps <- c("SACCLI", "HOPEOD", "AFZEXY", "TETRNU")

# expand.grid to make 100 values each for twi and williams_dec

newdata <- expand.grid(
    twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 10),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 10),
    Species = sps
)

newdata <- expand.grid(
    twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

# i<-2
for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(isocline_model_twi, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
    # pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
    pred[[i]] <- pred_sum

    # new preds
    pred_new <- predict(fit, newdata = newdata, re_formula = NULL)
    new_preds[[i]] <- cbind(newdata, pred_new)
    new_preds[[i]]$yr <- rep(yrs[i], nrow(new_preds[[i]]))
}


saveRDS(fits, "results/models/non_negative/fits_isoclines.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines.RDS")

saveRDS(fits, "results/models/non_negative/fits_isoclines_nore.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_nore.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_nore.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_nore.RDS")

new_preds_df <- do.call(rbind, new_preds)
coefs_df <- do.call(rbind, coefs)

coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

coefs_df

library(ggplot2)
# plot for four species
iso_plot <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = twi, fill = Estimate)
) +
    geom_tile() +
    geom_contour(aes(z = Estimate), colour = "black") +
    facet_grid(yr ~ Species) +
    theme_bw() +
    scale_fill_gradient2()

png("results/plots/non_negative/isocline_trial.png", width = 32, height = 4, units = "in", res = 300)
iso_plot
dev.off()


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
    scale_fill_gradient2()

png("results/plots/non_negative/isocline_trial_nore.png", width = 8, height = 4, units = "in", res = 300)
iso_plot
dev.off()


# model with size and deciduousness------------------------

library(brms)
isocline_model_dbh <- bf(sens.prop ~ 1 + calcDBH_min1 + williams_dec + calcDBH_min1:williams_dec)

coefs <- list()
pred <- list()
fits <- list()
new_preds <- list()

# make new data
sps <- unique(tree.time$Species)
sps <- setdiff(sps, c("ALPHVE", "MACASI"))
# sps <- c("SACCLI", "HOPEOD", "AFZEXY", "TETRNU")

# expand.grid to make 100 values each for twi and williams_dec

# newdata <- expand.grid(
#     twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 10),
#     williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 10),
#     Species = sps
# )

newdata <- expand.grid(
    calcDBH_min1 = seq(min(tree.time$calcDBH_min1, na.rm = T), max(tree.time$calcDBH_min1, na.rm = T), length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(isocline_model_dbh, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
    # pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
    pred[[i]] <- pred_sum

    # new preds
    pred_new <- predict(fit, newdata = newdata, re_formula = NULL)
    new_preds[[i]] <- cbind(newdata, pred_new)
    new_preds[[i]]$yr <- rep(yrs[i], nrow(new_preds[[i]]))
}


saveRDS(fits, "results/models/non_negative/fits_isoclines_dbh.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_dbh.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_dbh.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_dbh.RDS")

new_preds_df <- do.call(rbind, new_preds)
coefs_df <- do.call(rbind, coefs)

coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

coefs_df

# library(ggplot2)
# # plot for four species
# iso_plot <- ggplot(
#     new_preds_df,
#     aes(x = williams_dec, y = twi, fill = Estimate)
# ) +
#     geom_tile() +
#     geom_contour(aes(z = Estimate), colour = "black") +
#     facet_grid(yr ~ Species) +
#     theme_bw() +
#     scale_fill_gradient2()

# png("results/plots/non_negative/isocline_trial.png", width = 32, height = 4, units = "in", res = 300)
# iso_plot
# dev.off()


iso_plot <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = calcDBH_min1, fill = Estimate)
) +
    geom_tile() +
    labs(x = "Deciduousness", y = "DBH", fill = "Sensitivity") +
    # geom_contour(aes(z = Estimate), colour = "black") +
    # facet_grid(yr ~ Species) +
    facet_wrap(~yr) +
    theme_bw() +
    scale_fill_gradient2()

png("results/plots/non_negative/isocline_dbh.png", width = 8, height = 4, units = "in", res = 300)
iso_plot
dev.off()

# cii model with deciduousness-----------------------
isocline_model_cii <- bf(sens.prop ~ 1 + cii_min1 + williams_dec + cii_min1:williams_dec)

coefs <- list()
pred <- list()
fits <- list()
new_preds <- list()

# make new data
sps <- unique(tree.time$Species)
sps <- setdiff(sps, c("ALPHVE", "MACASI"))
# sps <- c("SACCLI", "HOPEOD", "AFZEXY", "TETRNU")

# expand.grid to make 100 values each for twi and williams_dec

newdata <- expand.grid(
    cii_min1 = seq(1, 5, length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(isocline_model_cii, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
    # pred[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sum)
    pred[[i]] <- pred_sum

    # new preds
    pred_new <- predict(fit, newdata = newdata, re_formula = NULL)
    new_preds[[i]] <- cbind(newdata, pred_new)
    new_preds[[i]]$yr <- rep(yrs[i], nrow(new_preds[[i]]))
}


saveRDS(fits, "results/models/non_negative/fits_isoclines_cii.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_cii.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_cii.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_cii.RDS")

new_preds_df <- do.call(rbind, new_preds)
coefs_df <- do.call(rbind, coefs)

coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

coefs_df

iso_plot <- ggplot(
    new_preds_df,
    aes(x = williams_dec, y = cii_min1, fill = Estimate)
) +
    geom_tile() +
    labs(x = "Deciduousness", y = "CII", fill = "Sensitivity") +
    # geom_contour(aes(z = Estimate), colour = "black") +
    # facet_grid(yr ~ Species) +
    facet_wrap(~yr) +
    theme_bw() +
    scale_fill_gradient2()

png("results/plots/non_negative/isocline_cii.png", width = 8, height = 4, units = "in", res = 300)
iso_plot
dev.off()
