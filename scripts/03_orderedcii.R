# models considering cii as an ordered factor

# Load required libraries---------------------
library(tidyverse)
library(brms)

# load data--------------------------
rm(list = ls())
tree.time <- read.csv("data/HKK-dendro/to_publish/sensitivity_dataset_for_review.csv")


# conditional dependencies---------------

# ## test conditional independences for each species-----------------------------

# # DBH | TWI
# # CII | TWI

# # 2015 data only

# cond_dep <- tree.time %>%
#     filter(yr == 2015) %>%
#     group_by(Species) %>%
#     dplyr::summarise(
#         DBH_TWI = cor.test(calcDBH_min1, twi)[4]$estimate,
#         DBH_TWI_p = cor.test(calcDBH_min1, twi)[3]$p.value,
#         CII_TWI = cor.test(cii_min1, twi)[4]$estimate,
#         CII_TWI_p = cor.test(cii_min1, twi)[3]$p.value
#     )

# # these variables are conditionally independent for the most part

# # plot the conditional independencies

# library(ggpubr)
# cond_dep_dbh_twi <- ggscatter(
#     data = tree.time %>% filter(Cno == 15),
#     x = "calcDBH_min1", y = "twi",
#     add = "reg.line", conf.int = TRUE,
#     cor.coef = TRUE, cor.method = "spearman",
#     xlab = "DBH", ylab = "TWI"
# ) +
#     facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

# png("doc/display/cond_dep_dbh_twi.png", width = 12, height = 12, units = "in", res = 300)
# cond_dep_dbh_twi
# dev.off()

# # CII | TWI
# cond_dep_cii_twi <- ggscatter(
#     data = tree.time %>% filter(Cno == 15),
#     x = "cii_min1", y = "twi",
#     add = "reg.line", conf.int = TRUE,
#     cor.coef = TRUE, cor.method = "pearson",
#     xlab = "CII", ylab = "TWI"
# ) +
#     facet_wrap(~ factor(Species, levels = names(sort(table(Species), decreasing = T))))

# png("doc/display/cond_dep_cii_twi.png", width = 12, height = 12, units = "in", res = 300)
# cond_dep_cii_twi
# dev.off()


# directory--------

plot_dir <- "results/plots/orderedcii"
models_dir <- "results/models/orderedcii"


# define and run the model

tree_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + mo(cii_min1) + twi_scaled + (1 | Species))

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

saveRDS(fits, file = paste0(models_dir, "/fits.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs.rds"))
saveRDS(pred, file = paste0(models_dir, "/pred.rds"))

# plot the results

# unlist coefs
coefs <- do.call(rbind, coefs)
coefs

summary(fits[[2]])


# model fitting both models----------------------------------------

tree_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + mo(cii_min1) + twi_scaled + (1 | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled, family = cumulative(link = "logit"))

# run the model

coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model + cii_model + set_rescor(FALSE),
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
    # pred_sum <- as.data.frame(t(apply(preds, 2, quantile, probs = c(.5, .05, .95))))
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_partialmed.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_partialmed.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_partialmed.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_partialmed.rds"))

# plot the results

# unlist coefs
coefs <- do.call(rbind, coefs)
coefs

summary(fits[[2]])


# pred vs obs plot

pred <- readRDS(paste0(models_dir, "/pred_sens_partialmed.rds"))
pred_df <- do.call(rbind, pred)

pred_plot <- ggplot(pred_df, aes(x = sens.prop, y = median)) +
    geom_point(alpha = 0.05) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, alpha = 0.05) +
    facet_wrap(~yr) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Predicted vs observed sensitivity",
        x = "Observed sensitivity",
        y = "Predicted sensitivity"
    )

# pred vs obs for cii
pred_cii <- readRDS(paste0(models_dir, "/pred_cii_partialmed.rds"))
pred_cii_df <- do.call(rbind, pred_cii)

pred_plot_cii <- ggplot(pred_cii_df, aes(x = cii_min1, y = median)) +
    geom_jitter(alpha = 0.05) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, alpha = 0.05) +
    facet_wrap(~yr) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Predicted vs observed CII",
        x = "Observed CII",
        y = "Predicted CII"
    )

png(paste0(plot_dir, "/pred_obs_partialmed.png"), width = 8, height = 8, units = "in", res = 300)
pred_plot / pred_plot_cii
dev.off()

# full mediation model------------------------------------------

tree_model <- bf(sens.prop ~ mo(cii_min1) + twi_scaled + (1 | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled, family = cumulative(link = "logit"))

# run the model

coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model + cii_model + set_rescor(FALSE),
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
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_fullmediation.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_fullmediation.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_fullmediation.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_fullmediation.rds"))

# pred vs obs plot
pred_sens <- readRDS(paste0(models_dir, "/pred_sens_fullmediation.rds"))
pred_sens_df <- do.call(rbind, pred_sens)

pred_plot <- ggplot(pred_sens_df, aes(x = sens.prop, y = median)) +
    geom_point(alpha = 0.05) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, alpha = 0.05) +
    facet_wrap(~yr) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Predicted vs observed sensitivity",
        x = "Observed sensitivity",
        y = "Predicted sensitivity"
    )

# pred vs obs for cii
pred_cii <- readRDS(paste0(models_dir, "/pred_cii_fullmediation.rds"))
pred_cii_df <- do.call(rbind, pred_cii)

pred_plot_cii <- ggplot(pred_cii_df, aes(x = cii_min1, y = median)) +
    geom_jitter(alpha = 0.05) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.1, alpha = 0.05) +
    facet_wrap(~yr) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Predicted vs observed CII",
        x = "Observed CII",
        y = "Predicted CII"
    )

library(patchwork)
png(paste0(plot_dir, "/pred_obs_fullmediation.png"), width = 8, height = 8, units = "in", res = 300)
pred_plot / pred_plot_cii
dev.off()



# comparing the model fits------------------------------------

parmed <- readRDS(paste0(models_dir, "/fits_partialmed.rds"))
fullmed <- readRDS(paste0(models_dir, "/fits_fullmediation.rds"))

LOO(parmed[[1]], fullmed[[1]])
LOO(parmed[[2]], fullmed[[2]])

# variance explained
performance::icc(parmed[[1]])

summary(parmed[[1]])


# ordered CII, partial mediation using species scaled variables-------------------------

# Model 3: models with variables scaled at species level
tree_model_rel <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + mo(cii_min1) + twi_scaled_sp + (1 | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled_sp, family = cumulative(link = "logit"))

# run the model

coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model_rel + cii_model + set_rescor(FALSE),
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
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_partialmed_rel.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_partialmed_rel.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_partialmed_rel.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_partialmed_rel.rds"))

coefs <- do.call(rbind, coefs)

# plot coefs

par.keep <- c("b_sensprop_calcDBH_min1_scaled_sp", "b_sensprop_twi_scaled_sp", "bsp_sensprop_mocii_min1")
par_names <- as_labeller(c(
    b_sensprop_calcDBH_min1_scaled_sp = "DBH",
    b_sensprop_twi_scaled_sp = "TWI",
    bsp_sensprop_mocii_min1 = "CII"
))

coefs_plot <- ggplot(coefs %>% filter(param %in% par.keep), aes(x = param, y = median, ymin = lwr, ymax = upr)) +
    geom_pointrange() +
    facet_wrap(~yr) +
    theme_bw() +
    scale_x_discrete(labels = par_names) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        title = "Model coefficients",
        x = "Parameter",
        y = "Coefficient"
    ) +
    coord_flip()

png(paste0(plot_dir, "/coefs_partialmed_rel.png"), width = 8, height = 8, units = "in", res = 300)
coefs_plot
dev.off()

# plot monotonic effects

fits <- readRDS(paste0(models_dir, "/fits_partialmed_rel.rds"))
cond1 <- conditional_effects(fits[[1]], "cii_min1", plot = F)

p1 <- plot(cond1, plot = F)[[1]] +
    geom_pointrange() + theme_bw()

p1

cond2 <- conditional_effects(fits[[2]], "cii_min1", plot = F)
p2 <- plot(cond2, plot = F)[[1]] +
    geom_pointrange() + theme_bw()

p2

# Model 4: varying slopes-----------------------------

tree_model_spre <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + mo(cii_min1) + twi_scaled_sp + (1 + calcDBH_min1_scaled_sp + twi_scaled_sp + mo(cii_min1) | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled_sp, family = cumulative(link = "logit"))

# run the model
coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model_spre + cii_model + set_rescor(FALSE),
        data = tree.time %>% filter(yr == yrs[i]), family = gaussian(),
        chains = 4, iter = 4000, warmup = 2000, cores = 4,
        control = list(adapt_delta = 0.9)
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
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_spre.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_spre.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_spre.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_spre.rds"))

# Model 5: varying slopes without species scaling------------------------------

tree_model_rel_spre <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + mo(cii_min1) + twi_scaled + (1 + calcDBH_min1_scaled + twi_scaled + mo(cii_min1) | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled, family = cumulative(link = "logit"))

# run the model
coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model_rel_spre + cii_model + set_rescor(FALSE),
        data = tree.time %>% filter(yr == yrs[i]), family = gaussian(),
        chains = 4, iter = 4000, warmup = 2000, cores = 4,
        control = list(adapt_delta = 0.9)
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
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_rel_spre.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_rel_spre.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_rel_spre.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_rel_spre.rds"))

# Model 2: species random effect on intercept only ---------------------------------

tree_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + mo(cii_min1) + twi_scaled + (1 | Species))
cii_model <- bf(ordered(cii_min1) ~ calcDBH_min1_scaled, family = cumulative(link = "logit"))

# run the model
coefs <- list()
pred_sens <- list()
pred_cii <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(tree_model + cii_model + set_rescor(FALSE),
        data = tree.time %>% filter(yr == yrs[i]), family = gaussian(),
        chains = 4, iter = 4000, warmup = 2000, cores = 4,
        control = list(adapt_delta = 0.9)
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
    pred_sens_df <- as.data.frame(t(apply(preds[, , 1], 2, quantile, probs = c(.5, .05, .95))))
    pred_cii_df <- as.data.frame(t(apply(preds[, , 2], 2, quantile, probs = c(.5, .05, .95))))
    colnames(pred_sens_df) <- colnames(pred_cii_df) <- c("median", "lwr", "upr")
    # add this to the observation data
    pred_sens[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_sens_df)
    pred_cii[[i]] <- cbind(tree.time %>% filter(yr == yrs[i]), pred_cii_df)
}

saveRDS(fits, file = paste0(models_dir, "/fits_tree.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_tree.rds"))
saveRDS(pred_sens, file = paste0(models_dir, "/pred_sens_tree.rds"))
saveRDS(pred_cii, file = paste0(models_dir, "/pred_cii_tree.rds"))
