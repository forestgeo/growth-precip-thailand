# models considering cii as an ordered factor and using independent causal queries

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

# directory--------

plot_dir <- "results/plots/causal"
models_dir <- "results/models/causal"


# define and run the models

dbh_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + (1 + calcDBH_min1_scaled | Species))
twi_model <- bf(sens.prop ~ 1 + twi_scaled + (1 + twi_scaled | Species))
cii_model <- bf(sens.prop ~ 1 + mo(cii_min1) + calcDBH_min1_scaled + (1 + mo(cii_min1) | Species))


# run the model for dbh----------------------

coefs <- list()
pred <- list()
fits <- list()
i <- 1
for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(dbh_model,
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

saveRDS(fits, file = paste0(models_dir, "/fits_dbh.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_dbh.rds"))
saveRDS(pred, file = paste0(models_dir, "/pred_dbh.rds"))

# run the twi models------------------------

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(twi_model,
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

saveRDS(fits, file = paste0(models_dir, "/fits_twi.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_twi.rds"))
saveRDS(pred, file = paste0(models_dir, "/pred_twi.rds"))

# run models for cii----------

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(cii_model,
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

saveRDS(fits, file = paste0(models_dir, "/fits_cii.rds"))
saveRDS(coefs, file = paste0(models_dir, "/coefs_cii.rds"))
saveRDS(pred, file = paste0(models_dir, "/pred_cii.rds"))

coefs_dbh <- readRDS(paste0(models_dir, "/coefs_dbh.rds"))
coefs_twi <- readRDS(paste0(models_dir, "/coefs_twi.rds"))
coefs_cii <- readRDS(paste0(models_dir, "/coefs_cii.rds"))

coefs_dbh_df <- do.call(rbind, coefs_dbh)
coefs_twi_df <- do.call(rbind, coefs_twi)
coefs_cii_df <- do.call(rbind, coefs_cii)

coefs_dbh_df$model <- "dbh"
coefs_twi_df$model <- "twi"
coefs_cii_df$model <- "cii"

coefs_df <- rbind(coefs_dbh_df, coefs_twi_df, coefs_cii_df[coefs_cii_df$param != "b_calcDBH_min1_scaled", ])
nrow(coefs_df)
head(coefs_df)

# par_names <- as_labeller(c("b_calcDBH_min1_scaled_sp" = "DBH effect", "b_cii_min1_scaled_sp" = "CII effect", "b_twi_scaled_sp" = "TWI effect"))

par_names <- as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "bsp_mocii_min1" = "exposure effect", "b_twi_scaled" = "wetness effect"))

coefs_all_sp <- ggplot(
    data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "bsp_mocii_min1", "b_twi_scaled")),
    aes(
        x = factor(param, levels = c("b_calcDBH_min1_scaled", "bsp_mocii_min1", "b_twi_scaled")),
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


png("doc/display/causal_inquiry.png", width = 8, height = 4, units = "in", res = 300)
coefs_all_sp
dev.off()
