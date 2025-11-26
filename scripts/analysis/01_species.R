# Load required libraries---------------------

required_pkgs <- c(
    "tidyverse", "brms", "patchwork"
)

# install any missing packages
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
}

# load packages quietly
invisible(lapply(required_pkgs, function(pkg) {
    suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}))

# load data--------------------------
rm(list = ls())
tree.time <- read.csv("data/dendro/sensitivity_dataset.csv")
median_incs <- read.csv("data/dendro/summaries_dataset.csv")

head(list.files(recursive = T))

# ensure output directory exists
if (!dir.exists("results/models/non_negative")) {
    dir.create("results/models/non_negative", recursive = TRUE)
    message("Created directory: results/models/non_negative")
}

if (!dir.exists("results/plots/non_negative")) {
    dir.create("results/plots/non_negative", recursive = TRUE)
    message("Created directory: results/plots/non_negative")
}

# plot the sensitivity values for all trees in 2010 and 2015
yrs <- c(2010, 2015, 2020)

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

head(sens_sp)
head(median_incs)

# make sp_vars data
sp_vars <- tree.time %>%
    group_by(Species, spfull) %>%
    dplyr::summarise(williams_dec = first(williams_dec))

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
    dplyr::mutate(mod = list(lm(sens_sp_mean ~ williams_dec, data = data))) %>%
    dplyr::reframe(broom::tidy(mod))

sens_sp_lms

# Model 0: intercept only model--------------------------------------------

intercept_model <- bf(sens.prop ~ 1 + (1 | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(intercept_model, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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

png("results/plots/non_negative/pp_intercept.png", width = 8, height = 4, res = 300, units = "in")
pp_check(fits[[1]]) + pp_check(fits[[2]])
dev.off()

# plot violin plots of observed sensitivity with a line for predictions
nrow(coefs_df)

# add intercept to ranef_df
# make intercept vector
intercepts <- coefs_df %>%
    filter(param %in% "Intercept") %>%
    select(median) %>%
    pull(median)

intercepts <- rep(intercepts, each = nrow(ranef_df) / 3)

ranef_df <- ranef_df %>%
    mutate(
        intercept = intercepts + median,
        lwr = intercepts + lwr,
        upr = intercepts + upr
    )

top_10_sp <- tree.time %>%
    filter(yr == 2015) %>%
    group_by(Species) %>%
    dplyr::summarise(n.tree = n()) %>%
    ungroup() %>%
    dplyr::arrange(desc(n.tree)) %>%
    head(10)

# plot intercepts against species characteristics

# join ranef_df with sp_vars
ranef_df <- merge(ranef_df, sp_vars, by = "Species", all.x = TRUE)

# Deciduousness model--------------------------------------------

# model with twi--------------------------------------------
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

# newdata <- expand.grid(
#     twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 10),
#     williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 10),
#     Species = sps
# )

newdata <- expand.grid(
    twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    # fit <- brm(isocline_model_twi, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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


saveRDS(fits, "results/models/non_negative/fits_isoclines_twi.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_twi.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_twi.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_twi.RDS")

# model with tpi--------------------------

isocline_model_tpi <- bf(sens.prop ~ 1 + tpi + williams_dec + tpi:williams_dec)

coefs <- list()
pred <- list()
fits <- list()
new_preds <- list()


newdata <- expand.grid(
    tpi = seq(min(tree.time$tpi, na.rm = T), max(tree.time$tpi, na.rm = T), length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(isocline_model_tpi, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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

saveRDS(fits, "results/models/non_negative/fits_isoclines_tpi.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_tpi.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_tpi.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_tpi.RDS")

# model with twi and tpi--------------------------
isocline_model_twi_tpi <- bf(sens.prop ~ 1 + twi + tpi + williams_dec + twi:williams_dec + tpi:williams_dec)
coefs <- list()
pred <- list()
fits <- list()
new_preds <- list()

newdata <- expand.grid(
    twi = seq(min(tree.time$twi, na.rm = T), max(tree.time$twi, na.rm = T), length.out = 15),
    tpi = seq(min(tree.time$tpi, na.rm = T), max(tree.time$tpi, na.rm = T), length.out = 15),
    williams_dec = seq(min(tree.time$williams_dec, na.rm = T), max(tree.time$williams_dec, na.rm = T), length.out = 15)
)

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(isocline_model_twi_tpi, data = tree.time %>% filter(yr == yrs[i]), family = gaussian(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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

saveRDS(fits, "results/models/non_negative/fits_isoclines_twi_tpi.RDS")
saveRDS(pred, "results/models/non_negative/pred_isoclines_twi_tpi.RDS")
saveRDS(coefs, "results/models/non_negative/coefs_isoclines_twi_tpi.RDS")
saveRDS(new_preds, "results/models/non_negative/new_preds_isoclines_twi_tpi.RDS")
