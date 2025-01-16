# models using the zero growth framework

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
        avg_inc_tree_scaled = scale(avg_inc_tree, center = TRUE, scale = TRUE),
        sens.bin = ifelse(sens.prop > 0, 1, 0)
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

# sensitivities for top 10 species using zero growth model

top_10_sp <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)


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


# model probablity of growth ---------------------------------

intercept_model <- bf(sens.bin ~ 1 + (1 | Species))

coefs <- list()
pred <- list()
fits <- list()

yrs <- c(2010, 2015)
for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(intercept_model, data = tree.time %>% filter(yr == yrs[i]), family = bernoulli(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
saveRDS(fits, "results/models/zero_growth/fits_interceptonly.RDS")


# save the coefs and predictions
coefs_df <- do.call(rbind, coefs)
pred_df <- do.call(rbind, pred)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/zero_growth/sensitivity_model_intercept.RData")

saveRDS(pred_df, "results/models/zero_growth/predictions_intercept.RData")


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

png("results/plots/zero_growth/ranefs_intercept.png", width = 8, height = 8, units = "in", res = 300)
ranef_plot
dev.off()

library(patchwork)
# png("results/plots/non_negative/pp_intercept.png", width = 8, height = 4, res = 300, units = "in")
# pp_check(fits[[1]]) + pp_check(fits[[2]])
# dev.off()

# plot predicted vs observed distributions within species

# make long df for plotting
pred_obs_sp <- pred_df %>%
    pivot_longer(c("sens.bin", "median")) %>%
    dplyr::mutate(name = ifelse(name == "sens.bin", "obs", "pred"))

# top10_predobs <- ggplot(data = pred_obs_sp %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = value, fill = Species, color = name)) +
#     geom_violin() + # geom_boxplot(width=0.1, fill="white") +
#     facet_wrap(~yr) +
#     # ylim(-5, 5)+
#     # make color scale viridis
#     scale_fill_viridis_d() +
#     # scale_color_manual(values=c())+
#     # make y axis percent
#     #scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
#     geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
#     labs(title = "Distribution of sensitivities for top 10 species", x = "Sensitivity", y = "Density") +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90))
# top10_predobs

# plot violin plots of observed sensitivity with a line for predictions

# add intercept to ranef_df
# make intercept vector
intercepts <- coefs_df %>%
    filter(param %in% "Intercept") %>%
    select(median) %>%
    pull(median)
intercepts <- rep(intercepts, each = nrow(ranef_df) / 2)
ranef_df$intercepts <- intercepts + ranef_df$median

top10_predobs <- ggplot(data = pred_df %>% filter(Species %in% top_10_sp$Species, yr %in% yrs), aes(x = Species, y = sens.bin, fill = Species)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~yr) +
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    # scale_y_continuous(labels = scales::percent, limits = c(-5, 5)) +
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


png("results/plots/zero_growth/ranefs_intercept_bysp_new.png", width = 12, height = 8, units = "in", res = 300)
top10_predobs
dev.off()



# M2 : tree is a tree model--------------------------------------

## define and run model -------------------------

tree_model <- bf(sens.bin ~ 1 + calcDBH_min1_scaled + cii_min1 + twi_scaled + (1 | Species))

# run the model

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    # yrs<-c(2010, 2015, 2020)
    fit <- brm(tree_model,
        # data = tree.time %>% filter(yr == yrs[i]), family = skew_normal(),
        data = tree.time %>% filter(yr == yrs[i]), family = bernoulli(),
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

saveRDS(fits, "results/models/zero_growth/fits_treeisatree.RDS")

summary(fits[[1]])

png("results/plots/zero_growth/fits_treeisatree.png")
plot(fits[[1]], variable = "^b", regex = T)
dev.off()

library(patchwork)
png("results/plots/zero_growth/pp_treeisatree.png", width = 8, height = 4, units = "in", res = 300)
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

saveRDS(coefs_df, "results/models/zero_growth/sensitivity_model_treeisatree.RData")

saveRDS(pred_df, "results/models/zero_growth/predictions_treeisatree.RData")

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
    facet_grid(~ factor(yr, levels = c(2010, 2015)), scales = "free") +
    labs(title = "Effect of parameters on growth occurrence", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

coefs_tree

# write these as pngs
png("results/plots/zero_growth/coefs_tree.png", width = 6, height = 4, units = "in", res = 300)
coefs_tree
dev.off()




# plot predicted vs observed
pred_obs <- pred_df %>%
    pivot_longer(c("sens.bin", "median")) %>%
    dplyr::mutate(name = ifelse(name == "sens.bin", "obs", "pred"))

pred_plot <- ggplot(data = pred_df, aes(x = sens.bin, y = median, col = Species)) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col = Species), width = 0.1) +
    scale_color_viridis_d() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~yr, scales = "free") +
    labs(title = "Predictions vs observations", x = "Observations", y = "Predictions") +
    guides(color = "none") +
    theme_bw()

pred_plot

# make a confusion matrix


# plot this as tiles
confusion_plot <- ggplot(data = pred_df, aes(x = sens.bin, y = pred, fill = accuracy)) +
    geom_tile() +
    geom_text(aes(label = round(accuracy, 2)), vjust = 1) +
    scale_fill_viridis_c() +
    labs(title = "Model accuracy", x = "Year", y = "") +
    theme_bw()
confusion_plot

png("results/plots/zero_growth/pred_vs_obs_treeisatree.png", width = 8, height = 4, units = "in", res = 300)
pred_plot
dev.off()


# M3 : tree is a tree model with species random effect----------------------

tree_model_rel <- bf(sens.bin ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 | Species))

coefs <- list()
pred <- list()
fits <- list()

yrs <- c(2010, 2015)

for (i in 1:length(yrs)) {
    fit <- brm(tree_model_rel, data = tree.time %>% filter(yr == yrs[i]), family = bernoulli(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
saveRDS(fits, "results/models/zero_growth/fits_treeisatreerel.RDS")

library(patchwork)
# plot posterior predictive checks
png("results/plots/zero_growth/pp_treeisatreerel.png", width = 8, height = 4, res = 300, units = "in")
pp_check(fits[[1]]) + pp_check(fits[[2]])
dev.off()

# plot chains

p1 <- plot(fits[[1]], variable = "^b", regex = T)
p2 <- plot(fits[[2]], variable = "^b", regex = T)
p <- cbind(p1, p2)
plot(p)

library(gtable)
library(gridExtra)

# png("doc/display/diagnostics_treeisatreerel.png", width = 8, height = 4, units = "in", res = 300)

# dev.off()



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

saveRDS(coefs_df, "results/models/zero_growth/sensitivity_model_treeisatreerel.RData")
saveRDS(pred_df, "results/models/zero_growth/predictions_treeisatree_rel.RData")

# M4: model with varying slopes---------------------------------------

all_slope_model <- bf(sens.bin ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp + (1 + cii_min1_scaled_sp + calcDBH_min1_scaled_sp + twi_scaled_sp | Species))

coefs <- list()
pred <- list()
fits <- list()

for (i in 1:length(yrs)) {
    fit <- brm(all_slope_model, data = tree.time %>% filter(yr == yrs[i]), family = bernoulli(), iter = 3000, warmup = 1000, chains = 4, cores = 4)
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
saveRDS(fits, "results/models/zero_growth/fits_spre.RDS")

# unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

pred_df <- do.call(rbind, pred)
head(pred_df)

# add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
    ifelse(coefs_df$lwr > 0, "pos", "neg")
)

saveRDS(coefs_df, "results/models/zero_growth/sensitivity_model_spre.RData")
saveRDS(pred_df, "results/models/zero_growth/predictions_spre.RData")


# Combined plots---------------------------------------
# make a plot with three different coefs together

coefs_df_tree <- readRDS("results/models/zero_growth/sensitivity_model_treeisatree.RData")
coefs_df_tree_rel <- readRDS("results/models/zero_growth/sensitivity_model_treeisatreerel.RData")
coefs_df_spre <- readRDS("results/models/zero_growth/sensitivity_model_spre.RData")

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
    labs(title = "Effect of parameters on growth occurrence", x = "", y = "coefficient") +
    guides(color = "none") +
    theme_bw() +
    coord_flip()

png("results/plots/zero_growth/coefs_all_plot.png", width = 8, height = 6, units = "in", res = 300)
coefs_all_plot
dev.off()
