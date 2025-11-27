# models considering cii as an ordered factor

# Load required libraries---------------------
required_pkgs <- c(
    "tidyverse", "brms", "patchwork", "performance"
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
yrs <- c(2010, 2015, 2020)

# for this model, we have to remove trees with NA values because TPI cannot be calculated for plot margins

tree.time <- tree.time %>%
    filter(!is.na(tpi_scaled))

# directories--------

plot_dir <- "results/plots/orderedcii_tpi"
models_dir <- "results/models/orderedcii_tpi"

# ensure output directory exists
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message(paste0("Created directory: ", plot_dir))
}

if (!dir.exists(models_dir)) {
    dir.create(models_dir, recursive = TRUE)
    message(paste0("Created directory: ", models_dir))
}


# define and run the model
# Model 1: varying slopes without species scaling------------------------------

tree_model_rel_spre <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + mo(cii_min1) + tpi_scaled + (1 + calcDBH_min1_scaled + tpi_scaled + mo(cii_min1) | Species))
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
