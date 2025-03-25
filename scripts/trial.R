library(brms)
library(tidyverse)
library(tidybayes)
library(patchwork) # to look at plots side-by-side

# example data from ?brms::mo

# generate some data
income_options <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(income_options, 100, TRUE),
    levels = income_options, ordered = TRUE
)
mean_ls <- c(30, 60, 70, 75)
ls <- mean_ls[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)

# add some city-level variation
dat$city <- rep(1:10, each = 10)
var_city <- rnorm(10, sd = 10)
dat$ls <- dat$ls + var_city[dat$city]


# plot raw data:
ggplot(dat, aes(income, ls)) +
    geom_point()

ggplot(dat, aes(income, ls, color = city)) +
    geom_point()

# fit a simple monotonic model
get_prior(ls ~ mo(income), data = dat)
fit1 <- brm(ls ~ mo(income), data = dat, chains = 4, cores = 4, backend = "cmdstanr")
summary(fit1)
p_ce <- plot(conditional_effects(fit1))

# now try manually with spread_draws:
get_variables(fit1) # what's available

p_manual_ce <- fit1 %>%
    spread_draws(
        b_Intercept, bsp_moincome,
        simo_moincome1[i]
    ) %>%
    # simplex parameters are grouped by variable `i`, pivot into same row
    pivot_wider(names_from = `i`, names_prefix = "simo_moincome1_", values_from = simo_moincome1) %>%
    mutate(
        # D is equal to number of categories minus 1
        D = 3,
        # for each level of income, get posterior mu
        mu_cat_1 = b_Intercept,
        mu_cat_2 = b_Intercept + (bsp_moincome * D * (simo_moincome1_1)),
        mu_cat_3 = b_Intercept + (bsp_moincome * D * (simo_moincome1_1 + simo_moincome1_2)),
        mu_cat_4 = b_Intercept + (bsp_moincome * D * (simo_moincome1_1 + simo_moincome1_2 + simo_moincome1_3))
    ) %>%
    # at this point or after pivoting again, you could apply the inverse link function on mu
    pivot_longer(starts_with("mu_cat"), names_to = "mu_cat", values_to = "post_mu") %>%
    ggplot(aes(y = post_mu, x = mu_cat)) +
    stat_halfeye()

(p_ce$income + scale_y_continuous(limits = c(20, 81))) +
    (p_manual_ce + scale_y_continuous(limits = c(20, 81)))

# looks equivalent

# fit a model with a random effect------------------------------------

fit6 <- brm(ls ~ mo(income) + (mo(income) | city), data = dat)

summary(fit6)
p_ce <- plot(conditional_effects(fit6, ))

# now try manually with spread_draws:
get_variables(fit1) # what's available
