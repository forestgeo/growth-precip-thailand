#models to assess sensitivity of tree growth to precipitation

# Load required libraries
library(tidyverse)
library(brms)

# load data
rm(list=ls())
datasets<-readRDS("data/HKK-dendro/sensitivity_data.RData")

tree.time<-datasets$tree.time
tree_vars<-datasets$tree_vars
sp_vars<-datasets$sp_vars

colnames(sp_vars)[1]<-"Species"

#merge tree_vars and sp_vars
tree_vars<-merge(tree_vars, sp_vars, by="Species", all.x=TRUE)
#merge tree_vars and tree.time
tree.time<-merge(tree.time, tree_vars, by="Tag", all.x=TRUE)

colnames(tree.time)

# test conditional independences for each species

# DBH | TWI
# CII | TWI

#2015 data only

cond_dep <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        DBH_TWI = cor(calcDBH_min1, twi),
        CII_TWI = cor(cii_min1, twi)
    )

# these variables are conditionally independent for the most part

# plot the conditional independencies

# DBH | TWI
cond_dep_dbh_twi <- ggplot() +
    geom_point(data = tree.time %>% filter (Cno == 15), 
    aes (x = calcDBH_min1, y = twi)) +
    geom_abline(intercept = 0, slope = 0) +
    facet_wrap(~Species) +
    geom_text(data = cond_dep, aes(label = paste0("cor = ", round(DBH_TWI,2))), x = 20, y = 10, hjust = 0, vjust = 0)

png("doc/display/cond_dep_dbh_twi.png", width = 8, height = 8, units="in", res=300)
cond_dep_dbh_twi
dev.off()

# CII | TWI
cond_dep_cii_twi <- ggplot() +
    geom_point(data = tree.time, aes (x = cii_min1, y = twi)) +
    geom_abline(intercept = 0, slope = 0) +
    facet_wrap(~Species) +
    geom_text(data = cond_dep, aes(label = paste0("cor = ", round(CII_TWI,2))), x = 2, y = 10, hjust = 0, vjust = 0)

png("doc/display/cond_dep_cii_twi.png", width = 8, height = 8, units="in", res=300)
cond_dep_cii_twi
dev.off()


# model from DAG

sp_model <- bf(sens.prop ~ 1 + calcDBH_min1 + cii_min1 + twi)

# run the model for the top 10 species for 2015

top_10_sp <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)

# make a function to run the model and return coefs

run_model <- function(data, model) {
    fit <- brm(model, data = data, family = gaussian(), chains = 2, cores = 2)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("mean", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    return(post_sum)
}

# run this function for the top 10 species

coefs <- list()

for (i in 1:nrow(top_10_sp)) {
    sp <- top_10_sp$Species[i]
    data <- tree.time %>% filter(Species == sp)
    post_sum <- run_model(data, sp_model)
    coefs[[i]] <- post_sum
}

#unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

#add species names by repeating each element of top_10_sp$Species 8 times
coefs_df$Species <- rep(top_10_sp$Species, each = 8)

#add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no", 
ifelse(coefs_df$lwr>0, "pos", "neg"))


#plot the mean and 95% credible interval for each parameter

coefs_sp<-ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1", "b_cii_min1", "b_twi")), 
aes(x = Species, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    #make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~param, scales = "free_y") +
    guides(color="none")+ theme_bw()+
    coord_flip()

coefs_sp

#write these as pngs
png("doc/display/coefs_sp.png", width = 8, height = 8, units="in", res=300)
coefs_sp
dev.off()

#plot these coefs against species variables

#merge the sp_vars with the coefs_df
coefs_df<-merge(coefs_df, sp_vars, by="Species", all.x=TRUE)

#plot maxDBH against DBH, exposure and TWI effects
maxDBH_bDBH_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1")), aes(x = maxDBH, y = mean, 
color=factor(signif, levels=c("neg", "no", "post")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color=factor(signif, levels=c("neg", "no", "pos"))), width=0.1) +
    scale_color_manual(values = c("red", "black", "blue")) +
    ggtitle("maxDBH vs DBH effect") +
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bDBH_plot

maxDBH_bCII_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_cii_min1")), aes(x = maxDBH, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("maxDBH vs CII effect") +
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bCII_plot

maxDBH_bTWI_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_twi")), aes(x = maxDBH, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("maxDBH vs TWI effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bTWI_plot

#plot deciduousness against DBH, CII and TWI effet

decid_bDBH_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "no", "pos"))), width=0.1) +
    scale_color_manual(values = c("red", "black", "blue")) +
    ggtitle("Deciduousness vs DBH effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bDBH_plot

decid_bTWI_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_twi")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("Deciduousness vs TWI effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bTWI_plot

decid_bCII_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_cii_min1")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("Deciduousness vs CII effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bCII_plot

#write these as pngs
library(gridExtra)
png("doc/display/coefs_sp_values.png", width = 12, height = 8, units="in", res=300)
grid.arrange(maxDBH_bDBH_plot, maxDBH_bCII_plot, maxDBH_bTWI_plot, 
decid_bDBH_plot, decid_bCII_plot, decid_bTWI_plot, ncol=3)
dev.off()