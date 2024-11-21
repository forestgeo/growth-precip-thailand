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
        DBH_TWI = cor.test(calcDBH_min1, twi)[4]$estimate,
        DBH_TWI_p= cor.test(calcDBH_min1, twi)[3]$p.value,
        CII_TWI = cor.test(cii_min1, twi)[4]$estimate,
        CII_TWI_p= cor.test(cii_min1, twi)[3]$p.value
    )

# these variables are conditionally independent for the most part

# plot the conditional independencies

library(ggpubr)
cond_dep_dbh_twi<-ggscatter(data = tree.time %>% filter (Cno == 15), 
         x = "calcDBH_min1", y = "twi", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "DBH", ylab = "TWI")+
          facet_wrap(~factor(Species, levels = names(sort(table(Species), decreasing=T)))) 
    
png("doc/display/cond_dep_dbh_twi.png", width = 12, height = 12, units="in", res=300)
cond_dep_dbh_twi
dev.off()

# CII | TWI
cond_dep_cii_twi<-ggscatter(data = tree.time %>% filter (Cno == 15), 
         x = "cii_min1", y = "twi", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CII", ylab = "TWI")+
          facet_wrap(~factor(Species, levels = names(sort(table(Species), decreasing=T))) )

png("doc/display/cond_dep_cii_twi.png", width = 12, height = 12, units="in", res=300)
cond_dep_cii_twi
dev.off()

#Species-specific models-----------------------------

#standardise the variables within species and across all

tree.time<-tree.time %>%
dplyr::mutate(calcDBH_min1_scaled = scale(calcDBH_min1, center = TRUE, scale = TRUE),
cii_min1_scaled = scale(cii_min1, center = TRUE, scale = TRUE),
twi_scaled = scale(twi, center = TRUE, scale = TRUE))%>%
group_by(Species) %>%
#scale while retaining original values
dplyr::mutate(calcDBH_min1_scaled_sp = scale(calcDBH_min1, center = TRUE, scale = TRUE),
cii_min1_scaled_sp = scale(cii_min1, center = TRUE, scale = TRUE),
twi_scaled_sp = scale(twi, center = TRUE, scale = TRUE))%>%
ungroup()


# model from DAG

sp_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled_sp + cii_min1_scaled_sp + twi_scaled_sp)

# run the model for the top 10 species for 2015

top_10_sp <- tree.time %>%
    filter(Cno == 15) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) %>%
    head(10)

#plot the distribution of sensitivities for the top 10 species
sens_top10<-ggplot(data = tree.time %>% filter(Species %in% top_10_sp$Species), aes(x = Species, y=sens.prop, fill=Species)) +
    geom_violin() + geom_boxplot(width=0.1, fill="white") +
    facet_wrap(~yr) +
    #ylim(-5, 5)+
    #make color scale viridis
    scale_fill_viridis_d() +
    #make y axis percent
    scale_y_continuous(labels = scales::percent, limits=c(-5, 5)) +
    geom_hline(yintercept=c(-1, 0, 1), lty=2)+
  labs(title = "Distribution of sensitivities for top 10 species", x = "Sensitivity", y = "Density") +
    theme_bw()+
  theme(axis.text.x = element_text(angle=90))

png("doc/display/sens_top10.png", width = 12, height = 8, units="in", res=300)
sens_top10
dev.off()

#pull out the species with rings
ring_sp<-tree.time %>% 
filter(Species %in% c("AFZEXY", "NEOLOB", "TOONCI", "CHUKTA", "MELIAZ"))

remove_sp<-tree.time %>%
dplyr::mutate(Species = "ALL")

ring_sp<-rbind(ring_sp, remove_sp)

library(viridis)

ring_sp$Species<-factor(ring_sp$Species, levels=c("AFZEXY", "NEOLOB", "TOONCI", "CHUKTA", "MELIAZ", "ALL"))

#plot the distribution of sensitivities for the species with rings
sens_ring<-ggplot(data = ring_sp, aes(x = Species, y=sens.prop, 
fill=Species)) +
    geom_violin() +
    geom_boxplot(width=0.1, fill="white") +
    facet_wrap(~yr) +
    #make color scale viridis and one grey
    scale_fill_manual(values = c(viridis(3), "grey80")) +
    #make y axis percent
    scale_y_continuous(labels = scales::percent, limits=c(-5, 5)) +
    geom_hline(yintercept=c(-1, 0, 1), lty=2)+
    #add n for the number of trees
    geom_text(
        #make n for each species
        data = ring_sp %>% group_by(Species, yr) %>% dplyr::summarise(n = n()),
        aes(label = paste0("n = ", n), x = Species, y = 5), size=5.5, hjust = 0.5, vjust = 0)+
  labs(title = "Distribution of sensitivities for species with rings", x = "Sensitivity", y = "Density") +
    theme_bw()+
  theme(axis.text.x = element_text(angle=90))

png("doc/display/sens_ring.png", width = 12, height = 8, units="in", res=300)
sens_ring
dev.off()

# make a function to run the model and return coefs

run_model <- function(data, model) {
    fit <- brm(model, data = data, family = gaussian(), chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("mean", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    return(post_sum)
}

# run this function for all 30 species in 2015

tree.time.2015<-tree.time %>% filter(Cno == 15)

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
    data <- tree.time.2015 %>% filter(Species == sp)
    post_sum <- run_model(data, sp_model)
    coefs[[i]] <- post_sum
}

fit <- brm(sp_model, data = data, family = gaussian(), chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("mean", "lwr", "upr")
    post_sum$param <- rownames(post_sum)

#traceplot for the model from brms
plot(fit)


#unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

#add species names by repeating each element of top_10_sp$Species 8 times
coefs_df$Species <- rep(all_sp$Species, each = 8)

#add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no", 
ifelse(coefs_df$lwr>0, "pos", "neg"))

saveRDS(list(coefs_df, all_sp), "data/HKK-dendro/sensitivity_model.RData")

#plot the mean and 95% credible interval for each parameter

#read the RDS file
coefs_df<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[1]]
all_sp<-readRDS("data/HKK-dendro/sensitivity_model.RData")[[2]]

#first make labels
par_names<-as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_sp<-ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled")), 
aes(x = factor(Species, levels=rev(all_sp$Species)), y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    #make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~param, scales = "free", labeller = par_names) +
    labs(title = "Effect of parameters on growth sensitivity", x = "Species", y = "Mean") +
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
maxDBH_bDBH_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled")), aes(x = maxDBH, y = mean, 
color=factor(signif, levels=c("neg", "no", "post")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color=factor(signif, levels=c("neg", "no", "pos"))), width=0.1) +
    scale_color_manual(values = c("red", "black", "blue")) +
    ggtitle("maxDBH vs DBH effect") +
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bDBH_plot

maxDBH_bCII_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_cii_min1_scaled")), aes(x = maxDBH, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("maxDBH vs CII effect") +
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bCII_plot

maxDBH_bTWI_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_twi_scaled")), aes(x = maxDBH, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    ggtitle("maxDBH vs TWI effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

maxDBH_bTWI_plot

#plot deciduousness against DBH, CII and TWI effet

decid_bDBH_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "no", "pos"))), width=0.1) +
    scale_color_manual(values = c("red", "black", "blue")) +
    xlab("deciduousness")+
    ggtitle("Deciduousness vs DBH effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bDBH_plot

decid_bTWI_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_twi_scaled")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab("deciduousness")+
    ggtitle("Deciduousness vs TWI effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bTWI_plot

decid_bCII_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_cii_min1_scaled")), aes(x = williams_dec, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    xlab("deciduousness")+
    ggtitle("Deciduousness vs CII effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

decid_bCII_plot

#plot growth rate against DBH, CII and TWI effect

gr_bDBH_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled")), aes(x = median_inc, y = mean, col=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "no", "pos"))), width=0.01) +
    scale_color_manual(values = c("red", "black", "blue")) +
    xlab("median of annualised increment")+
    ggtitle("Growth rate vs DBH effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

gr_bDBH_plot

gr_bCII_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_cii_min1_scaled")), aes(x = median_inc, y = mean, col=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "no", "pos"))), width=0.01) +
    scale_color_manual(values = c("red", "black", "blue")) +
    xlab("median of annualised increment")+
    ggtitle("Growth rate vs CII effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

gr_bCII_plot

gr_bTWI_plot <- ggplot(data = coefs_df %>% filter(param %in% c("b_twi_scaled")), aes(x = median_inc, y = mean, col=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "no", "pos"))), width=0.01) +
    scale_color_manual(values = c("red", "black", "blue")) +
    xlab("median of annualised increment")+
    ggtitle("Growth rate vs TWI effect") + guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

gr_bTWI_plot

#write these as pngs
library(gridExtra)
png("doc/display/coefs_sp_values.png", width = 12, height = 12, units="in", res=300)
grid.arrange(maxDBH_bDBH_plot, maxDBH_bCII_plot, maxDBH_bTWI_plot, 
decid_bDBH_plot, decid_bCII_plot, decid_bTWI_plot, 
gr_bDBH_plot, gr_bCII_plot, gr_bTWI_plot, ncol=3)
dev.off()

#run model for 2010 data

tree.time.2010<-tree.time %>% filter(Cno == 5)

all_sp <- tree.time %>%
    filter(Cno == 5) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n)) 

coefs <- list()

for (i in 1:nrow(all_sp)) {
    print(paste0("Running model for species ", i, " : ", all_sp$Species[i]))
    sp <- all_sp$Species[i]
    data <- tree.time.2010 %>% filter(Species == sp)
    post_sum <- run_model(data, sp_model)
    coefs[[i]] <- post_sum
}

#unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

#add species names by repeating each element of top_10_sp$Species 8 times
coefs_df$Species <- rep(all_sp$Species, each = 8)

#add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no", 
ifelse(coefs_df$lwr>0, "pos", "neg"))

saveRDS(list(coefs_df, all_sp), "data/HKK-dendro/sensitivity_model_2010.RData")

#first make labels
par_names<-as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_sp<-ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled")), 
aes(x = factor(Species, levels=rev(all_sp$Species)), y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    #make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~param, scales = "free", labeller = par_names) +
    labs(title = "Effect of parameters on growth sensitivity", x = "Species", y = "Mean") +
    guides(color="none")+ theme_bw()+
    coord_flip()

coefs_sp

#write these as pngs
png("doc/display/coefs_sp_2010.png", width = 8, height = 8, units="in", res=300)
coefs_sp
dev.off()

#mapply to make the plots of coefs against species variables

#merge with sp_vars
coefs_df<-merge(coefs_df, sp_vars, by="Species", all.x=TRUE)

#function to make one plot

coef_var_plot<- function (coefdf, par, parname, var, varname) {
    ggplot(data = coefdf %>% dplyr::filter(param == par), aes_string(x = var, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() #+
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "black", "blue")) +
    ggtitle(paste0(parname, " effect with ", varname)) +
    xlab(paste0(varname)) + ylab ("mean effect")+
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()
}

p1<-coef_var_plot(coefdf=coefs_df, par="b_calcDBH_min1_scaled", parname="DBH", var="maxDBH", varname="maxDBH")

# TODO: fix this function and plot species vars

#run model for 2020

tree.time.2020<-tree.time %>% filter(Cno == 25)

all_sp <- tree.time %>%
    filter(Cno == 25) %>%
    group_by(Species) %>%
    dplyr::summarise(
        n = n()
    ) %>%
    arrange(desc(n))

coefs <- list()

for (i in 1:nrow(all_sp)) {
    print(paste0("Running model for species ", i, " : ", all_sp$Species[i]))
    sp <- all_sp$Species[i]
    data <- tree.time.2020 %>% filter(Species == sp)
    post_sum <- run_model(data, sp_model)
    coefs[[i]] <- post_sum
    pct_complete <- i/nrow(all_sp) * 100
    print(paste0("Progress:", rep("=", pct_complete %/% 10), " ", pct_complete, "%"))
}

#unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

#add species names by repeating each element of top_10_sp$Species 8 times
coefs_df$Species <- rep(all_sp$Species, each = 8)

#add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
ifelse(coefs_df$lwr>0, "pos", "neg"))

saveRDS(list(coefs_df, all_sp), "data/HKK-dendro/sensitivity_model_2020.RData")

#first make labels
par_names<-as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_sp<-ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled")),
aes(x = factor(Species, levels=rev(all_sp$Species)), y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    #make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "black")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~param, scales = "free", labeller = par_names) +
    labs(title = "Effect of parameters on growth sensitivity", x = "Species", y = "Mean") +
    guides(color="none")+ theme_bw()+
    coord_flip()

coefs_sp

#write these as pngs
png("doc/display/coefs_sp_2020.png", width = 8, height = 8, units="in", res=300)
coefs_sp
dev.off()


# Tree is a tree model (model with all trees together)

#first test conditional dependencies DBH | TWI and CII | TWI
#make a long dataframe for plotting conditional independencies

cond_dep_all <- tree.time %>%
    #filter(Cno == 15) %>%
    group_by(Cno) %>%
    pivot_longer(cols = c("calcDBH_min1", "cii_min1"), names_to = "var", values_to = "value")%>%
    dplyr::mutate(varnames=ifelse(var=="calcDBH_min1", "DBH", "CII"))
    
# these variables are conditionally independent for the most part

# plot the conditional independencies

library(ggpubr)
cond_dep_all_plot<-ggscatter(data = cond_dep_all, 
         x = "value", y = "twi", 
          add = "reg.line", conf.int = TRUE, alpha=0.3,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "variable", ylab = "TWI")+
          facet_wrap(var~yr, scales="free")
    
png("doc/display/cond_dep_alltrees.png", width = 8, height = 8, units="in", res=300)
cond_dep_all_plot
dev.off()



#first find sensitivity relative to the mean growth across all trees of that species

# #make a new column for the mean growth rate of each species
# tree.time<-tree.time %>% group_by(Species, yr) %>% 
# dplyr::mutate(sp_mean_inc = mean(inc_annual, na.rm = TRUE))%>%
# ungroup()%>%
# dplyr::mutate(sens_prop_sp = (inc_annual - sp_mean_inc)/sp_mean_inc)


#tree is a tree model

tree_model <- bf(sens.prop ~ 1 + calcDBH_min1_scaled + cii_min1_scaled + twi_scaled + (1|Species))

#run the model

coefs <- list()

for (i in 1:3){    
    yrs<-c(2010, 2015, 2020)
    fit <- brm(tree_model, data = tree.time%>%filter(yr==yr[i]), family = gaussian(), chains = 4, cores = 4)
    post <- posterior_samples(fit)
    post_sum <- as.data.frame(t(apply(post, 2, quantile, probs = c(.5, .05, .95))))
    colnames(post_sum) <- c("mean", "lwr", "upr")
    post_sum$param <- rownames(post_sum)
    coefs[[i]] <- post_sum
}

#unlist coefs, make a df and add species names
coefs_df <- do.call(rbind, coefs)
head(coefs_df)

#add species names by repeating each element of top_10_sp$Species 8 times
coefs_df$yr <- rep(yrs, each = 39)

#add a column for significance
coefs_df$signif <- ifelse(coefs_df$lwr < 0 & coefs_df$upr > 0, "no",
ifelse(coefs_df$lwr>0, "pos", "neg"))

saveRDS(coefs_df, "data/HKK-dendro/sensitivity_model_treeisatree.RData")

#first make labels
par_names<-as_labeller(c("b_calcDBH_min1_scaled" = "DBH effect", "b_cii_min1_scaled" = "CII effect", "b_twi_scaled" = "TWI effect"))

coefs_tree<-ggplot(data = coefs_df %>% filter(param %in% c("b_calcDBH_min1_scaled", "b_cii_min1_scaled", "b_twi_scaled")),
aes(x = param , y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    scale_x_discrete(labels=par_names)+
    #make error bars with narrow heads
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop=FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    #facet_grid(~param, scales = "free", labeller = par_names) +
    facet_grid(~factor(yr, levels=c(2010, 2015, 2020)), scales="free")+
    labs(title = "Effect of parameters on growth sensitivity", x = "", y = "coefficient") +
    guides(color="none")+ theme_bw()+
    coord_flip()

coefs_tree

#write these as pngs
png("doc/display/coefs_tree_allyrs.png", width = 6, height = 4, units="in", res=300)
coefs_tree
dev.off()

#make plots for the species random effects

#first subset the coefs_df for the random effects
ranef_df<-coefs_df %>% filter(grepl("r_Species", param))
#add a column with the species names by removing the "r_Species[" and ",Intercept]" from the param column
ranef_df$Species<-gsub("r_Species\\[|\\,Intercept\\]", "", ranef_df$param)

#plot the random effects
ranef_plot<-ggplot(data = ranef_df, aes(x = Species, y = mean, col=factor(signif, levels=c("neg", "pos", "no")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, col=factor(signif, levels=c("neg", "pos", "no"))), width=0.1) +
    scale_color_manual(values = c("red", "blue", "grey40"), drop=F) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(~factor(yr, levels=c(2010, 2015, 2020)), scales="free")+
    labs(title = "Species random effects", x = "Species", y = "coefficient") +
    guides(color="none")+ theme_bw()+
    coord_flip()

png("doc/display/ranefs_tree_allyrs.png", width = 8, height = 8, units="in", res=300)
ranef_plot
dev.off()

# plotting the random effects against species characteristics

#merge the sp_vars with the coefs_df
ranef_df<-merge(ranef_df, sp_vars, by="Species", all.x=TRUE)

head(ranef_df)

#plot ranefs against maxDBH, deciduousness and growth rate

#first make long dataframe for plotting
ranef_long_2015<-ranef_df %>% filter(yr==2015) %>%
pivot_longer(cols = c("maxDBH", "williams_dec", "median_inc"), names_to = "sp_var", values_to = "value")

#make labels
sp_var_names<-as_labeller(c("maxDBH" = "maxDBH", "williams_dec" = "deciduousness", "median_inc" = "growth rate"))

ranef_sp_plot <- ggplot(data = ranef_long_2015, aes(x = value, y = mean, 
color=factor(signif, levels=c("neg", "no", "pos")))) +
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color=factor(signif, levels=c("neg", "no", "pos"))), width=0.1) +
    scale_color_manual(values = c("red", "grey40", "blue")) +
    ggtitle("Species variables and random effect") + ylab("species random effect")+
    facet_wrap(~sp_var, scales="free", labeller = sp_var_names)+
    guides(color="none") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw()

png("doc/display/ranef_sp_plot.png", width = 10, height = 4, units="in", res=300)
ranef_sp_plot
dev.off()