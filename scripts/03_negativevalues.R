#assessing and understanding the negative values in the dataset

# Load required libraries---------------------
library(tidyverse)
library(brms)

# load data--------------------------
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

#standardise the variables within species and across all--------------------

tree.time<-tree.time %>%
dplyr::mutate(calcDBH_min1_scaled = scale(calcDBH_min1, center = TRUE, scale = TRUE),
cii_min1_scaled = scale(cii_min1, center = TRUE, scale = TRUE),
cii_min1_scaled = ifelse(cii_min1_scaled=="NaN", 0, cii_min1_scaled),
twi_scaled = scale(twi, center = TRUE, scale = TRUE),
median_inc_scaled = scale(median_inc, center = TRUE, scale=TRUE))%>%
group_by(Species) %>%
#scale while retaining original values
dplyr::mutate(calcDBH_min1_scaled_sp = scale(calcDBH_min1, center = TRUE, scale = TRUE),
cii_min1_scaled_sp = scale(cii_min1, center = TRUE, scale = TRUE),
cii_min1_scaled_sp = ifelse(cii_min1_scaled_sp=="NaN", 0, cii_min1_scaled_sp),
twi_scaled_sp = scale(twi, center = TRUE, scale = TRUE))%>%
ungroup()%>%
#remove large outliers for each year
group_by(yr)%>%
#find sens.prop values that are 3 sds from the mean
dplyr::mutate(sens.prop = ifelse(sens.prop > mean(sens.prop, na.rm=TRUE) + 3*sd(sens.prop, na.rm=TRUE), NA, sens.prop)) %>%
filter(!is.na(sens.prop) & !is.na(cii_min1) & !is.na(calcDBH_min1) & !is.na(twi))%>%
ungroup()

# how many trees have negative values for inc_annual each year?------------------------

# make a new column for negative or not
tree.time<-tree.time%>%
dplyr::mutate(neg_inc_annual = ifelse(inc_annual < 0, 1, 0))

#plot these for each year
neg_plot<-ggplot(tree.time, aes(x=yr, fill=factor(neg_inc_annual))) +
geom_bar() +
scale_fill_manual(values=c("grey40", "red")) +
theme_minimal() +
labs(title="Number of trees with negative increment values each year",
x="Year",
y="Number of trees") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("doc/display/neg_inc_annual.png", width=8, height=6, units="in", res=300)
neg_plot
dev.off()

#what proportion of trees in each species has negative values for inc_annual in 2010 and 2015?------------------------


neg_species<-tree.time%>%
filter(yr %in% c(2010, 2015))%>%
group_by(Species, yr)%>%
dplyr::summarise(inc_annual = sum(neg_inc_annual)/n())

#plot these
neg_species_plot<-ggplot(neg_species, aes(x=Species, y=inc_annual, fill=factor(yr))) +
geom_bar(stat="identity", position="dodge") +
scale_fill_manual(values=c("goldenrod", "seagreen")) +
theme_minimal() +
labs(title="Proportion of trees with negative increment values in 2010 and 2015",
x="Species",
y="Proportion of trees") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))


png("doc/display/neg_species.png", width=8, height=6, units="in", res=300)
neg_species_plot
dev.off()


#does TWI predict negative values for inc_annual?------------------------

#plot the relationship between twi and neg_inc_annual
twi_plot<-ggplot(tree.time%>%filter(yr %in% c(2010, 2015)), aes(x=twi, y=neg_inc_annual)) +
geom_point() +
theme_minimal() +
facet_wrap(~yr) +
labs(title="Relationship between TWI and negative increment values",
x="TWI",
y="Negative increment values") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("doc/display/twi_plot.png", width=8, height=6, units="in", res=300)
twi_plot
dev.off()

#model the relationship
twi_model<-brm(neg_inc_annual ~ twi, data=tree.time%>%filter(yr==2010), family=bernoulli(), chains=4, iter=2000)

#run a glm
twi_glm<-glm(neg_inc_annual ~ twi+yr, data=tree.time%>%filter(yr%in%c(2010, 2015)), family=binomial())
summary(twi_glm)

size_glm<-glm(neg_inc_annual ~ calcDBH_min1, data=tree.time%>%filter(yr%in%c(2010, 2015)), family=binomial())
summary(size_glm)

#plot the relationship between size and neg_inc_annual
size_plot<-ggplot(tree.time%>%filter(yr %in% c(2010, 2015)), aes(x=calcDBH_min1, y=neg_inc_annual)) +
geom_point() +
theme_minimal() +
facet_wrap(~yr) +
labs(title="Relationship between size and negative increment values",
x="Size",
y="Negative increment values") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("doc/display/size_plot.png", width=8, height=6, units="in", res=300)
size_plot
dev.off()

#model the relationship with size, exposure and twi
shrink_model<-glm(neg_inc_annual ~ calcDBH_min1 + twi + cii_min1 + yr, data=tree.time%>%filter(yr%in% c(2010, 2015)), family=binomial())
summary(shrink_model)
