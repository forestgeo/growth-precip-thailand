# trying sensitivity with tree rings

# load libraries
# load the packages used by this script
# library(plyr)
# library(tidyverse)
library(dplR)
# library(reshape2)
# library(janitor)
# library(lubridate)
# library(httr)

# global variables
start_date <- 1960
end_date <- 2022


# load tree ring data

#-------core data-------------------------------------------------------------------------
# Read in the tree ring data [PRIVATE]
#--------------------------------------------------------------------------------=

# #reading the GitHub PAT (stored privately and locally)
# pat<-as.character(read.table("data/HKK/git_pat.txt")[1,1])


# remotes:::download("data/HKK/ax.txt", "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/TreeRing/Raw/Vlam_HKKcores/RWL_HKK_2013_IN/Ax.IN.txt",
#                    auth_token = pat, quiet=F)

# remotes:::download("data/HKK/ct.txt", "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/TreeRing/Raw/Vlam_HKKcores/RWL_HKK_2013_IN/Ct.IN2.txt",
#                    auth_token = pat, quiet=F)

# remotes:::download("data/HKK/ma.txt", "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/TreeRing/Raw/Vlam_HKKcores/RWL_HKK_2013_IN/Ma.IN.txt",
#                    auth_token = pat, quiet=F)

# remotes:::download("data/HKK/tc.txt", "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/TreeRing/Raw/Vlam_HKKcores/RWL_HKK_2013_IN/Tc.IN.TXT",
#                    auth_token = pat, quiet=F)


ax <- read.rwl("data/HKK-dendro/ax.txt")
tc <- read.rwl("data/HKK-dendro/tc.txt")
ma <- read.rwl("data/HKK-dendro/ma.txt")
ct <- read.rwl("data/HKK-dendro/ct.txt")


rw.vlam <- data.frame(combine.rwl(combine.rwl(combine.rwl(ax, tc), ma), ct))
rw.vlam$year <- as.numeric(row.names(rw.vlam))

rw.hkk <- reshape2::melt(rw.vlam, id.vars = "year", na.rm = T) # using reshape2 package to "melt" the dataframe into 3 columns: year, tree, and increment measurement
colnames(rw.hkk) <- c("year", "coreID", "core_measure") # rename the columns
head(rw.hkk) # view the first few records

set.seed(123)

hkk.cores <- rw.hkk %>%
    group_by(coreID) %>%
    dplyr::summarise(
        treeID = paste0("HKK.", str_sub(first(coreID), start = 1, end = -2)),
        nmeas = length(coreID),
        start.yr = min(year),
        end.yr = max(year)
    ) %>%
    ungroup() %>%
    group_by(treeID) %>%
    dplyr::mutate(keep = ifelse(coreID %in% sample(coreID, 1), 1, 0))

# see issue #7 https://github.com/EcoClimLab/bayesian-data-fusion/issues/7

rw.hkk.long <- rw.hkk %>%
    dplyr::mutate(
        Year = year,
        diainc = (core_measure) * 2 / 10,
        sp = str_sub(coreID, start = 1, end = 2),
        treeID = paste0("HKK.", str_sub(coreID, start = 1, end = -2)),
        sp = ifelse(sp == "ax", "AFZEXY",
            ifelse(sp == "tc", "TOONCI",
                ifelse(sp == "ma", "MELIAZ", "CHUKTA")
            )
        )
    ) %>% # convert to diameter and cm
    # group_by(sp, treeID, year) %>%
    # dplyr::summarise(
    #  diainc=mean(diainc)
    # )%>%
    filter(coreID %in% hkk.cores[hkk.cores$keep == 1, ]$coreID) %>%
    filter(start_date <= year & year <= end_date) %>% # start to end date
    filter(!is.na(diainc)) %>%
    select(sp, treeID, year, diainc) %>%
    ungroup()


# making a metadata file
rw.meta <- rw.hkk.long %>%
    group_by(treeID) %>%
    dplyr::summarise(
        sp = first(sp),
        mean.inc = mean(diainc, na.rm = T)
    )

rw.hkk.2010 <- rw.hkk.long %>%
    filter(year == 2010)

rw.hkk.2010 <- merge(rw.hkk.2010, select(rw.meta, -sp), by = "treeID", all.x = T)

rw.hkk.2010 <- rw.hkk.2010 %>%
    dplyr::mutate(sens.prop = (diainc - mean.inc) / mean.inc)

# plot the distribution of these across species
sens.all.tr <- ggplot(data = rw.hkk.2010, aes(x = sens.prop)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed") +
    # facet_wrap(~yr) +
    # xlim(-5, 5)+
    labs(title = "Distribution of drought sensitivities \nfor all trees", x = "Sensitivity", y = "Density") +
    theme_bw()

png("doc/display/sens_all_tr.png", width = 4, height = 4, units = "in", res = 300)
sens.all.tr
dev.off()

library(tidyverse)

# load data--------------------------
datasets <- readRDS("data/HKK-dendro/sensitivity_data.RData")

load("data/HKK-dendro/DendroByBand.RData")

tree.time <- datasets$tree.time
tree_vars <- datasets$tree_vars
sp_vars <- datasets$sp_vars

colnames(sp_vars)[1] <- "Species"

# merge tree_vars and sp_vars
# tree_vars <- merge(tree_vars, sp_vars, by = "Species", all.x = TRUE)
# merge tree_vars and tree.time
# tree.time <- merge(tree.time, tree_vars, by = "Tag")

tree.time <- merge(tree.time, select(trees, c(Tag, SPCODE.UPDATE)), by = "Tag", all.x = T)

tree.time <- tree.time %>%
    dplyr::rename(Species = SPCODE.UPDATE)

# standardise the variables within species and across all--------------------

tree.time <- tree.time %>%
    dplyr::mutate(
        calcDBH_min1_scaled = scale(calcDBH_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled = scale(cii_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled = ifelse(cii_min1_scaled == "NaN", 0, cii_min1_scaled),
        # twi_scaled = scale(twi, center = TRUE, scale = TRUE),
        # median_inc_scaled = scale(median_inc, center = TRUE, scale = TRUE),
        # avg_inc_tree_scaled = scale(avg_inc_tree, center = TRUE, scale = TRUE)
    ) %>%
    group_by(Species) %>%
    # scale while retaining original values
    dplyr::mutate(
        calcDBH_min1_scaled_sp = scale(calcDBH_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled_sp = scale(cii_min1, center = TRUE, scale = TRUE),
        cii_min1_scaled_sp = ifelse(cii_min1_scaled_sp == "NaN", 0, cii_min1_scaled_sp),
        # twi_scaled_sp = scale(twi, center = TRUE, scale = TRUE)
    ) %>%
    ungroup() %>%
    # remove large outliers for each year
    group_by(yr) %>%
    # find sens.prop values that are 3 sds from the mean
    #    dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 3 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 3 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
    dplyr::mutate(sens.prop = ifelse(mean(sens.prop, na.rm = TRUE) + 4 * sd(sens.prop, na.rm = TRUE) > sens.prop & sens.prop > mean(sens.prop, na.rm = TRUE) - 4 * sd(sens.prop, na.rm = TRUE), sens.prop, NA)) %>%
    # filter(!is.na(sens.prop) & !is.na(cii_min1) & !is.na(calcDBH_min1) & !is.na(twi)) %>%
    filter(!is.na(sens.prop))
ungroup()

sens_dendro_tr_2010 <- tree.time %>%
    filter(yr == 2010) %>%
    dplyr::rename(sp = Species, year = yr, diainc = inc_annual) %>%
    select(sp, treeID, year, diainc, sens.prop) %>%
    dplyr::mutate(source = "dendroband")

sens_dendro_tr_2010 <- sens_dendro_tr_2010[which(sens_dendro_tr_2010$sp %in% unique(rw.hkk.2010$sp)), ]

rw.hkk.2010 <- rw.hkk.2010 %>%
    select(-mean.inc) %>%
    dplyr::mutate(source = "tree.rings")

both.2010 <- rbind(rw.hkk.2010, sens_dendro_tr_2010)

sens_trsp <- ggplot(data = both.2010, aes(x = sp, y = sens.prop, fill = sp)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~source) +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(title = "Distribution of drought sensitivities from two sources in 2010", x = "Species", y = "Sensitivity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

png("doc/display/tr_dendro_2010.png", width = 6, height = 4, units = "in", res = 300)
sens_trsp
dev.off()


# combine with tree metadata -------------

# read metadata file
meta <- readxl::read_excel("data/HKK-dendro/coremeta.xlsx", sheet = 2)

# canopy position for each tree is given by the column Dawk
# ref : https://github.com/forestgeo/HKK-tree-growth-data/blob/main/data/TreeRing/Raw/Vlam_HKKcores/metadata%20and%20mapping/2020_7_7_email%20regarding%20canopy%20positions.pdf

# ref for meaning of codes : https://github.com/forestgeo/HKK-tree-growth-data/blob/main/data/TreeRing/Raw/Vlam_HKKcores/metadata%20and%20mapping/Classification-of-crown-position-and-crown-shape-as-proposed-by-Dawkins-adapted-from.png

# size of the tree - DSH should be used and not DBH.
# DSH is the height at coring

head(rw.meta)

meta$treeID <- paste0("HKK.", meta$Tree)
head(meta$treeID)
head(rw.hkk.2010)

# join Dawk and DSH to the rw.hkk.2010 data
rw.hkk.2010 <- merge(rw.hkk.2010, select(meta, c("treeID", "Dawk", "DSH")), by = "treeID", all.x = T)

rw.hkk.2010$DSH <- as.numeric(rw.hkk.2010$DSH)

# remove NAs
rw.hkk.2010 <- rw.hkk.2010 %>%
    filter(!is.na(Dawk) & !is.na(DSH))

# distribution of Dawk and DSH across all individuals and species
dbh_cii_all <- ggplot(rw.hkk.2010, aes(
    x = DSH, y = factor(Dawk),
    fill = stat(x)
)) +
    ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    viridis::scale_fill_viridis(name = "growth", option = "C", direction = -1) +
    labs(title = "DSH in 2010 (cm) from tree rings", y = "canopy illumination") +
    # theme_ipsum() +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
    )

# make this by species
dbh_cii_sp <- ggplot(rw.hkk.2010, aes(
    x = DSH, y = factor(Dawk),
    fill = stat(x)
)) +
    ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    viridis::scale_fill_viridis(name = "growth", option = "C", direction = -1) +
    facet_wrap(~sp) +
    labs(title = "DSH in 2010 (cm) from tree rings", y = "canopy illumination") +
    # theme_ipsum() +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
    )

library(patchwork)
png("doc/display/explore/tree_rings_dbh_cii.png", width = 8, height = 4, units = "in", res = 300)
dbh_cii_all + dbh_cii_sp
dev.off()

# plot sensitivity distibution by canopy illumination for all and by species

sens_tr_cii <- ggplot(data = rw.hkk.2010, aes(x = Dawk, y = sens.prop)) +
    geom_violin(fill = "grey40") +
    geom_boxplot(width = 0.1, fill = "white") +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(title = "Distribution of drought sensitivities \nfrom tree rings in 2010", x = "Crown illumination index", y = "Sensitivity") +
    theme_bw()

sens_tr_ciisp <- ggplot(data = rw.hkk.2010, aes(x = Dawk, y = sens.prop, fill = sp)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") +
    facet_wrap(~sp) +
    # ylim(-5, 5)+
    # make color scale viridis
    scale_fill_viridis_d() +
    # make y axis percent
    scale_y_continuous(labels = scales::percent) +
    geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
    labs(x = "Crown illumination index", y = "Sensitivity") +
    theme_bw()

# save plots

png("doc/display/explore/tree_ring_cii_sens.png", width = 8, height = 4, units = "in", res = 300)
sens_tr_cii + sens_tr_ciisp
dev.off()
