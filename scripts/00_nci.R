#script to calculate NCI from the last census before dendrobands were installed

# Load libraries
library(spatstat)


library(spatstat)
library(dplyr)
library(sf)
library(tidyverse)

#load data
load("data/HKK-dendro/stem4.RData")
load("data/HKK-dendro/DendroByBand.RData")

hkk_fordist<-hkk.stem4 %>%
filter(tag==StemTag) %>%
filter(status=="A")%>%
filter(!is.na(gx) & !is.na(gy))


#function to calculate NCI for each tree-------------
calculate_nci_tree <- function(df, tag, max_radius = 15, alpha = 1, beta = 2) {
  # Convert dataframe to sf object
  df_sf <- st_as_sf(df, coords = c("gx", "gy"))

  # Calculate NCI for the tree
    focal_tree <- df_sf[df_sf$tag == tag, ]

    # Find neighbors within max_radius
    neighbors <- st_is_within_distance(focal_tree, df_sf, dist = max_radius)

    # Check if neighbors list is empty
    if (length(neighbors) == 0 || length(neighbors[[1]]) == 0) {
      return(0)  # No neighbors found
    }

    neighbors <- neighbors[[1]]
    neighbors <- neighbors[neighbors != which(df_sf$tag==tag)]  # Remove self from neighbors

    if (length(neighbors) == 0) return(0)

    neighbor_trees <- df_sf[neighbors, ]
    distances <- st_distance(focal_tree, neighbor_trees)
    distances <- as.numeric(distances)

    nci<-sum((df$dbh[neighbors]^alpha / distances^beta))

  return(nci)
}

#apply this to the df-----------------------
nci <- mapply(function(x)calculate_nci_tree(hkk_fordist, tag=trees$Tag[x]), x=1:nrow(trees))

nci_df<-data.frame(tag = trees$Tag, nci_15 = nci)
write.csv(nci_df, "data/HKK-dendro/nci.csv") #write csv
