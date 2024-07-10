
# load libraries---------------------------
library(tidyverse)
library(lubridate)
library(httr)
library(janitor)

#reading the GitHub PAT (stored privately and locally)
#pat<-as.character(read.table("../data/HKK-dendro/git_pat.txt")[1,1])

# remotes:::download("DendroByBand.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-dendrobands/main/data/modified_data/DendroByBand.RData",
#                    auth_token = pat)

load("growth-precip-thailand/data/HKK-dendro/DendroByBand.RData")


dendrobyband<-merge(dendrobyband, trees, by="Tag")

hkk_OBSID$year<-year(hkk_OBSID$Date2)
hkk_OBSID$month<-month(hkk_OBSID$Date2)
hkk_OBSID$day<-day(hkk_OBSID$Date2)

#correcting two sets of wrong dates

hkk_OBSID$year[hkk_OBSID$year==2209]<-2009
hkk_OBSID$year[hkk_OBSID$month==1 & hkk_OBSID$Cno==5]<-2011

dendrobyband<-merge(dendrobyband, hkk_OBSID, by=c("OBSID", "Cno", "Date2", "Tag"))

dendrobyband$Date2<-as.Date(paste0(dendrobyband$year,"-", dendrobyband$month, "-", dendrobyband$day))


# calculate dbh by stem-------------------------------

#remove all data with decimal errors, misplaced decimal

dendro.data<-dendrobyband %>%
  dplyr::rename(sp=SPCODE.UPDATE)%>%
  mutate(year=year(Date2)) %>%
  filter(is.na(DBH1_flag), DBH1_error=="none",
         win.decimal.error=="none", is.na(win.flag))%>%
  group_by(Cno, sp, Tag) %>%
  dplyr::summarise(calcDBH=calcDBH[which(n.dendro==min(n.dendro))],
                   Date2=Date2[which(n.dendro==min(n.dendro))],
                   n.dendro=min(n.dendro), #which dendroband is it from?
                   Cii=max(Cii, na.rm=T))

dbh.data<-dendrobyband %>%
  dplyr::rename(sp=SPCODE.UPDATE)%>%
  mutate(year=year(Date2)) %>%
  filter(is.na(dbh.flag), dbh.decimal.error=="none")%>%
  group_by(Cno, sp, Tag) %>%
  dplyr::summarise(dbh=DbhDendro[which(n.dendro==min(n.dendro))],
                   Date2=Date2[which(n.dendro==min(n.dendro))],
                   n.dendro=min(n.dendro),
                   Cii=max(Cii, na.rm=T))

#calculate long-term growth rates-----------------------------

# remotes:::download("../data/HKK-dendro/stem1.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem1.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem2.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem2.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem3.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem3.rdata",
#                    auth_token = pat)
# remotes:::download("../data/HKK-dendro/stem4.RData",
#                    "https://raw.githubusercontent.com/forestgeo/HKK-tree-growth-data/main/data/Census/raw_data/hkk.stem4.rdata",
#                    auth_token = pat)

load("growth-precip-thailand/data/HKK-dendro/stem1.RData")
load("growth-precip-thailand/data/HKK-dendro/stem2.RData")
load("growth-precip-thailand/data/HKK-dendro/stem3.RData")
load("growth-precip-thailand/data/HKK-dendro/stem4.RData")

hkk_census<-bind_rows(hkk.stem1, hkk.stem2, hkk.stem3, hkk.stem4)

hkk_census <- hkk_census %>%
  filter(status=="A")%>%
  filter(tag %in% dendro.data$Tag)%>%
  group_by(tag, StemTag)%>%
  dplyr::mutate(dbhmin1 = dplyr::lag(dbh, n=1, default=NA, order_by=StemTag),
                inc_annual = (dbh-dbhmin1)/5)%>%
  filter(!is.na(dbhmin1))%>%
  #average increment in cm
  dplyr::summarise(avg_inc = mean(inc_annual/10, na.rm=T),
                   sd_inc = sd(inc_annual/10, na.rm=T))
                  

#since dendrobands are installed on trees with tag==StemTag
hkk_census<-hkk_census %>%
filter(tag==StemTag)

#make the colnames compatible
hkk_census<-hkk_census %>%
rename(Tag=tag)  

#calculate increments-----------------------------

#add long-term increment to dendro data
dendro.data<-merge(dendro.data, hkk_census, by="Tag", all.x=T)
dbh.data<-merge(dbh.data, hkk_census, by="Tag", all.x=T)

#making increments for each year

alt_census<-c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)

dendro_inc<- dendro.data %>%
  filter(Cno %in% alt_census)%>%
  arrange(Cno)%>%
  group_by(Tag) %>% #grouping by treeID to calculate increment
  dplyr::mutate(calcDBH_min1=dplyr::lag(calcDBH, n=1, default=NA, order_by=Tag),
  n.dendro_diff=n.dendro-dplyr::lag(n.dendro, n=1, default=NA, order_by=Tag),
         cii_min1=dplyr::lag(Cii, n=1, default=NA, order_by=Tag),
         tdif=difftime(Date2, dplyr::lag(Date2, n=1, default=NA, order_by=Tag), units = "days"),
         inc=calcDBH-calcDBH_min1,
         inc_annual=(inc/as.numeric(tdif))*365,
         inc.time=Cno-dplyr::lag(Cno, n=1, default=NA, order_by=Tag),
         #removing increment value if the previous measure was before the last year
         inc_annual=ifelse(inc.time==2, inc_annual, NA),
         treeID=paste0("HKK.", Tag))%>%
  filter(!is.na(inc))%>%
  filter(!is.na(inc_annual))%>%
  ungroup()%>%
  dplyr::mutate(
    diainc.scaled = scale(inc_annual),
    out.inc = ifelse(abs(diainc.scaled)>3, 1, 0))%>%
  filter(out.inc==0)


dbh_inc<- dbh.data %>%
  filter(Cno %in% alt_census)%>%
  arrange(Cno)%>%
  group_by(Tag) %>% #grouping by treeID to calculate increment
  dplyr::mutate(dbh_min1=dplyr::lag(dbh, n=1, default=NA, order_by=Tag),
  n.dendro_diff=n.dendro-dplyr::lag(n.dendro, n=1, default=NA, order_by=Tag),
         cii_min1=dplyr::lag(Cii, n=1, default=NA, order_by=Tag),
         tdif=difftime(Date2, dplyr::lag(Date2, n=1, default=NA, order_by=Tag), units = "days"),
         inc=dbh-dbh_min1,
         inc_annual=(inc/as.numeric(tdif))*365,
         inc.time=Cno-dplyr::lag(Cno, n=1, default=NA, order_by=Tag),
         #removing increment value if the previous measure was before the last year
         inc_annual=ifelse(inc.time==2, inc_annual, NA),
         treeID=paste0("HKK.", Tag))%>%
  filter(!is.na(inc))%>%
  filter(!is.na(inc_annual))%>%
  ungroup()%>%
  dplyr::mutate(
    diainc.scaled = scale(inc_annual),
    out.inc = ifelse(abs(diainc.scaled)>3, 1, 0))%>%
  filter(out.inc==0)

#identiy outliers------------------------------

#merge dbh and dendro data by Tag and Cno
inc_all<-merge(dendro_inc, dbh_inc, by=c("sp", "Tag", "Cno"), all=T)

#Perpendicular distance from point 'a' to a line with 'slope' and 'intercept'
dist_point_line <- function(a, slope, intercept) {
    b = c(1, intercept+slope)
    c = c(-intercept/slope,0)       
    v1 <- b - c
    v2 <- a - b
    m <- cbind(v1,v2)
    return(abs(det(m))/sqrt(sum(v1*v1)))
}

inc_all<-inc_all %>%
  dplyr::rename(inc_annual.dendro=inc_annual.x, 
  inc_annual.dbh=inc_annual.y,
  band.diff.dendro=n.dendro_diff.x,
  band.diff.dbh=n.dendro_diff.y) %>%
  #split into three size classes
  dplyr::mutate(size_class= ifelse(dbh<35, "small",
  ifelse(dbh<70, "medium", "large")),
  size_class=factor(size_class, levels=c("small", "medium", "large"))
  #find the deviation of the increment from the 1:1 line
  )%>%
  select(sp, Tag, Cno, inc_annual.dendro, inc_annual.dbh, dbh, 
  size_class, band.diff.dbh, band.diff.dendro)

#using the distance to the line
dev.1<-mapply(function(x)dist_point_line(c(inc_all$inc_annual.dendro[x], inc_all$inc_annual.dbh[x]), 1, 0), x=1:nrow(inc_all))
inc_all$dev.1<-dev.1

#remove outliers-----------------------------------

#outlier dataframe
outliers<-inc_all %>% filter(dev.1>=0.5)
write.csv(outliers, "growth-precip-thailand/data/HKK-dendro/inc_outliers.csv", row.names=F)

#clean dataframe

dendro_inc<-merge(dendro_inc, inc_all %>% select(Tag, Cno, sp, dev.1), 
by=c("Tag", "Cno", "sp"), all.x=T)

dendro_inc_clean<-dendro_inc %>% filter(dev.1<0.5)
