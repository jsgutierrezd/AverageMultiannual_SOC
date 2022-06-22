#===============================================================
# Proc02 Data preparation - Point data set
#===============================================================
rm(list = ls())


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AverageMultiannual_SOC")



# 2) Libraries ------------------------------------------------------------


pckg <- c('DescTools',     
          'magrittr',
          'tidyr',
          'readr',
          'dplyr',
          'trend'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)




# 3) Point data set preparation -------------------------------------------


# 3.1) Data loading -------------------------------------------------------

data <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/DataJoined.csv",
                   delim = ";",
                   col_types = "cfdddcddddfcdddddfd" ) %>% 
  select(PointID,ID,X,Y,Depth,Year,OrgC,BDwhole) %>% 
  dplyr::mutate(
    PointID = as.factor(formatC(PointID,flag=0,width=3))
  )
head(data)

# 3.2) Data reshaping to wide format --------------------------------------

coord.oc.wide <- data %>%
  pivot_wider(names_from = Year, values_from = c(OrgC,BDwhole))


# 3.4) Calculating the number of OC values per point ----------------------
names(coord.oc.wide)
coord.oc.wide <- coord.oc.wide %>%  
  rowwise %>%  
  mutate(Values_per_row = sum(!is.na(c_across(OrgC_1986:OrgC_2019))))


# 3.5) Filtering the data set and keeping points with four OC values and depth interval 0-25 cm --------

coord.oc1.wide <- coord.oc.wide %>%
  filter(Depth==1) %>% filter(Values_per_row==4) %>% 
  select(-c(BDwhole_1986,BDwhole_1997,BDwhole_2009))


# 3.5.1) SOC stock calculation with bulk density and thickness ------------

coord.oc1.wide <- coord.oc1.wide %>% 
  mutate(SOC_1986 = OrgC_1986*BDwhole_2019*25,
         SOC_1997 = OrgC_1997*BDwhole_2019*25,
         SOC_2009 = OrgC_2009*BDwhole_2019*25,
         SOC_2019 = OrgC_2019*BDwhole_2019*25)


# 3.5.2) OC temporal trend calculation --------------------------------------
coord.oc1.wide <- coord.oc1.wide %>% na.omit


devtools::install_github("nxskok/mkac")
library(mkac)

names(coord.oc1.wide)
slopes <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- theil_sen_slope(y = c(unlist(coord.oc1.wide[i,11:14])), x = c(1986,1997,2009,2019)) 
  slopes <- c(slopes,tmp)
}

pval <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- mk.test(c(unlist(coord.oc1.wide[i,11:14])))$p.value 
  pval <- c(pval,tmp)
}

coord.oc1.wide <- coord.oc1.wide %>% bind_cols(Slope=slopes,Trendpvalue=pval)

coord.oc1.wide <- coord.oc1.wide %>% mutate(SlopeClass=ifelse(Slope>=0,"Gain","Loss"))

# 3.5.3) Area under the curve calculation -----------------------------------

auc <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- AUC(x = c(0,11,23,33),y = c(unlist(coord.oc1.wide[i,11:14])),method = "trapezoid") 
  auc <- c(auc,tmp)
}
coord.oc1.wide <- coord.oc1.wide %>% bind_cols(AUC=auc)

#Calculate the relative accumulation rate = ((AUC-BaseArea)/BaseArea)*100
coord.oc1.wide <- coord.oc1.wide %>% mutate(BaseArea= SOC_1986*(2019-1986),
                                            AUCupdated=AUC-BaseArea,
                                            RAR=((AUC-BaseArea)/BaseArea)*100)

coord.oc1.wide <- coord.oc1.wide %>% mutate(AUCupdClass=ifelse(AUCupdated>=0,"Gain","Loss"))


# 3.5.4) SOC average and range --------------------------------------------

names(coord.oc1.wide)

mean <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- mean(c(unlist(coord.oc1.wide[i,11:14]))) 
  mean <- c(mean,tmp)
}

std <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- sd(c(unlist(coord.oc1.wide[i,11:14]))) 
  std <- c(std,tmp)
}

diff <- c()
for (i in 1:dim(coord.oc1.wide)[1]) {
  tmp <- max(c(unlist(coord.oc1.wide[i,11:14]))) -  min(c(unlist(coord.oc1.wide[i,11:14]))) 
  diff <- c(diff,tmp)
}


coord.oc1.wide <- coord.oc1.wide %>% 
  bind_cols(SOCAverage=mean,SOCsd=std,SOCdiff=diff) %>% 
  mutate(SOCcv=(SOCsd/SOCAverage)*100)



# 3.5.4) Final data set exporting -------------------------------------------

write_csv(coord.oc1.wide,"C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC00_25cm.csv")

# 3.6) Filtering the data set and keeping points with four OC values Depth interval 25-50 cm --------

coord.oc2.wide <- coord.oc.wide %>%
  filter(Depth==2) %>% filter(Values_per_row==4)

# 3.6.1) SOC stock calculation with bilk density and thickness ------------

bd2 <- bd %>% dplyr::select(Point,Depth,BDwhole) %>%  filter(Depth==2) %>% dplyr::select(-Depth)
bd2

coord.oc2.wide <- coord.oc2.wide %>% inner_join(bd2,by = "Point")
coord.oc2.wide


# 3.6.2) OC temporal trend calculation --------------------------------------

coord.oc2.wide <- coord.oc2.wide %>% 
  mutate(SOC1986 = Year1986*BDwhole*25,
         SOC1997 = Year1997*BDwhole*25,
         SOC2009 = Year2009*BDwhole*25,
         SOC2019 = Year2019*BDwhole*25)

coord.oc2.wide <- coord.oc2.wide %>% na.omit
#devtools::install_github("nxskok/mkac")
library(mkac)

slopes <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- theil_sen_slope(y = c(unlist(coord.oc2.wide[i,12:15])), x = c(1986,1997,2009,2019)) 
  slopes <- c(slopes,tmp)
}

pval <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- mk.test(c(unlist(coord.oc2.wide[i,12:15])))$p.value 
  pval <- c(pval,tmp)
}

coord.oc2.wide <- coord.oc2.wide %>% bind_cols(Slope=slopes,Trendpvalue=pval)

coord.oc2.wide <- coord.oc2.wide %>% mutate(SlopeClass=ifelse(Slope>=0,"Gain","Loss"))


# 3.6.3) Area under the curve calculation -----------------------------------

auc <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- AUC(x = c(0,11,23,33),y = c(unlist(coord.oc2.wide[i,12:15])),method = "trapezoid") 
  auc <- c(auc,tmp)
}
coord.oc2.wide <- coord.oc2.wide %>% bind_cols(AUC=auc)

#Calculate the relative accumulation rate = ((AUC-BaseArea)/BaseArea)*100
coord.oc2.wide <- coord.oc2.wide %>% mutate(BaseArea= SOC1986*(2019-1986),
                                            AUCupdated=AUC-BaseArea,
                                            RAR=((AUC-BaseArea)/BaseArea)*100)
coord.oc2.wide <- coord.oc2.wide %>% mutate(AUCupdClass=ifelse(AUCupdated>=0,"Gain","Loss"))

# 3.6.4) SOC average and range --------------------------------------------

names(coord.oc2.wide)

mean <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- mean(c(unlist(coord.oc2.wide[i,12:15]))) 
  mean <- c(mean,tmp)
}

std <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- sd(c(unlist(coord.oc2.wide[i,12:15]))) 
  std <- c(std,tmp)
}

diff <- c()
for (i in 1:dim(coord.oc2.wide)[1]) {
  tmp <- max(c(unlist(coord.oc2.wide[i,12:15]))) -  min(c(unlist(coord.oc2.wide[i,12:15]))) 
  diff <- c(diff,tmp)
}


coord.oc2.wide <- coord.oc2.wide %>% 
  bind_cols(SOCAverage=mean,SOCsd=std,SOCdiff=diff) %>% 
  mutate(SOCcv=(SOCsd/SOCAverage)*100)


# 3.6.4) Final data set exporting -------------------------------------------

write_csv(coord.oc2.wide,"C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC25_50cm.csv")
