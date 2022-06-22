#===============================================================
# Proc03 Data preparation - Regression Matrix
#===============================================================
rm(list = ls())


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AverageMultiannual_SOC")



# 2) Libraries ------------------------------------------------------------


pckg <- c('terra',     
          'magrittr',
          'tidyr',
          'readr',
          'dplyr'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)


# 3) Loading environmental layers -----------------------------------------

covs <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif")
names(covs) <- readRDS("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")


# 4) Loading point data set -----------------------------------------------


# 4.1) Point data set 0-25 cm ---------------------------------------------

data <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC00_25cm.csv",
                 delim=",")
data <- data %>% dplyr::select(X,Y,SOC_1986:SOCcv)
names(data)

data_sp <- vect(data, geom=c("X", "Y"),crs="epsg:25832")# %>%  project("epsg:25832")#Data frame as spatial points data frame
data_sp


# 4.1.2) Multi-point extraction 0-25 cm -----------------------------------

data <- cbind(data,terra::extract(covs,data_sp))
data$ID <- NULL
data <- data %>% na.omit

# 4.1.3) Exporting regression matrix 0-25 cm ------------------------------

write_csv(data,"RegMat00_25cm.csv")

# 4.2) Point data set 25-50 cm ---------------------------------------------

data <- read_delim("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/Points/PointsOC25_50cm.csv",
                   delim=",")
data <- data %>% dplyr::select(x,y,SOC1986:SOCcv)
names(data)

data_sp <- vect(data, geom=c("x", "y"),crs="epsg:25832")# %>%  project("epsg:25832")#Data frame as spatial points data frame
data_sp


# 4.2.2) Multi-point extraction 25-50 cm -----------------------------------

data <- cbind(data,terra::extract(covs,data_sp))
data$ID <- NULL

# 4.3.3) Exporting regression matrix 25-50 cm -----------------------------

write_csv(data,"RegMat25_50cm.csv")

