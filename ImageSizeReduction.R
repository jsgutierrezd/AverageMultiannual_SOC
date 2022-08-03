# 1984-1986 ---------------------------------------------------------------

setwd("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/Multiannual86_19")
#install.packages("terra")
library(terra)
library(parallel)

mask <- vect("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/CoastLine/Kyst.shp")

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/Multiannual1986_2019/Multiannual86_19/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)
all <- rast(all)
names(all)

all <- round(all, digits=2)

extended <- extend(all
                   , y = mask)

cropped <- crop(extended, mask,mask=T)

# rst <- cropped*mask

names <- paste0(getwd(),"/",names(all),".tif")

terra::writeRaster(cropped, filename=names,
                   datatype='FLT4S',overwrite=T)
