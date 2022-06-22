#===============================================================
# Proc01 Data preparation - Environmental layers
#===============================================================
rm(list = ls())


# 1) Working directory ----------------------------------------------------


setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AverageMultiannual_SOC")



# 2) Libraries ------------------------------------------------------------


pckg <- c('terra',     
          'magrittr',
          'magrittr',
          'devtools',
          'raster',
          'parallel'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)



# 3) Environmental layers harmonization -----------------------------------

# This section aims to harmonize environmental layers geometry using the coast line shapefile as the reference layer

coast <- vect("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/CoastLine/Kyst.shp")
#data.frame(coast)
#plot(coast)

# 3.1) Bioclimatic layers -----------------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Climate/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

climate <- rast(all) %>% crop(coast,mask=T)


# 3.2) Climatic variables -------------------------------------------------

all <- list.files("O:/Tech_AGRO/Jord/Sebastian/BIOLAYERS_CHELSA1985_2021_30M/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all[7:39],function(x){
  stack(x)
})




beginCluster(detectCores()-5)
biostack <- all[[1]][[1]]
for (i in 2:length(all)) {
  temp <- all[[i]][[1]]
  biostack <- stack(biostack,temp)
  bio1 <- clusterR(biostack,
                  calc, 
                  args=list(fun=mean)
                  )
}
endCluster()
Sys.time()-start


biostack <- all[[1]][[2]]
for (i in 2:length(all)) {
  temp <- all[[i]][[2]]
  biostack <- c(biostack,temp)
  bio2 <- app(biostack,fun="mean",cores=10)
}

biostack <- all[[1]][[3]]
for (i in 2:length(all)) {
  temp <- all[[i]][[3]]
  biostack <- c(biostack,temp)
  bio3 <- app(biostack,fun="mean",cores=10)
}

biostack <- all[[1]][[4]]
for (i in 2:length(all)) {
  temp <- all[[i]][[4]]
  biostack <- c(biostack,temp)
  bio4 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[5]]
for (i in 2:length(all)) {
  temp <- all[[i]][[5]]
  biostack <- c(biostack,temp)
  bio5 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[6]]
for (i in 2:length(all)) {
  temp <- all[[i]][[6]]
  biostack <- c(biostack,temp)
  bio6 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[7]]
for (i in 2:length(all)) {
  temp <- all[[i]][[7]]
  biostack <- c(biostack,temp)
  bio7 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[8]]
for (i in 2:length(all)) {
  temp <- all[[i]][[8]]
  biostack <- c(biostack,temp)
  bio8 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[9]]
for (i in 2:length(all)) {
  temp <- all[[i]][[9]]
  biostack <- c(biostack,temp)
  bio9 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[10]]
for (i in 2:length(all)) {
  temp <- all[[i]][[10]]
  biostack <- c(biostack,temp)
  bio10 <- app(biostack,fun="mean",cores=10)
}



biostack <- all[[1]][[11]]
for (i in 2:length(all)) {
  temp <- all[[i]][[11]]
  biostack <- c(biostack,temp)
  bio11 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[12]]
for (i in 2:length(all)) {
  temp <- all[[i]][[12]]
  biostack <- c(biostack,temp)
  bio12 <- app(biostack,fun="mean",cores=10)
}



biostack <- all[[1]][[13]]
for (i in 2:length(all)) {
  temp <- all[[i]][[13]]
  biostack <- c(biostack,temp)
  bio13 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[14]]
for (i in 2:length(all)) {
  temp <- all[[i]][[14]]
  biostack <- c(biostack,temp)
  bio14 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[15]]
for (i in 2:length(all)) {
  temp <- all[[i]][[15]]
  biostack <- c(biostack,temp)
  bio15 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[16]]
for (i in 2:length(all)) {
  temp <- all[[i]][[16]]
  biostack <- c(biostack,temp)
  bio16 <- app(biostack,fun="mean",cores=10)
}



biostack <- all[[1]][[17]]
for (i in 2:length(all)) {
  temp <- all[[i]][[17]]
  biostack <- c(biostack,temp)
  bio17 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[18]]
for (i in 2:length(all)) {
  temp <- all[[i]][[18]]
  biostack <- c(biostack,temp)
  bio18 <- app(biostack,fun="mean",cores=10)
}


biostack <- all[[1]][[19]]
for (i in 2:length(all)) {
  temp <- all[[i]][[19]]
  biostack <- c(biostack,temp)
  bio19 <- app(biostack,fun="mean",cores=10)
}

climate <- rast(all) %>% crop(coast,mask=T)


# 3.3) Geomorphology layers -----------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Geomorphology/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

geomorphology <- rast(all) %>% crop(coast,mask=T)


# 3.4) Geology ------------------------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Geology/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

geology <- rast(all) %>% crop(coast,mask=T)


# 3.5) Soil layers --------------------------------------------------------


# 3.5.1) Clay layers for different depth intervals ------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Soil/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

soil <- rast(all) %>% crop(coast,mask=T)

# 3.6) Organimsms layers --------------------------------------------------

# 3.6.1) Existing covariates ----------------------------------------------

all <- list.files("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/",
                  pattern = "tiff$|tif$",
                  full.names = TRUE
)

all <- lapply(all,function(x){
  rast(x)
})

organisms <- rast(all) %>% crop(coast,mask=T)


# 3.6.2) New covariates ---------------------------------------------------

# Remote sensing data --> Landsat 7

#Function to calculate spectral indices

indicesL7 <- function(bandstack){
  ndvi <- (bandstack[[4]]-bandstack[[3]])/(bandstack[[4]]+bandstack[[3]])
  
  evi <- ((bandstack[[4]]-bandstack[[3]])/(bandstack[[4]] + 6*bandstack[[3]] - 7.5*bandstack[[1]] + 1))*2.5
  
  bsi <- ((bandstack[[6]] + bandstack[[3]]) - (bandstack[[3]]+bandstack[[1]]))/((bandstack[[6]] + bandstack[[3]]) + (bandstack[[3]]+bandstack[[1]]))
  
  savi <- (bandstack[[4]]-bandstack[[3]])/((bandstack[[4]]+bandstack[[3]] + 0.5)*1.5)
  
  str <- ((1-bandstack[[6]])^2)/(2*bandstack[[6]])
  
  brightness <- bandstack[[1]]*0.2043 + bandstack[[2]]*0.4158 + bandstack[[3]]*0.5524 + 
    bandstack[[4]]*0.5741 + bandstack[[5]]*0.3124 + bandstack[[6]]*0.2303
  
  greenness <- bandstack[[1]]*(-0.1603) + bandstack[[2]]*(-0.2819) + bandstack[[3]]*(-0.4934) + 
    bandstack[[4]]*0.7940 + bandstack[[5]]*0.0002 + bandstack[[6]]*(-0.1446) 
  
  wetness <- bandstack[[1]]*0.0315 + bandstack[[2]]*0.2021 + bandstack[[3]]*0.3102 + 
    bandstack[[4]]*0.1594 + bandstack[[5]]*0.6806 + bandstack[[6]]*(-0.6109)
  indices <- c(ndvi,evi,savi,bsi,str,brightness,greenness,wetness)
  return(indices)
}

# Load images
l7med <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/Landsat7Bands/L7Med86_19.tif") %>% 
  crop(coast,mask=T)
l7p10 <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/Landsat7Bands/L7P1086_19.tif") %>% 
  crop(coast,mask=T)
l7p90 <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Organisms/Landsat7Bands/L7P9086_19.tif") %>% 
  crop(coast,mask=T)


start <- Sys.time()
l7medInd <- app(l7med,   fun=indicesL7, cores =16)  %>% 
  resample(organisms)
names(l7medInd) <- c("NDVIMed","EVIMed","SAVIMed","BSIMed","STRMed","BRIGHTMed","GREENMed","WETMed")
l7p10Ind <- app(l7p10,   fun=indicesL7, cores =16)  %>% 
  resample(organisms)
names(l7p10Ind) <- c("NDVIP10","EVIP10","SAVIP10","BSIP10","STRP10","BRIGHTP10","GREENP10","WETP10") 
l7p90Ind <- app(l7p90,   fun=indicesL7, cores =16)  %>% 
  resample(organisms)
names(l7p90Ind) <- c("NDVIP90","EVIP90","SAVIP90","BSIP90","STRP90","BRIGHTP90","GREENP90","WETP90")
Sys.time()-start

organisms <- c(organisms,l7medInd,l7p10Ind,l7p90Ind)

# Land use time series 1990/2005-2019

# Categories
# 11. Settlement
# 12. Permanent wetland
# 13. Woodland
# 14. Trees crop
# 15. Cropland
# 16. Grassland
# 17. Periodic wetland
# 18. Unclassified
# 19. Other land uses
# 20. Water

files <- list.files(path = "O:/Tech_AGRO/Jord/Sebastian/LandcoverBaseMap1990-2018/LULC_1990_2020",
                    pattern = "tiff$|tif$",
                    full.names = TRUE
)

all <- rast(files)
all

LU <- resample(all,geology,method="near") %>% crop(coast,mask=T)

# i=1
LU_dummy <- list(0)
for (i in 1:nlyr(LU)) {
  LU_dummy[[i]] <- segregate(LU[[i]])
}

LU11 <- c(LU_dummy[[1]][[1]],
          LU_dummy[[2]][[1]],
          LU_dummy[[3]][[1]],
          LU_dummy[[4]][[1]],
          LU_dummy[[5]][[1]],
          LU_dummy[[6]][[1]],
          LU_dummy[[7]][[1]],
          LU_dummy[[8]][[1]],
          LU_dummy[[9]][[1]],
          LU_dummy[[10]][[1]],
          LU_dummy[[11]][[1]])

LU11 <- app(LU11,fun="sum") 

LU12 <- c(LU_dummy[[1]][[2]],
          LU_dummy[[2]][[2]],
          LU_dummy[[3]][[2]],
          LU_dummy[[4]][[2]],
          LU_dummy[[5]][[2]],
          LU_dummy[[6]][[2]],
          LU_dummy[[7]][[2]],
          LU_dummy[[8]][[2]],
          LU_dummy[[9]][[2]],
          LU_dummy[[10]][[2]],
          LU_dummy[[11]][[2]])
LU12 <- app(LU12,fun="sum")

LU13 <- c(LU_dummy[[1]][[3]],
          LU_dummy[[2]][[3]],
          LU_dummy[[3]][[3]],
          LU_dummy[[4]][[3]],
          LU_dummy[[5]][[3]],
          LU_dummy[[6]][[3]],
          LU_dummy[[7]][[3]],
          LU_dummy[[8]][[3]],
          LU_dummy[[9]][[3]],
          LU_dummy[[10]][[3]],
          LU_dummy[[11]][[3]])
LU13 <- app(LU13,fun="sum") 

LU14 <- c(LU_dummy[[1]][[4]],
          LU_dummy[[2]][[4]],
          LU_dummy[[3]][[4]],
          LU_dummy[[4]][[4]],
          LU_dummy[[5]][[4]],
          LU_dummy[[6]][[4]],
          LU_dummy[[7]][[4]],
          LU_dummy[[8]][[4]],
          LU_dummy[[9]][[4]],
          LU_dummy[[10]][[4]],
          LU_dummy[[11]][[4]])
LU14 <- app(LU14,fun="sum") 

LU15 <- c(LU_dummy[[1]][[5]],
          LU_dummy[[2]][[5]],
          LU_dummy[[3]][[5]],
          LU_dummy[[4]][[5]],
          LU_dummy[[5]][[5]],
          LU_dummy[[6]][[5]],
          LU_dummy[[7]][[5]],
          LU_dummy[[8]][[5]],
          LU_dummy[[9]][[5]],
          LU_dummy[[10]][[5]],
          LU_dummy[[11]][[5]])
LU15 <- app(LU15,fun="sum") 
plot(LU15)

LU16 <- c(LU_dummy[[1]][[6]],
          LU_dummy[[2]][[6]],
          LU_dummy[[3]][[6]],
          LU_dummy[[4]][[6]],
          LU_dummy[[5]][[6]],
          LU_dummy[[6]][[6]],
          LU_dummy[[7]][[6]],
          LU_dummy[[8]][[6]],
          LU_dummy[[9]][[6]],
          LU_dummy[[10]][[6]],
          LU_dummy[[11]][[6]])
LU16 <- app(LU16,fun="sum") 

LU17 <- c(LU_dummy[[1]][[7]],
          LU_dummy[[2]][[7]],
          LU_dummy[[3]][[7]],
          LU_dummy[[4]][[7]],
          LU_dummy[[5]][[7]],
          LU_dummy[[6]][[7]],
          LU_dummy[[7]][[7]],
          LU_dummy[[8]][[7]],
          LU_dummy[[9]][[7]],
          LU_dummy[[10]][[7]],
          LU_dummy[[11]][[7]])
LU17 <- app(LU17,fun="sum") 

LU18 <- c(LU_dummy[[1]][[8]],
          LU_dummy[[2]][[8]],
          LU_dummy[[3]][[8]],
          LU_dummy[[4]][[8]],
          LU_dummy[[5]][[8]],
          LU_dummy[[6]][[8]],
          LU_dummy[[7]][[8]],
          LU_dummy[[8]][[8]],
          LU_dummy[[9]][[8]],
          LU_dummy[[10]][[8]],
          LU_dummy[[11]][[8]])
LU18 <- app(LU18,fun="sum") 

LU19 <- c(LU_dummy[[1]][[9]],
          LU_dummy[[2]][[9]],
          LU_dummy[[3]][[9]],
          LU_dummy[[4]][[9]],
          LU_dummy[[5]][[9]],
          LU_dummy[[6]][[9]],
          LU_dummy[[7]][[9]],
          LU_dummy[[8]][[9]],
          LU_dummy[[9]][[9]],
          LU_dummy[[10]][[9]],
          LU_dummy[[11]][[9]])
LU19 <- app(LU19,fun="sum") 


LU20 <- c(LU_dummy[[1]][[10]],
          LU_dummy[[2]][[10]],
          LU_dummy[[3]][[10]],
          LU_dummy[[4]][[10]],
          LU_dummy[[5]][[10]],
          LU_dummy[[6]][[10]],
          LU_dummy[[7]][[10]],
          LU_dummy[[8]][[10]],
          LU_dummy[[9]][[10]],
          LU_dummy[[10]][[10]],
          LU_dummy[[11]][[10]])
LU20 <- app(LU20,fun="sum") 


LUsum <- LU12+LU13+LU14+LU15+LU16+LU17
plot(LUsum)
LU12p <- LU12/LUsum
LU13p <- LU13/LUsum 
LU14p <- LU14/LUsum  
LU15p <- LU15/LUsum 
LU16p <- LU16/LUsum  
LU17p <- LU17/LUsum 

LU12p[LU12p==0] <- NA
LU13p[LU13p==0] <- NA
LU14p[LU14p==0] <- NA
LU15p[LU15p==0] <- NA
LU16p[LU16p==0] <- NA
LU17p[LU17p==0] <- NA

LU12ln <- app(LU12p,fun="log")
LU13ln <- app(LU13p,fun="log")
LU14ln <- app(LU14p,fun="log")
LU15ln <- app(LU15p,fun="log")
LU16ln <- app(LU16p,fun="log")
LU17ln <- app(LU17p,fun="log")

LU12pln <- LU12p*LU12ln
LU13pln <- LU13p*LU13ln
LU14pln <- LU14p*LU14ln
LU15pln <- LU15p*LU15ln
LU16pln <- LU16p*LU16ln
LU17pln <- LU17p*LU17ln

Index <- c(abs(LU12pln),abs(LU13pln),abs(LU14pln),
           abs(LU15pln),abs(LU16pln),abs(LU17pln))

plot(Index)

DivIndex <- app(Index, fun="sum",na.rm=T)
DivIndex[DivIndex>1] <- NA
DivIndex <- round(DivIndex,1)

organisms <- c(organisms,DivIndex,LU12,LU13,LU14,LU15,LU16,LU17)


# covs <- rast("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif")
# names(covs) <- readRDS("C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")
# names(covs)
# 
# covs <- c(covs[[1:31]],DivIndex,LU12,LU13,LU14,LU15,LU16,LU17,covs[[32:75]])
# names(covs)[32:38] <- c("DivIndLU","YearsPermWetland","YearsWoodland",
#                        "YearsTreeCrop","YearsCropland","YearsGrassland",
#                        "YearsPerdWetland")
# names(covs)
# 3.7) Spatial information ------------------------------------------------
coord <- as.data.frame(organisms[[1]], xy=TRUE)

x <- data.frame(coord$x,coord$y,coord$x) 
x <- rast(x, type="xyz")
crs(x)  <- "epsg:25832"
x <- x %>% crop(coast,mask=T) %>% resample(organisms[[1]])
names(x) <- "CoordX"

y <- data.frame(coord$x,coord$y,coord$y)
y <- rast(y, type="xyz")   
crs(y)  <- "epsg:25832"
y <- y %>% crop(coast,mask=T) %>% resample(organisms[[1]])
names(y) <- "CoordY"

# 3.8) Whole set of environmental layers set ------------------------------
rm(cov)
covs <- c(x,y,organisms,soil,geology,geomorphology,climate)
names(covs)

# 4) Exporting environmental layers ---------------------------------------
writeRaster(cov,"C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif",overwrite=T)
names(cov)
saveRDS(names(cov),"C:/Users/au704633/OneDrive - Aarhus Universitet/Documents/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")












