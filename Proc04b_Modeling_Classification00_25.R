#===============================================================
# Proc05 Modeling - Classification
#===============================================================
rm(list = ls())


# 1) Working directory ----------------------------------------------------

setwd("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/AverageMultiannual_SOC")



# 2) Libraries ------------------------------------------------------------

pckg <- c('caret',     
          'magrittr',
          'tidyr',
          'readr',
          'dplyr',
          'MASS',
          'parallel',
          'doParallel',
          'ranger',
          'nnet',
          'terra',
          'raster',
          'snow',
          'data.table',
          'Boruta'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)

# 3) Loading data ---------------------------------------------------------
data <- read_delim("RegMat00_25cm.csv",
                   delim=",") %>% data.frame

# data$wetland <- as.factor(data$wetland)
# data$geology <- as.factor(data$geology)
# data$georeg <- as.factor(data$georeg)
data$SlopeClass <- as.factor(data$SlopeClass)
data$AUCupdClass <- as.factor(data$AUCupdClass)
# str(data)


# 4) Slope class ----------------------------------------------------------

# 4.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(9,19:100)]
names(data1)

set.seed(17)
inTrain <- createDataPartition(y = data1$SlopeClass,
                               times = 100,
                               p = .75, 
                               list = T)
train_data_list <- list()
for (i in 1:length(inTrain)) {
  train_data_list[[i]] <- data1[ inTrain[[i]],]
}

test_data_list <- list()
for (i in 1:length(inTrain)) {
  test_data_list[[i]] <- data1[ -inTrain[[i]],]
}
#train_data <- data1[ inTrain[2],]
#test_data <- data1[-inTrain,]

# 4.2) Recursive feature elimination --------------------------------------
names(data1)
{
  start <- Sys.time()
  cl <- parallel::makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5,allowParallel = TRUE)
  rfe_list <- list()
  for (i in 1:length(inTrain)) {
    set.seed(i)
    rfe_list[[i]] <- rfe(x=train_data_list[[i]][,c(2:83)], 
                         y=train_data_list[[i]][,1], 
                         sizes=c(1:10), 
                         rfeControl=control2) %>% predictors()
    
  }
  print(Sys.time() - start)
}

lapply(rfe_list,length)
#rfe_list[[1]]
#plot(rfe, type=c("g", "o"))
#predictors(rfe)


# 4.3) Boruta -------------------------------------------------------------

names(data1)
#i=1
start <- Sys.time()
bor_list <- list()
for (i in 1:100) {
  bor <- Boruta(SlopeClass ~ ., data = train_data_list[[i]], doTrace = 0, ntree = 500,maxRuns=500)
  bor_list[[i]] <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])
  }
print(Sys.time() - start)


bor_list
# rfe_list[[1]]
# rfe_list[[2]]
# bor <- Boruta(SlopeClass ~ ., data = train_data_list[[i]], doTrace = 0, ntree = 500,maxRuns=500)
# plot(bor, xlab = "", xaxt = "n")
# lz<-lapply(1:ncol(bor$ImpHistory),function(i)
#   bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
# names(lz) <- colnames(bor$ImpHistory)
# Labels <- sort(sapply(lz,median))
# axis(side = 1,las=2,labels = names(Labels),
#      at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
# print(Sys.time() - start)
# 
# print(bor)
# names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])
# 
# final.bor <- TentativeRoughFix(bor)
# print(final.bor)
# # 
# # getSelectedAttributes(bor, withTentative = F)
# # boruta.df <- attStats(final.bor)


# 4.4) Models -------------------------------------------------------------


# 4.4.1) RFE models -------------------------------------------------------

fm_list <- list()
for (i in 1:length(rfe_list)) {
  fm_list[[i]] <- as.formula(paste0("SlopeClass~",paste0(as.character(rfe_list[[i]]),collapse = "+")))
}



# a) Random Forest ----------------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       returnResamp = "all",
                       search = "grid",
                       classProbs = TRUE
)

grid <- expand.grid(mtry = c(1,2,3,4,5),
                    splitrule = c("gini", "extratrees"),
                    min.node.size = c(1,2,3,4,5)
)

# set.seed(37)
model_rf_list <- list()
start <- Sys.time()
for (i in 1:length(train_data_list)) {
  set.seed(i)
  model_rf_list[[i]] <- train(fm_list[[i]],
                              data=train_data_list[[i]],
                              method = "ranger",
                              trControl = rctrlG,
                              tuneGrid = grid,
                              num.trees = 500,
                              importance = "impurity"
  ) 
}

print(Sys.time() - start)
var_imp_list <- list()
for (i in 1:length(model_rf_list)) {
  var_imp_list[[i]] <- model_rf_list[[i]][[11]][[6]]
}

pred_rf_list <- list()
for (i in 1:length(model_rf_list)) {
  pred_rf_list[[i]] <- predict(model_rf_list[[i]], newdata = test_data_list[[i]])
}

confMat_dataframe <- data.frame()
for (i in 1:length(pred_rf_list)) {
  tmp <- c(confusionMatrix(pred_rf_list[[i]] ,test_data_list[[i]][,1])$overall)
  confMat_dataframe <- rbind(confMat_dataframe,tmp)
  names(confMat_dataframe) <- c("Accuracy","Kappa","AccuracyLower",
                                "AccuracyUpper","AccuracyNull","AccuracyPValue",
                                "McNemarPValue")
  
}
confMat_dataframe

boxplot(confMat_dataframe$Kappa)
boxplot(confMat_dataframe$Accuracy)
summary(confMat_dataframe$Kappa)


library(ROCR)
data(ROCR.simple)
pred <- prediction(ROCR.simple$predictions,ROCR.simple$labels)
pred


x11()
hist(confMat_dataframe$Kappa,main="Histogram Kappa coefficient (n=100)",xlab="Kappa coefficient")
x11()
hist(confMat_dataframe$Accuracy,main="Histogram Accuracy coefficient (n=100)",xlab="Accuracy coefficient")

# 4.5) SlopeClass Model generalization ------------------------------------

cov <- stack("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif")
names(cov) <- readRDS("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")

# sp_pred <- predict(cov[[rfe_list[[1]]]],
#                 model_rf_list[[1]])
# ,
#                 filename = "SALIDAS/LAYERS/RANGER_ORDEN_03082021.tif",
#                 format = "GTiff", overwrite = T)

start <- Sys.time()
no_cores <- detectCores() - 2
beginCluster(no_cores)
sp_prob <- raster::stack()
for (i in 1:2) {
  i=1
  sp_prob <- clusterR(cov, 
                           fun = predict,  
                           args=list(model=model_rf_list[[95]], 
                                     type="prob"))
  
}
writeRaster(sp_prob,"ProbMapGains.tif")
plot(sp_prob)
endCluster()

print(Sys.time() - start)


# 5) AUC class ----------------------------------------------------------

# 5.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(14,19:100)]
names(data1)

set.seed(60)
inTrain <- createDataPartition(y = data1$AUCupdClass,
                               times = 100,
                               p = .75, 
                               list = T)
train_data_list <- list()
for (i in 1:length(inTrain)) {
  train_data_list[[i]] <- data1[ inTrain[[i]],]
}

test_data_list <- list()
for (i in 1:length(inTrain)) {
  test_data_list[[i]] <- data1[ -inTrain[[i]],]
}
#train_data <- data1[ inTrain[2],]
#test_data <- data1[-inTrain,]

# 5.2) Recursive feature elimination --------------------------------------
names(data1)
{
  start <- Sys.time()
  cl <- parallel::makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5,allowParallel = TRUE)
  rfe_list <- list()
  for (i in 1:length(inTrain)) {
    set.seed(i)
    rfe_list[[i]] <- rfe(x=train_data_list[[i]][,c(2:83)], 
                         y=train_data_list[[i]][,1], 
                         sizes=c(1:10), 
                         rfeControl=control2) %>% predictors()
    
  }
  print(Sys.time() - start)
}

lapply(rfe_list,length)
#rfe_list[[1]]
#plot(rfe, type=c("g", "o"))
#predictors(rfe)


# 5.3) Boruta -------------------------------------------------------------

names(data1)
#i=1
start <- Sys.time()
bor_list <- list()
for (i in 1:100) {
  bor <- Boruta(AUCupdClass ~ ., data = train_data_list[[i]], doTrace = 0, ntree = 500,maxRuns=500)
  bor_list[[i]] <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])
}
print(Sys.time() - start)


bor_list
# rfe_list[[1]]
# rfe_list[[2]]
# bor <- Boruta(SlopeClass ~ ., data = train_data_list[[i]], doTrace = 0, ntree = 500,maxRuns=500)
# plot(bor, xlab = "", xaxt = "n")
# lz<-lapply(1:ncol(bor$ImpHistory),function(i)
#   bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
# names(lz) <- colnames(bor$ImpHistory)
# Labels <- sort(sapply(lz,median))
# axis(side = 1,las=2,labels = names(Labels),
#      at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
# print(Sys.time() - start)
# 
# print(bor)
# names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])
# 
# final.bor <- TentativeRoughFix(bor)
# print(final.bor)
# # 
# # getSelectedAttributes(bor, withTentative = F)
# # boruta.df <- attStats(final.bor)


# 5.4) Models -------------------------------------------------------------


# 5.4.1) RFE models -------------------------------------------------------

fm_list <- list()
for (i in 1:length(rfe_list)) {
  fm_list[[i]] <- as.formula(paste0("AUCupdClass~",paste0(as.character(rfe_list[[i]]),collapse = "+")))
}



# a) Random Forest ----------------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       returnResamp = "all",
                       search = "grid",
                       classProbs = TRUE
)

grid <- expand.grid(mtry = c(1,2,3,4,5),
                    splitrule = c("gini", "extratrees"),
                    min.node.size = c(1,2,3,4,5)
)

# set.seed(37)
model_rf_list <- list()
start <- Sys.time()
for (i in 1:length(train_data_list)) {
  set.seed(i)
  model_rf_list[[i]] <- train(fm_list[[i]],
                              data=train_data_list[[i]],
                              method = "ranger",
                              trControl = rctrlG,
                              tuneGrid = grid,
                              num.trees = 500,
                              importance = "impurity"
  ) 
}

print(Sys.time() - start)
var_imp_list <- list()
for (i in 1:length(model_rf_list)) {
  var_imp_list[[i]] <- model_rf_list[[i]][[11]][[6]]
}

var_imp_list[1:10]

pred_rf_list <- list()
for (i in 1:length(model_rf_list)) {
  pred_rf_list[[i]] <- predict(model_rf_list[[i]], newdata = test_data_list[[i]])
}

confMat_dataframe <- data.frame()
for (i in 1:length(pred_rf_list)) {
  tmp <- c(confusionMatrix(pred_rf_list[[i]] ,test_data_list[[i]][,1])$overall)
  confMat_dataframe <- rbind(confMat_dataframe,tmp)
  names(confMat_dataframe) <- c("Accuracy","Kappa","AccuracyLower",
                                "AccuracyUpper","AccuracyNull","AccuracyPValue",
                                "McNemarPValue")
  
}
confMat_dataframe

boxplot(confMat_dataframe$Kappa)
boxplot(confMat_dataframe$Accuracy)
summary(confMat_dataframe$Kappa)



x11()
hist(confMat_dataframe$Kappa,main="Histogram Kappa coefficient (n=100)",xlab="Kappa coefficient")
x11()
hist(confMat_dataframe$Accuracy,main="Histogram Accuracy coefficient (n=100)",xlab="Accuracy coefficient")

# 5.5) SlopeClass Model generalization ------------------------------------

cov <- stack("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif")
names(cov) <- readRDS("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")

# sp_pred <- predict(cov[[rfe_list[[1]]]],
#                 model_rf_list[[1]])
# ,
#                 filename = "SALIDAS/LAYERS/RANGER_ORDEN_03082021.tif",
#                 format = "GTiff", overwrite = T)

start <- Sys.time()
no_cores <- detectCores() - 2
beginCluster(no_cores)
#sp_prob <- raster::stack()
#for (i in 1:2) {
#  i=1
  sp_prob <- clusterR(cov, 
                      fun = predict,  
                      args=list(model=model_rf_list[[1]], 
                                type="prob"))
  
#}
writeRaster(sp_prob,"ProbMapGains.tif",overwrite=T)
plot(sp_prob)
endCluster()

print(Sys.time() - start)
