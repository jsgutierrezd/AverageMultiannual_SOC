#===============================================================
# Proc06 Models ensemble - Classification
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
          'Boruta',
          'caretEnsemble'
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

data$wetland <- as.factor(data$wetland)
data$geology <- as.factor(data$geology)
data$georeg <- as.factor(data$georeg)
data$SlopeClass <- as.factor(data$SlopeClass)

data$AUCupdClass <- ifelse(data$AUCupdClass=="Gain","X1","X0")
data$AUCupdClass <- as.factor(data$AUCupdClass)
str(data)



# 4) AUC class ----------------------------------------------------------

# 4.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(14,31,59,60,64)]
names(data1)

names(getModelInfo())
set.seed(107)
inTrain <- createDataPartition(y = data1$AUCupdClass, p = .75, list = FALSE)
training <- data1[ inTrain,]
testing <- data1[-inTrain,]
my_control <- trainControl(
  method="boot",
  number=25,
  savePredictions="final",
  classProbs=TRUE,
  index=createResample(training$AUCupdClass, 25),
  summaryFunction=twoClassSummary
)

model_list <- caretList(
  AUCupdClass~., data=training,
  trControl=my_control,
  methodList=c("ranger", "nnet","svmLinear")
)

p <- as.data.frame(predict(model_list, newdata=testing))
print(p)

xyplot(resamples(model_list))
modelCor(resamples(model_list))


greedy_ensemble <- caretEnsemble(
  model_list, 
  metric="ROC",
  trControl=trainControl(
    number=5,
    summaryFunction=twoClassSummary,
    classProbs=TRUE
  ))
summary(greedy_ensemble)

library("caTools")
model_preds <- lapply(model_list, predict, newdata=testing, type="prob")
model_preds <- lapply(model_preds, function(x) x[,"X1"])
model_preds <- data.frame(model_preds)
(ens_preds <- predict(greedy_ensemble, newdata=testing))##, type="prob"))
complete_preds <- data.frame(model_preds,ens_preds,testing=testing$AUCupdClass)
complete_preds$ranger <- as.factor(ifelse(complete_preds$ranger>=0.5,"X1","X0"))
complete_preds$svmLinear <- as.factor(ifelse(complete_preds$svmLinear>=0.5,"X1","X0"))
complete_preds$nnet  <- as.factor(ifelse(complete_preds$nnet>=0.5,"X1","X0"))
complete_preds


confusionMatrix(complete_preds$ranger ,complete_preds$testing)
confusionMatrix(complete_preds$svmLinear ,complete_preds$testing)
confusionMatrix(complete_preds$ens_preds ,complete_preds$testing)
confusionMatrix(complete_preds$nnet ,complete_preds$testing)
model_preds$ensemble <- ens_preds
caTools::colAUC(model_preds, testing$AUCupdClass)

varImp(greedy_ensemble)


cov <- stack("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/Cov8619.tif")
names(cov) <- readRDS("~/AARHUS_PhD/DSMactivities/MonitoringGridData/SOCModeling/EnvironmLayers_30m/NamesCov8619.rds")

start <- Sys.time()
no_cores <- detectCores() - 2
beginCluster(no_cores)
#sp_prob <- raster::stack()
#for (i in 1:2) {
#  i=1
sp_prob <- clusterR(cov[[names(training[,-1])]], 
                    fun = predict,  
                    args=list(model=greedy_ensemble))

sp_prob <- clusterR(cov[[names(training[,-1])]], 
                    fun = predict,  
                    args=list(model=greedy_ensemble))

summary(greedy_ensemble)

plot(sp_prob)
writeRaster(sp_prob,"ProbMapGainEnsemble.tif")

#}