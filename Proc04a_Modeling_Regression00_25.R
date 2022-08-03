#===============================================================
# Proc04 Modeling - Regression
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
          'e1071',
          'ranger',
          'hydroGOF',
          'Boruta',
          'quantregForest',
          'prospectr'
)

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
lapply(pckg,usePackage)

goof <- function(observed,predicted, plot.it = FALSE, type = "DSM"){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  # Standard error of prediction ˆ2 or MSE
  SEP2 <- mean((observed - predicted)^2)
  # Standard error of prediction or RMSE
  SEP <- sqrt(SEP2)
  #Bias
  bias <- mean(predicted) - mean(observed)
  # residual variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[3] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  if (plot.it==TRUE){
    eqscplot(observed, predicted)
    abline(a = 0, b = 1, col = "brown4")
  }
  if (type == "DSM"){
    gf <- data.frame(R2 = R2,
                     concordance = ccc,
                     MSE = SEP2,
                     RMSE = SEP,
                     bias = bias,
                     row.names = NULL
    )
  }
  else if (type == "spec"){
    gf <- data.frame(R2 = R2,
                     concordance = ccc,
                     MSE = SEP2,
                     RMSE = SEP,
                     bias = bias,
                     MSEc = SEP2c,
                     RMSEc = SEPc,
                     RPD = RPD,
                     RPIQ = RPIQ,
                     row.names = NULL
    )
  }
  else {
    stop("ERROR: Revise the type of output you require. Select from either DSM or spec")
  }
  return(gf)
}

# 3) Loading data ---------------------------------------------------------
data <- read_delim("RegMat00_25cm.csv",
                   delim=",") %>% data.frame

# data$wetland <- as.factor(data$wetland)
# data$geology <- as.factor(data$geology)
# data$georeg <- as.factor(data$georeg)
# data$SlopeClass <- as.factor(data$SlopeClass)
# data$AUCupdClass <- as.factor(data$AUCupdClass)
# str(data)



# 4) Temporal trend slope -------------------------------------------------


# 4.1) Data splitting -----------------------------------------------------
names(data)
data <- data[,-c(1:6,8:18,51:57)]
names(data)
set.seed(1912)
# inTrain <- createDataPartition(y =data$Slope, p = .70, list = FALSE) #Random
inTrain <- kenStone(data, k = nrow(data)*0.70, metric = "mahal") # Kennard Stone

# train_data <- data[ inTrain,] %>% data.frame #Random
train_data <- data[ inTrain$model,] %>% data.frame #Kennard Stone
y_train <- train_data[,1]
x_train <- train_data[,c(2,3,34:93)]# without including the spectral indices
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(Slope=y_train,x_train)

# 4.2) PCA on spectral indices ----------------------------------------------
names(train_data)
pcaMed<-prcomp(train_data[,c(4:13)], scale=TRUE) 
summary(pcaMed)
(corvar <- pcaMed$rotation %*% diag(pcaMed$sdev))
Pred.pcs<-predict(pcaMed,train_data[,c(4:13)])
x_train$PCA1Med=Pred.pcs[,1] 
x_train$PCA2Med=Pred.pcs[,2]
x_train$PCA3Med=Pred.pcs[,3] 

pcaP10<-prcomp(train_data[,c(14:23)], scale=TRUE) 
summary(pcaP10)
(corvar <- pcaP10$rotation %*% diag(pcaP10$sdev))
Pred.pcs<-predict(pcaP10,train_data[,c(14:23)])
x_train$PCA1P10=Pred.pcs[,1] 
x_train$PCA2P10=Pred.pcs[,2]
x_train$PCA3P10=Pred.pcs[,3]

pcaP90<-prcomp(train_data[,c(24:33)], scale=TRUE) 
summary(pcaP90)
(corvar <- pcaP90$rotation %*% diag(pcaP90$sdev))
Pred.pcs<-predict(pcaP90,train_data[,c(24:33)])
x_train$PCA1P90=Pred.pcs[,1] 
x_train$PCA2P90=Pred.pcs[,2]
x_train$PCA3P90=Pred.pcs[,3]

x_train

# test_data <- data[-inTrain,] #Random
test_data <- data[inTrain$test,] # Kennard Stone
# test_data <- data[-indx$index_samples,] # CLHS


y_test <- test_data[,1]
x_test <- test_data[c(2,3,34:93)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(Slope=y_test,x_test)


Pred.pcs<-predict(pcaMed,test_data[,c(4:13)])
x_test$PCA1Med=Pred.pcs[,1] 
x_test$PCA2Med=Pred.pcs[,2]
x_test$PCA3Med=Pred.pcs[,3] 

Pred.pcs<-predict(pcaP10,test_data[,c(14:23)])
x_test$PCA1P10=Pred.pcs[,1] 
x_test$PCA2P10=Pred.pcs[,2]
x_test$PCA3P10=Pred.pcs[,3]

Pred.pcs<-predict(pcaP90,test_data[,c(24:33)])
x_test$PCA1P90=Pred.pcs[,1] 
x_test$PCA2P90=Pred.pcs[,2]
x_test$PCA3P90=Pred.pcs[,3]

x_train$geology9 <- NULL
x_train$geology11 <- NULL

x_test$geology9 <- NULL
x_test$geology11 <- NULL

train_data <- x_train
test_data <- x_test

summary(x_train)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]

hist(train_data$Slope+3)


# 6.2) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(1941)
  (bor <- Boruta(x = train_data[,c(2:70)],
                 y = train_data[,1], 
                 #data = train_data, 
                 doTrace = 0, 
                 ntree = 500,
                 maxRuns=500))
  plot(bor, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(bor$ImpHistory),function(i)
    bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
  names(lz) <- colnames(bor$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
  print(Sys.time() - start)
}

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

bor <- TentativeRoughFix(bor)
print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

preds <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

saveRDS(preds,"Outputs/NamesPreds/SlopeNum/PredictorsSlope_02082022.rds")

# 4.3) Models -------------------------------------------------------------

names(train_data)
fm <- as.formula(paste0("Slope~",paste0(preds,collapse = "+")))
fm

# 4.3.1) MLR --------------------------------------------------------------
model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data)

pred_lm <- predict(model_lm_up, newdata = test_data[,preds])
(lm.goof <- goof(observed = test_data$Slope, predicted = exp(pred_lm)))

# 4.3.2) Cubist -----------------------------------------------------------

ctrl <- trainControl(method = "boot",
                     summaryFunction = defaultSummary,
                     selectionFunction = "best"
)
grid <- expand.grid(committees = c(1, 10, 50), 
                    neighbors = c(1, 5, 9))
set.seed(49)
model_cubist <- train(fm,
                      data=train_data,
                      method = "cubist",
                      tuneGrid = grid,
                      trControl=ctrl
)

model_cubist
summary(model_cubist)

pred_cub <- predict(model_cubist, newdata = test_data[,preds])
(cub.goof <- goof(observed = test_data$Slope, predicted = pred_cub))


# 4.3.3) Random Forest ----------------------------------------------------
rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)
grid <- expand.grid(mtry = c(2,3,4),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4)
)

set.seed(37)
model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 500,
                  importance = "impurity"
)
model_rf

pred_rf <- predict(model_rf, newdata = test_data[,preds])
(rf.goof <- goof(observed = test_data$Slope, predicted = pred_rf))


# 4.3.4) Support vector machines ------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.2,0.02),
                                 cost = c(5,7,15,20)))

model_svm <- tuneResult$best.model
print(model_svm)

pred_svm <- predict(model_svm, newdata = test_data[,preds])
(svm.goof <- goof(observed = test_data$Slope, predicted = pred_svm))

# 4.3.5) Quantile random forest --------------------------------------------

set.seed(835)
model_qrf <- quantregForest(y = train_data[,"Slope"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = as.numeric(model_rf$bestTune))
model_qrf
importance(model_qrf,type = 2)


pred_qrf <- predict(model_qrf, newdata = test_data[,preds])
(qrf.goof <- goof(observed = test_data$Slope, predicted = pred_qrf[,2]))


#saveRDS(model_qrf,"Outputs/Models/.........")


# 5) Temporal trend slope -------------------------------------------------


# 3) Loading data ---------------------------------------------------------
data <- read_delim("RegMat00_25cm.csv",
                   delim=",") %>% data.frame

# data$wetland <- as.factor(data$wetland)
# data$geology <- as.factor(data$geology)
# data$georeg <- as.factor(data$georeg)
# data$SlopeClass <- as.factor(data$SlopeClass)
# data$AUCupdClass <- as.factor(data$AUCupdClass)
# str(data)



# 4) Temporal trend slope -------------------------------------------------


# 4.1) Data splitting -----------------------------------------------------
names(data)
data <- data[,c(12,19:50,58:117)]
names(data)
set.seed(905)
# inTrain <- createDataPartition(y =data$Slope, p = .70, list = FALSE) #Random
inTrain <- kenStone(data, k = nrow(data)*0.70, metric = "mahal") # Kennard Stone

# train_data <- data[ inTrain,] %>% data.frame #Random
train_data <- data[ inTrain$model,] %>% data.frame #Kennard Stone
y_train <- train_data[,1]
x_train <- train_data[,c(2,3,34:93)]# without including the spectral indices
max_train <- apply(x_train, 2, max)
min_train <- apply(x_train, 2, min)
x_train <- scale(x_train, center = min_train, scale = max_train-min_train)
x_train <- data.frame(AUC=y_train,x_train)

# 4.2) PCA on spectral indices ----------------------------------------------
names(train_data)
pcaMed<-prcomp(train_data[,c(4:13)], scale=TRUE) 
summary(pcaMed)
(corvar <- pcaMed$rotation %*% diag(pcaMed$sdev))
Pred.pcs<-predict(pcaMed,train_data[,c(4:13)])
x_train$PCA1Med=Pred.pcs[,1] 
x_train$PCA2Med=Pred.pcs[,2]
x_train$PCA3Med=Pred.pcs[,3] 

pcaP10<-prcomp(train_data[,c(14:23)], scale=TRUE) 
summary(pcaP10)
(corvar <- pcaP10$rotation %*% diag(pcaP10$sdev))
Pred.pcs<-predict(pcaP10,train_data[,c(14:23)])
x_train$PCA1P10=Pred.pcs[,1] 
x_train$PCA2P10=Pred.pcs[,2]
x_train$PCA3P10=Pred.pcs[,3]

pcaP90<-prcomp(train_data[,c(24:33)], scale=TRUE) 
summary(pcaP90)
(corvar <- pcaP90$rotation %*% diag(pcaP90$sdev))
Pred.pcs<-predict(pcaP90,train_data[,c(24:33)])
x_train$PCA1P90=Pred.pcs[,1] 
x_train$PCA2P90=Pred.pcs[,2]
x_train$PCA3P90=Pred.pcs[,3]

x_train

# test_data <- data[-inTrain,] #Random
test_data <- data[inTrain$test,] # Kennard Stone
# test_data <- data[-indx$index_samples,] # CLHS


y_test <- test_data[,1]
x_test <- test_data[c(2,3,34:93)]
x_test <- scale(x_test, center = min_train, scale = max_train-min_train)
x_test <- data.frame(AUC=y_test,x_test)


Pred.pcs<-predict(pcaMed,test_data[,c(4:13)])
x_test$PCA1Med=Pred.pcs[,1] 
x_test$PCA2Med=Pred.pcs[,2]
x_test$PCA3Med=Pred.pcs[,3] 

Pred.pcs<-predict(pcaP10,test_data[,c(14:23)])
x_test$PCA1P10=Pred.pcs[,1] 
x_test$PCA2P10=Pred.pcs[,2]
x_test$PCA3P10=Pred.pcs[,3]

Pred.pcs<-predict(pcaP90,test_data[,c(24:33)])
x_test$PCA1P90=Pred.pcs[,1] 
x_test$PCA2P90=Pred.pcs[,2]
x_test$PCA3P90=Pred.pcs[,3]

x_train$geology9 <- NULL
x_train$geology11 <- NULL

x_test$geology9 <- NULL
x_test$geology11 <- NULL

train_data <- x_train
test_data <- x_test

summary(x_train)
train_data <- train_data[complete.cases(train_data),]
summary(test_data)
test_data <- test_data[complete.cases(test_data),]




# 6.2) Boruta algorithm ---------------------------------------------------
names(train_data)
train_data <- train_data %>% na.omit %>% data.frame
{
  start <- Sys.time()
  set.seed(907)
  (bor <- Boruta(x = train_data[,c(2:70)],
                 y = train_data[,1], 
                 #data = train_data, 
                 doTrace = 0, 
                 ntree = 500,
                 maxRuns=500))
  plot(bor, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(bor$ImpHistory),function(i)
    bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
  names(lz) <- colnames(bor$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
  print(Sys.time() - start)
}

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

bor <- TentativeRoughFix(bor)
print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

preds <- names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

#saveRDS(preds,"Outputs/NamesPreds/SlopeNum/PredictorsSlope_02082022.rds")

# 4.3) Models -------------------------------------------------------------

names(train_data)
fm <- as.formula(paste0("AUC~",paste0(preds,collapse = "+")))
fm

# 4.3.1) MLR --------------------------------------------------------------
model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data[,preds])
(lm.goof <- goof(observed = test_data$AUC, predicted = pred_lm))

# 4.3.2) Cubist -----------------------------------------------------------

ctrl <- trainControl(method = "boot",
                     summaryFunction = defaultSummary,
                     selectionFunction = "best"
)
grid <- expand.grid(committees = c(1, 10, 50), 
                    neighbors = c(1, 5, 9))
set.seed(49)
model_cubist <- train(fm,
                      data=train_data,
                      method = "cubist",
                      tuneGrid = grid,
                      trControl=ctrl
)

model_cubist
summary(model_cubist)

pred_cub <- predict(model_cubist, newdata = test_data[,preds])
(cub.goof <- goof(observed = test_data$AUC, predicted = pred_cub))


# 4.3.3) Random Forest ----------------------------------------------------
rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)
grid <- expand.grid(mtry = c(2,3,4),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4)
)

set.seed(37)
model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 500,
                  importance = "impurity"
)
model_rf

pred_rf <- predict(model_rf, newdata = test_data[,preds])
(rf.goof <- goof(observed = test_data$AUC, predicted = pred_rf))


# 4.3.4) Support vector machines ------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.2,0.02),
                                 cost = c(5,7,15,20)))

model_svm <- tuneResult$best.model
print(model_svm)

pred_svm <- predict(model_svm, newdata = test_data[,preds])
(svm.goof <- goof(observed = test_data$AUC, predicted = pred_svm))

# 4.3.5) Quantile random forest --------------------------------------------

set.seed(835)
model_qrf <- quantregForest(y = train_data[,"AUC"],
                            x = train_data[,preds],
                            data=train_data,
                            keep.inbag=TRUE,
                            mtry = as.numeric(model_rf$bestTune))
model_qrf
importance(model_qrf,type = 2)


pred_qrf <- predict(model_qrf, newdata = test_data[,preds])
(qrf.goof <- goof(observed = test_data$AUC, predicted = pred_qrf[,2]))


#saveRDS(model_qrf,"Outputs/Models/.........")

