preds <- as.data.frame(do.call(cbind, rfe_list)) 

uniqvars <- data.frame(x=1:78)
for (i in 1:ncol(preds)) {
  tmp <- unique(preds[,i])
  tmp <- c(tmp,rep(NA,times=78-length(tmp)))
  uniqvars <- cbind(uniqvars,tmp)
}

uniqvars[,1] <- NULL
uniqvars <- transpose(uniqvars)

find_mode <- function(x) {
  u <- unique(na.omit(x))
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

imppred <- c()
for (i in 1:78) {
  tmp <- find_mode(uniqvars[,i])
  imppred <- unique(c(imppred,tmp))
}

# Main possible predictors according to repeated RFE
imppred




data(mdrr)
mdrrDescr <- mdrrDescr[,-nearZeroVar(mdrrDescr)]
mdrrDescr <- mdrrDescr[, -findCorrelation(cor(mdrrDescr), .8)]

set.seed(1)
inTrain <- createDataPartition(mdrrClass, p = .75, list = FALSE)[,1]

train <- mdrrDescr[ inTrain, ]
test  <- mdrrDescr[-inTrain, ]
trainClass <- mdrrClass[ inTrain]
testClass  <- mdrrClass[-inTrain]

set.seed(2)

ldaProfile <- rfe(train_data_list[[1]][,c(2:)], train_data_list[[1]][,1],
                  sizes = c(1:30),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))
plot(ldaProfile, type = c("o", "g"))

postResample(predict(ldaProfile, test), testClass)

ldaProfile$optVariables



set.seed(2)
ldaProfile <- rfe(train, trainClass,
                  sizes = c(1:10, 15, 30),
                  rfeControl = rfeControl(functions = ldaFuncs, method = "cv"))
plot(ldaProfile, type = c("o", "g"))

postResample(predict(ldaProfile, test), testClass)



names(data1)
names(train_data_list[[i]])
install.packages("doMC")
library(doMC)
{
  i=1
  start <- Sys.time()
  cl <- parallel::makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  control2 <- rfeControl(functions=caretFuncs, method="cv",allowParallel = TRUE)
  rfe_vars <- list()
  for (i in 1:length(inTrain)) {
    set.seed(i)
    rfe_vars[[i]] <- rfe(x=train_data_list[[i]][,c(2:7,9:32,36:39,42:78)], 
                         y=train_data_list[[i]][,1], 
                         sizes=c(1:30), 
                         rfeControl=control2,method="svmRadial")# %>% predictors()
    
  }
  print(Sys.time() - start)
}




##Seleccion de variables --> algoritmo Boruta
start <- Sys.time()
(bor <- Boruta(pH.0_30_Sum.Pond ~ ., data = data, doTrace = 0, ntree = 500,maxRuns=500))
plot(bor, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(bor$ImpHistory),function(i)
  bor$ImpHistory[is.finite(bor$ImpHistory[,i]),i])
names(lz) <- colnames(bor$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(bor$ImpHistory), cex.axis = 0.7)
print(Sys.time() - start)

print(bor)
names(bor$finalDecision[bor$finalDecision %in% c("Confirmed")])

# final.bor <- TentativeRoughFix(bor)
# print(final.bor)
# 
# getSelectedAttributes(bor, withTentative = F)
# boruta.df <- attStats(final.bor)
plot(tmp, type = c("o", "g"))
tmp







# 4.3) Models -------------------------------------------------------------

fm_list <- list()
for (i in 1:length(rfe_list)) {
  fm_list[[i]] <- as.formula(paste0("SlopeClass~",paste0(as.character(rfe_list[[i]]),collapse = "+")))
}

fm_list
fm <- as.formula(paste0("SlopeClass~",paste0(as.character(imppred),collapse = "+")))
fm

# 4.3.1) MLP --------------------------------------------------------------
# start <- Sys.time()
# model_mlp_list <- list()
# for (i in 1:length(fm_list)) {
#   model_mlp_list[[i]] <- nnet(fm_list[[i]],
#                               data=train_data_list[[i]],
#                               size = 2,
#                               decay = 1e-5,
#                               maxit = 100
#   ) 
# }
# print(Sys.time() - start)
# 
# pred_mlp_list <- list()
# for (i in 1:length(model_mlp_list)) {
#   pred_mlp_list[[i]] <- predict(model_mlp_list[[i]], newdata = test_data_list[[i]], type = "class")
# }
# 
# confMat_dataframe <- data.frame()
# for (i in 1:length(pred_mlp_list)) {
#   tmp <- c(confusionMatrix(as.factor(pred_mlp_list[[i]]) ,test_data_list[[i]][,1])$overall)
#   confMat_dataframe <- rbind(confMat_dataframe,tmp)
#   names(confMat_dataframe) <- c("Accuracy","Kappa","AccuracyLower",
#                                 "AccuracyUpper","AccuracyNull","AccuracyPValue",
#                                 "McNemarPValue")
#   
# }
# confMat_dataframe
# 4.3.2) Random Forest ----------------------------------------------------

rctrlG <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       returnResamp = "all",
                       search = "grid",
                       classProbs = TRUE
)

grid <- expand.grid(mtry = c(2,3,4,5),
                    splitrule = c("gini", "extratrees"),
                    min.node.size = c(1,2,3,4,5)
)

# set.seed(37)
model_rf_list <- list()
start <- Sys.time()
for (i in 1:length(train_data_list)) {
  set.seed(i)
  model_rf_list[[i]] <- train(fm,
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
hist(confMat_dataframe$Kappa)
hist(confMat_dataframe$Accuracy)
# 4.3.3) Support vector machines ------------------------------------------



# 4.4) SlopeClass Model generalization ------------------------------------

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
  sp_prob[[i]] <- clusterR(cov, 
                           fun = predict,  
                           args=list(model=model_rf_list[[i]], 
                                     type="prob"))
  
}
sp_prob
endCluster()
print(Sys.time() - start)

# 5) AUCupdClass ----------------------------------------------------------

# 5.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(9,19:93)]

set.seed(17)
inTrain <- createDataPartition(y = data1$SlopeClass,
                               times = 100,
                               p = .70, 
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
  control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, repeats=10,allowParallel = TRUE)
  rfe_list <- list()
  for (i in 1:length(inTrain)) {
    set.seed(i)
    rfe_list[[i]] <- rfe(x=train_data_list[[i]][,c(2:76)], 
                         y=train_data_list[[i]][,1], 
                         sizes=c(1:70), 
                         rfeControl=control2) %>% predictors()
  }
  print(Sys.time() - start)
}
rfe_list

#plot(rfe, type=c("g", "o"))
#predictors(rfe)

# 5.3) Models -------------------------------------------------------------

fm_list <- list()
for (i in 1:length(bor_list)) {
  fm_list[[i]] <- as.formula(paste0("SlopeClass ~",paste0(as.character(bor_list[[i]]),collapse = "+")))
}

fm_list
# fm <- as.formula(paste0("AUCupdClass~",paste0(as.character(rfe_list[[i]]),collapse = "+")))
# fm

# 5.3.1) MLP --------------------------------------------------------------
# start <- Sys.time()
# model_mlp_list <- list()
# for (i in 1:length(fm_list)) {
#   model_mlp_list[[i]] <- nnet(fm_list[[i]],
#                               data=train_data_list[[i]],
#                               size = 2,
#                               decay = 1e-5,
#                               maxit = 100
#   ) 
# }
# print(Sys.time() - start)
# 
# pred_mlp_list <- list()
# for (i in 1:length(model_mlp_list)) {
#   pred_mlp_list[[i]] <- predict(model_mlp_list[[i]], newdata = test_data_list[[i]], type = "class")
# }
# 
# confMat_dataframe <- data.frame()
# for (i in 1:length(pred_mlp_list)) {
#   tmp <- c(confusionMatrix(as.factor(pred_mlp_list[[i]]) ,test_data_list[[i]][,1])$overall)
#   confMat_dataframe <- rbind(confMat_dataframe,tmp)
#   names(confMat_dataframe) <- c("Accuracy","Kappa","AccuracyLower",
#                                 "AccuracyUpper","AccuracyNull","AccuracyPValue",
#                                 "McNemarPValue")
#   
# }
confMat_dataframe
# 5.3.2) Random Forest ----------------------------------------------------

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
for (i in 1:length(fm_list)) {
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
hist(confMat_dataframe$Kappa)
hist(confMat_dataframe$Accuracy)
