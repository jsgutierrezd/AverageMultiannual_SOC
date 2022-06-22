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
          'hydroGOF'
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
data$AUCupdClass <- as.factor(data$AUCupdClass)
str(data)



# 4) Temporal trend slope -------------------------------------------------


# 4.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(7,19:100)]

inTrain <- createDataPartition(y = data1$Slope, p = .70, list = FALSE)

train_data <- data1[ inTrain,]
test_data <- data1[-inTrain,]


# 4.2) Recursive feature elimination --------------------------------------
names(train_data)
{
  set.seed(2622)
  start <- Sys.time()
  cl <- parallel::makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=20, repeats=20,allowParallel = TRUE)
  (rfe <- rfe(x=train_data[,c(2:83)], y=train_data[,1], sizes=c(1:70), rfeControl=control2))
  print(Sys.time() - start)
}

plot(rfe, type=c("g", "o"))
predictors(rfe)


# 4.3) Models -------------------------------------------------------------

names(train_data)
fm <- as.formula(paste0("Slope~",paste0(as.character(predictors(rfe)),collapse = "+")))
fm

# 4.3.1) MLR --------------------------------------------------------------
model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data)

lm.goof <- gof(sim = pred_lm,obs = test_data$Slope)
lm.goof

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

pred_cubist <- predict(model_cubist, newdata = test_data)

cubist.goof <- gof(sim = pred_cubist,obs = test_data$Slope)
cubist.goof


# 4.3.3) Random Forest ----------------------------------------------------
rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)
grid <- expand.grid(mtry = c(2,3,4,5,6),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4,5,10)
)

set.seed(37)
model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 1000,
                  importance = "impurity"
)
model_rf

pred_rf <- predict(model_rf, newdata = test_data)

rf.goof <- gof(sim = pred_rf,obs = test_data$Slope)

rf.goof

# 4.3.4) Support vector machines ------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.2,0.02),
                                 cost = c(5,7,15,20)))

model_svm <- tuneResult$best.model
print(model_svm)

pred_svm <- predict(model_svm, newdata = test_data)
svm.goof <- gof(sim = pred_svm,obs = test_data$Slope)
svm.goof


# 5) Temporal trend slope -------------------------------------------------


# 5.1) Data splitting -----------------------------------------------------
names(data)
data1 <- data[,c(12,19:100)]

inTrain <- createDataPartition(y = data1$AUCupdated, p = .70, list = FALSE)

train_data <- data1[ inTrain,]
test_data <- data1[-inTrain,]


# 5.2) Recursive feature elimination --------------------------------------
names(train_data)
{
  set.seed(2622)
  start <- Sys.time()
  cl <- parallel::makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  control2 <- rfeControl(functions=rfFuncs, method="repeatedcv", number=20, repeats=20,allowParallel = TRUE)
  (rfe <- rfe(x=train_data[,c(2:83)], y=train_data[,1], sizes=c(1:20), rfeControl=control2))
  print(Sys.time() - start)
}

plot(rfe, type=c("g", "o"))
predictors(rfe)


# 5.3) Models -------------------------------------------------------------

names(train_data)
fm <- as.formula(paste0("AUCupdated~",paste0(as.character(predictors(rfe)),collapse = "+")))
fm

# 5.3.1) MLR --------------------------------------------------------------
model_lm <- lm(fm,
               data=train_data)
summary(model_lm)       
model_lm_up <- stepAIC(model_lm,direction="both")

summary(model_lm_up)

pred_lm <- predict(model_lm_up, newdata = test_data)

lm.goof <- gof(sim = pred_lm,obs = test_data$AUCupdated)
lm.goof

# 5.3.2) Cubist -----------------------------------------------------------

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

pred_cubist <- predict(model_cubist, newdata = test_data)

cubist.goof <- gof(sim = pred_cubist,obs = test_data$AUCupdated)
cubist.goof


# 5.3.3) Random Forest ----------------------------------------------------
rctrlG <- trainControl(method = "repeatedcv",
                       number = 20,
                       repeats = 20,
                       returnResamp = "all",
                       search = "grid"
)
grid <- expand.grid(mtry = c(2,3,4,5,6),
                    splitrule = c("variance", "extratrees"),
                    min.node.size = c(2,3,4,5,10)
)

set.seed(37)
model_rf <- train(fm,
                  data=train_data,
                  method = "ranger",
                  trControl = rctrlG,
                  tuneGrid = grid,
                  num.trees = 1000,
                  importance = "impurity"
)
model_rf

pred_rf <- predict(model_rf, newdata = test_data)

rf.goof <- gof(sim = pred_rf,obs = test_data$AUCupdated)

rf.goof

# 5.3.4) Support vector machines ------------------------------------------

tuneResult <- tune(svm, fm, data = train_data,
                   ranges = list(epsilon = seq(0.1,0.2,0.02),
                                 cost = c(5,7,15,20)))

model_svm <- tuneResult$best.model
print(model_svm)

pred_svm <- predict(model_svm, newdata = test_data)
svm.goof <- gof(sim = pred_svm,obs = test_data$AUCupdated)
svm.goof





