###################################################
### code Read and partition data 
###################################################
# read data with T2
rpart_inputdata_wT2 <- function(subdata, typefeat) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "newT2wUpdatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                             INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                             INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)
                             INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesionfields =  names(lesionsQuery)
  lesioninfo = lesionsQuery[c(1:26)]
  dynfeatures = lesionsQuery[c(29:62)]
  morphofeatures = lesionsQuery[c(65:83)]
  texfeatures = lesionsQuery[c(86:129)]
  stage1features = lesionsQuery[c(132:231)]
  T2info = lesionsQuery[c(253:265)]
  T2features = lesionsQuery[c(263:264,266:287,232:251)]
  
  
  # combine  features according to typefeat
  if(typefeat=="both"){
    allfeatures = cbind(lesioninfo[c(1:2)], dynfeatures, morphofeatures, texfeatures, stage1features, T2features)   
  }
  if(typefeat=="T1w"){
    allfeatures = cbind(lesioninfo[c(1:2)], dynfeatures, morphofeatures, texfeatures, stage1features)   
  }
  if(typefeat=="T2w"){
    allfeatures = cbind(lesioninfo[c(1:2)], T2features)   
  }
  
  # subdata accordingly
  if(subdata=="stage2"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "NC", "C") -> M$lesion_label
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "NC", "C") -> N$lesion_label
    Fo<-subset(allfeatures, lesion_label=="fociB" | lesion_label=="fociM")
    ifelse( Fo$lesion_label == "fociB", "NC", "C") -> Fo$lesion_label
    allfeatures = data.frame(rbind(M,N,Fo)) 
  }
  if(subdata=="stage1"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "mass", "mass") -> M$lesion_label
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "nonmass", "nonmass") -> N$lesion_label
    Fo<-subset(allfeatures, lesion_label=="fociB" | lesion_label=="fociM")
    ifelse( Fo$lesion_label == "fociB", "foci", "foci") -> Fo$lesion_label
    allfeatures = data.frame(rbind(M,N,Fo)) 
  }
  if(subdata=="oneshot"){
    # organized the data by subdata
    M<-subset(allfeatures, lesion_label=="massB" | lesion_label=="massM")
    ifelse( M$lesion_label == "massB", "NC", "C") -> M$lesion_label
    N<-subset(allfeatures, lesion_label=="nonmassB" | lesion_label=="nonmassM")
    ifelse( N$lesion_label == "nonmassB", "NC", "C") -> N$lesion_label
    Fo<-subset(allfeatures, lesion_label=="fociB" | lesion_label=="fociM")
    ifelse( Fo$lesion_label == "fociB", "NC", "C") -> Fo$lesion_label
    allfeatures = data.frame(rbind(M,N,Fo)) 
  }
  # procees data
  allfeatures$lesion_label <- as.factor(allfeatures$lesion_label)
  # allfeatures$peakCr_inside <- as.factor(allfeatures$peakCr_inside)
  # allfeatures$peakVr_inside <- as.factor(allfeatures$peakVr_inside)
  # allfeatures$peakCr_countor <- as.factor(allfeatures$peakCr_countor)
  # allfeatures$peakVr_countor <- as.factor(allfeatures$peakVr_countor)
  # allfeatures$k_Max_Margin_Grad <- as.factor(allfeatures$k_Max_Margin_Grad)
  # allfeatures$max_RGH_mean_k <- as.factor(allfeatures$max_RGH_mean_k)
  # allfeatures$max_RGH_var_k <- as.factor(allfeatures$max_RGH_var_k)
  
  # 2) Query info
  lesionsQueryinfo <- dbGetQuery(conn, "SELECT *
                                 FROM  lesion
                                 INNER JOIN f_T2 ON (f_T2.lesion_id = lesion.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesionsfields = names(lesionsQueryinfo)
  lesionsQueryinfo = lesionsQueryinfo[c(1:26,29:39)]
  
  output <- list(features=allfeatures, info=lesionsQueryinfo)
  return(output)
}

###################################################
### code to create a cross-validation set up: 
### cvfoldk = number of cv folds typically 5 or 10
### out: particvfoldK = all cv-K ids
###################################################
library(caret)

cvfold_partition <- function(dat, cvfoldK){
  ndat = nrow(dat)
  outcomesetDi  <- dat$lesion_label
  #For multiple k-fold cross-validation, completely independent folds are created.
  #when y is a factor in an attempt to balance the class distributions within the splits.
  #The names of the list objects will denote the fold membership using the pattern 
  #"Foldi.Repj" meaning the ith section (of k) of the jth cross-validation set (of times).
  partitionsetDi <- createFolds(y = outcomesetDi, ## the outcome data are needed
                                k = cvfoldK, ## The percentage of data in the training set
                                list = TRUE) ## The format of the results. 
  return(partitionsetDi)
}

###################################################
### code to sample kparti from a cross-validation set up: 
### kparti = k fold to exclude
### outs: cvTrainsetD, cvTestsetD
###################################################
kparti_sample <- function(dat, particvfoldK, cvfoldK, kparti){
  allparti = 1:cvfoldK
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, particvfoldK[[kadd]])
  }
  # partition data
  cvTrainsetD <-  dat[ cvfoldadd ,]
  cvTestsetD <-   dat[-cvfoldadd ,]
  
  output <- list(cvTrainsetD=cvTrainsetD, cvTestsetD=cvTestsetD)
  return(output)
}



### code to perform RRF featsel for lop
### parameters, featTrain, type
RRF_featsel <- function(featTrain, type){
  library(RRF)
  ### 
  Xfs = na.omit(featTrain)
  Ys = Xfs$lesion_label
  Xfs = Xfs[,-1]
  
  set.seed(1)
  bestmtry = tuneRRF(Xfs, Ys, mtryStart = 1, ntreeTry=1000, doBest=FALSE, plot=FALSE, trace=FALSE)
  mtryind = which.min(as.data.frame(bestmtry)$OOBError)
  fs_rrf = RRF(Xfs, Ys, mtry=bestmtry[mtryind], flagReg = 1, 
               ntree=2000, 
               localImp=TRUE,
               proximity=TRUE)
  #print(fs_rrf)
  # overall feature importance
  #varImp_rrf = data.frame(varImpPlot(fs_rrf, sort=TRUE))
  
  varImp_rrf = data.frame(fs_rrf$importance)
  # sort features by MeanDecAccuracy
  varImp = varImp_rrf[ order(-varImp_rrf$MeanDecreaseGini), ] 
  # pick only non-zero variables
  varImp = unique(varImp)
  df = data.frame(selfeat=rownames(varImp))
  df$MeanDecreaseGini = varImp$MeanDecreaseGini 
  df$SelectedFeatureGroup = type
  # adject
  df$selfeat = as.character(df$selfeat)
  print(cat("Selected features for group: MeanDecreaseGini",type,"\n========="))
  print(df$selfeat)
  
  return(df)
}


###################################################
### code to Train and compare models
###################################################
library(caret)
library(dplyr)         # Used by caret
library(kernlab)       # support vector machine 

train_tune_Models <- function(trainData, testData, typeSampling){
  set.seed(1)
  # Setup for cross validation
  ctrl <- trainControl(method="repeatedcv",   # cross validation
                       number=10, # 10fold
                       repeats=1,		    # do 5 repititions of cv
                       summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                       classProbs=TRUE)
  
  ######################
  if(typeSampling=="oversmin"){
    # add class balance: oversample minclassn classn to match majclassn
    majclassn = max(summary(trainData$lesion_label))
    minclassf = trainData[trainData$lesion_label=="C",]
    setD = minclassf[sample(1:nrow(minclassf), majclassn, replace=TRUE),]
    setD = rbind(setD, trainData[trainData$lesion_label=="NC",])
    trainX = setD[,c(3:ncol(trainData))]
    trainy = setD$lesion_label
    pander(summary(trainy))
  }
  if(typeSampling=="undersmaj"){
    # add class balance: downsample majority classn to match minclassn
    minclassn = min(summary(trainData$lesion_label))
    maxclassf = trainData[trainData$lesion_label=="NC",]
    setD = maxclassf[sample(1:nrow(maxclassf), minclassn, replace=FALSE),]
    setD = rbind(setD, trainData[trainData$lesion_label=="C",])
    trainX = setD[,c(3:ncol(trainData))]
    trainy = setD$lesion_label
    pander(summary(trainy))
  }
  if(typeSampling=="regular"){
    trainX = trainData[,c(3:ncol(trainData))]
    trainy = trainData$lesion_label
    pander(summary(trainy))
  }
  
  ######################
  #Train and Tune the models
  # the same random number seed is set prior to each tunning.
  # This ensures that the same resampling sets are used, 
  # when we compare the resampling profiles between models.
  ######################
  #Train and Tune the SVM
  set.seed(1)
  svmRadial <- train(x=trainX,
                     y=trainy,
                     method = "svmRadial",   # Radial kernel
                     tuneLength = 10,					# 9 values of the cost function
                     preProc = c("center","scale"),  # Center and scale data
                     metric="ROC",
                     trControl=ctrl)
  print(svmRadial)
  print("Best performing parameters")
  svmRadial$results[svmRadial$results$ROC==max(svmRadial$results$ROC),]
  
  ## second pass
  # Use the expand.grid to specify the search space	
  sigma = svmRadial$bestTune$sigma
  C = svmRadial$bestTune$C
  tryC = c()
  for(i in (C-4):(C+4)){
    if(i%% 2 == 0 && i>0){
      tryC = c(tryC, i)}}
  
  grid <- expand.grid(sigma = c(sigma-0.01, sigma-0.005, sigma, sigma+0.005, sigma+0.01),
                      C = tryC)
  
  # Train and Tune the SVM 2
  set.seed(1)
  svmRadial2 <- train(x=trainX,
                      y=trainy,
                      method = "svmRadial",
                      preProc = c("center","scale"),
                      metric="ROC",
                      tuneGrid = grid,
                      trControl=ctrl,
                      verboseIter=FALSE)
  
  print("Best performing parameters")
  svmRadial2$results[svmRadial2$results$ROC==max(svmRadial2$results$ROC),]
  
  ###########################
  # 2) Train and Tune a RF
  ########################### 
  set.seed(1)
  rf <- train(x=trainX,
              y=trainy,
              method = "rf",   # Random Forest model
              tuneLength = 9,					# 3 values of mtry
              preProc = c("center","scale"),  # Center and scale data
              metric="ROC",
              trControl=ctrl)
  print(rf)
  print("10fold cross validation RF Best performing parameters")
  rf$results[rf$results$ROC == max(rf$results$ROC),]
  
  ###########################
  # 3) Train and Tune an obliqueRF
  ###########################
  set.seed(1)
  ORF <- train(x=trainX,
               y=trainy,
               method = "ORFlog",   # Oblique Random Forest model
               tuneLength = 3,					# 3 values of mtry
               preProc = c("center","scale"),  # Center and scale data
               metric="ROC",
               trControl=ctrl,
               verbose=FALSE)
  print(ORF)
  print("10fold cross validation obliqueRF Best performing parameters")
  ORF$results[ORF$results$ROC == max(ORF$results$ROC),]
  
  
  ###########################
  # 4) Train ensembles of Boosted Classification Trees
  ########################### 
  set.seed(1)
  ada_btrees <- train(x=trainX,
                      y=trainy,
                      method = "ada",   # Ensembles of Generalized Lienar Models
                      tuneLength = 3,					# 3 values of iter, maxdepth, nu
                      preProc = c("center","scale"),  # Center and scale data
                      metric="ROC",
                      trControl=ctrl)
  print(ada_btrees)
  print("10fold cross validation ada_btrees Best performing parameters")
  ada_btrees$results[ada_btrees$results$ROC == max(ada_btrees$results$ROC),]
  
  output = list(RBFsigma=svmRadial2, rf=rf, ORF=ORF, ada_btrees=ada_btrees)
  return(output)  
}


plot_train_tune_Models <- function(outputs, testx, testy, title){
  ########## plot comparison
  # resample imgT1w
  outputs[[1]]$resample = na.omit(outputs[[1]]$resample)
  res_1 = data.frame(AUC = mean(outputs[[1]]$resample$ROC),
                          std = sd(outputs[[1]]$resample$ROC),
                          group = "cv",
                          classifier=names(outputs[1]))
  #append test resulst
  proby_1_test = predict(outputs[[1]], newdata=testX, type="prob")
  ROCF_1 <- roc(testy, proby_1_test$C, plot=FALSE)
  res_1 = rbind(res_1, 
                     data.frame(AUC = ROCF_1$auc,
                                std = 0,
                                group = "test",
                                classifier=names(outputs[1])))
  ############
  # resample withT2w
  outputs[[2]]$resample = na.omit(outputs[[2]]$resample)
  res_2 = data.frame(AUC = mean(outputs[[2]]$resample$ROC),
                     std = sd(outputs[[2]]$resample$ROC),
                     group = "cv",
                     classifier=names(outputs[2]))
  #append test resulst
  proby_2_test = predict(outputs[[2]], newdata=testX, type="prob")
  ROCF_2 <- roc(testy, proby_2_test$C, plot=FALSE)
  res_2 = rbind(res_2, 
                data.frame(AUC = ROCF_2$auc,
                           std = 0,
                           group = "test",
                           classifier=names(outputs[2])))
  ############
  # resample withLMSIR
  outputs[[3]]$resample = na.omit(outputs[[3]]$resample)
  res_3 = data.frame(AUC = mean(outputs[[3]]$resample$ROC),
                     std = sd(outputs[[3]]$resample$ROC),
                     group = "cv",
                     classifier=names(outputs[3]))
  #append test resulst
  proby_3_test = predict(outputs[[3]], newdata=testX, type="prob")
  ROCF_3 <- roc(testy, proby_3_test$C, plot=FALSE)
  res_3 = rbind(res_3, 
                data.frame(AUC = ROCF_3$auc,
                           std = 0,
                           group = "test",
                           classifier=names(outputs[3])))
  
  
  # append all 5 classifier resutls
  res_outputs = rbind(res_1, res_2, res_3)
  print(res_outputs)
  
  # bar plot
  p <- ggplot(data=res_outputs, aes(x=classifier, y=AUC, fill=group, group=group)) +
    geom_bar(stat="identity",color="black", position=position_dodge()) +
    geom_errorbar(aes(ymin=AUC-std, ymax=AUC+std), width=.2, position=position_dodge(0.9)) +
    geom_text(data=data.frame(res_outputs), 
              aes(label=formatC(AUC,digits=2, format="f")),
              position=position_dodge(0.9), vjust=1.5,
              size=4)
  
  # Finished bar plot
  print(p + labs(title=paste0("Comparison Models ",title)) +
          scale_y_continuous(limits = c(0, 1), breaks=0:10/10))
  
  return(res_outputs)
}


###########################
# 4) Train ensembles of Boosted Classification Trees (selective features)
########################### 
library(caret) 
library(dplyr)         # Used by caret
library(kernlab)       # support vector machine
library(adabag)
require(ggplot2)
library(pROC)

train_selectiveAdaboost <- function(dat, testData){
  
  ###########################
  # Train ensembles of Boosted Classification Trees
  ########################### 
  # Setup for cross validation
  ctrl <- trainControl(method="repeatedcv",   # cross validation
                       number=10, # 10fold
                       repeats=1,		    # do 5 repititions of cv
                       summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                       classProbs=TRUE)
  
  grid = expand.grid(iter=c(100,150),
                     maxdepth=c(3,5,10),
                     nu=c(0.01,0.05,0.1))
  
  ########################### 
  set.seed(312)
  adacv <- train(x=dat[,2:ncol(dat)],
                      y=dat$lesion_label,
                      method = "ada",   # Ensembles of Generalized Lienar Models
                      metric="ROC",
                      trControl=ctrl,
                      tuneGrid = grid,
                      verbose=TRUE)
  print(adacv)
  print("10fold cross validation adacv Best performing parameters")
  beastperf = adacv$results[adacv$results$ROC == max(adacv$results$ROC),]
  print(beastperf)
  
  # train final ensemble
  adaf <- boosting(lesion_label ~ ., 
                        data = dat,
                        boos=TRUE, 
                        mfinal=beastperf$iter, 
                        control = rpart.control(maxdepth = beastperf$maxdepth),
                        coeflearn='Freund')
  
  ######### plot errors
  errTrain = errorevol(adaf,dat)$error
  errTest = errorevol(adaf,testData)$error
  
  errorTrees = rbind(data.frame(bTrees=1:length(errTrain), errorRate=errTrain, type="train"),
                     data.frame(bTrees=1:length(errTest), errorRate=errTest, type="test"))
  #plot
  print(ggplot(data=errorTrees, aes(x=bTrees, y=errorRate, group=type, color=type)) +
          geom_line() +
          geom_point() )
  
  ######### find min errorTest ensembles
  minTest_ntrees = (1:length(errTest))[errTest==min(errTest)]
  bestada = list()
  bestada$trees = adaf$trees[1:minTest_ntrees]
  bestada$weights = adaf$weights[1:minTest_ntrees]
  
  ntrees = length(bestada$trees)
  sum_wpred = 0
  for (t in 1:ntrees){
    # Calcultate posterior Probabilities on grid points
    temp <- predict(bestada$trees[[t]], newdata = testData) 
    tempw <- temp*bestada$weights[[t]]
    # weight each tree by its alpha coefficient
    sum_wpred = sum_wpred+tempw
  }
  testsel = sum_wpred/ntrees
  
  ######### predict AUC
  testprob_adacv = predict(adacv, newdata=testData, type="prob")
  ROCF_testadacv <- roc(testData$lesion_label, testprob_adacv$C, plot=FALSE)
  print(ROCF_testadacv)
  
  testprob_adaf = predict(adaf, newdata=testData, type="prob")
  ROCF_testadaf <- roc(testData$lesion_label, testprob_adaf$prob[,1], plot=FALSE)
  print(ROCF_testadaf)
  
  ROCF_testsel <- roc(testData$lesion_label, testsel[,1], plot=FALSE)
  print(ROCF_testsel)
  
  # append perfm for heldout ROC
  output <- list(tunning = list(beastperf = beastperf, 
                                resample = adacv$resample),
                 
                 ensemble_heldoperf = list( testprob_adacv=testprob_adacv,
                                            testprob_adaf=testprob_adaf,
                                            testsel=testsel,
                                            adacv = adacv,
                                            adaf = adaf,
                                            bestada = bestada))
  return(output)
}


create_selectiveAdaboost_ensemble <- function(dat_T1w, dat_both, particvfoldK, cvK, ntrees, typeSampling){
  #inint dat= dat_both
  tunning = list()
  ensemble_heldoperf=list()
  
  for(r in 1:cvK){
    ## pick one of cvfold for held-out test, train on the rest
    print(paste("cvfold for held-out test = ",r))
    kparti_setT1w = kparti_sample(dat_T1w, particvfoldK, cvK, r)
    kparti_setboth = kparti_sample(dat_both, particvfoldK, cvK, r)
    
    # # Boruta on $cvTrainsetD
    selfeatures_kfold = RRF_featsel(kparti_setT1w$cvTrainsetD[,2:ncol(kparti_setT1w$cvTrainsetD)], type="")

    # Subset
    ###################### 
    trainData <-  na.omit(kparti_setboth$cvTrainsetD)
    testData <-  na.omit(kparti_setboth$cvTestsetD)
    
    # Sampling
    ######################
    if(typeSampling=="oversmin"){
      # add class balance: oversample minclassn classn to match majclassn
      majclassn = max(summary(trainData$lesion_label))
      minclassf = trainData[trainData$lesion_label=="C",]
      setD = minclassf[sample(1:nrow(minclassf), majclassn, replace=TRUE),]
      setD = rbind(setD, trainData[trainData$lesion_label=="NC",])
      pander(summary(setD$lesion_label))
    }
    if(typeSampling=="undersmaj"){
      # add class balance: downsample majority classn to match minclassn
      minclassn = min(summary(trainData$lesion_label))
      maxclassf = trainData[trainData$lesion_label=="NC",]
      setD = maxclassf[sample(1:nrow(maxclassf), minclassn, replace=FALSE),]
      setD = rbind(setD, trainData[trainData$lesion_label=="C",])
      pander(summary(setD$lesion_label))
      
    }
    if(typeSampling=="regular"){
      setD = trainData
      pander(summary(setD$lesion_label))
    }

    ######################
    # Train and Tune the models
    ######################
    # Setup for cross validation
    errorTrees = data.frame()
    
    # set datasets
    data_onlyT1w = setD[,c("lesion_label",selfeatures_kfold$selfeat)]
    data_withLMSIR = cbind(data_onlyT1w, setD[,c(200)])
    data_withT2w = cbind(data_onlyT1w, setD[,c(200:243)])
    
    ###########################
    # Train ensembles of Boosted Classification Trees
    ########################### 
    # Setup for cross validation
    ctrl <- trainControl(method="repeatedcv",   # cross validation
                         number=10, # 10fold
                         repeats=1,		    # do 5 repititions of cv
                         summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                         classProbs=TRUE)
    
    grid = expand.grid(iter=c(25,50,100),
                       maxdepth=c(1,3,5),
                       nu=c(0.01,0.1))
    
    ###########################  for data_onlyT1w
    set.seed(312)
    adacv_wT2w <- train(x=data_withT2w[,2:ncol(data_withT2w)],
                        y=data_withT2w$lesion_label,
                        method = "ada",   # Ensembles of Generalized Lienar Models
                        metric="ROC",
                        trControl=ctrl,
                        tuneGrid = grid,
                        verbose=TRUE)
    print(adacv_wT2w)
    print("10fold cross validation adacv_wT2w Best performing parameters")
    beastperf = adacv_wT2w$results[adacv_wT2w$results$ROC == max(adacv_wT2w$results$ROC),]
    print(beastperf)
    
    # train final ensemble
    adaf_wT2w <- boosting(lesion_label ~ ., 
                            data = data_withT2w,
                            boos=TRUE, 
                            mfinal=beastperf$iter, 
                            control = rpart.control(maxdepth = beastperf$maxdepth),
                            coeflearn='Freund')
    
    ######### plot errors
    errTrain = errorevol(adaf_wT2w,trainData)$error
    errTest = errorevol(adaf_wT2w,testData)$error
    
    errorTrees = rbind(data.frame(bTrees=1:length(errTrain), errorRate=errTrain, type="train"),
                       data.frame(bTrees=1:length(errTest), errorRate=errTest, type="test"))
    #plot
    print(ggplot(data=errorTrees, aes(x=bTrees, y=errorRate, group=type, color=type)) +
      geom_line() +
      geom_point() )
    
    ######### find min errorTest ensembles
    minTest_ntrees = (1:length(errTest))[errTest==min(errTest)]
    bestada_wT2w = list()
    bestada_wT2w$trees = adaf_wT2w$trees[1:minTest_ntrees]
    bestada_wT2w$weights = adaf_wT2w$weights[1:minTest_ntrees]
    
    ntrees = length(bestada_wT2w$trees)
    sum_wpred = 0
    for (t in 1:ntrees){
      # Calcultate posterior Probabilities on grid points
      temp <- predict(bestada_wT2w$trees[[t]], newdata = testData) 
      tempw <- temp*bestada_wT2w$weights[[t]]
      # weight each tree by its alpha coefficient
      sum_wpred = sum_wpred+tempw
    }
    testsel = sum_wpred/ntrees
    
    ######### predict AUC
    testprob_adacv = predict(adacv_wT2w, newdata=testData, type="prob")
    ROCF_testadacv <- roc(testData$lesion_label, testprob_adacv$C, plot=FALSE)
    print(ROCF_testadacv)
    
    testprob_adaf = predict(adaf_wT2w, newdata=testData, type="prob")
    ROCF_testadaf <- roc(testData$lesion_label, testprob_adaf$prob[,1], plot=FALSE)
    print(ROCF_testadaf)
    
    ROCF_testsel <- roc(testData$lesion_label, testsel[,1], plot=TRUE)
    print(ROCF_testsel)
    
    
    # append perfm for heldout ROC
    tunning = append(tunning, list(tunning = list(beastperf = beastperf,
                                                  resample = RF_trees$resample)))
    
    ensemble_heldoperf = append(ensemble_heldoperf, 
                                list(ensemble_heldoperf = list( 
                                  adacv_wT2w = adacv_wT2w,
                                  adaf_wT2w = adaf_wT2w,
                                  bestada_wT2w = bestada_wT2w)))
  }
  
  output <- list(tunning = tunning, ensemble_heldoperf=ensemble_heldoperf)
  return(output)
}    


###################################################
# Build a classifier with internal cv of parameters
###################################################
boosting_Train_wcv_wperf <- function(TrainsetD, TestsetD, typeSampling){
  # TrainsetD=imgT2pLMSIR
  # TestsetD=T1T2test[c(names(imgT2pLMSIR))]
  # typeSampling="oversmin"
  library(pROC)
  # to eval
  test_roc <- function(model, testdata) {
    # create call
    testpred <- predict(model, newdata=testdata, type="prob")
    roc_obj <- roc(testdata$lesion_label, testpred$prob[,1])
    output <- list(auc=roc_obj$auc, ci=ci(roc_obj), testpred=testpred$prob)
    return(output)
  }
  
  ## Sampling
  ######################
  if(typeSampling=="oversmin"){
    # add class balance: oversample minclassn classn to match majclassn
    majclassn = max(summary(TrainsetD$lesion_label))
    minclassf = TrainsetD[TrainsetD$lesion_label=="C",]
    setD = minclassf[sample(1:nrow(minclassf), majclassn, replace=TRUE),]
    setD = rbind(setD, TrainsetD[TrainsetD$lesion_label=="NC",])
    pander(summary(setD$lesion_label))
  }
  if(typeSampling=="undersmaj"){
    # add class balance: downsample majority classn to match minclassn
    minclassn = min(summary(TrainsetD$lesion_label))
    maxclassf = TrainsetD[TrainsetD$lesion_label=="NC",]
    setD = maxclassf[sample(1:nrow(maxclassf), minclassn, replace=FALSE),]
    setD = rbind(setD, TrainsetD[TrainsetD$lesion_label=="C",])
    pander(summary(setD$lesion_label))
    
  }
  if(typeSampling=="regular"){
    setD = TrainsetD
    pander(summary(setD$lesion_label))
  }
  
  # train
  ###################################################
  # create grid of evaluation points
  gmaxdepth = c(3,5,15) 
  gmfinal = c(50,75,100,250,750)
  grd <- expand.grid(maxdepth = gmaxdepth, mfinal = gmfinal)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    maxdepth=grd$maxdepth[k]
    mfinal=grd$mfinal[k]
    
    # Build in 
    cat("maxdepth: ", maxdepth, "mfinal: ", mfinal, "\n")
    rpart_adaboost <- boosting(lesion_label ~ ., data = setD, 
                               mfinal = mfinal,
                               boos = TRUE,
                               coeflearn = 'Freund',
                               control = rpart.control(maxdepth = maxdepth))
    
    varImps = rpart_adaboost$imp[order(rpart_adaboost$imp, decreasing = TRUE)]
    #mostImp = varImps[varImps>1]
    #print(mostImp)
    
    #trainACU = print( sum(rpart_adaboost$class == TrainsetD$lesion_label) / length(TrainsetD$lesion_label) )
    ## in testing
    test_predboosting <- predict.boosting(rpart_adaboost, newdata = TestsetD) 
    #testAUC = print( sum(test_predboosting$class == TestsetD$lesion_label) / length(TestsetD$lesion_label) )
    
    ## to show the error evolution usefulness
    # evol.test <- errorevol(rpart_adaboost, TestsetD) 
    # evol.train <- errorevol(rpart_adaboost, TrainsetD) 
    # plot(evol.test$error, type = "l", ylim = c(0, 1), 
    #      main = "Boosting error versus number of trees", 
    #      xlab = "Iterations", 
    #      ylab = "Error", col = "red", lwd = 2) 
    # lines(evol.train$error, cex = .5, col = "blue", lty = 2, lwd = 2) 
    # legend("topright", c("test", "train"), col = c("red", "blue"), lty = 1:2, lwd = 2)
    # 
    # collect results
    # for TrainsetD
    trainperf <- test_roc(rpart_adaboost, setD)
    trainperf$labelstest = setD$lesion_label
    testperf <- test_roc(rpart_adaboost, TestsetD)
    testperf$labelstest = TestsetD$lesion_label
    
    # for ROC
    grdperf$rocTrain[k] <- trainperf$auc
    grdperf$rocTest[k] <- testperf$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(maxdepth=maxdepth,
                              mfinal=mfinal,
                              trainperf=trainperf,
                              testperf=testperf,
                              forest=rpart_adaboost,
                              varImps=varImps)))
  }
  
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestuned = M[index]$M
  return(bestuned)
}


stage1_boosting_Train_wcv_wperf <- function(TrainsetD, TestsetD, typeSampling){
  ## TrainsetD=stage1_data
  ## TestsetD=T1T2test
  ## typeSampling="oversmin"
  library(pROC)
  # to eval
  test_roc <- function(model, testdata) {
    # create call
    testpred <- predict(model, newdata=testdata, type="prob")
    roc_obj <- roc(testdata$lesion_label, testpred$prob[,1])
    output <- list(auc=roc_obj$auc, ci=ci(roc_obj), testpred=testpred$prob)
    return(output)
  }
  
  ## Sampling
  ######################
  if(typeSampling=="oversmin"){
    # add class balance: oversample minclassn classn to match majclassn
    majclassn = max(summary(TrainsetD$find_t2_signal_int))
    minclassf = TrainsetD[TrainsetD$find_t2_signal_int=="reported",]
    setD = minclassf[sample(1:nrow(minclassf), majclassn, replace=TRUE),]
    setD = rbind(setD, TrainsetD[TrainsetD$find_t2_signal_int=="notReported",])
    pander(summary(setD$find_t2_signal_int))
    setD = setD[,-c(1)]
  }
  if(typeSampling=="undersmaj"){
    # add class balance: downsample majority classn to match minclassn
    minclassn = min(summary(TrainsetD$find_t2_signal_int))
    maxclassf = TrainsetD[TrainsetD$find_t2_signal_int=="notReported",]
    setD = maxclassf[sample(1:nrow(maxclassf), minclassn, replace=FALSE),]
    setD = rbind(setD, TrainsetD[TrainsetD$find_t2_signal_int=="reported",])
    pander(summary(setD$find_t2_signal_int))
    setD = setD[,-c(1)]
  }
  if(typeSampling=="regular"){
    setD = TrainsetD
    pander(summary(setD$find_t2_signal_int))
    setD = setD[,-c(1)]
  }
  
  # train
  ###################################################
  # create grid of evaluation points
  gmaxdepth = c(1,3,5,15) 
  gmfinal = c(50,100,200)
  grd <- expand.grid(maxdepth = gmaxdepth, mfinal = gmfinal)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    maxdepth=grd$maxdepth[k]
    mfinal=grd$mfinal[k]
    
    # Build in 
    cat("maxdepth: ", maxdepth, "mfinal: ", mfinal, "\n")
    rpart_adaboost <- boosting(find_t2_signal_int ~ ., data = setD, 
                               mfinal = mfinal,
                               boos = FALSE,
                               control = rpart.control(maxdepth = maxdepth))
    
    varImps = rpart_adaboost$imp[order(rpart_adaboost$imp, decreasing = TRUE)]
    #mostImp = varImps[varImps>1]
    #print(mostImp)
    
    #trainACU = print( sum(rpart_adaboost$class == TrainsetD$lesion_label) / length(TrainsetD$lesion_label) )
    ## in testing
    test_predboosting <- predict.boosting(rpart_adaboost, newdata = TestsetD) 
    #testAUC = print( sum(test_predboosting$class == TestsetD$lesion_label) / length(TestsetD$lesion_label) )
    
    ## to show the error evolution usefulness
    # evol.test <- errorevol(rpart_adaboost, TestsetD) 
    # evol.train <- errorevol(rpart_adaboost, TrainsetD) 
    # plot(evol.test$error, type = "l", ylim = c(0, 1), 
    #      main = "Boosting error versus number of trees", 
    #      xlab = "Iterations", 
    #      ylab = "Error", col = "red", lwd = 2) 
    # lines(evol.train$error, cex = .5, col = "blue", lty = 2, lwd = 2) 
    # legend("topright", c("test", "train"), col = c("red", "blue"), lty = 1:2, lwd = 2)
    # 
    # collect results
    # for TrainsetD
    trainperf <- test_roc(rpart_adaboost, TrainsetD)
    trainperf$labelstest = setD$find_t2_signal_int
    testperf <- test_roc(rpart_adaboost, TestsetD)
    testperf$labelstest = TestsetD$find_t2_signal_int
    
    # for ROC
    grdperf$rocTrain[k] <- trainperf$auc
    grdperf$rocTest[k] <- testperf$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(maxdepth=maxdepth,
                              mfinal=mfinal,
                              trainperf=trainperf,
                              testperf=testperf,
                              forest=rpart_adaboost,
                              varImps=varImps)))
  }
  
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestuned = M[index]$M
  return(bestuned)
}


###################################################
# Build a nnet classifier
###################################################
nnet_Train_wcv_wperf <- function(TrainsetD, TestsetD){
  # TrainsetD=imgT2pLMSIR
  # TestsetD=T1T2test

  library(pROC)
  library(RSNNS)
  
  # datasets
  ###################################################
  # equal vars in train and test
  setD = na.omit(TrainsetD)
  TestsetD = TestsetD[,c(names(setD))]
  
  # split in vars/labels
  trainValues <- setD[,2:ncol(setD)]
  trainTargets <- decodeClassLabels(setD[,c(1)])
  
  setDtest = na.omit(TestsetD)
  testValues <- setDtest[,2:ncol(setD)]
  testTargets <- decodeClassLabels(setDtest[,c(1)])
  
  # format inputs to NNET
  data <- list(inputsTrain = trainValues,
               targetsTrain = trainTargets,
               inputsTest = testValues,
               targetsTest = testTargets)
  
  ## Data normalization:  The input matrix is column-wise normalized.
  ## norm: (def) the data is normalized to mean zero, variance one
  # to get normalized params: attr(dataInput$inputsTest, "normParams")
  dataInput <- normTrainingAndTestSet(data)
  
  # train
  ###################################################
  # create grid of evaluation points
  gsize1 = c(100,250,500) 
  gsize2 = c(12,24,36)
  gsize3 = c(10)
  grd <- expand.grid(size1 = gsize1, size2=gsize2, size3=gsize3)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    size1=grd$size1[k]
    size2=grd$size2[k]
    size3=grd$size3[k]
    
    # Build in 
    # cat("size1: ", size1, "size2: ", size2, "size3: ", size3)
    nnetmodel <- mlp(dataInput$inputsTrain, dataInput$targetsTrain, 
                     size=c(size1,size2,size3), learnFuncParams=c(0.1),
                     maxit=200, 
                     inputsTest=dataInput$inputsTest, 
                     targetsTest=dataInput$targetsTest)

    #par(mfrow=c(2,2))
    #plotIterativeError(nnetmodel)
    predictions <- predict(nnetmodel, dataInput$inputsTest)
    
    #plotRegressionError(predictions[,2], dataInput$targetsTest[,2])
    #confusionMatrix(dataInput$targetsTrain,fitted.values(nnetmodel))
    #confusionMatrix(dataInput$targetsTest,predictions)
    #plotROC(fitted.values(nnetmodel)[,2], dataInput$targetsTrain[,2])
    #plotROC(predictions[,2], dataInput$targetsTest[,2])

    # collect results
    # for TrainsetD
    roc_train <- roc(dataInput$targetsTrain[,2], fitted.values(nnetmodel)[,2])
    trainperf <- list(auc=roc_train$auc, ci=ci(roc_train), testpred=fitted.values(nnetmodel))
    trainperf$labelstest = setD$lesion_label
    
    roc_test <- roc(dataInput$targetsTest[,2], predictions[,2])
    testperf <- list(auc=roc_test$auc, ci=ci(roc_test), testpred=predictions)
    testperf$labelstest = setDtest$lesion_label
    
    # for ROC
    grdperf$rocTrain[k] <- trainperf$auc
    grdperf$rocTest[k] <- testperf$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(size1=size1,
                              size2=size2,
                              size3=size3,
                              trainperf=trainperf,
                              testperf=testperf,
                              nnet=nnetmodel)))
  }
  
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestuned = M[index]$M
  return(bestuned)
}

plotnnetImp <-function(out.var,mod.in,bar.plot=T,struct=NULL,x.lab=NULL,
                    y.lab=NULL, wts.only = F){
    
    require(ggplot2)
    
    # function works with neural networks from neuralnet, nnet, and RSNNS package
    # manual input vector of weights also okay
    
    #sanity checks
    if('numeric' %in% class(mod.in)){
      if(is.null(struct)) stop('Three-element vector required for struct')
      if(length(mod.in) != ((struct[1]*struct[2]+struct[2]*struct[3])+(struct[3]+struct[2])))
        stop('Incorrect length of weight matrix for given network structure')
      if(substr(out.var,1,1) != 'Y' | 
         class(as.numeric(gsub('^[A-Z]','', out.var))) != 'numeric')
        stop('out.var must be of form "Y1", "Y2", etc.')
    }
    if('train' %in% class(mod.in)){
      if('nnet' %in% class(mod.in$finalModel)){
        mod.in<-mod.in$finalModel
        warning('Using best nnet model from train output')
      }
      else stop('Only nnet method can be used with train object')
    }
    
    #gets weights for neural network, output is list
    #if rescaled argument is true, weights are returned but rescaled based on abs value
    nnet.vals<-function(mod.in,nid,rel.rsc,struct.out=struct){
      
      require(scales)
      require(reshape)
      
      if('numeric' %in% class(mod.in)){
        struct.out<-struct
        wts<-mod.in
      }
      
      #neuralnet package
      if('nn' %in% class(mod.in)){
        struct.out<-unlist(lapply(mod.in$weights[[1]],ncol))
        struct.out<-struct.out[-length(struct.out)]
        struct.out<-c(
          length(mod.in$model.list$variables),
          struct.out,
          length(mod.in$model.list$response)
        )      	
        wts<-unlist(mod.in$weights[[1]])   
      }
      
      #nnet package
      if('nnet' %in% class(mod.in)){
        struct.out<-mod.in$n
        wts<-mod.in$wts
      }
      
      #RSNNS package
      if('mlp' %in% class(mod.in)){
        struct.out<-c(mod.in$nInputs,mod.in$archParams$size,mod.in$nOutputs)
        hid.num<-length(struct.out)-2
        wts<-mod.in$snnsObject$getCompleteWeightMatrix()
        
        #get all input-hidden and hidden-hidden wts
        inps<-wts[grep('Input',row.names(wts)),grep('Hidden_2',colnames(wts)),drop=F]
        inps<-melt(rbind(rep(NA,ncol(inps)),inps))$value
        uni.hids<-paste0('Hidden_',1+seq(1,hid.num))
        for(i in 1:length(uni.hids)){
          if(is.na(uni.hids[i+1])) break
          tmp<-wts[grep(uni.hids[i],rownames(wts)),grep(uni.hids[i+1],colnames(wts)),drop=F]
          inps<-c(inps,melt(rbind(rep(NA,ncol(tmp)),tmp))$value)
        }
        
        #get connections from last hidden to output layers
        outs<-wts[grep(paste0('Hidden_',hid.num+1),row.names(wts)),grep('Output',colnames(wts)),drop=F]
        outs<-rbind(rep(NA,ncol(outs)),outs)
        
        #weight vector for all
        wts<-c(inps,melt(outs)$value)
        assign('bias',F,envir=environment(nnet.vals))
      }
      
      if(nid) wts<-rescale(abs(wts),c(1,rel.rsc))
      
      #convert wts to list with appropriate names 
      hid.struct<-struct.out[-c(length(struct.out))]
      row.nms<-NULL
      for(i in 1:length(hid.struct)){
        if(is.na(hid.struct[i+1])) break
        row.nms<-c(row.nms,rep(paste('hidden',i,seq(1:hid.struct[i+1])),each=1+hid.struct[i]))
      }
      row.nms<-c(
        row.nms,
        rep(paste('out',seq(1:struct.out[length(struct.out)])),each=1+struct.out[length(struct.out)-1])
      )
      out.ls<-data.frame(wts,row.nms)
      out.ls$row.nms<-factor(row.nms,levels=unique(row.nms),labels=unique(row.nms))
      out.ls<-split(out.ls$wts,f=out.ls$row.nms)
      
      assign('struct',struct.out,envir=environment(nnet.vals))
      
      out.ls
      
    }
    
    # get model weights
    best.wts<-nnet.vals(mod.in,nid=F,rel.rsc=5,struct.out=NULL)
    
    # weights only if T
    if(wts.only) return(best.wts)
    
    # get column index value for response variable to measure
    if('numeric' %in% class(mod.in)){
      out.ind <-  as.numeric(gsub('^[A-Z]','',out.var))
    } else {
      out.ind<-which(out.var==colnames(eval(mod.in$call$y)))
    }
    
    #get variable names from mod.in object
    #change to user input if supplied
    if('numeric' %in% class(mod.in)){
      x.names<-paste0(rep('X',struct[1]),seq(1:struct[1]))
      y.names<-paste0(rep('Y',struct[3]),seq(1:struct[3]))
    }
    if('mlp' %in% class(mod.in)){
      all.names<-mod.in$snnsObject$getUnitDefinitions()
      x.names<-all.names[grep('Input',all.names$unitName),'unitName']
      y.names<-all.names[grep('Output',all.names$unitName),'unitName']
    }
    if('nn' %in% class(mod.in)){
      x.names<-mod.in$model.list$variables
      y.names<-mod.in$model.list$response
    }
    if('xNames' %in% names(mod.in)){
      x.names<-mod.in$xNames
      y.names<-attr(terms(mod.in),'factor')
      y.names<-row.names(y.names)[!row.names(y.names) %in% x.names]
    }
    if(!'xNames' %in% names(mod.in) & 'nnet' %in% class(mod.in)){
      if(is.null(mod.in$call$formula)){
        x.names<-colnames(eval(mod.in$call$x))
        y.names<-colnames(eval(mod.in$call$y))
      }
      else{
        forms<-eval(mod.in$call$formula)
        x.names<-mod.in$coefnames
        facts<-attr(terms(mod.in),'factors')
        y.check<-mod.in$fitted
        if(ncol(y.check)>1) y.names<-colnames(y.check)
        else y.names<-as.character(forms)[2]
      } 
    }
    #change variables names to user sub 
    if(!is.null(x.lab)){
      if(length(x.names) != length(x.lab)) stop('x.lab length not equal to number of input variables')
      else x.names<-x.lab
    }
    if(!is.null(y.lab)){
      if(length(y.names) != length(y.lab)) stop('y.lab length not equal to number of output variables')
      else y.names<-y.lab
    }
    
    #get input-hidden weights and hidden-output weights, remove bias
    inp.hid<-data.frame(
      do.call('cbind',best.wts[grep('hidden',names(best.wts))])[-1,],
      row.names=x.names#c(colnames(eval(mod.in$call$x)))
    )
    hid.out<-best.wts[[grep(paste('out',out.ind),names(best.wts))]][-1]
    
    #multiply hidden-output connection for each input-hidden weight
    mult.dat<-data.frame(
      sapply(1:ncol(inp.hid),function(x) inp.hid[,x]*hid.out[x]),
      row.names=rownames(inp.hid)
    )    
    names(mult.dat)<-colnames(inp.hid)
    
    #get relative contribution of each input variable to each hidden node, sum values for each input
    #inp.cont<-rowSums(apply(mult.dat,2,function(x) abs(x)/sum(abs(x))))
    inp.cont<-rowSums(mult.dat)
    
    #get relative contribution
    #inp.cont/sum(inp.cont)
    
    rel.imp<-{
      signs<-sign(inp.cont)
      signs*rescale(abs(inp.cont),c(0,1))
    }
    
    if(!bar.plot){
      return(list(
        mult.wts=mult.dat,
        inp.cont=inp.cont,
        rel.imp=rel.imp
      ))
    }
    
    to_plo <- data.frame(rel.imp,x.names)[order(rel.imp),,drop = F]
    to_plo$x.names <- factor(x.names[order(rel.imp)], levels = x.names[order(rel.imp)])
    out_plo <- ggplot(to_plo, aes(x = x.names, y = rel.imp, fill = rel.imp,
                                  colour = rel.imp)) + 
      geom_bar(stat = 'identity') + 
      scale_x_discrete(element_blank()) +
      scale_y_continuous(y.names)
    
    return(out_plo)
    
  }

plot_nnet_relImpVars <- function(bestuned){
  #import 'gar.fun' from Github
  library(devtools)
  source_url('https://gist.githubusercontent.com/fawda123/6206737/raw/d6f365c283a8cae23fb20892dc223bc5764d50c7/gar_fun.r')#import 'gar.fun' from Github
  # The function is very simple to implement and has the following arguments:
  #   
  #   out.var	character string indicating name of response variable in the neural network object to be evaluated, only one input is allowed for models with multivariate response
  # mod.in	model object for input created from nnet function
  # bar.plot	logical value indicating if a figure is also created in the output, default T
  # x.names	character string indicating alternative names to be used for explanatory variables in the figure, default is taken from mod.in
  
  # The new arguments for gar.fun are as follows:
  #   
  #   out.var	character string indicating name of response variable in the neural network object to be evaluated, only one input is allowed for models with multivariate response, must be of form ‘Y1’, ‘Y2’, etc. if using numeric values as weight inputs for mod.in
  # mod.in	model object for neural network created using the nnet, RSNNS, or neuralnet packages, alternatively, a numeric vector specifying model weights in specific order, see example
  # bar.plot	logical value indicating if a figure or relative importance values are returned, default T
  # struct	numeric vector of length three indicating structure of the neural network, e.g., 2, 2, 1 for two inputs, two hidden, and one response, only required if mod.in is a vector input
  # x.lab	character string indicating alternative names to be used for explanatory variables in the figure, default is taken from mod.in
  # y.lab	character string indicating alternative names to be used for response variable in the figure, default is taken from out.var
  # wts.only	logical indicating of only model weights should be returned
  
  # The function returns a list with three elements, the most important of which is the last element named rel.imp. 
  # This element indicates the relative importance of each input variable for the named response variable as a value 
  # from -1 to 1. From these data, we can get an idea of what the neural network is telling us about the specific 
  # importance of each explanatory for the response variable. Here’s the function in action:
  
  nnetm = bestuned$nnet
  num.vars = nnetm$nInputs
  
  #create a pretty color vector for the bar plot
  cols<-colorRampPalette(c('lightgreen','lightblue'))(num.vars)
  
  #use the function on the model created above
  par(mar=c(3,4,1,1),family='serif')
  x.names=names(trainValues)
  
  vals.only <- gar.fun('lesion_label', nnetm, struct = c(nnetm$archParams$size), wts.only = F)
  print(vals.only)
  
  p1 <- gar.fun('C',vals.only, struct = c(nnetm$archParams$size) )
  p1
  
}