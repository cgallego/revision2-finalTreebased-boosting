# Rule extraction - experiments new paper



# Incorporating T2w breast MRI in CADx of breast lesions:
## Overall aim:
To investigate whether we can obtain an increase in discrimination ability for the classification task of cancerous and non-cancerous lesions using T2w imaging derived features in addition to more conventional T1w CE-MRI based-features

After running 10 folds of cross-validation (cv), below is the AUC distributions achieved on held-out test sets:

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7600  0.7899  0.8083  0.8088  0.8177  0.8879 
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7200  0.7871  0.8188  0.8182  0.8500  0.9010 
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7829  0.7913  0.8103  0.8247  0.8448  0.9394 
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## Performance of combined T1w+T2w vs. only T1w classifiers
### in all lesions (pooled data across cv-folds and plot pooled results)

```
[1] "Results for T1w-only features classifier:"
```

```
Area under the curve: 0.8014
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4706139 0.6041    0.6653  0.7224 0.7696    0.8089  0.8482
```

```
[1] "Results for T2w measured LMSIR features classifier:"
```

```
Area under the curve: 0.8149
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4040565 0.7469       0.8   0.849 0.6518    0.6963  0.7408
```

```
[1] "Results for T2w predicted LMSIR features classifier:"
```

```
Area under the curve: 0.828
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.409482 0.6776    0.7347  0.7918 0.7304    0.7723  0.8141
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### Significance t-test and Effect size:

Power analysis is an important aspect of experimental design. It allows us to determine the sample size required to detect an effect of a given size with a given degree of confidence. Conversely, it allows us to determine the probability of detecting an effect of a given size with a given level of confidence, under sample size constraints. If the probability is unacceptably low, we would be wise to alter or abandon the experiment.

The following four quantities have an intimate relationship:
1.  sample size
2.  effect size
3.  significance level = P(Type I error) = probability of finding an effect that is not there
4.  power = 1 - P(Type II error) = probability of finding an effect that is there

Given any three, we can determine the fourth.

Cohen (1988) hesitantly defined effect sizes as 
"small,d = .2,"
"medium,d = .5," and 
"large,d =.8", 
stating that "there is a certain risk in inherent in offering conventional operational definitions for those terms for use in power analysis in as diverse a field of inquiry as behavioral science" (p. 25).

Effect sizes can also be thought of as the average percentile standing of the average treated (or experimental) participant relative to the average untreated (or control) participant. An ES of 0.0 indicates that the mean of thetreated group is at the 50 thpercentile of the untreated group. An ES of 0.8 indicates that the mean of the treated group is at the 79th percentile of the untreated group. An effect size of 1.7 indicates that the mean of thetreated group is at the 95.5 percentile of the untreated group.


```

	Bootstrap test for two correlated ROC curves

data:  p2mLMSIR$ROC and p1all$ROC
D = 0.92937, boot.n = 2000, boot.stratified = 1, p-value = 0.1763
alternative hypothesis: true difference in AUC is greater than 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8149054   0.8013730 
```

```

	Bootstrap test for two correlated ROC curves

data:  p2pLMSIR$ROC and p1all$ROC
D = 2.0863, boot.n = 2000, boot.stratified = 1, p-value = 0.01848
alternative hypothesis: true difference in AUC is greater than 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8280158   0.8013730 
```

```

	Bootstrap test for two correlated ROC curves

data:  p2mLMSIR$ROC and p2pLMSIR$ROC
D = -1.0196, boot.n = 2000, boot.stratified = 1, p-value = 0.3079
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8149054   0.8280158 
```

```

     Two-sample t test power calculation 

              n = 245
              d = 0.1617076
      sig.level = 0.1805556
          power = 0.6741573
    alternative = two.sided

NOTE: n is number in *each* group
```
#### Aim Experiment 4: Rule extraction and likelihood ratios of T1w vs. T2w diagnostic predictions

```r
# read k-fold-out
allrulesres = c()

for (k in 1:10) {
    # 1:10f cv
    load(paste0("Z:/Cristina/Section2/finalTreebased-boosting/Outputs/T2SIvspredLMSIR_boost_addeddiagvalue_cv", 
        k, ".RData"))
    source("Z:/Cristina/Section2/finalTreebased-boosting/FunctionsRules.R")
    
    ################## Define datasets
    imgT1train = T1train[, -c(ncol(T1train))]
    imgT2train = T1T2train[, -c(ncol(T1T2train))]  #'T1w+T2wSI'
    imgT2mLMSIR = cbind(imgT2train[, -c(201, 223:242)], trainpredT2LMSIR[2])  #'T1w -T2wSI + measureLMSIR
    imgT2pLMSIR = cbind(imgT2train[, -c(201, 223:242)], trainpredT2LMSIR[1])  #'T1w+T2wSI_predictedLMSIR'
    
    # test set
    imgT1test = T1test[, -c(ncol(T1test))]
    T1T2test = allfT1T2[[2]][-c(201, 223:ncol(allfT1T2[[2]]))]
    imgT2mLMSIRtest = cbind(T1T2test, testpredT2LMSIR[2])
    imgT2pLMSIRtest = cbind(T1T2test, testpredT2LMSIR[1])
    classes = levels(imgT2pLMSIRtest[, "lesion_label"])
    
    # form feature dictionary
    fdict = feature_dictionary(imgT2pLMSIR)
    fnames = fdict$fnnames
    
    ############ Build a train trees check performance using 10-fold cross-validation
    cat("\n============ rule gbm trees imgT2pLMSIR \n")
    # train treees create grid of evaluation points
    gbmControl <- trainControl(method = "cv", number = 10, returnResamp = "final", savePredictions = "final", 
        summaryFunction = twoClassSummary, classProbs = TRUE)
    
    # for tunned parameters
    gbmGrid = expand.grid(.n.trees = c(1000, 1500, 2000), .interaction.depth = c(2, 3, 4, 6), .shrinkage = c(0.01), 
        .n.minobsinnode = c(10))
    
    # caret train cross-validated gbmModel, finalmodel has optimal parameters on training set
    imgT2pLMSIR = na.omit(imgT2pLMSIR)
    gbmModel <- train(lesion_label ~ ., data = imgT2pLMSIR, method = "gbm", distribution = "adaboost", 
        trControl = gbmControl, tuneGrid = gbmGrid, verbose = FALSE, metric = "ROC")
    
    # predictions <- predict(object=gbmModel$, X, type='prob') summarize the model # plot the effect of
    # parameters on accuracy
    print(gbmModel)
    plot(gbmModel)
    print(gbmModel$results)
    
    # transform rf object to an inTrees' format
    X <- imgT2pLMSIR[, 2:(ncol(imgT2pLMSIR))]
    target <- imgT2pLMSIR[, "lesion_label"]
    treeList <- myGBM2List(gbmModel$finalModel, X)
    
    # extract rules
    myruleExec = myextractRules(treeList, X, ntree = gbmModel$finalModel$n.trees)
    ruleExec = myruleExec$ruleExec
    rules = myruleExec$rules
    rulesProb = myruleExec$rulesProb
    
    # how good are these rules or what classification we obtain by applying them.  This can be done with
    # the handy function getRuleMetric. For the first 10 extracted rules:
    ruleMetric <- mygetRuleMetric(rules, rulesProb, X, target, gbmModel)  # get rule metrics
    ruleMetric[1:10, ]
    
    # a matrix including the conditions, prediction, and metrics, ordered by priority.
    alearner <- ruleMetric
    restest = c()
    allselRulesIx = list()
    for (i in 1:nrow(imgT2pLMSIRtest)) {
        X = imgT2pLMSIRtest[i, 2:ncol(imgT2pLMSIRtest)]
        y = imgT2pLMSIRtest[i, "lesion_label"]
        rulesoutput = myapplyLearner(alearner, X, y, minerr = 0.1, minfrq = 10/627, classes, gbmModel)
        
        # collect top scoring rules
        allselRulesIx[[i]] = rulesoutput[[1]]
        restest = rbind(restest, rulesoutput[[2]])
    }
    
    # compare
    confusionMatrix(perf_T2wpLMSIR$pred, perf_T2wpLMSIR$obs)
    dfrestest = as.data.frame(restest, stringsAsFactors = FALSE)
    dfrestest$predRules = ifelse(dfrestest$probModel >= 0.5, "C", "NC")
    confusionMatrix(dfrestest$predRules, restest[, "obsR"])
    
    allrulesres = rbind(allrulesres, dfrestest)
    pander(allrulesres)
    
    # plot iter plot ROCs each pass individually in l-o-p heldout test cases
    par(mfrow = c(1, 1))
    n = 15
    colors = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
    # plot 1/4
    p1 = calcAUC_plot(perf_imgT1$obs, perf_imgT1$C, xptext = 0.35, yptext = 0.85, 1, colors[2], atitle = "")
    par(new = TRUE)
    p2 = calcAUC_plot(perf_T2wmLMSIR$obs, perf_T2wmLMSIR$C, xptext = 0.45, yptext = 0.75, 2, colors[9], 
        atitle = "")
    par(new = TRUE)
    p3 = calcAUC_plot(perf_T2wpLMSIR$obs, perf_T2wpLMSIR$C, xptext = 0.55, yptext = 0.65, 3, colors[11], 
        atitle = "")
    
    par(new = TRUE)
    p4 = calcAUC_plot(perf_T2wpLMSIR$obs, as.numeric(restest[, "probRules"]), xptext = 0.65, yptext = 0.55, 
        4, colors[14], atitle = "")
    par(new = TRUE)
    p5 = calcAUC_plot(perf_T2wpLMSIR$obs, as.numeric(restest[, "probModel"]), xptext = 0.75, yptext = 0.45, 
        5, colors[15], atitle = paste0("ROCs 10f-patient out cv test k-fold= ", k))
    
    legend("bottomright", legend = c(paste0("T1w"), paste0("T2w+mLMSIR"), paste0("T1w+predLMSIR"), paste0("probRules"), 
        paste0("probModel")), col = c(colors[2], colors[9], colors[11], colors[14], colors[15]), lty = c(1, 
        2, 3, 4, 5), lwd = 2)
    
    # save current state k patient out
    save.image(paste0("Outputs/gmbRule_noT2SIpredLMSIR_boost_addeddiagvalue_cv", k, ".RData"))
}
```

```

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

549 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 494, 494, 494, 495, 494, 494, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8123492  0.5963203  0.8401070
  2                  1500     0.8143420  0.6054113  0.8370766
  2                  2000     0.8149590  0.6192641  0.8429590
  3                  1000     0.8146532  0.6008658  0.8521390
  3                  1500     0.8165826  0.6194805  0.8369875
  3                  2000     0.8197803  0.6333333  0.8369875
  4                  1000     0.8166144  0.5969697  0.8612299
  4                  1500     0.8211901  0.6426407  0.8430481
  4                  2000     0.8236230  0.6290043  0.8429590
  6                  1000     0.8215544  0.6103896  0.8581105
  6                  1500     0.8215125  0.6194805  0.8520499
  6                  2000     0.8230377  0.6238095  0.8551693

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 4, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8123492 0.5963203 0.8401070 0.08250904
4       0.01                 3             10    1000 0.8146532 0.6008658 0.8521390 0.08141301
7       0.01                 4             10    1000 0.8166144 0.5969697 0.8612299 0.08700230
10      0.01                 6             10    1000 0.8215544 0.6103896 0.8581105 0.08304086
2       0.01                 2             10    1500 0.8143420 0.6054113 0.8370766 0.07866763
5       0.01                 3             10    1500 0.8165826 0.6194805 0.8369875 0.07895563
8       0.01                 4             10    1500 0.8211901 0.6426407 0.8430481 0.08592654
11      0.01                 6             10    1500 0.8215125 0.6194805 0.8520499 0.08057237
3       0.01                 2             10    2000 0.8149590 0.6192641 0.8429590 0.07351966
6       0.01                 3             10    2000 0.8197803 0.6333333 0.8369875 0.07738703
9       0.01                 4             10    2000 0.8236230 0.6290043 0.8429590 0.07866589
12      0.01                 6             10    2000 0.8230377 0.6238095 0.8551693 0.07896270
       SensSD     SpecSD
1  0.10693921 0.06766866
4  0.10975921 0.05539319
7  0.09746297 0.06658237
10 0.10641096 0.07399926
2  0.11230432 0.07229550
5  0.11460194 0.06834386
8  0.12139429 0.04617496
11 0.10410533 0.05936421
3  0.10721243 0.06172839
6  0.11641990 0.05870705
9  0.11650036 0.05471884
12 0.11511408 0.05952938
10000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
9679   4 0.069 0.158
                                                                                                                             condition
9679 V12<=4.15872171672308 & dce3SE15>0.464838174729024 & lateSE17>0.858933269740971 & T2texture_correlation_nondir>0.0841295312637462
     pred                                       prob sumtempProb
9679   NC 0.6143961, 0.5155817, 0.5006583, 0.5675936   0.5495574
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4959   3 0.086 0.106
                                                                               condition pred
4959 irregularity<=0.909859616438443 & V13>5.24330958981904 & T2RGH_var>345.831042045511   NC
                                prob sumtempProb
4959 0.5741232, 0.6082200, 0.5286033   0.5703155
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 27 out of 112"
     len  freq   err
4815   3 0.056 0.097
                                                                                                condition
4815 texture_variance_nondir_post1>192.620387925557 & V5>3.51255119080561 & T2skew_F_r_i>1.22986842858167
     pred                            prob sumtempProb
4815    C 0.4011533, 0.4281161, 0.6413092   0.4901929
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
2014   3 0.264 0.152
                                                                                                    condition
2014 texture_variance_nondir_post1<=153.171430752191 & lateSE9>1.2281372126798 & T2_lesionSI>53.3142075787255
     pred                            prob sumtempProb
2014   NC 0.4799925, 0.3713933, 0.5005795   0.4506551
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 13 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 20 out of 112"
[1] "test complies with # rules: 15 out of 112"
     len  freq   err
8608   3 0.062 0.088
                                                                                   condition pred
8608 dce2SE14>0.432624096558318 & dce3SE13<=0.660588526384995 & T2_lesionSI>43.1284281668948   NC
                                prob sumtempProb
8608 0.4585595, 0.4432739, 0.4748801   0.4589045
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4241   3 0.155 0.118
                                                                                                             condition
4241 earlySE8<=0.956225206128341 & T2grad_margin_var<=1523.4375605976 & T2texture_sumvariance_nondir<=972.655393130357
     pred                            prob sumtempProb
4241   NC 0.3919588, 0.5853889, 0.5176618   0.4983365
[1] "test complies with # rules: 18 out of 112"
[1] "test complies with # rules: 19 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 18 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 14 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
1713   3 0.217 0.118
                                                                                           condition
1713 SER_inside<=0.6680364195439 & T2kurt_F_r_i<=4.46086165558844 & LMSIR_predicted>1.97039398616693
     pred                            prob sumtempProb
1713   NC 0.5304901, 0.4858707, 0.5115495   0.5093034
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
6859   4 0.066 0.139
                                                                                                                                                condition
6859 Kpeak_inside>-0.00173418493839903 & Vr_increasingRate_inside<=0.158266757432716 & V3>3.59029918667748 & T2texture_sumaverage_nondir>41.1799811180647
     pred                                       prob sumtempProb
6859   NC 0.4450344, 0.4982416, 0.4757835, 0.5864314   0.5013727
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 16 out of 112"
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 7 out of 112"
     len  freq   err
4815   3 0.056 0.097
                                                                                                condition
4815 texture_variance_nondir_post1>192.620387925557 & V5>3.51255119080561 & T2skew_F_r_i>1.22986842858167
     pred                            prob sumtempProb
4815    C 0.4011533, 0.4281161, 0.6413092   0.4901929
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
8259   3 0.097 0.132
                                                                                                         condition
8259 circularity<=0.855726514911507 & texture_sumaverage_nondir_post1>125.817629568624 & T2RGH_var>345.89795954425
     pred                            prob sumtempProb
8259   NC 0.6441597, 0.5996009, 0.4813968   0.5750525
[1] "test complies with # rules: 2000 out of 112"
     len  freq  err
1042   3 0.078 0.14
                                                                                                                               condition
1042 texture_diffvariance_nondir_post1<=124.098515348777 & texture_energy_nondir_post4<=0.00567725672571473 & T2RGH_var>706.065771712568
     pred                            prob sumtempProb
1042   NC 0.4661382, 0.4903879, 0.3816503   0.4460588
[1] "test complies with # rules: 7 out of 112"
     len  freq   err
5637   3 0.259 0.092
                                                                                        condition
5637 SER_inside<=0.684074475459813 & skew_F_r_i<=0.670444046675986 & T2_lesionSI>47.8757379930724
     pred                            prob sumtempProb
5637   NC 0.3961950, 0.4296275, 0.4587135   0.4281787
[1] "test complies with # rules: 6 out of 112"
    len  freq   err
884   4 0.168 0.098
                                                                                                                  condition
884 irregularity<=0.920362096693739 & V2>5.72871946106871 & V5>4.48442653681277 & T2texture_entropy_nondir>2.63686154430888
    pred                                       prob sumtempProb
884   NC 0.4926052, 0.4893595, 0.5431961, 0.5185153   0.5109191
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4241   3 0.155 0.118
                                                                                                             condition
4241 earlySE8<=0.956225206128341 & T2grad_margin_var<=1523.4375605976 & T2texture_sumvariance_nondir<=972.655393130357
     pred                            prob sumtempProb
4241   NC 0.3919588, 0.5853889, 0.5176618   0.4983365
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len freq   err
7014   3 0.12 0.167
                                                                                                   condition
7014 mean_F_r_i>760.47219110698 & dce3SE4>1.22621084356217 & T2texture_diffvariance_nondir<=519.719345448646
     pred                            prob sumtempProb
7014    C 0.5166803, 0.6039917, 0.6390194   0.5865638
[1] "test complies with # rules: 2000 out of 112"
     len freq   err
7014   3 0.12 0.167
                                                                                                   condition
7014 mean_F_r_i>760.47219110698 & dce3SE4>1.22621084356217 & T2texture_diffvariance_nondir<=519.719345448646
     pred                            prob sumtempProb
7014    C 0.5166803, 0.6039917, 0.6390194   0.5865638
[1] "test complies with # rules: 8 out of 112"
     len  freq   err
5637   3 0.259 0.092
                                                                                        condition
5637 SER_inside<=0.684074475459813 & skew_F_r_i<=0.670444046675986 & T2_lesionSI>47.8757379930724
     pred                            prob sumtempProb
5637   NC 0.3961950, 0.4296275, 0.4587135   0.4281787
[1] "test complies with # rules: 15 out of 112"
[1] "test complies with # rules: 9 out of 112"
     len freq   err
7086   3 0.04 0.091
                                                                                               condition
7086 max_RGH_var<=0.0830499307748255 & V9<=3.5187690383568 & T2texture_variance_nondir<=325.568685326362
     pred                            prob sumtempProb
7086   NC 0.5486425, 0.6723696, 0.4721203   0.5643775
[1] "test complies with # rules: 11 out of 112"
     len  freq   err
8608   3 0.062 0.088
                                                                                   condition pred
8608 dce2SE14>0.432624096558318 & dce3SE13<=0.660588526384995 & T2_lesionSI>43.1284281668948   NC
                                prob sumtempProb
8608 0.4585595, 0.4432739, 0.4748801   0.4589045
```

```
Area under the curve: 0.76
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4500341   0.48      0.68    0.84 0.6286    0.7714  0.9143
```

```
Area under the curve: 0.72
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3870426   0.64       0.8    0.96 0.5421    0.6857  0.8286
```

```
Area under the curve: 0.7829
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3974667  0.599      0.76    0.92    0.6    0.7429  0.8857
```

```
Area under the curve: 0.6629
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5340025    0.2      0.36    0.56 0.8286    0.9143       1
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```
Area under the curve: 0.7543
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.2795234   0.64       0.8    0.96 0.6286    0.7714  0.9143

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

540 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 485, 486, 486, 487, 485, 486, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8020456  0.5772727  0.8626894
  2                  1500     0.8023084  0.5727273  0.8410038
  2                  2000     0.8047172  0.5727273  0.8442235
  3                  1000     0.8035604  0.5913420  0.8533144
  3                  1500     0.8064800  0.5913420  0.8440341
  3                  2000     0.8084129  0.6056277  0.8501894
  4                  1000     0.8063980  0.5729437  0.8500000
  4                  1500     0.8062873  0.5961039  0.8531250
  4                  2000     0.8059434  0.6058442  0.8530303
  6                  1000     0.8136692  0.5865801  0.8531250
  6                  1500     0.8154233  0.5961039  0.8562500
  6                  2000     0.8158451  0.6101732  0.8592803

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 6, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8020456 0.5772727 0.8626894 0.05368740
4       0.01                 3             10    1000 0.8035604 0.5913420 0.8533144 0.05650952
7       0.01                 4             10    1000 0.8063980 0.5729437 0.8500000 0.05456490
10      0.01                 6             10    1000 0.8136692 0.5865801 0.8531250 0.05376567
2       0.01                 2             10    1500 0.8023084 0.5727273 0.8410038 0.05087710
5       0.01                 3             10    1500 0.8064800 0.5913420 0.8440341 0.05354417
8       0.01                 4             10    1500 0.8062873 0.5961039 0.8531250 0.05212182
11      0.01                 6             10    1500 0.8154233 0.5961039 0.8562500 0.05213514
3       0.01                 2             10    2000 0.8047172 0.5727273 0.8442235 0.04986438
6       0.01                 3             10    2000 0.8084129 0.6056277 0.8501894 0.05340471
9       0.01                 4             10    2000 0.8059434 0.6058442 0.8530303 0.05540713
12      0.01                 6             10    2000 0.8158451 0.6101732 0.8592803 0.05368424
       SensSD     SpecSD
1  0.07504336 0.05886635
4  0.06436983 0.05873325
7  0.07054528 0.04704067
10 0.07130822 0.04761384
2  0.08758555 0.04734733
5  0.08439276 0.05048936
8  0.06828920 0.04281413
11 0.07858200 0.04060281
3  0.08758555 0.04101167
6  0.06285592 0.04138064
9  0.07629132 0.04343187
12 0.05995363 0.04592602
14000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq   err
2266   3 0.174 0.138
                                                                                           condition
2266 SER_countor>0.626636929139987 & irregularity>0.968408057848163 & T2kurt_F_r_i<=4.66181507751613
     pred                            prob sumtempProb
2266    C 0.7980711, 0.5773883, 0.5266421   0.6340339
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 45 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 2000 out of 194"
    len  freq   err
834   3 0.131 0.113
                                                                                                      condition
834 iiMin_change_Variance_uptake<=0.232882294446878 & earlySE8<=0.887942203214057 & T2RGH_var<=348.454579378676
    pred                            prob sumtempProb
834   NC 0.4963651, 0.3699555, 0.5985809   0.4883005
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 10 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 27 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 8 out of 194"
     len  freq   err
7286   4 0.159 0.093
                                                                                                                                                                                 condition
7286 irregularity>0.925762714811578 & texture_variance_nondir_post1>188.263791397188 & texture_diffvariance_nondir_post4<=352.460877999769 & T2texture_diffentropy_nondir>1.24666771645525
     pred                                       prob sumtempProb
7286    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq err
6656   5 0.074 0.1
                                                                                                                                                        condition
6656 max_RGH_mean_k<=3.5 & texture_sumvariance_nondir_post1>415.025384577296 & dce2SE2>1.69971466829934 & dce2SE4>1.34623918866759 & T2_lesionSI>52.7200517318938
     pred                                                  prob sumtempProb
6656    C 0.4868222, 0.5391032, 0.6519079, 0.7495536, 0.5621248   0.5979023
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 21 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 2000 out of 194"
    len  freq   err
834   3 0.131 0.113
                                                                                                      condition
834 iiMin_change_Variance_uptake<=0.232882294446878 & earlySE8<=0.887942203214057 & T2RGH_var<=348.454579378676
    pred                            prob sumtempProb
834   NC 0.4963651, 0.3699555, 0.5985809   0.4883005
[1] "test complies with # rules: 21 out of 194"
     len  freq                err
1560   4 0.059 0.0620000000000001
                                                                                                                                  condition
1560 washoutRate_inside>0.00337769286510059 & max_RGH_var>0.0706830683129879 & T2min_F_r_i<=2.5 & T2texture_entropy_nondir>3.44135312373849
     pred                                       prob sumtempProb
1560    C 0.5070210, 0.6912999, 0.7587984, 0.6857119   0.6607078
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 23 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len freq   err
3455   4 0.03 0.125
                                                                                                                                                                condition
3455 irregularity>0.976454144711226 & texture_sumvariance_nondir_post3<=332.960432019186 & T2kurt_F_r_i<=4.66181507751613 & T2texture_diffentropy_nondir>1.24354568156965
     pred                                       prob sumtempProb
3455    C 0.5264760, 0.7506956, 0.5101984, 0.5449949   0.5830912
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 30 out of 194"
[1] "test complies with # rules: 26 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 10 out of 194"
      len  freq   err
13991   4 0.106 0.088
                                                                                                             condition
13991 mean_F_r_i>781.648597277964 & V19>4.23672233467332 & dce3SE5<=0.913188263359108 & T2skew_F_r_i<=1.38531917646111
      pred                                       prob sumtempProb
13991   NC 0.4503873, 0.5643563, 0.5306730, 0.2458475   0.4478161
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 16 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 32 out of 194"
[1] "test complies with # rules: 9 out of 194"
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 19 out of 194"
      len  freq                err
10069   3 0.063 0.0590000000000001
                                                                                                                               condition
10069 earlySE0<=0.384462022449602 & T2texture_inversediffmoment_nondir<=0.200298613838423 & T2texture_sumaverage_nondir>23.6389015999873
      pred                            prob sumtempProb
10069   NC 0.5055803, 0.4846939, 0.5071630   0.4991458
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 15 out of 194"
     len  freq   err
7286   4 0.159 0.093
                                                                                                                                                                                 condition
7286 irregularity>0.925762714811578 & texture_variance_nondir_post1>188.263791397188 & texture_diffvariance_nondir_post4<=352.460877999769 & T2texture_diffentropy_nondir>1.24666771645525
     pred                                       prob sumtempProb
7286    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq   err
8026   3 0.061 0.121
                                                                                              condition
8026 V12<=4.26731764525766 & dce3SE8>0.889522814281942 & T2texture_correlation_nondir>0.124732725643069
     pred                            prob sumtempProb
8026   NC 0.5140166, 0.4364616, 0.4869693   0.4791492
[1] "test complies with # rules: 6 out of 194"
     len  freq   err
6375   4 0.202 0.092
                                                                                                           condition
6375 irregularity<=0.925751210761359 & V5>3.55748428286012 & V10>5.29508712685388 & LMSIR_predicted>1.98623476838631
     pred                                       prob sumtempProb
6375   NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 7 out of 194"
```

```
Area under the curve: 0.7976
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5053268 0.4516    0.6129  0.7742  0.925     0.975       1
```

```
Area under the curve: 0.8855
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4771523 0.5806    0.7419   0.871  0.775     0.875   0.975
```

```
Area under the curve: 0.8508
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4292579 0.5161    0.6774  0.8387  0.775     0.875   0.975
```

```
Area under the curve: 0.8081
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5114321 0.5806    0.7419   0.871  0.625      0.75   0.875
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```
Area under the curve: 0.8548
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5067642 0.5161    0.6774  0.8387    0.8       0.9   0.975

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

553 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 497, 498, 499, 498, 498, 497, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8111291  0.5800395  0.8601326
  2                  1500     0.8133534  0.6025692  0.8448864
  2                  2000     0.8155170  0.6160079  0.8297348
  3                  1000     0.8132271  0.5936759  0.8540720
  3                  1500     0.8148337  0.6160079  0.8540720
  3                  2000     0.8182917  0.6160079  0.8448864
  4                  1000     0.8215404  0.5843874  0.8571023
  4                  1500     0.8225283  0.6023715  0.8632576
  4                  2000     0.8208642  0.6023715  0.8541667
  6                  1000     0.8205565  0.5843874  0.8661932
  6                  1500     0.8233598  0.5891304  0.8723485
  6                  2000     0.8235335  0.6069170  0.8632576

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 6, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8111291 0.5800395 0.8601326 0.06816381
4       0.01                 3             10    1000 0.8132271 0.5936759 0.8540720 0.06996264
7       0.01                 4             10    1000 0.8215404 0.5843874 0.8571023 0.06889478
10      0.01                 6             10    1000 0.8205565 0.5843874 0.8661932 0.07016599
2       0.01                 2             10    1500 0.8133534 0.6025692 0.8448864 0.07468625
5       0.01                 3             10    1500 0.8148337 0.6160079 0.8540720 0.07144251
8       0.01                 4             10    1500 0.8225283 0.6023715 0.8632576 0.07572169
11      0.01                 6             10    1500 0.8233598 0.5891304 0.8723485 0.07855032
3       0.01                 2             10    2000 0.8155170 0.6160079 0.8297348 0.07528202
6       0.01                 3             10    2000 0.8182917 0.6160079 0.8448864 0.07592558
9       0.01                 4             10    2000 0.8208642 0.6023715 0.8541667 0.07641870
12      0.01                 6             10    2000 0.8235335 0.6069170 0.8632576 0.07910714
      SensSD     SpecSD
1  0.1039143 0.07595729
4  0.1119282 0.05685754
7  0.1035704 0.08701874
10 0.1141162 0.07863847
2  0.1179623 0.07372294
5  0.1283116 0.06033988
8  0.1329191 0.08353301
11 0.1413913 0.06666978
3  0.1283116 0.08117445
6  0.1386315 0.06166511
9  0.1392367 0.07653744
12 0.1479545 0.06879749
14000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 17 out of 170"
     len  freq                err
2355   5 0.166 0.0649999999999999
                                                                                                                                                                                             condition
2355 maxVr_inside>0.132933970351532 & irregularity<=0.926301266037192 & texture_inversediffmoment_nondir_post4<=0.181096333105491 & T2max_F_r_i>237.5 & T2texture_sumvariance_nondir<=1892.48143275042
     pred                                                  prob sumtempProb
2355   NC 0.5671989, 0.5243535, 0.2820157, 0.4817899, 0.5310127   0.4772741
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 18 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
1023   4 0.146 0.136
                                                                                                                                            condition
1023 washoutRate_countor<=0.00237263537187332 & var_F_r_i<=60724.8179962257 & dce2SE4<=0.86655964324336 & T2texture_variance_nondir<=380.057119029906
     pred                                       prob sumtempProb
1023   NC 0.4613493, 0.5235084, 0.4794116, 0.3752399   0.4598773
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 10 out of 170"
      len  freq   err
10557   3 0.103 0.088
                                                                                     condition pred
10557 SER_countor<=0.713066433977772 & lateSE3<=1.00154434150194 & T2RGH_var<=348.899261925448   NC
                                 prob sumtempProb
10557 0.2830341, 0.4749445, 0.5770729   0.4450172
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 20 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 7 out of 170"
     len  freq   err
7487   5 0.128 0.099
                                                                                                                                                                                  condition
7487 iAUC1_inside<=1867.20963919127 & max_F_r_i<=1997 & texture_variance_nondir_post1>210.870898937361 & texture_diffvariance_nondir_post2<=324.174200101235 & T2_lesionSI>51.1411574636996
     pred                                                  prob sumtempProb
7487    C 0.4257563, 0.4854569, 0.4236016, 0.4099355, 0.4474628   0.4384426
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 19 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 9 out of 170"
[1] "test complies with # rules: 11 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
9620   4 0.069 0.105
                                                                                                                                               condition
9620 texture_sumaverage_nondir_post4<=168.539244885183 & V11<=51.2306200640228 & earlySE8<=0.974715284142456 & T2texture_entropy_nondir>2.47386447490301
     pred                                       prob sumtempProb
9620   NC 0.5130526, 0.4818617, 0.5684597, 0.5257756   0.5222874
[1] "test complies with # rules: 21 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
1812   4 0.146 0.111
                                                                                                                               condition
1812 texture_variance_nondir_post1>162.468635587987 & lateSE4>0.90412425359233 & lateSE19>1.04549581201775 & T2RGH_var<=605.680930822889
     pred                                       prob sumtempProb
1812    C 0.5779666, 0.4707040, 0.6166417, 0.5418299   0.5517855
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
5968   5 0.354 0.168
                                                                                                                                                                                             condition
5968 var_F_r_i<=73622.3679819383 & texture_sumvariance_nondir_post1<=262.013409085089 & texture_energy_nondir_post4<=0.00584870387417049 & dce3SE17>0.522817659910163 & T2kurt_F_r_i>-0.55567242234486
     pred                                                  prob sumtempProb
5968   NC 0.4572754, 0.3762719, 0.4043077, 0.4322990, 0.4773361    0.429498
[1] "test complies with # rules: 26 out of 170"
[1] "test complies with # rules: 24 out of 170"
      len  freq   err
13111   5 0.074 0.073
                                                                                                                                                                            condition
13111 texture_sumentropy_nondir_post1>1.94006995925894 & dce2SE7>0.524160803062216 & lateSE5>0.945952590906062 & lateSE8>1.13819699152219 & T2texture_entropy_nondir>2.33506839869461
      pred                                                  prob sumtempProb
13111    C 0.4301397, 0.4668065, 0.4695489, 0.5370841, 0.4446505    0.469646
[1] "test complies with # rules: 17 out of 170"
     len  freq   err
7487   5 0.128 0.099
                                                                                                                                                                                  condition
7487 iAUC1_inside<=1867.20963919127 & max_F_r_i<=1997 & texture_variance_nondir_post1>210.870898937361 & texture_diffvariance_nondir_post2<=324.174200101235 & T2_lesionSI>51.1411574636996
     pred                                                  prob sumtempProb
7487    C 0.4257563, 0.4854569, 0.4236016, 0.4099355, 0.4474628   0.4384426
[1] "test complies with # rules: 12 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq err                                                                     condition
7061   3 0.036 0.1 A_countor>0.878898973669725 & max_F_r_i>2132 & T2kurt_F_r_i<=5.44474459334273
     pred                            prob sumtempProb
7061    C 0.5297514, 0.5577010, 0.5004046   0.5292857
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 24 out of 170"
[1] "test complies with # rules: 9 out of 170"
```

```
Area under the curve: 0.8121
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4912968 0.4737    0.6842  0.8947 0.6944    0.8333  0.9444
```

```
Area under the curve: 0.8129
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5665428 0.3684    0.5789  0.7895 0.8611    0.9444       1
```

```
Area under the curve: 0.8026
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3918007 0.5789    0.7895  0.9474 0.6667    0.8056  0.9167
```

```
Area under the curve: 0.7164
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5252716 0.3684    0.5789  0.7895 0.6667    0.8056  0.9167
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```
Area under the curve: 0.7734
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6399194 0.3158    0.5263  0.7368 0.8056    0.9167       1

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

555 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 500, 499, 500, 500, 499, 499, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8168825  0.5833333  0.8691622
  2                  1500     0.8189087  0.6108225  0.8602496
  2                  2000     0.8179746  0.6203463  0.8573975
  3                  1000     0.8185296  0.6017316  0.8691622
  3                  1500     0.8170660  0.6019481  0.8631016
  3                  2000     0.8184000  0.6067100  0.8661319
  4                  1000     0.8224084  0.5969697  0.8691622
  4                  1500     0.8227794  0.6203463  0.8749554
  4                  2000     0.8228176  0.6155844  0.8690731
  6                  1000     0.8229437  0.6108225  0.8721925
  6                  1500     0.8226810  0.6153680  0.8661319
  6                  2000     0.8239397  0.6292208  0.8690731

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 6, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8168825 0.5833333 0.8691622 0.03598241
4       0.01                 3             10    1000 0.8185296 0.6017316 0.8691622 0.03326726
7       0.01                 4             10    1000 0.8224084 0.5969697 0.8691622 0.03445598
10      0.01                 6             10    1000 0.8229437 0.6108225 0.8721925 0.03165785
2       0.01                 2             10    1500 0.8189087 0.6108225 0.8602496 0.03392892
5       0.01                 3             10    1500 0.8170660 0.6019481 0.8631016 0.03539762
8       0.01                 4             10    1500 0.8227794 0.6203463 0.8749554 0.03263147
11      0.01                 6             10    1500 0.8226810 0.6153680 0.8661319 0.03095617
3       0.01                 2             10    2000 0.8179746 0.6203463 0.8573975 0.03131150
6       0.01                 3             10    2000 0.8184000 0.6067100 0.8661319 0.03676419
9       0.01                 4             10    2000 0.8228176 0.6155844 0.8690731 0.03457480
12      0.01                 6             10    2000 0.8239397 0.6292208 0.8690731 0.03089254
      SensSD     SpecSD
1  0.1210542 0.05975956
4  0.1131492 0.05450550
7  0.1276090 0.06623782
10 0.1177696 0.05592403
2  0.1271431 0.06236142
5  0.1029471 0.06613434
8  0.1125458 0.04422939
11 0.1253463 0.04696029
3  0.1185072 0.07615324
6  0.1184934 0.05795405
9  0.1256829 0.05627370
12 0.1227656 0.04896727
14000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
3745   4 0.031 0.118
                                                                                                                                                                                       condition
3745 Vr_increasingRate_inside>0.119387287517364 & texture_variance_nondir_post1>158.282915050042 & texture_inversediffmoment_nondir_post2>0.144212576138197 & T2grad_margin_var>3296.95677657013
     pred                                       prob sumtempProb
3745    C 0.4794050, 0.3560235, 0.3728296, 0.5111716   0.4298574
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 16 out of 187"
[1] "test complies with # rules: 17 out of 187"
     len  freq   err
3603   5 0.124 0.087
                                                                                                                                condition
3603 irregularity<=0.92495165286538 & V0>7.25458872783576 & V6>5.44452206416059 & lateSE6>0.535319521686134 & T2RGH_mean>21.8641883788852
     pred                                                  prob sumtempProb
3603   NC 0.5119074, 0.4878885, 0.5521666, 0.4338532, 0.4158351   0.4803301
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 15 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err
12196   3 0.124 0.101
                                                                                                                 condition
12196 SER_inside<=0.656166577515658 & texture_sumaverage_nondir_post4<=207.067310504744 & LMSIR_predicted>1.72917461715247
      pred                            prob sumtempProb
12196   NC 0.4866332, 0.4729611, 0.5208541   0.4934828
[1] "test complies with # rules: 10 out of 187"
     len  freq   err
3603   5 0.124 0.087
                                                                                                                                condition
3603 irregularity<=0.92495165286538 & V0>7.25458872783576 & V6>5.44452206416059 & lateSE6>0.535319521686134 & T2RGH_mean>21.8641883788852
     pred                                                  prob sumtempProb
3603   NC 0.5119074, 0.4878885, 0.5521666, 0.4338532, 0.4158351   0.4803301
[1] "test complies with # rules: 21 out of 187"
     len  freq                err
3190   4 0.029 0.0620000000000001
                                                                                                                                                          condition
3190 skew_F_r_i>-0.0687768048959314 & irregularity>0.925751210761359 & texture_inversediffmoment_nondir_post4<=0.0923255421666012 & T2kurt_F_r_i>-0.349415185608738
     pred                                       prob sumtempProb
3190   NC 0.5586323, 0.4173910, 0.4476405, 0.4538664   0.4693826
[1] "test complies with # rules: 24 out of 187"
[1] "test complies with # rules: 25 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 33 out of 187"
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len freq   err
6273   5 0.05 0.107
                                                                                                                                   condition
6273 A_inside<=1.61072218817226 & circularity<=0.855734481103919 & dce2SE11<=0.904584163031592 & lateSE0<=1.15150229077305 & T2min_F_r_i<=23
     pred                                                  prob sumtempProb
6273   NC 0.4461920, 0.6252745, 0.7221661, 0.6488880, 0.5635662   0.6012174
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 7 out of 187"
     len  freq   err
7128   5 0.137 0.092
                                                                                                                                                                                                      condition
7128 UptakeRate_inside<=1.00418756143313 & texture_energy_nondir_post1<=0.0130423566941265 & texture_variance_nondir_post1<=53.5902830040147 & earlySE10>0.34987914202482 & T2grad_margin_var<=6402.45036359788
     pred                                                  prob sumtempProb
7128   NC 0.4445401, 0.4908341, 0.4556805, 0.5115279, 0.4697661   0.4744697
[1] "test complies with # rules: 10 out of 187"
      len  freq   err
12569   5 0.061 0.088
                                                                                                                                                                                                                         condition
12569 var_F_r_i<=71220.2954156667 & texture_inversediffmoment_nondir_post2>0.200974618005213 & earlySE8<=0.86237226360737 & T2texture_inversediffmoment_nondir>0.0996454758962476 & T2texture_diffentropy_nondir<=1.52849412755374
      pred                                                  prob sumtempProb
12569   NC 0.4320504, 0.3620002, 0.5822682, 0.4883600, 0.4268523   0.4583062
[1] "test complies with # rules: 26 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
7574   4 0.079 0.114
                                                                                                                   condition
7574 SER_inside>0.742111684078729 & Vr_post_1_inside>0.171768754142706 & V13>25.2275266554215 & T2_lesionSI>69.3437696866711
     pred                                       prob sumtempProb
7574    C 0.5765578, 0.4723026, 0.4088882, 0.6560729   0.5284554
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 12 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 12 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 21 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 7 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq  err
2610   4 0.045 0.12
                                                                                                                                                                               condition
2610 circularity>0.855735583523288 & texture_sumaverage_nondir_post1>128.049051372998 & texture_energy_nondir_post4<=0.00765036336922699 & T2texture_diffentropy_nondir>1.75785784833312
     pred                                       prob sumtempProb
2610    C 0.5030484, 0.4860747, 0.5198417, 0.6825833    0.547887
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err
12138   5 0.047 0.115
                                                                                                                                               condition
12138 A_inside>2.41863937468322 & SER_inside>0.717590940320822 & irregularity>0.879980421830873 & lateSE19>1.10405688473339 & T2RGH_var>354.155127809895
      pred                                                  prob sumtempProb
12138    C 0.5931262, 0.5791689, 0.5281811, 0.7094207, 0.4091477   0.5638089
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
1427   4 0.032 0.111
                                                                                                                                                                              condition
1427 irregularity<=0.984668990886858 & texture_variance_nondir_post1>191.786000082219 & T2texture_correlation_nondir>0.263972379118135 & T2texture_correlation_nondir>0.361597155528048
     pred                                       prob sumtempProb
1427    C 0.4937985, 0.4742500, 0.4468187, 0.5846211   0.4998721
[1] "test complies with # rules: 19 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err                                                      condition pred
13176   2 0.032 0.111 max_RGH_mean<=0.495991412376266 & T2_lesionSI>113.683860342556   NC
                      prob sumtempProb
13176 0.5332697, 0.3199769   0.4266233
[1] "test complies with # rules: 9 out of 187"
[1] "test complies with # rules: 2000 out of 187"
```

```
Area under the curve: 0.8167
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4064521  0.625    0.7917  0.9583 0.5333       0.7  0.8667
```

```
Area under the curve: 0.85
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3199375  0.875    0.9583       1    0.5    0.6667  0.8333
```

```
Area under the curve: 0.8181
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5001467 0.3333    0.5417    0.75 0.8333    0.9333  1.0000
  0.3260255 0.7500    0.8750    1.00 0.4333    0.6000  0.7667
```

```
Area under the curve: 0.5861
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5068191  0.375    0.5833    0.75    0.5    0.6667  0.8333
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-4.png)<!-- -->

```
Area under the curve: 0.8236
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.1913509    0.5    0.7083   0.875 0.6333       0.8  0.9333

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

550 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 495, 494, 495, 496, 494, 495, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8117956  0.5660173  0.8679144
  2                  1500     0.8082400  0.5889610  0.8559715
  2                  2000     0.8053979  0.6073593  0.8588235
  3                  1000     0.8153983  0.5935065  0.8588235
  3                  1500     0.8116854  0.6121212  0.8498217
  3                  2000     0.8114499  0.6025974  0.8590018
  4                  1000     0.8128471  0.5935065  0.8557041
  4                  1500     0.8114896  0.6073593  0.8617647
  4                  2000     0.8094669  0.6162338  0.8468806
  6                  1000     0.8138885  0.5841991  0.8589127
  6                  1500     0.8101564  0.5935065  0.8589127
  6                  2000     0.8090940  0.6073593  0.8588235

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 1000, interaction.depth = 3, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8117956 0.5660173 0.8679144 0.06865282
4       0.01                 3             10    1000 0.8153983 0.5935065 0.8588235 0.06969911
7       0.01                 4             10    1000 0.8128471 0.5935065 0.8557041 0.07246357
10      0.01                 6             10    1000 0.8138885 0.5841991 0.8589127 0.07933798
2       0.01                 2             10    1500 0.8082400 0.5889610 0.8559715 0.07192912
5       0.01                 3             10    1500 0.8116854 0.6121212 0.8498217 0.07457000
8       0.01                 4             10    1500 0.8114896 0.6073593 0.8617647 0.07440672
11      0.01                 6             10    1500 0.8101564 0.5935065 0.8589127 0.08120705
3       0.01                 2             10    2000 0.8053979 0.6073593 0.8588235 0.07501765
6       0.01                 3             10    2000 0.8114499 0.6025974 0.8590018 0.07849322
9       0.01                 4             10    2000 0.8094669 0.6162338 0.8468806 0.07555994
12      0.01                 6             10    2000 0.8090940 0.6073593 0.8588235 0.07782917
      SensSD     SpecSD
1  0.1105581 0.05549245
4  0.1131804 0.06426348
7  0.1111336 0.07451233
10 0.1227037 0.06089525
2  0.1047963 0.06005040
5  0.1181756 0.06738902
8  0.1246640 0.06472774
11 0.1155789 0.06089525
3  0.1246640 0.06426348
6  0.1009493 0.04785395
9  0.1110598 0.06936467
12 0.1077005 0.06274999
4000 rules (length<=6) were extracted from the first 1000 trees.
[1] "test complies with # rules: 1000 out of 25"
     len  freq  err
3815   3 0.195 0.14
                                                                                                                           condition
3815 irregularity<=0.929432771288307 & texture_sumaverage_nondir_post4>170.690447525835 & T2texture_variance_nondir>134.392530089492
     pred                            prob sumtempProb
3815   NC 0.5956767, 0.5393015, 0.4608842   0.5319541
[1] "test complies with # rules: 1000 out of 25"
     len  freq  err
3815   3 0.195 0.14
                                                                                                                           condition
3815 irregularity<=0.929432771288307 & texture_sumaverage_nondir_post4>170.690447525835 & T2texture_variance_nondir>134.392530089492
     pred                            prob sumtempProb
3815   NC 0.5956767, 0.5393015, 0.4608842   0.5319541
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                               condition pred
2967   2 0.049 0.185 earlySE3<=1.43562762520037 & T2RGH_var>915.990383661592   NC
                     prob sumtempProb
2967 0.5473018, 0.4849801   0.5161409
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                    condition pred
2043   2 0.202 0.189 Tpeak_countor>12.0617409688816 & T2RGH_var<=370.226486471733   NC
                     prob sumtempProb
2043 0.3730531, 0.5089100   0.4409815
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                        condition pred
1560   2 0.093 0.157 irregularity>0.984788961473386 & T2kurt_F_r_i>-0.785487357074651    C
                     prob sumtempProb
1560 0.4785715, 0.5020898   0.4903307
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
2283   3 0.058 0.156
                                                                                            condition
2283 Vr_post_1_inside<=1.9595435001623 & earlySE8>1.84525127427783 & LMSIR_predicted>1.95958273598223
     pred                            prob sumtempProb
2283    C 0.4964885, 0.5170830, 0.5275262   0.5136992
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err
486   3 0.391 0.205
                                                                                                              condition
486 washoutRate_inside<=0.0290387586166104 & var_F_r_i<=25298.2264620836 & T2texture_sumaverage_nondir>24.8365868891366
    pred                            prob sumtempProb
486   NC 0.4555950, 0.3438633, 0.3683227   0.3892603
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err
933   3 0.176 0.175
                                                                                                 condition
933 Vr_decreasingRate_inside<=0.0123396939400547 & dce2SE7<=1.32152828351672 & T2RGH_var<=354.155127809895
    pred                            prob sumtempProb
933   NC 0.4807406, 0.4020555, 0.5836313   0.4888091
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len freq   err
1210   3 0.06 0.152
                                                                                                          condition
1210 texture_sumaverage_nondir_post1>178.317115217582 & V19<=4.69302425134638 & T2grad_margin_var<=5974.84492572428
     pred                            prob sumtempProb
1210   NC 0.6718742, 0.5561653, 0.5057248   0.5779214
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 978 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
2734   3 0.175 0.115
                                                                                                            condition
2734 iiMin_change_Variance_uptake>0.318089063969771 & T2kurt_F_r_i<=4.55628566003179 & T2kurt_F_r_i<=2.86325409387345
     pred                            prob sumtempProb
2734   NC 0.4691266, 0.4615858, 0.5104885   0.4804003
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                        condition pred
1560   2 0.093 0.157 irregularity>0.984788961473386 & T2kurt_F_r_i>-0.785487357074651    C
                     prob sumtempProb
1560 0.4785715, 0.5020898   0.4903307
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
2366   2 0.229 0.143
                                                                           condition pred
2366 irregularity<=0.928021468790062 & T2texture_diffentropy_nondir>1.45569379526772   NC
                     prob sumtempProb
2366 0.5544469, 0.4850301   0.5197385
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1995   3 0.031 0.176
                                                                                                          condition
1995 Vr_increasingRate_countor>0.0473468061227342 & max_RGH_mean<=0.535957127734903 & T2kurt_F_r_i>4.79798690590228
     pred                            prob sumtempProb
1995   NC 0.5666348, 0.4564168, 0.7178750   0.5803089
[1] "test complies with # rules: 5 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 7 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                       condition pred
1827   2 0.169 0.183 irregularity>0.977785790197412 & T2kurt_F_r_i<=4.65160033948678    C
                     prob sumtempProb
1827 0.5495094, 0.5000662   0.5247878
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                                    condition
1605   2 0.462 0.232 SER_countor<=0.705775175050863 & T2texture_variance_nondir<=648.653091721668
     pred                 prob sumtempProb
1605   NC 0.4425944, 0.4801895    0.461392
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                    condition pred
2043   2 0.202 0.189 Tpeak_countor>12.0617409688816 & T2RGH_var<=370.226486471733   NC
                     prob sumtempProb
2043 0.3730531, 0.5089100   0.4409815
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
3718   3 0.135 0.189
                                                                                                                           condition
3718 Vr_increasingRate_countor<=0.113820610943595 & texture_sumaverage_nondir_post2<=223.35047590633 & T2skew_F_r_i>1.21484263542727
     pred                            prob sumtempProb
3718   NC 0.5761680, 0.6450753, 0.5418709   0.5877047
[1] "test complies with # rules: 1000 out of 25"
```

```
Area under the curve: 0.8044
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4215601 0.7037    0.8519   0.963 0.4062    0.5938    0.75
```

```
Area under the curve: 0.7616
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4408423 0.4444    0.6296  0.8148  0.625    0.7812  0.9062
```

```
Area under the curve: 0.7894
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4635855 0.4815    0.6667  0.8148   0.75     0.875  0.9688
```

```
Area under the curve: 0.7755
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4964488 0.7778    0.8889       1 0.5312    0.6875  0.8438
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-5.png)<!-- -->

```
Area under the curve: 0.7512
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3403943 0.4444    0.6296  0.8148  0.625    0.7812  0.9062

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

546 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 491, 491, 492, 490, 492, 492, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8203492  0.6003953  0.8543561
  2                  1500     0.8224179  0.6233202  0.8421402
  2                  2000     0.8258472  0.6187747  0.8357955
  3                  1000     0.8272177  0.6007905  0.8420455
  3                  1500     0.8284405  0.6053360  0.8482955
  3                  2000     0.8292103  0.6189723  0.8420455
  4                  1000     0.8251213  0.6142292  0.8606061
  4                  1500     0.8254535  0.6096838  0.8482008
  4                  2000     0.8248304  0.6187747  0.8389205
  6                  1000     0.8248392  0.5871542  0.8420455
  6                  1500     0.8240648  0.6142292  0.8327652
  6                  2000     0.8262342  0.6144269  0.8357008

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 3, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8203492 0.6003953 0.8543561 0.07031925
4       0.01                 3             10    1000 0.8272177 0.6007905 0.8420455 0.06030782
7       0.01                 4             10    1000 0.8251213 0.6142292 0.8606061 0.06184086
10      0.01                 6             10    1000 0.8248392 0.5871542 0.8420455 0.06727245
2       0.01                 2             10    1500 0.8224179 0.6233202 0.8421402 0.06574981
5       0.01                 3             10    1500 0.8284405 0.6053360 0.8482955 0.06352300
8       0.01                 4             10    1500 0.8254535 0.6096838 0.8482008 0.06436847
11      0.01                 6             10    1500 0.8240648 0.6142292 0.8327652 0.06457492
3       0.01                 2             10    2000 0.8258472 0.6187747 0.8357955 0.06586353
6       0.01                 3             10    2000 0.8292103 0.6189723 0.8420455 0.06281742
9       0.01                 4             10    2000 0.8248304 0.6187747 0.8389205 0.05937237
12      0.01                 6             10    2000 0.8262342 0.6144269 0.8357008 0.06140174
       SensSD     SpecSD
1  0.09172554 0.08311759
4  0.08808211 0.07117030
7  0.10159376 0.07943364
10 0.10296276 0.08602646
2  0.11098597 0.07071982
5  0.11055443 0.07547533
8  0.09202052 0.07281349
11 0.10725398 0.08060833
3  0.10924106 0.08206245
6  0.12067677 0.07267892
9  0.11928660 0.07195094
12 0.09287112 0.07968594
8000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                                 condition pred
4163   3 0.084 0.109 V5>3.51255119080561 & V11<=5.03787177277192 & T2_lesionSI>52.876161736488   NC
                                prob sumtempProb
4163 0.4877236, 0.4626612, 0.5053751   0.4852533
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
4779   3 0.092 0.12
                                                                                                      condition
4779 Tpeak_inside>11.755177337987 & mean_F_r_i<=776.259473090249 & T2texture_sumaverage_nondir>39.9589002202509
     pred                            prob sumtempProb
4779   NC 0.4713799, 0.4513984, 0.5841410   0.5023064
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
1959   3 0.073 0.15
                                                                                                                                                       condition
1959 iiMin_change_Variance_uptake>0.331545153209991 & texture_inversediffmoment_nondir_post2>0.179893974878372 & T2texture_correlation_nondir<=0.370601717909582
     pred                            prob sumtempProb
1959   NC 0.6521837, 0.4705059, 0.5648111   0.5625002
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
5631   3 0.125 0.162
                                                                                                condition
5631 min_F_r_i>386.5 & texture_sumvariance_nondir_post2<=1048.75233234583 & T2kurt_F_r_i>1.14883296665627
     pred                            prob sumtempProb
5631   NC 0.5786830, 0.4887894, 0.4385449   0.5020058
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
    len  freq   err
999   3 0.057 0.129
                                                                                 condition pred
999 earlySE10>0.398275630549396 & dce3SE11<=0.736907120014966 & T2RGH_var>289.255846492511   NC
                               prob sumtempProb
999 0.5753032, 0.5538037, 0.4956976   0.5416015
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
4215   3 0.205 0.143
                                                                                   condition pred
4215 irregularity<=0.907737082350515 & V4>2.45154174506599 & T2kurt_F_r_i>-0.521339102394693   NC
                                prob sumtempProb
4215 0.4551003, 0.5220241, 0.5059258   0.4943501
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
3746   3 0.178 0.175
                                                                                    condition pred
3746 mean_F_r_i<=1117.32137367602 & lateSE5<=0.952831756279501 & T2_lesionSI>52.9125674273991   NC
                                prob sumtempProb
3746 0.4856048, 0.4487932, 0.5044204   0.4796061
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
4735   3 0.027 0.133
                                                                                           condition
4735 max_RGH_mean>0.531482089317235 & earlySE12>0.86558442166171 & LMSIR_predicted<=1.96220151440786
     pred                            prob sumtempProb
4735    C 0.4969997, 0.4322224, 0.5708283   0.5000168
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
6925   3 0.082 0.133
                                                                                                condition
6925 max_F_r_i<=833.625 & texture_energy_nondir_post4<=0.00231572627771039 & T2RGH_mean<=45.7207332374131
     pred                            prob sumtempProb
6925   NC 0.4568386, 0.5386332, 0.4971891   0.4975536
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1060   3 0.114 0.177
                                                                             condition pred
1060 max_F_r_i>1454.125 & dce2SE9>0.644619457119457 & LMSIR_predicted>2.91344151959955    C
                                prob sumtempProb
1060 0.4692342, 0.6110102, 0.7437256     0.60799
[1] "test complies with # rules: 6 out of 43"
     len  freq   err
6322   3 0.038 0.095
                                                                                          condition
6322 Kpeak_countor>-0.00631091014572267 & dce2SE6<=0.621422517899093 & T2RGH_mean<=28.5072354996839
     pred                            prob sumtempProb
6322   NC 0.4934857, 0.3293870, 0.5828066   0.4685598
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
6488   3 0.081 0.205
                                                                                                               condition
6488 SER_inside>0.824019225011374 & Vr_post_1_countor>0.225922745750899 & T2texture_correlation_nondir>0.220883358585512
     pred                            prob sumtempProb
6488    C 0.5431932, 0.7044468, 0.5513300   0.5996567
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 1857 out of 43"
    len  freq   err
270   3 0.033 0.167
                                                                                                            condition
270 max_RGH_mean<=0.531406259436042 & texture_variance_nondir_post1<=189.286036058763 & T2kurt_F_r_i>4.73267665300798
    pred                            prob sumtempProb
270   NC 0.4224090, 0.5101070, 0.6310541     0.52119
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1075   3 0.115 0.175
                                                                                               condition
1075 max_F_r_i<=1472.625 & irregularity>0.976213715508256 & T2texture_sumentropy_nondir>1.61393912493361
     pred                            prob sumtempProb
1075    C 0.4853057, 0.4139845, 0.5153649   0.4715517
[1] "test complies with # rules: 2000 out of 43"
     len freq   err                                                                       condition
5725   2 0.19 0.135 circularity<=0.855726888231947 & T2texture_correlation_nondir<=0.32564507505023
     pred                 prob sumtempProb
5725   NC 0.4705375, 0.2400161   0.3552768
[1] "test complies with # rules: 2000 out of 43"
     len freq   err                                                                       condition
5725   2 0.19 0.135 circularity<=0.855726888231947 & T2texture_correlation_nondir<=0.32564507505023
     pred                 prob sumtempProb
5725   NC 0.4705375, 0.2400161   0.3552768
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1060   3 0.114 0.177
                                                                             condition pred
1060 max_F_r_i>1454.125 & dce2SE9>0.644619457119457 & LMSIR_predicted>2.91344151959955    C
                                prob sumtempProb
1060 0.4692342, 0.6110102, 0.7437256     0.60799
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                          condition pred
5422   2 0.062 0.147 beta_countor>0.0714608053643474 & T2_lesionSIstd<=54.2682625055716   NC
                     prob sumtempProb
5422 0.5066059, 0.6093130   0.5579594
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                 condition pred
5545   2 0.081 0.159 dce2SE10<=0.735088866532529 & T2RGH_var<=346.429594497573   NC
                     prob sumtempProb
5545 0.4999955, 0.6814744    0.590735
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
2730   2 0.121 0.121
                                                                                          condition
2730 texture_variance_nondir_post1<=48.1430868084606 & T2texture_sumaverage_nondir>44.8353718717041
     pred                 prob sumtempProb
2730   NC 0.4030437, 0.4857393   0.4443915
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
2499   3 0.042 0.13
                                                                                                   condition
2499 texture_sumaverage_nondir_post1>131.395713357342 & V10<=3.29645166883552 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
2499   NC 0.4583520, 0.4879422, 0.4728997   0.4730646
[1] "test complies with # rules: 1798 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1664 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
5631   3 0.125 0.162
                                                                                                condition
5631 min_F_r_i>386.5 & texture_sumvariance_nondir_post2<=1048.75233234583 & T2kurt_F_r_i>1.14883296665627
     pred                            prob sumtempProb
5631   NC 0.5786830, 0.4887894, 0.4385449   0.5020058
```

```
Area under the curve: 0.8181
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4406966    0.7      0.85       1 0.5957    0.7234  0.8511
```

```
Area under the curve: 0.85
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3014082   0.85      0.95       1 0.5314    0.6596  0.7872
```

```
Area under the curve: 0.8266
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4370225   0.55      0.75    0.95  0.617    0.7447  0.8511
```

```
Area under the curve: 0.5553
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5177135   0.35      0.55    0.75 0.4468    0.5957  0.7447
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-6.png)<!-- -->

```
Area under the curve: 0.8064
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5713629   0.45      0.65    0.85 0.7234    0.8298  0.9362

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

548 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 494, 493, 493, 493, 493, 493, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8063008  0.6034632  0.8298295
  2                  1500     0.8077918  0.6034632  0.8237689
  2                  2000     0.8058265  0.6216450  0.8207386
  3                  1000     0.8097808  0.5805195  0.8329545
  3                  1500     0.8118295  0.5941558  0.8177083
  3                  2000     0.8099231  0.6127706  0.8086174
  4                  1000     0.8102400  0.5896104  0.8451705
  4                  1500     0.8114967  0.6125541  0.8358902
  4                  2000     0.8115358  0.6309524  0.8177083
  6                  1000     0.8104167  0.5761905  0.8391098
  6                  1500     0.8136368  0.6082251  0.8422348
  6                  2000     0.8148982  0.6127706  0.8330492

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 6, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8063008 0.6034632 0.8298295 0.05059575
4       0.01                 3             10    1000 0.8097808 0.5805195 0.8329545 0.05630040
7       0.01                 4             10    1000 0.8102400 0.5896104 0.8451705 0.05680557
10      0.01                 6             10    1000 0.8104167 0.5761905 0.8391098 0.06382138
2       0.01                 2             10    1500 0.8077918 0.6034632 0.8237689 0.05332938
5       0.01                 3             10    1500 0.8118295 0.5941558 0.8177083 0.05996379
8       0.01                 4             10    1500 0.8114967 0.6125541 0.8358902 0.05645014
11      0.01                 6             10    1500 0.8136368 0.6082251 0.8422348 0.06583967
3       0.01                 2             10    2000 0.8058265 0.6216450 0.8207386 0.05275029
6       0.01                 3             10    2000 0.8099231 0.6127706 0.8086174 0.06201096
9       0.01                 4             10    2000 0.8115358 0.6309524 0.8177083 0.05820047
12      0.01                 6             10    2000 0.8148982 0.6127706 0.8330492 0.06679959
       SensSD     SpecSD
1  0.10755347 0.05367285
4  0.09320940 0.05358209
7  0.10316300 0.06268664
10 0.11144699 0.06675588
2  0.09629156 0.06497076
5  0.08478938 0.06976974
8  0.09448079 0.05740447
11 0.12556600 0.07902290
3  0.09651835 0.06595946
6  0.10586385 0.06524011
9  0.09506987 0.06202832
12 0.12198463 0.06853754
14000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
2358   3 0.104 0.105
                                                                                                         condition
2358 iiMin_change_Variance_uptake>0.283856875469795 & irregularity<=0.970416148881159 & T2RGH_var>353.975029351905
     pred                            prob sumtempProb
2358   NC 0.5610004, 0.2744690, 0.4693122   0.4349272
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 22 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 17 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 14 out of 160"
     len  freq   err
8602   6 0.128 0.071
                                                                                                                                                                                               condition
8602 SER_inside>0.838500730200154 & irregularity>0.940117684730477 & V5>3.49086708315126 & dce2SE7>0.625778063035859 & lateSE2>0.806322857503407 & T2texture_inversediffmoment_nondir<=0.201898180625362
     pred                                                             prob sumtempProb
8602    C 0.6746109, 0.4492241, 0.4756217, 0.4966826, 0.4517067, 0.6258111   0.5289428
[1] "test complies with # rules: 1893 out of 160"
[1] "test complies with # rules: 15 out of 160"
[1] "test complies with # rules: 30 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
4965   5 0.084 0.109
                                                                                                                                                                                                                                   condition
4965 UptakeRate_inside<=0.331212825563122 & texture_variance_nondir_post1<=188.263791397188 & texture_inversediffmoment_nondir_post1>0.208643298985366 & dce2SE17<=0.837805700721211 & T2texture_inversediffmoment_nondir<=0.221325574807973
     pred                                                  prob sumtempProb
4965   NC 0.3571228, 0.4603757, 0.3649369, 0.5094863, 0.4878860   0.4359615
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 5 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
3267   5 0.186 0.137
                                                                                                                                                                   condition
3267 SER_inside<=0.794461756516995 & texture_sumaverage_nondir_post2<=250.578758618558 & lateSE9>0.770890728378905 & lateSE9>1.24938940369975 & T2_lesionSI>54.1854160477181
     pred                                                  prob sumtempProb
3267   NC 0.3634066, 0.4839432, 0.4559588, 0.4990479, 0.5065191   0.4617751
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
      len  freq   err
11190   3 0.086 0.128
                                                                                                           condition
11190 texture_sumvariance_nondir_post1<=415.770720390637 & earlySE14>1.17045467703765 & T2_lesionSI>53.2834535991432
      pred                            prob sumtempProb
11190   NC 0.4248698, 0.3256649, 0.4487710   0.3997685
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
2924   4 0.133 0.137
                                                                                                                                                                    condition
2924 texture_variance_nondir_post2<=180.876506362501 & texture_diffvariance_nondir_post3<=338.678697077251 & T2grad_margin_var>4935.73456149471 & T2RGH_mean>23.9432421868704
     pred                                       prob sumtempProb
2924   NC 0.5609790, 0.5372390, 0.3844625, 0.4686015   0.4878205
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
5857   3 0.035 0.105
                                                                                                                        condition
5857 skew_F_r_i>-0.0668608838335453 & texture_inversediffmoment_nondir_post4<=0.0913900931574419 & T2skew_F_r_i>0.134325109470664
     pred                            prob sumtempProb
5857   NC 0.4848748, 0.3959785, 0.4346656   0.4385063
[1] "test complies with # rules: 13 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 7 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 18 out of 160"
[1] "test complies with # rules: 20 out of 160"
     len  freq                err
3964   4 0.091 0.0600000000000001
                                                                                                                                        condition
3964 texture_energy_nondir_post3<=0.00710250947915509 & V12<=4.03021189922947 & earlySE10>0.412776483160783 & T2grad_margin_var<=13087.2498193737
     pred                                       prob sumtempProb
3964   NC 0.5148103, 0.4888280, 0.5169547, 0.4973049   0.5044745
[1] "test complies with # rules: 23 out of 160"
     len  freq                err
6132   4 0.058 0.0620000000000001
                                                                                                                  condition
6132 mean_F_r_i>811.726530585355 & earlySE4>0.944321360871241 & earlySE8>0.962940140610659 & T2grad_margin>31.6517431192661
     pred                                       prob sumtempProb
6132    C 0.4419462, 0.6174297, 0.4525536, 0.5221466    0.508519
[1] "test complies with # rules: 10 out of 160"
     len  freq   err
5658   3 0.177 0.093
                                                                                                                                 condition
5658 SER_countor<=0.675838603601509 & texture_diffvariance_nondir_post1<=52.6220679340441 & T2texture_correlation_nondir>0.107935297719362
     pred                            prob sumtempProb
5658   NC 0.5505793, 0.3818076, 0.5017071   0.4780313
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 16 out of 160"
```

```
Area under the curve: 0.7874
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5251699 0.3333    0.5417    0.75 0.8718    0.9487       1
```

```
Area under the curve: 0.7863
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3986131 0.5833      0.75  0.9167  0.641    0.7692  0.8974
```

```
Area under the curve: 0.797
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3684782 0.7083    0.8333  0.9583 0.5128    0.6667  0.8205
```

```
Area under the curve: 0.8056
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5238473 0.4583    0.6667  0.8333 0.8205    0.9231       1
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-7.png)<!-- -->

```
Area under the curve: 0.7671
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.0171481 0.6667    0.8333  0.9583 0.4615    0.6154  0.7692

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

544 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 491, 490, 490, 489, 491, 489, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8127091  0.6049784  0.8623106
  2                  1500     0.8131961  0.5956710  0.8714015
  2                  2000     0.8130989  0.6047619  0.8560606
  3                  1000     0.8149019  0.6004329  0.8652462
  3                  1500     0.8151835  0.5822511  0.8590909
  3                  2000     0.8175048  0.5958874  0.8468750
  4                  1000     0.8159337  0.5822511  0.8714015
  4                  1500     0.8185151  0.5774892  0.8681818
  4                  2000     0.8181541  0.5958874  0.8469697
  6                  1000     0.8215952  0.5774892  0.8743371
  6                  1500     0.8212755  0.5820346  0.8621212
  6                  2000     0.8239583  0.5863636  0.8590909

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 6, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8127091 0.6049784 0.8623106 0.05337614
4       0.01                 3             10    1000 0.8149019 0.6004329 0.8652462 0.04693579
7       0.01                 4             10    1000 0.8159337 0.5822511 0.8714015 0.05394636
10      0.01                 6             10    1000 0.8215952 0.5774892 0.8743371 0.05268320
2       0.01                 2             10    1500 0.8131961 0.5956710 0.8714015 0.05069888
5       0.01                 3             10    1500 0.8151835 0.5822511 0.8590909 0.04970707
8       0.01                 4             10    1500 0.8185151 0.5774892 0.8681818 0.05367539
11      0.01                 6             10    1500 0.8212755 0.5820346 0.8621212 0.05307959
3       0.01                 2             10    2000 0.8130989 0.6047619 0.8560606 0.05112734
6       0.01                 3             10    2000 0.8175048 0.5958874 0.8468750 0.04919247
9       0.01                 4             10    2000 0.8181541 0.5958874 0.8469697 0.05061290
12      0.01                 6             10    2000 0.8239583 0.5863636 0.8590909 0.04651912
       SensSD     SpecSD
1  0.12346836 0.08126940
4  0.11922982 0.05368955
7  0.11930315 0.06025784
10 0.11899645 0.06351689
2  0.13589443 0.08456863
5  0.12121212 0.06121426
8  0.12279425 0.05221956
11 0.13136644 0.06302105
3  0.12822917 0.06725041
6  0.12049300 0.06054528
9  0.09963531 0.05707399
12 0.12626654 0.05778464
14000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 2000 out of 168"
    len  freq err
505   3 0.055 0.1
                                                                                                       condition
505 SER_countor<=0.707326058134249 & var_F_r_i<=19659.5152094893 & T2texture_energy_nondir<=0.000632018983719888
    pred                            prob sumtempProb
505   NC 0.4811075, 0.5600435, 0.3258341   0.4556617
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 8 out of 168"
      len  freq   err
13497   4 0.059 0.094
                                                                                                                                                         condition
13497 texture_sumaverage_nondir_post2<=226.942106218264 & lateSE8<=0.920961171375684 & T2RGH_var<=521.991169084647 & T2texture_energy_nondir<=0.000898111923071435
      pred                                       prob sumtempProb
13497   NC 0.4324431, 0.5009374, 0.2695222, 0.3900091   0.3982279
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 23 out of 168"
[1] "test complies with # rules: 27 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
4377   3 0.053 0.103
                                                                                                          condition
4377 texture_variance_nondir_post3<=223.487783416407 & dce2SE1>0.868363061127592 & T2kurt_F_r_i<=-0.520345577574295
     pred                            prob sumtempProb
4377   NC 0.5659013, 0.6675191, 0.4893800   0.5742668
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 2000 out of 168"
      len  freq   err
12729   5 0.142 0.104
                                                                                                                                                                  condition
12729 maxVr_inside>0.10451017124894 & maxVr_countor>0.0568363361355432 & texture_sumaverage_nondir_post3<=265.972128929092 & earlySE8<=0.967348276913911 & T2max_F_r_i<=406
      pred                                                  prob sumtempProb
12729   NC 0.5068127, 0.4770168, 0.4879522, 0.4684545, 0.3920593   0.4664591
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 16 out of 168"
     len  freq   err
7178   6 0.086 0.085
                                                                                                                                                                                         condition
7178 mean_F_r_i<=1117.30129009231 & irregularity<=0.974279294868243 & texture_sumaverage_nondir_post1>126.361051313361 & V0<=11.82706363796 & V12<=4.67723478611074 & T2_lesionSI>43.6236326923077
     pred                                                             prob sumtempProb
7178   NC 0.5334060, 0.5257874, 0.4944856, 0.4224559, 0.5227059, 0.5457399   0.5074301
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 22 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq  err
6133   3 0.138 0.12
                                                                            condition pred
6133 A_inside<=226.058984902916 & V2<=5.23315348512354 & T2RGH_mean<=54.5119132719293   NC
                                prob sumtempProb
6133 0.5390625, 0.5116299, 0.5106702   0.5204542
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 14 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition pred
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077    C
                                prob sumtempProb
8099 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 7 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition pred
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077    C
                                prob sumtempProb
8099 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq err
6891   3 0.055 0.1
                                                                                         condition
6891 ivVariance>0.0331643163273542 & circularity<=0.855734293227423 & T2RGH_mean<=54.5391564448446
     pred                            prob sumtempProb
6891   NC 0.4852294, 0.5920242, 0.5623309   0.5465282
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 24 out of 168"
[1] "test complies with # rules: 2000 out of 168"
      len  freq   err
12729   5 0.142 0.104
                                                                                                                                                                  condition
12729 maxVr_inside>0.10451017124894 & maxVr_countor>0.0568363361355432 & texture_sumaverage_nondir_post3<=265.972128929092 & earlySE8<=0.967348276913911 & T2max_F_r_i<=406
      pred                                                  prob sumtempProb
12729   NC 0.5068127, 0.4770168, 0.4879522, 0.4684545, 0.3920593   0.4664591
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 15 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
1833   3 0.088 0.104
                                                                                              condition
1833 irregularity<=0.984466492362561 & T2kurt_F_r_i>4.34329574043788 & LMSIR_predicted>1.98323389316333
     pred                            prob sumtempProb
1833   NC 0.4854764, 0.4002918, 0.4276336   0.4378006
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 8 out of 168"
     len  freq   err
1934   3 0.156 0.082
                                                                                   condition pred
1934 dce2SE12>0.417337428870872 & dce2SE15<=0.905346458673219 & T2RGH_mean<=27.5207599962328   NC
                                prob sumtempProb
1934 0.3281155, 0.4799466, 0.2881014   0.3653878
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
1273   5 0.105 0.105
                                                                                                                                                                                 condition
1273 max_RGH_var>0.0733863679263791 & texture_sumvariance_nondir_post1>261.264914441484 & V16>7.87099690126964 & lateSE2>0.817789299079648 & T2texture_energy_nondir<=0.000543922788822772
     pred                                                  prob sumtempProb
1273    C 0.6443449, 0.4816719, 0.6141464, 0.6863575, 0.7402640   0.6333569
[1] "test complies with # rules: 12 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition pred
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077    C
                                prob sumtempProb
8099 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 10 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition pred
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077    C
                                prob sumtempProb
8099 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 11 out of 168"
     len  freq   err
5842   3 0.066 0.083
                                                                                                                       condition
5842 ivVariance<=0.0375174805963147 & T2texture_entropy_nondir<=3.48878658358655 & T2texture_diffentropy_nondir>1.70774926189513
     pred                            prob sumtempProb
5842   NC 0.4992823, 0.5380998, 0.4868618   0.5080813
[1] "test complies with # rules: 7 out of 168"
     len  freq   err
8358   5 0.039 0.095
                                                                                                                                                                                 condition
8358 texture_sumaverage_nondir_post3>225.1975459895 & lateSE0>0.871147809637137 & lateSE0>0.987704612189311 & T2_lesionSI>98.472330941509 & T2texture_correlation_nondir>0.263195249129851
     pred                                                  prob sumtempProb
8358    C 0.4192649, 0.5621379, 0.5517737, 0.6701021, 0.4812520   0.5369061
[1] "test complies with # rules: 13 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 26 out of 168"
[1] "test complies with # rules: 11 out of 168"
     len  freq   err
5842   3 0.066 0.083
                                                                                                                       condition
5842 ivVariance<=0.0375174805963147 & T2texture_entropy_nondir<=3.48878658358655 & T2texture_diffentropy_nondir>1.70774926189513
     pred                            prob sumtempProb
5842   NC 0.4992823, 0.5380998, 0.4868618   0.5080813
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 6 out of 168"
```

```
Area under the curve: 0.8305
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4711738   0.52      0.72    0.88 0.7143    0.8333  0.9286
```

```
Area under the curve: 0.8248
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4113339    0.8      0.92       1 0.4524     0.619  0.7619
```

```
Area under the curve: 0.8543
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4115153   0.68      0.84    0.96 0.5952    0.7381   0.881
```

```
Area under the curve: 0.6905
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4938438   0.52      0.72    0.88 0.4762     0.619  0.7619
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-8.png)<!-- -->

```
Area under the curve: 0.7505
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5408004   0.32      0.52    0.72  0.881    0.9524       1

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

541 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 487, 487, 487, 488, 487, 487, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.8193228  0.6339921  0.8550403
  2                  1500     0.8227584  0.6517787  0.8517137
  2                  2000     0.8280975  0.6654150  0.8454637
  3                  1000     0.8283940  0.6385375  0.8613911
  3                  1500     0.8324830  0.6519763  0.8675403
  3                  2000     0.8345553  0.6519763  0.8612903
  4                  1000     0.8309932  0.6298419  0.8610887
  4                  1500     0.8349484  0.6385375  0.8643145
  4                  2000     0.8409960  0.6383399  0.8643145
  6                  1000     0.8290135  0.6296443  0.8737903
  6                  1500     0.8328231  0.6298419  0.8706653
  6                  2000     0.8334825  0.6341897  0.8644153

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 2000, interaction.depth = 4, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.8193228 0.6339921 0.8550403 0.07752424
4       0.01                 3             10    1000 0.8283940 0.6385375 0.8613911 0.07918200
7       0.01                 4             10    1000 0.8309932 0.6298419 0.8610887 0.08288330
10      0.01                 6             10    1000 0.8290135 0.6296443 0.8737903 0.08956652
2       0.01                 2             10    1500 0.8227584 0.6517787 0.8517137 0.08068377
5       0.01                 3             10    1500 0.8324830 0.6519763 0.8675403 0.08107106
8       0.01                 4             10    1500 0.8349484 0.6385375 0.8643145 0.08377263
11      0.01                 6             10    1500 0.8328231 0.6298419 0.8706653 0.08947922
3       0.01                 2             10    2000 0.8280975 0.6654150 0.8454637 0.08218668
6       0.01                 3             10    2000 0.8345553 0.6519763 0.8612903 0.08630293
9       0.01                 4             10    2000 0.8409960 0.6383399 0.8643145 0.08370573
12      0.01                 6             10    2000 0.8334825 0.6341897 0.8644153 0.08950412
      SensSD     SpecSD
1  0.1149866 0.05789291
4  0.1305660 0.06268625
7  0.1259357 0.05983410
10 0.1244665 0.07694203
2  0.1285034 0.06111611
5  0.1323847 0.05892978
8  0.1239643 0.06987664
11 0.1180895 0.07466796
3  0.1296206 0.06834061
6  0.1307885 0.06802788
9  0.1154243 0.06669870
12 0.1125594 0.06785014
10000 rules (length<=6) were extracted from the first 2000 trees.
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
6458   4 0.044 0.167
                                                                                                                              condition
6458 maxVr_inside<=0.244543536024623 & irregularity<=0.975937406084662 & T2grad_margin_var>6224.16338048594 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
6458   NC 0.4425387, 0.5276371, 0.6897943, 0.5563826   0.5540882
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
5961   4 0.181 0.122
                                                                                                         condition
5961 earlySE8<=0.954301342179907 & lateSE4<=1.71632227630424 & T2max_F_r_i<=408.5 & T2kurt_F_r_i<=5.01194223263781
     pred                                       prob sumtempProb
5961   NC 0.5749753, 0.4856098, 0.4904649, 0.5319248   0.5207437
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4063   4 0.128 0.116
                                                                                                     condition
4063 min_F_r_i>389 & mean_F_r_i<=996.57033595173 & lateSE18<=2.01081577015833 & T2kurt_F_r_i>0.588717905867016
     pred                                       prob sumtempProb
4063   NC 0.5037287, 0.4757571, 0.4159564, 0.3816225   0.4442662
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9774   4 0.165 0.146
                                                                                                                                                               condition
9774 Tpeak_inside>2.36694862113409 & texture_diffvariance_nondir_post1<=120.496206781367 & texture_contrast_nondir_post2>194.701661779236 & T2_lesionSI>44.0377680747215
     pred                                       prob sumtempProb
9774   NC 0.5050244, 0.5244491, 0.4491952, 0.5204948   0.4997909
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err                                                                  condition
9172   3 0.065 0.114 dce2SE8>0.526488212923419 & lateSE13<=0.782897303276147 & T2min_F_r_i<=8.5
     pred                            prob sumtempProb
9172   NC 0.5711305, 0.5520445, 0.4931413   0.5387721
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9914   4 0.028 0.133
                                                                                                                         condition
9914 max_RGH_var>0.0719122415437442 & V13>12.797357398796 & V18<=6.75866036727473 & T2texture_diffvariance_nondir>28.0656417109182
     pred                                       prob sumtempProb
9914    C 0.4789980, 0.5192834, 0.6233938, 0.5384807    0.540039
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4063   4 0.128 0.116
                                                                                                     condition
4063 min_F_r_i>389 & mean_F_r_i<=996.57033595173 & lateSE18<=2.01081577015833 & T2kurt_F_r_i>0.588717905867016
     pred                                       prob sumtempProb
4063   NC 0.5037287, 0.4757571, 0.4159564, 0.3816225   0.4442662
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4063   4 0.128 0.116
                                                                                                     condition
4063 min_F_r_i>389 & mean_F_r_i<=996.57033595173 & lateSE18<=2.01081577015833 & T2kurt_F_r_i>0.588717905867016
     pred                                       prob sumtempProb
4063   NC 0.5037287, 0.4757571, 0.4159564, 0.3816225   0.4442662
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
7549   3 0.083 0.156
                                                                                                   condition
7549 iMax_Variance_uptake<=6.70552353229436 & lateSE6>0.586002782390577 & T2grad_margin_var>6372.15933273777
     pred                            prob sumtempProb
7549   NC 0.3006064, 0.4980637, 0.4829544   0.4272082
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9774   4 0.165 0.146
                                                                                                                                                               condition
9774 Tpeak_inside>2.36694862113409 & texture_diffvariance_nondir_post1<=120.496206781367 & texture_contrast_nondir_post2>194.701661779236 & T2_lesionSI>44.0377680747215
     pred                                       prob sumtempProb
9774   NC 0.5050244, 0.5244491, 0.4491952, 0.5204948   0.4997909
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9828   4 0.091 0.143
                                                                                                                                             condition
9828 iMax_Variance_uptake<=12.5030776866371 & texture_diffvariance_nondir_post2>204.302110598256 & V9<=51.1881043457442 & T2_lesionSI>44.6703853824138
     pred                                       prob sumtempProb
9828    C 0.6681251, 0.4155613, 0.4415098, 0.4555614   0.4951894
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 7 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
2249   4 0.102 0.145
                                                                                                                                                                                         condition
2249 Vr_increasingRate_countor>0.0954987823275664 & Vr_decreasingRate_countor<=0.0636349491408924 & texture_sumvariance_nondir_post1>381.737085837675 & T2texture_contrast_nondir>338.784043932143
     pred                                       prob sumtempProb
2249    C 0.5529671, 0.5149765, 0.4829129, 0.6764732   0.5568325
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4635   3 0.052 0.107
                                                                                       condition
4635 edge_sharp_mean>0.671502660560278 & dce3SE1>1.22681788388568 & T2_lesionSI>50.6292782580445
     pred                            prob sumtempProb
4635   NC 0.4823546, 0.2729908, 0.4978068   0.4177174
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
5111   3 0.166 0.133
                                                                                   condition pred
5111 SER_inside<=0.682123129359999 & lateSE9<=1.19076743320431 & T2RGH_var<=575.314918315023   NC
                                prob sumtempProb
5111 0.5941872, 0.5023692, 0.3768535   0.4911366
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 5 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
5111   3 0.166 0.133
                                                                                   condition pred
5111 SER_inside<=0.682123129359999 & lateSE9<=1.19076743320431 & T2RGH_var<=575.314918315023   NC
                                prob sumtempProb
5111 0.5941872, 0.5023692, 0.3768535   0.4911366
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4635   3 0.052 0.107
                                                                                       condition
4635 edge_sharp_mean>0.671502660560278 & dce3SE1>1.22681788388568 & T2_lesionSI>50.6292782580445
     pred                            prob sumtempProb
4635   NC 0.4823546, 0.2729908, 0.4978068   0.4177174
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
1153   3 0.207 0.125
                                                                                                          condition
1153 Tpeak_inside>7.90722296680836 & texture_energy_nondir_post4<=0.00618966000712137 & T2RGH_var<=349.491222787327
     pred                            prob sumtempProb
1153   NC 0.4795554, 0.3563718, 0.3132517   0.3830596
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 2000 out of 47"
```

```
Area under the curve: 0.7729
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5067831    0.5       0.7     0.9 0.6875    0.8125  0.9167
```

```
Area under the curve: 0.7896
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4316145   0.55      0.75    0.95 0.6458    0.7708   0.875
```

```
Area under the curve: 0.7865
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3746715    0.6       0.8    0.95 0.5833    0.7083  0.8333
```

```
Area under the curve: 0.6865
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5035049   0.55      0.75     0.9 0.5833    0.7083  0.8333
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-9.png)<!-- -->

```
Area under the curve: 0.7344
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
 0.05547199   0.75       0.9       1 0.3333    0.4792   0.625

============ rule gbm trees imgT2pLMSIR 
Stochastic Gradient Boosting 

546 samples
221 predictors
  2 classes: 'C', 'NC' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 491, 491, 492, 492, 492, 492, ... 
Resampling results across tuning parameters:

  interaction.depth  n.trees  ROC        Sens       Spec     
  2                  1000     0.7985474  0.5586580  0.8557932
  2                  1500     0.8041488  0.5725108  0.8557041
  2                  2000     0.8040003  0.5915584  0.8467023
  3                  1000     0.8037888  0.5636364  0.8679144
  3                  1500     0.8050563  0.5822511  0.8616756
  3                  2000     0.8045717  0.5961039  0.8496435
  4                  1000     0.8131664  0.5634199  0.8738859
  4                  1500     0.8192131  0.5915584  0.8769162
  4                  2000     0.8175714  0.6097403  0.8617647
  6                  1000     0.8042125  0.5588745  0.8769162
  6                  1500     0.8082919  0.5727273  0.8527629
  6                  2000     0.8121606  0.6008658  0.8437611

Tuning parameter 'shrinkage' was held constant at a value of 0.01
Tuning parameter
 'n.minobsinnode' was held constant at a value of 10
ROC was used to select the optimal model using  the largest value.
The final values used for the model were n.trees = 1500, interaction.depth = 4, shrinkage =
 0.01 and n.minobsinnode = 10. 
   shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD
1       0.01                 2             10    1000 0.7985474 0.5586580 0.8557932 0.03744263
4       0.01                 3             10    1000 0.8037888 0.5636364 0.8679144 0.04092713
7       0.01                 4             10    1000 0.8131664 0.5634199 0.8738859 0.04361742
10      0.01                 6             10    1000 0.8042125 0.5588745 0.8769162 0.04527763
2       0.01                 2             10    1500 0.8041488 0.5725108 0.8557041 0.04443425
5       0.01                 3             10    1500 0.8050563 0.5822511 0.8616756 0.04127416
8       0.01                 4             10    1500 0.8192131 0.5915584 0.8769162 0.04668390
11      0.01                 6             10    1500 0.8082919 0.5727273 0.8527629 0.04823514
3       0.01                 2             10    2000 0.8040003 0.5915584 0.8467023 0.04418395
6       0.01                 3             10    2000 0.8045717 0.5961039 0.8496435 0.04260630
9       0.01                 4             10    2000 0.8175714 0.6097403 0.8617647 0.04952344
12      0.01                 6             10    2000 0.8121606 0.6008658 0.8437611 0.04932091
       SensSD     SpecSD
1  0.08298482 0.05987023
4  0.12023848 0.05110924
7  0.11230432 0.05206302
10 0.10991845 0.05871607
2  0.08832566 0.06146280
5  0.10436103 0.06703387
8  0.11131730 0.05695188
11 0.12471824 0.05759657
3  0.10669456 0.06442370
6  0.10568545 0.07001218
9  0.11364976 0.06237791
12 0.10958359 0.05483847
7500 rules (length<=6) were extracted from the first 1500 trees.
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5012   3 0.176 0.167
                                                                                                                                  condition
5012 texture_diffvariance_nondir_post1<=60.9797692220328 & lateSE9>1.16541692526887 & T2texture_inversediffmoment_nondir<=0.246115364253849
     pred                            prob sumtempProb
5012   NC 0.4782870, 0.3576450, 0.4992456   0.4450592
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5637   4 0.104 0.123
                                                                                                                          condition
5637 A_inside<=1.5528673159683 & V9<=16.7536361430883 & earlySE6<=0.819922077850664 & T2texture_sumvariance_nondir>185.022460304943
     pred                                       prob sumtempProb
5637   NC 0.4832756, 0.5288236, 0.2834297, 0.5989886   0.4736294
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
2001   3 0.165 0.111
                                                                                                           condition
2001 SER_inside<=1.00080283237816 & texture_sumaverage_nondir_post2<=249.582981542509 & T2RGH_mean<=22.5838715481632
     pred                            prob sumtempProb
2001   NC 0.4448360, 0.3843497, 0.3452982   0.3914946
[1] "test complies with # rules: 8 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len freq   err
813   3 0.11 0.133
                                                                                                                                    condition
813 SER_inside<=0.806001166108312 & texture_inversediffmoment_nondir_post3>0.147244327606865 & T2texture_sumvariance_nondir<=381.387338187838
    pred                            prob sumtempProb
813   NC 0.6043602, 0.4880339, 0.4552047   0.5158663
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3657   3 0.209 0.123
                                                                                          condition
3657 SER_inside<=0.683454887275891 & skew_F_r_i<=0.724263044011434 & T2grad_margin>30.5173903603995
     pred                            prob sumtempProb
3657   NC 0.5138256, 0.2956922, 0.4628214   0.4241131
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq  err
5742   3 0.046 0.16
                                                                                                                      condition
5742 texture_diffvariance_nondir_post2<=24.0599300089026 & V15<=44.4612422149372 & T2texture_sumaverage_nondir>22.9708102937136
     pred                            prob sumtempProb
5742   NC 0.4570783, 0.5057757, 0.4847688   0.4825409
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5551   4 0.273 0.154
                                                                                                                                                condition
5551 A_inside<=312.396375500685 & edge_sharp_std<=0.991727213591119 & T2RGH_var<=348.526546275196 & T2texture_inversediffmoment_nondir<=0.254036393444923
     pred                                       prob sumtempProb
5551   NC 0.4330791, 0.2267610, 0.4129190, 0.4568792   0.3824096
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5778   4 0.068 0.135
                                                                                                                                                                           condition
5778 Vr_post_1_inside>0.114510008249947 & texture_sumaverage_nondir_post1<=194.065626177225 & texture_inversediffmoment_nondir_post3<=0.104179970009203 & T2RGH_mean>21.889655177071
     pred                                       prob sumtempProb
5778   NC 0.4829825, 0.5148890, 0.4299427, 0.4467408   0.4686388
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err                                                                      condition
6440   2 0.082 0.133 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052
     pred                 prob sumtempProb
6440    C 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6543   4 0.207 0.159
                                                                                                                             condition
6543 alpha_inside<=1.57751330223065 & circularity<=0.855732884884821 & V5>4.51415360177622 & T2texture_entropy_nondir>2.90485367714144
     pred                                       prob sumtempProb
6543   NC 0.4572355, 0.4567487, 0.4936095, 0.4366377   0.4610578
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err                                                                      condition
6440   2 0.082 0.133 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052
     pred                 prob sumtempProb
6440    C 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err                                                                      condition
6440   2 0.082 0.133 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052
     pred                 prob sumtempProb
6440    C 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3408   3 0.051 0.143
                                                                                                      condition
3408 A_inside<=1.23329492688817 & T2RGH_mean<=52.7500988159803 & T2texture_correlation_nondir>0.255954782848448
     pred                            prob sumtempProb
3408   NC 0.7046357, 0.5281514, 0.5599449   0.5975773
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1024   3 0.313 0.164
                                                                                                                                         condition
1024 texture_variance_nondir_post2<=179.048343064544 & T2texture_energy_nondir>0.000478831110531278 & T2texture_sumaverage_nondir>40.2571158312255
     pred                            prob sumtempProb
1024   NC 0.4499608, 0.2332170, 0.3358754   0.3396844
[1] "test complies with # rules: 9 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len freq   err
813   3 0.11 0.133
                                                                                                                                    condition
813 SER_inside<=0.806001166108312 & texture_inversediffmoment_nondir_post3>0.147244327606865 & T2texture_sumvariance_nondir<=381.387338187838
    pred                            prob sumtempProb
813   NC 0.6043602, 0.4880339, 0.4552047   0.5158663
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3675   3 0.059 0.125
                                                                                    condition pred
3675 max_RGH_mean>0.524313219791607 & lateSE4>0.953138290625636 & T2RGH_mean>47.5488747094981   NC
                                prob sumtempProb
3675 0.5694536, 0.5043591, 0.5279330   0.5339152
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len  freq   err
550   3 0.148 0.198
                                                                                                                                     condition
550 texture_variance_nondir_post1>167.08364168864 & texture_inversediffmoment_nondir_post4>0.11508889696368 & LMSIR_predicted>1.76548498403474
    pred                            prob sumtempProb
550    C 0.4777159, 0.6997331, 0.5019843   0.5598111
[1] "test complies with # rules: 7 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6277   4 0.134 0.137
                                                                                                                condition
6277 V17<=38.6058204388241 & dce2SE4>0.633227549675858 & dce3SE4<=0.950788958274215 & T2grad_margin_var<=6402.45036359788
     pred                                       prob sumtempProb
6277   NC 0.5104606, 0.5866031, 0.5544229, 0.5395111   0.5477494
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 10 out of 51"
     len  freq   err
1758   3 0.115 0.095
                                                                                                            condition
1758 circularity<=0.855726514911507 & texture_sumaverage_nondir_post2<=249.928310833517 & T2RGH_mean>26.4497935000351
     pred                            prob sumtempProb
1758   NC 0.5552296, 0.4969470, 0.4644038   0.5055268
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1976   3 0.156 0.129
                                                                                  condition pred
1976 earlySE11<=0.882260181441451 & lateSE3<=1.10244673946561 & T2RGH_var<=348.899261925448   NC
                                prob sumtempProb
1976 0.5758153, 0.4910063, 0.5130602   0.5266273
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
```

```
Area under the curve: 0.8879
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.427195 0.6667       0.8  0.9333 0.7879    0.9091       1
```

```
Area under the curve: 0.901
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3753419 0.7333    0.8667  0.9667 0.7576    0.8788  0.9697
```

```
Area under the curve: 0.9394
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3675162 0.8333    0.9333       1 0.7576    0.8788  0.9697
```

```
Area under the curve: 0.7273
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5011773 0.6333       0.8  0.9333 0.4242    0.6061  0.7583
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-4-10.png)<!-- -->

```
Area under the curve: 0.898
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.1736443 0.7333    0.8667  0.9667  0.697    0.8485  0.9697
```

```r
# save current state k patient out
save.image("Outputs/STEL_noT2SIpredLMSIR_boost_addeddiagvalue.RData")
```

```r
confusionMatrix(allrulesres[, "predR"], allrulesres[, "obsR"])
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  119  40
        NC 126 342
                                          
               Accuracy : 0.7352          
                 95% CI : (0.6989, 0.7694)
    No Information Rate : 0.6093          
    P-Value [Acc > NIR] : 2.175e-11       
                                          
                  Kappa : 0.4066          
 Mcnemar's Test P-Value : 4.188e-11       
                                          
            Sensitivity : 0.4857          
            Specificity : 0.8953          
         Pos Pred Value : 0.7484          
         Neg Pred Value : 0.7308          
             Prevalence : 0.3907          
         Detection Rate : 0.1898          
   Detection Prevalence : 0.2536          
      Balanced Accuracy : 0.6905          
                                          
       'Positive' Class : C               
                                          
```

```r
par(mfrow = c(1, 1))
pModel = calcAUC_plot(allrulesres$obsR, as.numeric(allrulesres[, "probRules"]), xptext = 0.75, yptext = 0.45, 
    1, colors[15], atitle = "")
```

```
Area under the curve: 0.7035
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5103875 0.5755    0.6367   0.698  0.644    0.6911  0.7356
```

```r
par(new = TRUE)
pRules = calcAUC_plot(allrulesres$obsR, as.numeric(allrulesres[, "probModel"]), xptext = 0.65, yptext = 0.55, 
    5, colors[2], atitle = "ROCs 10f-patient out cv test all folds ")
```

```
Area under the curve: 0.785
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.1910424 0.6653    0.7224  0.7796  0.678    0.7225  0.7644
```

```r
legend("bottomright", legend = c(paste0("probRules"), paste0("probModel")), col = c(colors[15], colors[2]), 
    lty = c(1, 5), lwd = 2)
```

![](Rulextraction_Section2_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
thresh = as.numeric(row.names(pModel$best_thr$sensitivity))
confusionMatrix(ifelse(allrulesres[, "probModel"] >= thresh, "C", "NC"), allrulesres[, "obsR"])
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  143  69
        NC 102 313
                                          
               Accuracy : 0.7273          
                 95% CI : (0.6906, 0.7618)
    No Information Rate : 0.6093          
    P-Value [Acc > NIR] : 3.675e-10       
                                          
                  Kappa : 0.413           
 Mcnemar's Test P-Value : 0.0144          
                                          
            Sensitivity : 0.5837          
            Specificity : 0.8194          
         Pos Pred Value : 0.6745          
         Neg Pred Value : 0.7542          
             Prevalence : 0.3907          
         Detection Rate : 0.2281          
   Detection Prevalence : 0.3381          
      Balanced Accuracy : 0.7015          
                                          
       'Positive' Class : C               
                                          
```

