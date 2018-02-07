# Boosting classification trees - T2SI vs predLMSIR with extra trees 
Cristina Gallego  
Sept 5, 2017  

* This uses boosted tree-ensembles of T1w, T1w+T2w and T1w+T2w (with predicted LMSIR) 
* This code analysis T2w added diagnostic value by comparing with ensembles of only T1w DCE-based features
* T2w discrimination ability (added AUC ROC value)



## Run boostingTrees with internal cross-validation of parameters

```r
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, "textureUpdatedFeatures.db")

# 2) all T1w features
lesionsQuery <- dbGetQuery(conn, "SELECT *
         FROM  stage1features
           INNER JOIN lesion ON (stage1features.lesion_id = lesion.lesion_id)
           INNER JOIN f_T2 ON (stage1features.lesion_id = f_T2.lesion_id)")
 
# For bootstrap resampling, createResample is used
# Randomization is done at the patient level unique ids
lesionsQuery = subset(lesionsQuery, lesion_label != "fociB" & lesion_label != "fociM" ) # exclude foci at this point
id_cad_pts = lesionsQuery$cad_pt_no_txt
uniq_cad = unique(lesionsQuery$cad_pt_no_txt)
npatients = 1:length(uniq_cad)
  
# when y is a factor in an attempt to balance the class distributions within the splits.
# The names of the list objects will denote the fold membership using the pattern 
# resamples." meaning the ith section (of k) of the jth cross-validation set (of times).
set.seed(1234)
npatients = length(uniq_cad)
kfcvpartitionsetD <- createFolds(y = 1:length(uniq_cad),## the outcome data are needed
                                k = 10, 
                                list = TRUE)

# perform k-fold-out
for(k in 1:10){  # 1:10f cv
  ## Create folds leave-one-patient-out
  allfT1 = read3Dtex_T1uniqcad_parti(id_cad_pts, uniq_cad, kfcvpartitionsetD, 10, k)
  allfT1T2 = read3Dtex_T1T2uniqcad_parti(id_cad_pts, uniq_cad, kfcvpartitionsetD, 10, k)
    
  ## formant
  T1train = allfT1[[1]]; T1traininfo = allfT1[[5]]; 
  T1test = allfT1[[2]]; T1testinfo = allfT1[[6]];  
  
  T1T2train = allfT1T2[[1]]; T1T2traininfo = allfT1T2[[5]]; 
  T1T2test = allfT1T2[[2]]; T1T2testinfo = allfT1T2[[6]];  
  
  # remove measured muscle-to-lesion SI / to add predicted
  # add predicted T2w features
  trainpredT2LMSIR = getid_predLMSIR(LMSIR_cv, T1T2traininfo$lesion_id)
  testpredT2LMSIR = getid_predLMSIR(LMSIR_cv, T1T2testinfo$lesion_id)
  # add to test set
  T1T2test = cbind(T1T2test, testpredT2LMSIR[1])
  
  ##################
  # Define datasets
  ##################
  imgT1train = T1train[,-c(ncol(T1train))]
  imgT1T2 = T1T2train[,-c(201,ncol(T1T2train))] #"T1w -T2wSI_measureLMSIR"
  imgT2pLMSIR = cbind(imgT1T2[,-c(199:200,222:241)], trainpredT2LMSIR[1]) #"T1w+T2wSI_predictedLMSIR"
  
  ##################
  # Get Test info data
  ##################   
  dfinfo = cbind(T1T2testinfo[,c(1,24:26)])
  pander(summary(as.factor(dfinfo$lesion_label)))
  
  ##################   
  # Build final classifiers
  ##################
  cat("\n============ boosting trees imgT1train \n")
  # train trees
  treedata_imgT1 <- boosting_Train_wcv_wperf(imgT1train, T1test, typeSampling="oversmin")

  ######## 
  cat("\n============ boosting trees imgT1T2 \n")
  # train trees
  treedata_imgT1T2 <- boosting_Train_wcv_wperf(imgT1T2, T1T2test[c(names(imgT1T2))], typeSampling="oversmin")

  #######  
  cat("\n============ boosting trees imgT2pLMSIR \n")
  # train trees
  treedata_T2wpLMSIR <- boosting_Train_wcv_wperf(imgT2pLMSIR, T1T2test[c(names(imgT2pLMSIR))], typeSampling="oversmin")

  ##################
  ### predict for each classifier
  ##################
  ## for treedata_imgT1
  perf_imgT1 = data.frame(C=treedata_imgT1$testperf$testpred[,1],
                    NC=treedata_imgT1$testperf$testpred[,2],
                    pred=ifelse(apply(treedata_imgT1$testperf$testpred, 1, which.max)==1,"C","NC"), 
                    obs=treedata_imgT1$testperf$labelstest,
                    id=T1testinfo$lesion_id)
  
  # for treedata_imgT1T2
  perf_imgT1T2 = data.frame(C=treedata_imgT1T2$testperf$testpred[,1],
                    NC=treedata_imgT1T2$testperf$testpred[,2],
                    pred=ifelse(apply(treedata_imgT1T2$testperf$testpred, 1, which.max)==1,"C","NC"),
                    obs=treedata_imgT1T2$testperf$labelstest,
                    id=T1T2testinfo$lesion_id) 
  
  # for treedata_T2wpLMSIR
  perf_T2wpLMSIR = data.frame(C=treedata_T2wpLMSIR$testperf$testpred[,1],
                    NC=treedata_T2wpLMSIR$testperf$testpred[,2],
                    pred=ifelse(apply(treedata_T2wpLMSIR$testperf$testpred, 1, which.max)==1,"C","NC"),
                    obs=treedata_T2wpLMSIR$testperf$labelstest,
                    id=T1T2testinfo$lesion_id) 
  
  # AUC
  rocperf_imgT1 = roc(perf_imgT1$obs, perf_imgT1$C)
  print(rocperf_imgT1)

  rocperf_imgT1T2 = roc(perf_imgT1T2$obs, perf_imgT1T2$C)
  print(rocperf_imgT1T2)
  
  rocperf_T2wpLMSIR = roc(perf_T2wpLMSIR$obs, perf_T2wpLMSIR$C)
  print(rocperf_T2wpLMSIR)
   
  # plot every 10 patients
  ## plot ROCs each pass individually in l-o-p heldout test cases
  par(mfrow=c(1,1))
  n=15
  colors = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
  # plot 1/4
  p1 = calcAUC_plot(perf_imgT1$obs, perf_imgT1$C, 
                             xptext=0.45, yptext=0.75 , 1, colors[2], atitle="")
  par(new=TRUE)
  p2 = calcAUC_plot(perf_imgT1T2$obs, perf_imgT1T2$C, 
                             xptext=0.55, yptext=0.65, 2, colors[9], atitle="")
  par(new=TRUE)
  p3 = calcAUC_plot(perf_T2wpLMSIR$obs, perf_T2wpLMSIR$C,
                             xptext=0.65, yptext=0.55, 3, colors[11], 
                    atitle=paste0("ROCs 10f-patient out cv test k-fold= ",k))
  
  legend("bottomright", 
         legend = c(paste0("T1w"),
                    paste0("T1w+T2wtex+T2wSI"),
                    paste0("T1w+T2wtex+predLMSIR")),
         col = c(colors[2],colors[9],colors[11]), lty=c(1,2,3), lwd = 2)

    
  # save current state k patient out
  save.image(paste0("Outputs/T1T2imgvsT2wtextpredLMSIR_boost_addeddiagvalue_extratrees_cv",k,".RData"))

}
```

```
   massB    massM nonmassB nonmassM 
     216      152      131       68 
   massB    massM nonmassB nonmassM 
      24       15       11       10 
   massB    massM nonmassB nonmassM 
     216      152      131       68 
   massB    massM nonmassB nonmassM 
      24       15       11       10 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.6674286
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7405714
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7005714
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7131429
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7131429
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.6914286
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7062857
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.6982857
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.6948571
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7074286
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain rocTest
11        5    250        1    0.72
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7371429
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7142857
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7348571
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7268571
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.6674286
2         5     50        1 0.7405714
3        15     50        1 0.7005714
4         3     75        1 0.7131429
5         5     75        1 0.7131429
6        15     75        1 0.6914286
7         3    100        1 0.7062857
8         5    100        1 0.6982857
9        15    100        1 0.6948571
10        3    250        1 0.7074286
11        5    250        1 0.7200000
12       15    250        1 0.7371429
13        3    750        1 0.7142857
14        5    750        1 0.7348571
15       15    750        1 0.7268571
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7405714

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.6205714
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6765714
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.6857143
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7405714
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7222857
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
6       15     75        1   0.776
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.6857143
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7337143
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7451429
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain rocTest
10        3    250        1   0.744
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7165714
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7588571
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7485714
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7714286
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7577143
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.6205714
2         5     50        1 0.6765714
3        15     50        1 0.6857143
4         3     75        1 0.7405714
5         5     75        1 0.7222857
6        15     75        1 0.7760000
7         3    100        1 0.6857143
8         5    100        1 0.7337143
9        15    100        1 0.7451429
10        3    250        1 0.7440000
11        5    250        1 0.7165714
12       15    250        1 0.7588571
13        3    750        1 0.7485714
14        5    750        1 0.7714286
15       15    750        1 0.7577143
  maxdepth mfinal rocTrain rocTest
6       15     75        1   0.776

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7097143
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7314286
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7211429
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
4        3     75        1   0.648
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7051429
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7691429
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7154286
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7371429
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7611429
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7177143
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7394286
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7371429
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7268571
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7622857
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7417143
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7097143
2         5     50        1 0.7314286
3        15     50        1 0.7211429
4         3     75        1 0.6480000
5         5     75        1 0.7051429
6        15     75        1 0.7691429
7         3    100        1 0.7154286
8         5    100        1 0.7371429
9        15    100        1 0.7611429
10        3    250        1 0.7177143
11        5    250        1 0.7394286
12       15    250        1 0.7371429
13        3    750        1 0.7268571
14        5    750        1 0.7622857
15       15    750        1 0.7417143
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7691429

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 25 controls (perf_imgT1$obs C) > 35 cases (perf_imgT1$obs NC).
Area under the curve: 0.7406

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 25 controls (perf_imgT1T2$obs C) > 35 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.776

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 25 controls (perf_T2wpLMSIR$obs C) > 35 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.7691
```

```
Area under the curve: 0.7406
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4981562   0.44      0.64     0.8 0.6857    0.8286  0.9429
```

```
Area under the curve: 0.776
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4592056   0.64       0.8    0.96 0.5429    0.7143  0.8571
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-1.png)<!-- -->

```
Area under the curve: 0.7691
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5006342   0.56      0.72    0.88 0.6857    0.8286  0.9429
   massB    massM nonmassB nonmassM 
     210      150      132       64 
   massB    massM nonmassB nonmassM 
      30       17       10       14 
   massB    massM nonmassB nonmassM 
     210      150      132       64 
   massB    massM nonmassB nonmassM 
      30       17       10       14 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7790323
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8241935
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8427419
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7919355
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8209677
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8379032
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8120968
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7870968
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8040323
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8217742
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8427419
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8241935
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8217742
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8435484
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8258065
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7790323
2         5     50        1 0.8241935
3        15     50        1 0.8427419
4         3     75        1 0.7919355
5         5     75        1 0.8209677
6        15     75        1 0.8379032
7         3    100        1 0.8120968
8         5    100        1 0.7870968
9        15    100        1 0.8040323
10        3    250        1 0.8217742
11        5    250        1 0.8427419
12       15    250        1 0.8241935
13        3    750        1 0.8217742
14        5    750        1 0.8435484
15       15    750        1 0.8258065
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8435484

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.8225806
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8193548
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7975806
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7983871
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8427419
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8387097
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8354839
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8596774
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8209677
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8596774
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8637097
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8645161
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8798387
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8717742
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8516129
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.8225806
2         5     50        1 0.8193548
3        15     50        1 0.7975806
4         3     75        1 0.7983871
5         5     75        1 0.8427419
6        15     75        1 0.8387097
7         3    100        1 0.8354839
8         5    100        1 0.8596774
9        15    100        1 0.8209677
10        3    250        1 0.8596774
11        5    250        1 0.8637097
12       15    250        1 0.8645161
13        3    750        1 0.8798387
14        5    750        1 0.8717742
15       15    750        1 0.8516129
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8798387

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7120968
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8048387
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7830645
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.8032258
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8322581
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8516129
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8177419
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7870968
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8403226
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8217742
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8129032
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8322581
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8225806
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8435484
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8282258
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7120968
2         5     50        1 0.8048387
3        15     50        1 0.7830645
4         3     75        1 0.8032258
5         5     75        1 0.8322581
6        15     75        1 0.8516129
7         3    100        1 0.8177419
8         5    100        1 0.7870968
9        15    100        1 0.8403226
10        3    250        1 0.8217742
11        5    250        1 0.8129032
12       15    250        1 0.8322581
13        3    750        1 0.8225806
14        5    750        1 0.8435484
15       15    750        1 0.8282258
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8516129

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 31 controls (perf_imgT1$obs C) > 40 cases (perf_imgT1$obs NC).
Area under the curve: 0.8435

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 31 controls (perf_imgT1T2$obs C) > 40 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8798

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 31 controls (perf_T2wpLMSIR$obs C) > 40 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8516
```

```
Area under the curve: 0.8435
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4944509 0.5484    0.7097   0.871  0.825     0.925       1
```

```
Area under the curve: 0.8798
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4488463 0.7742    0.9032       1  0.575     0.725    0.85
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-2.png)<!-- -->

```
Area under the curve: 0.8516
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5006417 0.5161    0.6774  0.8387  0.925     0.975       1
   massB    massM nonmassB nonmassM 
     217      154      129       72 
   massB    massM nonmassB nonmassM 
      23       13       13        6 
   massB    massM nonmassB nonmassM 
     217      154      129       72 
   massB    massM nonmassB nonmassM 
      23       13       13        6 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7046784
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6739766
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.6988304
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.6944444
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7280702
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7339181
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.6432749
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.6988304
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.6666667
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7280702
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
11        5    250        1 0.747076
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7076023
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7134503
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.6900585
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.6988304
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7046784
2         5     50        1 0.6739766
3        15     50        1 0.6988304
4         3     75        1 0.6944444
5         5     75        1 0.7280702
6        15     75        1 0.7339181
7         3    100        1 0.6432749
8         5    100        1 0.6988304
9        15    100        1 0.6666667
10        3    250        1 0.7280702
11        5    250        1 0.7470760
12       15    250        1 0.7076023
13        3    750        1 0.7134503
14        5    750        1 0.6900585
15       15    750        1 0.6988304
   maxdepth mfinal rocTrain  rocTest
11        5    250        1 0.747076

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7076023
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6622807
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7353801
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.6798246
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7763158
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7690058
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7807018
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7412281
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7046784
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7368421
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7222222
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7309942
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7573099
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7192982
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7368421
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7076023
2         5     50        1 0.6622807
3        15     50        1 0.7353801
4         3     75        1 0.6798246
5         5     75        1 0.7763158
6        15     75        1 0.7690058
7         3    100        1 0.7807018
8         5    100        1 0.7412281
9        15    100        1 0.7046784
10        3    250        1 0.7368421
11        5    250        1 0.7222222
12       15    250        1 0.7309942
13        3    750        1 0.7573099
14        5    750        1 0.7192982
15       15    750        1 0.7368421
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7807018

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.6988304
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7178363
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.6666667
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7207602
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7967836
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7573099
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7222222
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7002924
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8011696
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
10        3    250        1 0.752924
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8245614
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
12       15    250        1 0.748538
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7719298
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7339181
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
15       15    750        1 0.755848
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.6988304
2         5     50        1 0.7178363
3        15     50        1 0.6666667
4         3     75        1 0.7207602
5         5     75        1 0.7967836
6        15     75        1 0.7573099
7         3    100        1 0.7222222
8         5    100        1 0.7002924
9        15    100        1 0.8011696
10        3    250        1 0.7529240
11        5    250        1 0.8245614
12       15    250        1 0.7485380
13        3    750        1 0.7719298
14        5    750        1 0.7339181
15       15    750        1 0.7558480
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8245614

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 19 controls (perf_imgT1$obs C) > 36 cases (perf_imgT1$obs NC).
Area under the curve: 0.7471

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 19 controls (perf_imgT1T2$obs C) > 36 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.7807

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 19 controls (perf_T2wpLMSIR$obs C) > 36 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8246
```

```
Area under the curve: 0.7471
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4452338 0.5263    0.7368  0.9474 0.5833    0.7222  0.8611
```

```
Area under the curve: 0.7807
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5383994 0.2632    0.4737  0.6842 0.9167    0.9722       1
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-3.png)<!-- -->

```
Area under the curve: 0.8246
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.491993 0.3684    0.5789  0.7895 0.8611    0.9444       1
   massB    massM nonmassB nonmassM 
     224      151      128       70 
   massB    massM nonmassB nonmassM 
      16       16       14        8 
   massB    massM nonmassB nonmassM 
     224      151      128       70 
   massB    massM nonmassB nonmassM 
      16       16       14        8 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7166667
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8097222
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7972222
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7486111
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7152778
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7611111
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7236111
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7972222
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7166667
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7458333
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8152778
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7888889
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7680556
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7958333
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7861111
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7166667
2         5     50        1 0.8097222
3        15     50        1 0.7972222
4         3     75        1 0.7486111
5         5     75        1 0.7152778
6        15     75        1 0.7611111
7         3    100        1 0.7236111
8         5    100        1 0.7972222
9        15    100        1 0.7166667
10        3    250        1 0.7458333
11        5    250        1 0.8152778
12       15    250        1 0.7888889
13        3    750        1 0.7680556
14        5    750        1 0.7958333
15       15    750        1 0.7861111
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8152778

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7944444
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7527778
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8361111
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7861111
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8694444
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7083333
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7930556
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7805556
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8013889
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8222222
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7888889
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7652778
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8111111
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8180556
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7958333
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7944444
2         5     50        1 0.7527778
3        15     50        1 0.8361111
4         3     75        1 0.7861111
5         5     75        1 0.8694444
6        15     75        1 0.7083333
7         3    100        1 0.7930556
8         5    100        1 0.7805556
9        15    100        1 0.8013889
10        3    250        1 0.8222222
11        5    250        1 0.7888889
12       15    250        1 0.7652778
13        3    750        1 0.8111111
14        5    750        1 0.8180556
15       15    750        1 0.7958333
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8694444

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.8166667
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7986111
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain rocTest
3       15     50        1  0.7125
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7833333
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8458333
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7902778
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7791667
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7944444
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8333333
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7805556
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain rocTest
11        5    250        1   0.825
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8041667
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8097222
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8055556
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain rocTest
15       15    750        1     0.8
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.8166667
2         5     50        1 0.7986111
3        15     50        1 0.7125000
4         3     75        1 0.7833333
5         5     75        1 0.8458333
6        15     75        1 0.7902778
7         3    100        1 0.7791667
8         5    100        1 0.7944444
9        15    100        1 0.8333333
10        3    250        1 0.7805556
11        5    250        1 0.8250000
12       15    250        1 0.8041667
13        3    750        1 0.8097222
14        5    750        1 0.8055556
15       15    750        1 0.8000000
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8458333

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 24 controls (perf_imgT1$obs C) > 30 cases (perf_imgT1$obs NC).
Area under the curve: 0.8153

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 24 controls (perf_imgT1T2$obs C) > 30 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8694

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 24 controls (perf_T2wpLMSIR$obs C) > 30 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8458
```

```
Area under the curve: 0.8153
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.419474  0.625    0.7917  0.9583 0.5667    0.7333  0.8667
```

```
Area under the curve: 0.8694
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4198863   0.75     0.875       1 0.6667       0.8  0.9333
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-4.png)<!-- -->

```
Area under the curve: 0.8458
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.456928 0.7083     0.875       1 0.6667       0.8  0.9333
   massB    massM nonmassB nonmassM 
     223      148      127       70 
   massB    massM nonmassB nonmassM 
      17       19       15        8 
   massB    massM nonmassB nonmassM 
     223      148      127       70 
   massB    massM nonmassB nonmassM 
      17       19       15        8 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7002315
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7337963
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain  rocTest
3       15     50        1 0.787037
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7453704
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7997685
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7581019
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7222222
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7893519
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain rocTest
9       15    100        1 0.78125
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7824074
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7662037
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7997685
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7962963
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
14        5    750        1 0.806713
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7939815
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7002315
2         5     50        1 0.7337963
3        15     50        1 0.7870370
4         3     75        1 0.7453704
5         5     75        1 0.7997685
6        15     75        1 0.7581019
7         3    100        1 0.7222222
8         5    100        1 0.7893519
9        15    100        1 0.7812500
10        3    250        1 0.7824074
11        5    250        1 0.7662037
12       15    250        1 0.7997685
13        3    750        1 0.7962963
14        5    750        1 0.8067130
15       15    750        1 0.7939815
   maxdepth mfinal rocTrain  rocTest
14        5    750        1 0.806713

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7858796
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7766204
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7847222
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.6990741
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7581019
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7650463
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7766204
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7337963
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8240741
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7835648
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7800926
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7314815
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8078704
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
14        5    750        1 0.787037
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8055556
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7858796
2         5     50        1 0.7766204
3        15     50        1 0.7847222
4         3     75        1 0.6990741
5         5     75        1 0.7581019
6        15     75        1 0.7650463
7         3    100        1 0.7766204
8         5    100        1 0.7337963
9        15    100        1 0.8240741
10        3    250        1 0.7835648
11        5    250        1 0.7800926
12       15    250        1 0.7314815
13        3    750        1 0.8078704
14        5    750        1 0.7870370
15       15    750        1 0.8055556
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8240741

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain  rocTest
1        3     50        1 0.744213
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7662037
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7361111
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7407407
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7152778
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7673611
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7430556
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7511574
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7604167
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7916667
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7766204
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7777778
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7662037
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7847222
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7708333
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7442130
2         5     50        1 0.7662037
3        15     50        1 0.7361111
4         3     75        1 0.7407407
5         5     75        1 0.7152778
6        15     75        1 0.7673611
7         3    100        1 0.7430556
8         5    100        1 0.7511574
9        15    100        1 0.7604167
10        3    250        1 0.7916667
11        5    250        1 0.7766204
12       15    250        1 0.7777778
13        3    750        1 0.7662037
14        5    750        1 0.7847222
15       15    750        1 0.7708333
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7916667

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 27 controls (perf_imgT1$obs C) > 32 cases (perf_imgT1$obs NC).
Area under the curve: 0.8067

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 27 controls (perf_imgT1T2$obs C) > 32 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8241

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 27 controls (perf_T2wpLMSIR$obs C) > 32 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.7917
```

```
Area under the curve: 0.8067
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4105882 0.7037    0.8519   0.963 0.4688    0.6562  0.8125
```

```
Area under the curve: 0.8241
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4138349 0.7037    0.8519   0.963 0.5625    0.7188   0.875
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-5.png)<!-- -->

```
Area under the curve: 0.7917
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4624332 0.5185    0.7037  0.8528  0.625    0.7812  0.9062
   massB    massM nonmassB nonmassM 
     209      153      126       72 
   massB    massM nonmassB nonmassM 
      31       14       16        6 
   massB    massM nonmassB nonmassM 
     209      153      126       72 
   massB    massM nonmassB nonmassM 
      31       14       16        6 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7989362
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8329787
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8425532
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7148936
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8287234
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7765957
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7680851
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8457447
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8117021
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8340426
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8521277
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8138298
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
13        3    750        1 0.837234
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8255319
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
15       15    750        1 0.837234
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7989362
2         5     50        1 0.8329787
3        15     50        1 0.8425532
4         3     75        1 0.7148936
5         5     75        1 0.8287234
6        15     75        1 0.7765957
7         3    100        1 0.7680851
8         5    100        1 0.8457447
9        15    100        1 0.8117021
10        3    250        1 0.8340426
11        5    250        1 0.8521277
12       15    250        1 0.8138298
13        3    750        1 0.8372340
14        5    750        1 0.8255319
15       15    750        1 0.8372340
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8521277

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7478723
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7691489
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7808511
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7340426
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7553191
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8095745
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8117021
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8297872
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8468085
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8053191
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8574468
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8489362
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8340426
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8202128
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8457447
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7478723
2         5     50        1 0.7691489
3        15     50        1 0.7808511
4         3     75        1 0.7340426
5         5     75        1 0.7553191
6        15     75        1 0.8095745
7         3    100        1 0.8117021
8         5    100        1 0.8297872
9        15    100        1 0.8468085
10        3    250        1 0.8053191
11        5    250        1 0.8574468
12       15    250        1 0.8489362
13        3    750        1 0.8340426
14        5    750        1 0.8202128
15       15    750        1 0.8457447
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8574468

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7223404
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7617021
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7574468
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7489362
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7776596
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7659574
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7968085
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7946809
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7702128
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7308511
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7595745
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8095745
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7882979
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7861702
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7914894
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7223404
2         5     50        1 0.7617021
3        15     50        1 0.7574468
4         3     75        1 0.7489362
5         5     75        1 0.7776596
6        15     75        1 0.7659574
7         3    100        1 0.7968085
8         5    100        1 0.7946809
9        15    100        1 0.7702128
10        3    250        1 0.7308511
11        5    250        1 0.7595745
12       15    250        1 0.8095745
13        3    750        1 0.7882979
14        5    750        1 0.7861702
15       15    750        1 0.7914894
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8095745

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 20 controls (perf_imgT1$obs C) > 47 cases (perf_imgT1$obs NC).
Area under the curve: 0.8521

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 20 controls (perf_imgT1T2$obs C) > 47 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8574

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 20 controls (perf_T2wpLMSIR$obs C) > 47 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8096
```

```
Area under the curve: 0.8521
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4521919    0.6       0.8    0.95 0.6596    0.7872  0.8936
```

```
Area under the curve: 0.8574
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4391238   0.85      0.95       1 0.5745    0.7021  0.8298
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-6.png)<!-- -->

```
Area under the curve: 0.8096
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3897094   0.75       0.9       1 0.4681     0.617  0.7447
   massB    massM nonmassB nonmassM 
     222      148      121       73 
   massB    massM nonmassB nonmassM 
      18       19       21        5 
   massB    massM nonmassB nonmassM 
     222      148      121       73 
   massB    massM nonmassB nonmassM 
      18       19       21        5 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7083333
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6816239
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7617521
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7403846
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7061966
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
6       15     75        1    0.75
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8023504
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7371795
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7446581
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7393162
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
11        5    250        1 0.741453
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7478632
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
13        3    750        1 0.732906
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7702991
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7574786
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7083333
2         5     50        1 0.6816239
3        15     50        1 0.7617521
4         3     75        1 0.7403846
5         5     75        1 0.7061966
6        15     75        1 0.7500000
7         3    100        1 0.8023504
8         5    100        1 0.7371795
9        15    100        1 0.7446581
10        3    250        1 0.7393162
11        5    250        1 0.7414530
12       15    250        1 0.7478632
13        3    750        1 0.7329060
14        5    750        1 0.7702991
15       15    750        1 0.7574786
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8023504

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7820513
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6741453
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8044872
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain  rocTest
4        3     75        1 0.767094
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.6987179
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7916667
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7425214
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7361111
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7382479
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7542735
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7596154
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7403846
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain rocTest
13        3    750        1    0.75
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7510684
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7905983
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7820513
2         5     50        1 0.6741453
3        15     50        1 0.8044872
4         3     75        1 0.7670940
5         5     75        1 0.6987179
6        15     75        1 0.7916667
7         3    100        1 0.7425214
8         5    100        1 0.7361111
9        15    100        1 0.7382479
10        3    250        1 0.7542735
11        5    250        1 0.7596154
12       15    250        1 0.7403846
13        3    750        1 0.7500000
14        5    750        1 0.7510684
15       15    750        1 0.7905983
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8044872

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7115385
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7457265
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7532051
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7638889
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7735043
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
6       15     75        1    0.75
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7094017
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7532051
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7457265
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8023504
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7991453
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
12       15    250        1 0.758547
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7371795
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7863248
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7767094
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7115385
2         5     50        1 0.7457265
3        15     50        1 0.7532051
4         3     75        1 0.7638889
5         5     75        1 0.7735043
6        15     75        1 0.7500000
7         3    100        1 0.7094017
8         5    100        1 0.7532051
9        15    100        1 0.7457265
10        3    250        1 0.8023504
11        5    250        1 0.7991453
12       15    250        1 0.7585470
13        3    750        1 0.7371795
14        5    750        1 0.7863248
15       15    750        1 0.7767094
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8023504

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 24 controls (perf_imgT1$obs C) > 39 cases (perf_imgT1$obs NC).
Area under the curve: 0.8024

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 24 controls (perf_imgT1T2$obs C) > 39 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8045

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 24 controls (perf_T2wpLMSIR$obs C) > 39 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8024
```

```
Area under the curve: 0.8024
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4567265 0.6667    0.8333  0.9583 0.5897    0.7436  0.8718
```

```
Area under the curve: 0.8045
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.437486 0.5417    0.7083   0.875 0.6923    0.8205  0.9231
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-7.png)<!-- -->

```
Area under the curve: 0.8024
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4756004    0.5    0.7083   0.875 0.6923    0.8205  0.9231
   massB    massM nonmassB nonmassM 
     212      143      128       77 
   massB    massM nonmassB nonmassM 
      28       24       14        1 
   massB    massM nonmassB nonmassM 
     212      143      128       77 
   massB    massM nonmassB nonmassM 
      28       24       14        1 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.8142857
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7380952
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8409524
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.8304762
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain  rocTest
5        5     75        1 0.767619
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.7752381
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7257143
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain  rocTest
8        5    100        1 0.832381
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7780952
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7542857
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7761905
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7961905
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7828571
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.8219048
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8038095
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.8142857
2         5     50        1 0.7380952
3        15     50        1 0.8409524
4         3     75        1 0.8304762
5         5     75        1 0.7676190
6        15     75        1 0.7752381
7         3    100        1 0.7257143
8         5    100        1 0.8323810
9        15    100        1 0.7780952
10        3    250        1 0.7542857
11        5    250        1 0.7761905
12       15    250        1 0.7961905
13        3    750        1 0.7828571
14        5    750        1 0.8219048
15       15    750        1 0.8038095
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8409524

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7114286
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.5885714
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain  rocTest
3       15     50        1 0.712381
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain  rocTest
4        3     75        1 0.792381
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.7380952
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain  rocTest
6       15     75        1 0.752381
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8038095
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain  rocTest
8        5    100        1 0.687619
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7104762
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain  rocTest
10        3    250        1 0.732381
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7419048
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7942857
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7733333
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7819048
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
15       15    750        1 0.752381
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7114286
2         5     50        1 0.5885714
3        15     50        1 0.7123810
4         3     75        1 0.7923810
5         5     75        1 0.7380952
6        15     75        1 0.7523810
7         3    100        1 0.8038095
8         5    100        1 0.6876190
9        15    100        1 0.7104762
10        3    250        1 0.7323810
11        5    250        1 0.7419048
12       15    250        1 0.7942857
13        3    750        1 0.7733333
14        5    750        1 0.7819048
15       15    750        1 0.7523810
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8038095

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7733333
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7371429
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7790476
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.7495238
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8009524
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8561905
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7561905
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7933333
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7495238
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8542857
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8104762
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8085714
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
13        3    750        1 0.852381
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7971429
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.8038095
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7733333
2         5     50        1 0.7371429
3        15     50        1 0.7790476
4         3     75        1 0.7495238
5         5     75        1 0.8009524
6        15     75        1 0.8561905
7         3    100        1 0.7561905
8         5    100        1 0.7933333
9        15    100        1 0.7495238
10        3    250        1 0.8542857
11        5    250        1 0.8104762
12       15    250        1 0.8085714
13        3    750        1 0.8523810
14        5    750        1 0.7971429
15       15    750        1 0.8038095
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8561905

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 25 controls (perf_imgT1$obs C) > 42 cases (perf_imgT1$obs NC).
Area under the curve: 0.841

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 25 controls (perf_imgT1T2$obs C) > 42 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.8038

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 25 controls (perf_T2wpLMSIR$obs C) > 42 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.8562
```

```
Area under the curve: 0.841
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3985963   0.72      0.88       1  0.619    0.7619   0.881
```

```
Area under the curve: 0.8038
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4320352    0.8      0.92       1 0.4762     0.619  0.7619
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-8.png)<!-- -->

```
Area under the curve: 0.8562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3802358   0.72      0.88       1 0.5714    0.7143  0.8339
   massB    massM nonmassB nonmassM 
     209      154      125       71 
   massB    massM nonmassB nonmassM 
      31       13       17        7 
   massB    massM nonmassB nonmassM 
     209      154      125       71 
   massB    massM nonmassB nonmassM 
      31       13       17        7 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain rocTest
1        3     50        1  0.6375
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.6895833
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.6385417
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.6291667
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.6697917
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
6       15     75        1 0.65625
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.6802083
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain  rocTest
8        5    100        1 0.653125
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain rocTest
9       15    100        1   0.625
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.6479167
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7041667
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7385417
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
13        3    750        1 0.665625
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain rocTest
14        5    750        1 0.70625
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7104167
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.6375000
2         5     50        1 0.6895833
3        15     50        1 0.6385417
4         3     75        1 0.6291667
5         5     75        1 0.6697917
6        15     75        1 0.6562500
7         3    100        1 0.6802083
8         5    100        1 0.6531250
9        15    100        1 0.6250000
10        3    250        1 0.6479167
11        5    250        1 0.7041667
12       15    250        1 0.7385417
13        3    750        1 0.6656250
14        5    750        1 0.7062500
15       15    750        1 0.7104167
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7385417

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.6854167
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain  rocTest
2        5     50        1 0.740625
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7260417
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
4        3     75        1 0.63125
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.6979167
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
6       15     75        1 0.73125
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.6541667
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.7416667
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.7177083
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.6760417
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7302083
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain rocTest
12       15    250        1 0.75625
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7166667
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain rocTest
14        5    750        1 0.74375
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7104167
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.6854167
2         5     50        1 0.7406250
3        15     50        1 0.7260417
4         3     75        1 0.6312500
5         5     75        1 0.6979167
6        15     75        1 0.7312500
7         3    100        1 0.6541667
8         5    100        1 0.7416667
9        15    100        1 0.7177083
10        3    250        1 0.6760417
11        5    250        1 0.7302083
12       15    250        1 0.7562500
13        3    750        1 0.7166667
14        5    750        1 0.7437500
15       15    750        1 0.7104167
   maxdepth mfinal rocTrain rocTest
12       15    250        1 0.75625

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.7135417
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.7604167
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.7333333
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain rocTest
4        3     75        1 0.75625
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.6677083
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.6885417
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.7354167
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain rocTest
8        5    100        1 0.73125
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.6947917
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.7322917
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.7729167
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.7572917
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.7208333
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.7510417
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7802083
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.7135417
2         5     50        1 0.7604167
3        15     50        1 0.7333333
4         3     75        1 0.7562500
5         5     75        1 0.6677083
6        15     75        1 0.6885417
7         3    100        1 0.7354167
8         5    100        1 0.7312500
9        15    100        1 0.6947917
10        3    250        1 0.7322917
11        5    250        1 0.7729167
12       15    250        1 0.7572917
13        3    750        1 0.7208333
14        5    750        1 0.7510417
15       15    750        1 0.7802083
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.7802083

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 20 controls (perf_imgT1$obs C) > 48 cases (perf_imgT1$obs NC).
Area under the curve: 0.7385

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 20 controls (perf_imgT1T2$obs C) > 48 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.7562

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 20 controls (perf_T2wpLMSIR$obs C) > 48 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.7802
```

```
Area under the curve: 0.7385
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4211708   0.65      0.85       1 0.4792     0.625    0.75
```

```
Area under the curve: 0.7562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4596946   0.55      0.75     0.9 0.6458    0.7708   0.875
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-9.png)<!-- -->

```
Area under the curve: 0.7802
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4135362   0.75       0.9       1 0.4375    0.5833  0.7292
   massB    massM nonmassB nonmassM 
     218      150      131       65 
   massB    massM nonmassB nonmassM 
      22       17       11       13 
   massB    massM nonmassB nonmassM 
     218      150      131       65 
   massB    massM nonmassB nonmassM 
      22       17       11       13 

============ boosting trees imgT1train 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain  rocTest
1        3     50        1 0.859596
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8878788
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.9131313
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.8727273
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8878788
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain  rocTest
6       15     75        1 0.910101
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8313131
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8727273
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8676768
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8737374
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.9080808
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.8959596
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8808081
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
14        5    750        1 0.9070707
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.9020202
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.8595960
2         5     50        1 0.8878788
3        15     50        1 0.9131313
4         3     75        1 0.8727273
5         5     75        1 0.8878788
6        15     75        1 0.9101010
7         3    100        1 0.8313131
8         5    100        1 0.8727273
9        15    100        1 0.8676768
10        3    250        1 0.8737374
11        5    250        1 0.9080808
12       15    250        1 0.8959596
13        3    750        1 0.8808081
14        5    750        1 0.9070707
15       15    750        1 0.9020202
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.9131313

============ boosting trees imgT1T2 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.9060606
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8989899
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
3       15     50        1 0.8656566
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.8575758
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8606061
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8666667
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8838384
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8989899
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8838384
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8969697
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.9020202
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.9040404
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.8858586
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain rocTest
14        5    750        1     0.9
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.9050505
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.9060606
2         5     50        1 0.8989899
3        15     50        1 0.8656566
4         3     75        1 0.8575758
5         5     75        1 0.8606061
6        15     75        1 0.8666667
7         3    100        1 0.8838384
8         5    100        1 0.8989899
9        15    100        1 0.8838384
10        3    250        1 0.8969697
11        5    250        1 0.9020202
12       15    250        1 0.9040404
13        3    750        1 0.8858586
14        5    750        1 0.9000000
15       15    750        1 0.9050505
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.9060606

============ boosting trees imgT2pLMSIR 
maxdepth:  3 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
1        3     50        1 0.8585859
maxdepth:  5 mfinal:  50 
  maxdepth mfinal rocTrain   rocTest
2        5     50        1 0.8767677
maxdepth:  15 mfinal:  50 
  maxdepth mfinal rocTrain rocTest
3       15     50        1     0.9
maxdepth:  3 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
4        3     75        1 0.8979798
maxdepth:  5 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
5        5     75        1 0.8585859
maxdepth:  15 mfinal:  75 
  maxdepth mfinal rocTrain   rocTest
6       15     75        1 0.8434343
maxdepth:  3 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
7        3    100        1 0.8707071
maxdepth:  5 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
8        5    100        1 0.8888889
maxdepth:  15 mfinal:  100 
  maxdepth mfinal rocTrain   rocTest
9       15    100        1 0.8929293
maxdepth:  3 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
10        3    250        1 0.8808081
maxdepth:  5 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
11        5    250        1 0.8828283
maxdepth:  15 mfinal:  250 
   maxdepth mfinal rocTrain   rocTest
12       15    250        1 0.9080808
maxdepth:  3 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
13        3    750        1 0.9050505
maxdepth:  5 mfinal:  750 
   maxdepth mfinal rocTrain  rocTest
14        5    750        1 0.889899
maxdepth:  15 mfinal:  750 
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.9262626
   maxdepth mfinal rocTrain   rocTest
1         3     50        1 0.8585859
2         5     50        1 0.8767677
3        15     50        1 0.9000000
4         3     75        1 0.8979798
5         5     75        1 0.8585859
6        15     75        1 0.8434343
7         3    100        1 0.8707071
8         5    100        1 0.8888889
9        15    100        1 0.8929293
10        3    250        1 0.8808081
11        5    250        1 0.8828283
12       15    250        1 0.9080808
13        3    750        1 0.9050505
14        5    750        1 0.8898990
15       15    750        1 0.9262626
   maxdepth mfinal rocTrain   rocTest
15       15    750        1 0.9262626

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 30 controls (perf_imgT1$obs C) > 33 cases (perf_imgT1$obs NC).
Area under the curve: 0.9131

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 30 controls (perf_imgT1T2$obs C) > 33 cases (perf_imgT1T2$obs NC).
Area under the curve: 0.9061

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 30 controls (perf_T2wpLMSIR$obs C) > 33 cases (perf_T2wpLMSIR$obs NC).
Area under the curve: 0.9263
```

```
Area under the curve: 0.9131
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4096503 0.6333       0.8  0.9333 0.7879    0.9091       1
```

```
Area under the curve: 0.9061
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4421673 0.7667       0.9       1 0.6667    0.8182  0.9394
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_extratrees_files/figure-html/run-boosting-10.png)<!-- -->

```
Area under the curve: 0.9263
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4587403 0.6333       0.8  0.9333 0.8485    0.9394       1
```


