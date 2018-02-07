# Boosting classification trees - T2SI vs predLMSIR 
Cristina Gallego  
August 28, 2017  

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
  save.image(paste0("Outputs/T1T2imgvsT2wtextpredLMSIR_boost_addeddiagvalue_cv",k,".RData"))

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
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9569841      0.6457143 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9569841 0.6457143
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.6788571 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.6788571
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7005714 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7005714
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9661570      0.6342857 
  maxdepth mfinal rocTrain   rocTest
4        1    100 0.966157 0.6342857
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.6845714 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.6845714
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7051429 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7051429
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9817705      0.6445714 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9817705 0.6445714
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.6834286 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.6834286
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7451429 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7451429
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9870192      0.6445714 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9870192 0.6445714
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7017143 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7017143
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7554286 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7554286
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6422857 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6422857
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7165714 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7165714
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7508571 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7508571
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9569841 0.6457143
2         3     75 1.0000000 0.6788571
3         6     75 1.0000000 0.7005714
4         1    100 0.9661570 0.6342857
5         3    100 1.0000000 0.6845714
6         6    100 1.0000000 0.7051429
7         1    150 0.9817705 0.6445714
8         3    150 1.0000000 0.6834286
9         6    150 1.0000000 0.7451429
10        1    200 0.9870192 0.6445714
11        3    200 1.0000000 0.7017143
12        6    200 1.0000000 0.7554286
13        1   1000 1.0000000 0.6422857
14        3   1000 1.0000000 0.7165714
15        6   1000 1.0000000 0.7508571
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7554286

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9344443      0.6342857 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9344443 0.6342857
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         3.000         75.000          1.000          0.736 
  maxdepth mfinal rocTrain rocTest
2        3     75        1   0.736
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7097143 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7097143
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9484133      0.6468571 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9484133 0.6468571
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         3.000        100.000          1.000          0.744 
  maxdepth mfinal rocTrain rocTest
5        3    100        1   0.744
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7268571 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7268571
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9660989      0.6468571 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9660989 0.6468571
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7177143 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7177143
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7291429 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7291429
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9756164      0.6445714 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9756164 0.6445714
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7177143 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7177143
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7165714 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7165714
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6514286 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6514286
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7348571 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7348571
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         6.000       1000.000          1.000          0.744 
   maxdepth mfinal rocTrain rocTest
15        6   1000        1   0.744
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9344443 0.6342857
2         3     75 1.0000000 0.7360000
3         6     75 1.0000000 0.7097143
4         1    100 0.9484133 0.6468571
5         3    100 1.0000000 0.7440000
6         6    100 1.0000000 0.7268571
7         1    150 0.9660989 0.6468571
8         3    150 1.0000000 0.7177143
9         6    150 1.0000000 0.7291429
10        1    200 0.9756164 0.6445714
11        3    200 1.0000000 0.7177143
12        6    200 1.0000000 0.7165714
13        1   1000 1.0000000 0.6514286
14        3   1000 1.0000000 0.7348571
15        6   1000 1.0000000 0.7440000
   maxdepth mfinal rocTrain rocTest
5         3    100        1   0.744
15        6   1000        1   0.744

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9537327      0.6697143 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9537327 0.6697143
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7074286 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7074286
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7394286 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7394286
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9674194      0.6777143 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9674194 0.6777143
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7051429 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7051429
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7211429 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7211429
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9820445      0.6845714 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9820445 0.6845714
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7131429 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7131429
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         6.000        150.000          1.000          0.728 
  maxdepth mfinal rocTrain rocTest
9        6    150        1   0.728
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9884643      0.6811429 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9884643 0.6811429
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         3.000        200.000          1.000          0.704 
   maxdepth mfinal rocTrain rocTest
11        3    200        1   0.704
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7462857 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7462857
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6708571 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6708571
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7234286 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7234286
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7405714 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7405714
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9537327 0.6697143
2         3     75 1.0000000 0.7074286
3         6     75 1.0000000 0.7394286
4         1    100 0.9674194 0.6777143
5         3    100 1.0000000 0.7051429
6         6    100 1.0000000 0.7211429
7         1    150 0.9820445 0.6845714
8         3    150 1.0000000 0.7131429
9         6    150 1.0000000 0.7280000
10        1    200 0.9884643 0.6811429
11        3    200 1.0000000 0.7040000
12        6    200 1.0000000 0.7462857
13        1   1000 1.0000000 0.6708571
14        3   1000 1.0000000 0.7234286
15        6   1000 1.0000000 0.7405714
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7462857

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 35 controls (perf_imgT1$obs 0) > 25 cases (perf_imgT1$obs 1).
Area under the curve: 0.7554

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 35 controls (perf_imgT1T2$obs 0) > 25 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.744

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 35 controls (perf_T2wpLMSIR$obs 0) > 25 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.7463
```

```
Area under the curve: 0.7554
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.613964 0.4857    0.6571     0.8   0.68      0.84    0.96
```

```
Area under the curve: 0.744
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5063767 0.6857    0.8286  0.9429   0.44      0.64    0.84
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-1.png)<!-- -->

```
Area under the curve: 0.7463
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5601932 0.6571       0.8  0.9143   0.44      0.64    0.84
   massB    massM nonmassB nonmassM 
     210      150      132       64 
   massB    massM nonmassB nonmassM 
      30       17       10       14 
   massB    massM nonmassB nonmassM 
     210      150      132       64 
   massB    massM nonmassB nonmassM 
      30       17       10       14 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9556915      0.7846774 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9556915 0.7846774
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8177419 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8177419
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7774194 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7774194
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9645318      0.7790323 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9645318 0.7790323
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7967742 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7967742
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7758065 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7758065
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9778992      0.7645161 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9778992 0.7645161
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7903226 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7903226
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7693548 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7693548
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9845936      0.7717742 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9845936 0.7717742
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7879032 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7879032
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7822581 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7822581
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7774194 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7774194
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7951613 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7951613
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8016129 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8016129
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9556915 0.7846774
2         3     75 1.0000000 0.8177419
3         6     75 1.0000000 0.7774194
4         1    100 0.9645318 0.7790323
5         3    100 1.0000000 0.7967742
6         6    100 1.0000000 0.7758065
7         1    150 0.9778992 0.7645161
8         3    150 1.0000000 0.7903226
9         6    150 1.0000000 0.7693548
10        1    200 0.9845936 0.7717742
11        3    200 1.0000000 0.7879032
12        6    200 1.0000000 0.7822581
13        1   1000 1.0000000 0.7774194
14        3   1000 1.0000000 0.7951613
15        6   1000 1.0000000 0.8016129
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8177419

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9489843      0.7330645 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9489843 0.7330645
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7790323 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7790323
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8104839 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8104839
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9543919      0.7274194 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9543919 0.7274194
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7790323 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7790323
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8080645 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8080645
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9680158      0.7290323 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9680158 0.7290323
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7798387 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7798387
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8080645 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8080645
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9772152      0.7274194 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9772152 0.7274194
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8145161 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8145161
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8145161 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8145161
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      0.9999829      0.7814516 
   maxdepth mfinal  rocTrain   rocTest
13        1   1000 0.9999829 0.7814516
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8129032 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8129032
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8306452 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8306452
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9489843 0.7330645
2         3     75 1.0000000 0.7790323
3         6     75 1.0000000 0.8104839
4         1    100 0.9543919 0.7274194
5         3    100 1.0000000 0.7790323
6         6    100 1.0000000 0.8080645
7         1    150 0.9680158 0.7290323
8         3    150 1.0000000 0.7798387
9         6    150 1.0000000 0.8080645
10        1    200 0.9772152 0.7274194
11        3    200 1.0000000 0.8145161
12        6    200 1.0000000 0.8145161
13        1   1000 0.9999829 0.7814516
14        3   1000 1.0000000 0.8129032
15        6   1000 1.0000000 0.8306452
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8306452

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9499077      0.7588710 
  maxdepth mfinal  rocTrain  rocTest
1        1     75 0.9499077 0.758871
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000      75.000000       1.000000       0.816129 
  maxdepth mfinal rocTrain  rocTest
2        3     75        1 0.816129
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7919355 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7919355
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9605605      0.7741935 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9605605 0.7741935
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8274194 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8274194
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8193548 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8193548
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9785148      0.7604839 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9785148 0.7604839
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8282258 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8282258
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8225806 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8225806
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9876458      0.7750000 
   maxdepth mfinal  rocTrain rocTest
10        1    200 0.9876458   0.775
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8266129 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8266129
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8209677 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8209677
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7572581 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7572581
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8274194 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8274194
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8467742 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8467742
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9499077 0.7588710
2         3     75 1.0000000 0.8161290
3         6     75 1.0000000 0.7919355
4         1    100 0.9605605 0.7741935
5         3    100 1.0000000 0.8274194
6         6    100 1.0000000 0.8193548
7         1    150 0.9785148 0.7604839
8         3    150 1.0000000 0.8282258
9         6    150 1.0000000 0.8225806
10        1    200 0.9876458 0.7750000
11        3    200 1.0000000 0.8266129
12        6    200 1.0000000 0.8209677
13        1   1000 1.0000000 0.7572581
14        3   1000 1.0000000 0.8274194
15        6   1000 1.0000000 0.8467742
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8467742

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 40 controls (perf_imgT1$obs 0) > 31 cases (perf_imgT1$obs 1).
Area under the curve: 0.8177

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 40 controls (perf_imgT1T2$obs 0) > 31 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.8306

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 40 controls (perf_T2wpLMSIR$obs 0) > 31 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8468
```

```
Area under the curve: 0.8177
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5129804  0.675       0.8   0.925 0.5806    0.7419   0.871
```

```
Area under the curve: 0.8306
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4915374  0.925     0.975       1 0.4194    0.5806  0.7419
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-2.png)<!-- -->

```
Area under the curve: 0.8468
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6506404  0.475     0.625   0.775 0.8387    0.9355       1
   massB    massM nonmassB nonmassM 
     217      154      129       72 
   massB    massM nonmassB nonmassM 
      23       13       13        6 
   massB    massM nonmassB nonmassM 
     217      154      129       72 
   massB    massM nonmassB nonmassM 
      23       13       13        6 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9342987      0.7894737 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9342987 0.7894737
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7119883 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7119883
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000      75.000000       1.000000       0.752924 
  maxdepth mfinal rocTrain  rocTest
3        6     75        1 0.752924
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9527757      0.7880117 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9527757 0.7880117
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.6988304 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.6988304
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7178363 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7178363
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9685506      0.7587719 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9685506 0.7587719
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.6988304 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.6988304
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7207602 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7207602
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9783404      0.7646199 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9783404 0.7646199
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7192982 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7192982
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7397661 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7397661
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      0.9997995      0.7368421 
   maxdepth mfinal  rocTrain   rocTest
13        1   1000 0.9997995 0.7368421
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7076023 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7076023
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7076023 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7076023
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9342987 0.7894737
2         3     75 1.0000000 0.7119883
3         6     75 1.0000000 0.7529240
4         1    100 0.9527757 0.7880117
5         3    100 1.0000000 0.6988304
6         6    100 1.0000000 0.7178363
7         1    150 0.9685506 0.7587719
8         3    150 1.0000000 0.6988304
9         6    150 1.0000000 0.7207602
10        1    200 0.9783404 0.7646199
11        3    200 1.0000000 0.7192982
12        6    200 1.0000000 0.7397661
13        1   1000 0.9997995 0.7368421
14        3   1000 1.0000000 0.7076023
15        6   1000 1.0000000 0.7076023
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9342987 0.7894737

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9381453      0.7719298 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9381453 0.7719298
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7997076 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7997076
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7383041 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7383041
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9545257      0.7690058 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9545257 0.7690058
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7967836 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7967836
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7090643 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7090643
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9661950      0.7763158 
  maxdepth mfinal rocTrain   rocTest
7        1    150 0.966195 0.7763158
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8055556 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8055556
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7163743 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7163743
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9772963      0.7923977 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9772963 0.7923977
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7953216 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7953216
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7178363 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7178363
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      0.9999916      0.7733918 
   maxdepth mfinal  rocTrain   rocTest
13        1   1000 0.9999916 0.7733918
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7967836 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7967836
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000    1000.000000       1.000000       0.751462 
   maxdepth mfinal rocTrain  rocTest
15        6   1000        1 0.751462
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9381453 0.7719298
2         3     75 1.0000000 0.7997076
3         6     75 1.0000000 0.7383041
4         1    100 0.9545257 0.7690058
5         3    100 1.0000000 0.7967836
6         6    100 1.0000000 0.7090643
7         1    150 0.9661950 0.7763158
8         3    150 1.0000000 0.8055556
9         6    150 1.0000000 0.7163743
10        1    200 0.9772963 0.7923977
11        3    200 1.0000000 0.7953216
12        6    200 1.0000000 0.7178363
13        1   1000 0.9999916 0.7733918
14        3   1000 1.0000000 0.7967836
15        6   1000 1.0000000 0.7514620
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8055556

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9497811      0.7953216 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9497811 0.7953216
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7894737 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7894737
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7982456 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7982456
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9622356      0.7923977 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9622356 0.7923977
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8143275 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8143275
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7909357 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7909357
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9793094      0.7923977 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9793094 0.7923977
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8187135 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8187135
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7733918 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7733918
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9851482      0.8114035 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9851482 0.8114035
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7982456 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7982456
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7704678 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7704678
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7967836 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7967836
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7923977 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7923977
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7850877 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7850877
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9497811 0.7953216
2         3     75 1.0000000 0.7894737
3         6     75 1.0000000 0.7982456
4         1    100 0.9622356 0.7923977
5         3    100 1.0000000 0.8143275
6         6    100 1.0000000 0.7909357
7         1    150 0.9793094 0.7923977
8         3    150 1.0000000 0.8187135
9         6    150 1.0000000 0.7733918
10        1    200 0.9851482 0.8114035
11        3    200 1.0000000 0.7982456
12        6    200 1.0000000 0.7704678
13        1   1000 1.0000000 0.7967836
14        3   1000 1.0000000 0.7923977
15        6   1000 1.0000000 0.7850877
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8187135

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 36 controls (perf_imgT1$obs 0) > 19 cases (perf_imgT1$obs 1).
Area under the curve: 0.7895

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 36 controls (perf_imgT1T2$obs 0) > 19 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.8056

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 36 controls (perf_T2wpLMSIR$obs 0) > 19 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8187
```

```
Area under the curve: 0.7895
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4839599 0.6944    0.8333  0.9444 0.4211    0.6316  0.8421
```

```
Area under the curve: 0.8056
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5121938   0.75    0.8611  0.9722 0.5263    0.7368  0.8947
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-3.png)<!-- -->

```
Area under the curve: 0.8187
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4627353      1         1       1 0.3158    0.5263  0.7368
   massB    massM nonmassB nonmassM 
     224      151      128       70 
   massB    massM nonmassB nonmassM 
      16       16       14        8 
   massB    massM nonmassB nonmassM 
     224      151      128       70 
   massB    massM nonmassB nonmassM 
      16       16       14        8 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9441140      0.7902778 
  maxdepth mfinal rocTrain   rocTest
1        1     75 0.944114 0.7902778
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8152778 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8152778
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7930556 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7930556
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9601385      0.8000000 
  maxdepth mfinal  rocTrain rocTest
4        1    100 0.9601385     0.8
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8027778 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8027778
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7986111 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7986111
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9731647      0.7916667 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9731647 0.7916667
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7986111 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7986111
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7958333 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7958333
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9824299      0.8083333 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9824299 0.8083333
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7861111 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7861111
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7958333 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7958333
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7958333 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7958333
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8083333 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8083333
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7972222 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7972222
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9441140 0.7902778
2         3     75 1.0000000 0.8152778
3         6     75 1.0000000 0.7930556
4         1    100 0.9601385 0.8000000
5         3    100 1.0000000 0.8027778
6         6    100 1.0000000 0.7986111
7         1    150 0.9731647 0.7916667
8         3    150 1.0000000 0.7986111
9         6    150 1.0000000 0.7958333
10        1    200 0.9824299 0.8083333
11        3    200 1.0000000 0.7861111
12        6    200 1.0000000 0.7958333
13        1   1000 1.0000000 0.7958333
14        3   1000 1.0000000 0.8083333
15        6   1000 1.0000000 0.7972222
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8152778

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9339771      0.8000000 
  maxdepth mfinal  rocTrain rocTest
1        1     75 0.9339771     0.8
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8388889 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8388889
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8361111 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8361111
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9481857      0.8180556 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9481857 0.8180556
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         3.000        100.000          1.000          0.825 
  maxdepth mfinal rocTrain rocTest
5        3    100        1   0.825
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8208333 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8208333
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9690405      0.8180556 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9690405 0.8180556
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7986111 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7986111
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8222222 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8222222
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9783461      0.8111111 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9783461 0.8111111
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
           3.0          200.0            1.0            0.8 
   maxdepth mfinal rocTrain rocTest
11        3    200        1     0.8
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8166667 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8166667
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8166667 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8166667
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
        3.0000      1000.0000         1.0000         0.7875 
   maxdepth mfinal rocTrain rocTest
14        3   1000        1  0.7875
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7847222 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7847222
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9339771 0.8000000
2         3     75 1.0000000 0.8388889
3         6     75 1.0000000 0.8361111
4         1    100 0.9481857 0.8180556
5         3    100 1.0000000 0.8250000
6         6    100 1.0000000 0.8208333
7         1    150 0.9690405 0.8180556
8         3    150 1.0000000 0.7986111
9         6    150 1.0000000 0.8222222
10        1    200 0.9783461 0.8111111
11        3    200 1.0000000 0.8000000
12        6    200 1.0000000 0.8166667
13        1   1000 1.0000000 0.8166667
14        3   1000 1.0000000 0.7875000
15        6   1000 1.0000000 0.7847222
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8388889

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9205917      0.7930556 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9205917 0.7930556
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7777778 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7777778
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7958333 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7958333
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9403651      0.7916667 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9403651 0.7916667
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8055556 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8055556
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8013889 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8013889
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9673941      0.7986111 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9673941 0.7986111
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8069444 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8069444
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8236111 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8236111
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9756828      0.8000000 
   maxdepth mfinal  rocTrain rocTest
10        1    200 0.9756828     0.8
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8083333 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8083333
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8222222 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8222222
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8208333 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8208333
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8152778 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8152778
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8069444 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8069444
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9205917 0.7930556
2         3     75 1.0000000 0.7777778
3         6     75 1.0000000 0.7958333
4         1    100 0.9403651 0.7916667
5         3    100 1.0000000 0.8055556
6         6    100 1.0000000 0.8013889
7         1    150 0.9673941 0.7986111
8         3    150 1.0000000 0.8069444
9         6    150 1.0000000 0.8236111
10        1    200 0.9756828 0.8000000
11        3    200 1.0000000 0.8083333
12        6    200 1.0000000 0.8222222
13        1   1000 1.0000000 0.8208333
14        3   1000 1.0000000 0.8152778
15        6   1000 1.0000000 0.8069444
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8236111

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 30 controls (perf_imgT1$obs 0) > 24 cases (perf_imgT1$obs 1).
Area under the curve: 0.8153

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 30 controls (perf_imgT1T2$obs 0) > 24 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.8389

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 30 controls (perf_T2wpLMSIR$obs 0) > 24 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8236
```

```
Area under the curve: 0.8153
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5306807    0.6    0.7667     0.9 0.5833      0.75  0.9167
```

```
Area under the curve: 0.8389
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5257692 0.7333    0.8667  0.9667    0.5    0.6667  0.8333
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-4.png)<!-- -->

```
Area under the curve: 0.8236
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.488011    0.9    0.9667       1 0.4167     0.625  0.8333
   massB    massM nonmassB nonmassM 
     223      148      127       70 
   massB    massM nonmassB nonmassM 
      17       19       15        8 
   massB    massM nonmassB nonmassM 
     223      148      127       70 
   massB    massM nonmassB nonmassM 
      17       19       15        8 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9450204      0.6666667 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9450204 0.6666667
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7372685 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7372685
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7418981 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7418981
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9581306      0.6597222 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9581306 0.6597222
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7395833 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7395833
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7534722 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7534722
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9759592      0.6620370 
  maxdepth mfinal  rocTrain  rocTest
7        1    150 0.9759592 0.662037
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7534722 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7534722
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7708333 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7708333
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9847102      0.6770833 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9847102 0.6770833
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7569444 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7569444
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7662037 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7662037
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6956019 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6956019
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7511574 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7511574
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7789352 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7789352
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9450204 0.6666667
2         3     75 1.0000000 0.7372685
3         6     75 1.0000000 0.7418981
4         1    100 0.9581306 0.6597222
5         3    100 1.0000000 0.7395833
6         6    100 1.0000000 0.7534722
7         1    150 0.9759592 0.6620370
8         3    150 1.0000000 0.7534722
9         6    150 1.0000000 0.7708333
10        1    200 0.9847102 0.6770833
11        3    200 1.0000000 0.7569444
12        6    200 1.0000000 0.7662037
13        1   1000 1.0000000 0.6956019
14        3   1000 1.0000000 0.7511574
15        6   1000 1.0000000 0.7789352
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7789352

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9430531      0.7025463 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9430531 0.7025463
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7465278 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7465278
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8252315 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8252315
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9485224      0.7071759 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9485224 0.7071759
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7650463 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7650463
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8333333 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8333333
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9666531      0.7152778 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9666531 0.7152778
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7592593 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7592593
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       6.00000      150.00000        1.00000        0.84375 
  maxdepth mfinal rocTrain rocTest
9        6    150        1 0.84375
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9775755      0.7199074 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9775755 0.7199074
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7592593 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7592593
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8483796 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8483796
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7118056 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7118056
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8020833 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8020833
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8032407 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8032407
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9430531 0.7025463
2         3     75 1.0000000 0.7465278
3         6     75 1.0000000 0.8252315
4         1    100 0.9485224 0.7071759
5         3    100 1.0000000 0.7650463
6         6    100 1.0000000 0.8333333
7         1    150 0.9666531 0.7152778
8         3    150 1.0000000 0.7592593
9         6    150 1.0000000 0.8437500
10        1    200 0.9775755 0.7199074
11        3    200 1.0000000 0.7592593
12        6    200 1.0000000 0.8483796
13        1   1000 1.0000000 0.7118056
14        3   1000 1.0000000 0.8020833
15        6   1000 1.0000000 0.8032407
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8483796

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9442286      0.7210648 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9442286 0.7210648
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7893519 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7893519
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8321759 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8321759
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9510939      0.7245370 
  maxdepth mfinal  rocTrain  rocTest
4        1    100 0.9510939 0.724537
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8113426 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8113426
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8194444 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8194444
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9707265      0.7164352 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9707265 0.7164352
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7997685 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7997685
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000     150.000000       1.000000       0.806713 
  maxdepth mfinal rocTrain  rocTest
9        6    150        1 0.806713
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9819429      0.7037037 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9819429 0.7037037
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7893519 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7893519
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8148148 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8148148
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7303241 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7303241
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8043981 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8043981
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8055556 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8055556
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9442286 0.7210648
2         3     75 1.0000000 0.7893519
3         6     75 1.0000000 0.8321759
4         1    100 0.9510939 0.7245370
5         3    100 1.0000000 0.8113426
6         6    100 1.0000000 0.8194444
7         1    150 0.9707265 0.7164352
8         3    150 1.0000000 0.7997685
9         6    150 1.0000000 0.8067130
10        1    200 0.9819429 0.7037037
11        3    200 1.0000000 0.7893519
12        6    200 1.0000000 0.8148148
13        1   1000 1.0000000 0.7303241
14        3   1000 1.0000000 0.8043981
15        6   1000 1.0000000 0.8055556
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8321759

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 32 controls (perf_imgT1$obs 0) > 27 cases (perf_imgT1$obs 1).
Area under the curve: 0.7789

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 32 controls (perf_imgT1T2$obs 0) > 27 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.8484

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 32 controls (perf_T2wpLMSIR$obs 0) > 27 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8322
```

```
Area under the curve: 0.7789
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4765379      1         1       1 0.2593    0.4444  0.6296
```

```
Area under the curve: 0.8484
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.564311 0.8438    0.9375       1 0.5185    0.7037  0.8519
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-5.png)<!-- -->

```
Area under the curve: 0.8322
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6376378    0.5    0.6562  0.8125 0.7778    0.8889       1
   massB    massM nonmassB nonmassM 
     209      153      126       72 
   massB    massM nonmassB nonmassM 
      31       14       16        6 
   massB    massM nonmassB nonmassM 
     209      153      126       72 
   massB    massM nonmassB nonmassM 
      31       14       16        6 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9422900      0.7617021 
  maxdepth mfinal rocTrain   rocTest
1        1     75  0.94229 0.7617021
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8021277 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8021277
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8031915 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8031915
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9517487      0.7659574 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9517487 0.7659574
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7904255 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7904255
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8148936 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8148936
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9641613      0.7606383 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9641613 0.7606383
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7861702 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7861702
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8255319 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8255319
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9726799      0.7585106 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9726799 0.7585106
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7904255 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7904255
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8255319 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8255319
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7861702 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7861702
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8074468 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8074468
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8297872 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8297872
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9422900 0.7617021
2         3     75 1.0000000 0.8021277
3         6     75 1.0000000 0.8031915
4         1    100 0.9517487 0.7659574
5         3    100 1.0000000 0.7904255
6         6    100 1.0000000 0.8148936
7         1    150 0.9641613 0.7606383
8         3    150 1.0000000 0.7861702
9         6    150 1.0000000 0.8255319
10        1    200 0.9726799 0.7585106
11        3    200 1.0000000 0.7904255
12        6    200 1.0000000 0.8255319
13        1   1000 1.0000000 0.7861702
14        3   1000 1.0000000 0.8074468
15        6   1000 1.0000000 0.8297872
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8297872

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9443172      0.7893617 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9443172 0.7893617
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7914894 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7914894
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8223404 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8223404
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9591268      0.7936170 
  maxdepth mfinal  rocTrain  rocTest
4        1    100 0.9591268 0.793617
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7787234 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7787234
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8265957 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8265957
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9763422      0.8021277 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9763422 0.8021277
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8287234 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8287234
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8478723 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8478723
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9883805      0.7861702 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9883805 0.7861702
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8202128 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8202128
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
          6.00         200.00           1.00           0.85 
   maxdepth mfinal rocTrain rocTest
12        6    200        1    0.85
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      1.000000    1000.000000       1.000000       0.812766 
   maxdepth mfinal rocTrain  rocTest
13        1   1000        1 0.812766
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8297872 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8297872
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
          6.00        1000.00           1.00           0.85 
   maxdepth mfinal rocTrain rocTest
15        6   1000        1    0.85
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9443172 0.7893617
2         3     75 1.0000000 0.7914894
3         6     75 1.0000000 0.8223404
4         1    100 0.9591268 0.7936170
5         3    100 1.0000000 0.7787234
6         6    100 1.0000000 0.8265957
7         1    150 0.9763422 0.8021277
8         3    150 1.0000000 0.8287234
9         6    150 1.0000000 0.8478723
10        1    200 0.9883805 0.7861702
11        3    200 1.0000000 0.8202128
12        6    200 1.0000000 0.8500000
13        1   1000 1.0000000 0.8127660
14        3   1000 1.0000000 0.8297872
15        6   1000 1.0000000 0.8500000
   maxdepth mfinal rocTrain rocTest
12        6    200        1    0.85
15        6   1000        1    0.85

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9365917      0.8021277 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9365917 0.8021277
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7765957 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7765957
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8553191 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8553191
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9519091      0.8031915 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9519091 0.8031915
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7648936 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7648936
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8510638 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8510638
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9668078      0.8031915 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9668078 0.8031915
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7712766 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7712766
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.8531915 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.8531915
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9794787      0.7978723 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9794787 0.7978723
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7797872 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7797872
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8478723 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8478723
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8117021 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8117021
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7978723 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7978723
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000    1000.000000       1.000000       0.843617 
   maxdepth mfinal rocTrain  rocTest
15        6   1000        1 0.843617
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9365917 0.8021277
2         3     75 1.0000000 0.7765957
3         6     75 1.0000000 0.8553191
4         1    100 0.9519091 0.8031915
5         3    100 1.0000000 0.7648936
6         6    100 1.0000000 0.8510638
7         1    150 0.9668078 0.8031915
8         3    150 1.0000000 0.7712766
9         6    150 1.0000000 0.8531915
10        1    200 0.9794787 0.7978723
11        3    200 1.0000000 0.7797872
12        6    200 1.0000000 0.8478723
13        1   1000 1.0000000 0.8117021
14        3   1000 1.0000000 0.7978723
15        6   1000 1.0000000 0.8436170
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8553191

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 47 controls (perf_imgT1$obs 0) > 20 cases (perf_imgT1$obs 1).
Area under the curve: 0.8298

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 47 controls (perf_imgT1T2$obs 0) > 20 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.85

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 47 controls (perf_T2wpLMSIR$obs 0) > 20 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8553
```

```
Area under the curve: 0.8298
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5750095 0.5532    0.7021  0.8298   0.75       0.9       1
```

```
Area under the curve: 0.85
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5757877 0.5745    0.7021  0.8298   0.75       0.9       1
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-6.png)<!-- -->

```
Area under the curve: 0.8553
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5155473 0.7447    0.8511  0.9367   0.55      0.75     0.9
   massB    massM nonmassB nonmassM 
     222      148      121       73 
   massB    massM nonmassB nonmassM 
      18       19       21        5 
   massB    massM nonmassB nonmassM 
     222      148      121       73 
   massB    massM nonmassB nonmassM 
      18       19       21        5 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9442239      0.7500000 
  maxdepth mfinal  rocTrain rocTest
1        1     75 0.9442239    0.75
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7606838 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7606838
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7286325 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7286325
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9578152      0.7585470 
  maxdepth mfinal  rocTrain  rocTest
4        1    100 0.9578152 0.758547
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7617521 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7617521
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7361111 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7361111
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9726985      0.7574786 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9726985 0.7574786
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7713675 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7713675
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7510684 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7510684
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9813173      0.7628205 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9813173 0.7628205
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7638889 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7638889
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7542735 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7542735
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      1.000000    1000.000000       1.000000       0.741453 
   maxdepth mfinal rocTrain  rocTest
13        1   1000        1 0.741453
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7596154 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7596154
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000    1000.000000       1.000000       0.767094 
   maxdepth mfinal rocTrain  rocTest
15        6   1000        1 0.767094
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9442239 0.7500000
2         3     75 1.0000000 0.7606838
3         6     75 1.0000000 0.7286325
4         1    100 0.9578152 0.7585470
5         3    100 1.0000000 0.7617521
6         6    100 1.0000000 0.7361111
7         1    150 0.9726985 0.7574786
8         3    150 1.0000000 0.7713675
9         6    150 1.0000000 0.7510684
10        1    200 0.9813173 0.7628205
11        3    200 1.0000000 0.7638889
12        6    200 1.0000000 0.7542735
13        1   1000 1.0000000 0.7414530
14        3   1000 1.0000000 0.7596154
15        6   1000 1.0000000 0.7670940
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7713675

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9325366      0.7414530 
  maxdepth mfinal  rocTrain  rocTest
1        1     75 0.9325366 0.741453
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.6794872 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.6794872
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7553419 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7553419
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9450314      0.7532051 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9450314 0.7532051
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7232906 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7232906
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7542735 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7542735
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9599997      0.7264957 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9599997 0.7264957
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7489316 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7489316
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7596154 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7596154
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9713300      0.7168803 
   maxdepth mfinal rocTrain   rocTest
10        1    200  0.97133 0.7168803
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7628205 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7628205
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7489316 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7489316
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6826923 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6826923
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7435897 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7435897
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7371795 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7371795
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9325366 0.7414530
2         3     75 1.0000000 0.6794872
3         6     75 1.0000000 0.7553419
4         1    100 0.9450314 0.7532051
5         3    100 1.0000000 0.7232906
6         6    100 1.0000000 0.7542735
7         1    150 0.9599997 0.7264957
8         3    150 1.0000000 0.7489316
9         6    150 1.0000000 0.7596154
10        1    200 0.9713300 0.7168803
11        3    200 1.0000000 0.7628205
12        6    200 1.0000000 0.7489316
13        1   1000 1.0000000 0.6826923
14        3   1000 1.0000000 0.7435897
15        6   1000 1.0000000 0.7371795
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7628205

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9485843      0.7692308 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9485843 0.7692308
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8215812 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8215812
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7425214 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7425214
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9594557      0.7702991 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9594557 0.7702991
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7991453 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7991453
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7553419 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7553419
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9747809      0.7692308 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9747809 0.7692308
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7991453 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7991453
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7799145 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7799145
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9853462      0.7724359 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9853462 0.7724359
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8023504 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8023504
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7991453 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7991453
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7489316 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7489316
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7820513 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7820513
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000    1000.000000       1.000000       0.792735 
   maxdepth mfinal rocTrain  rocTest
15        6   1000        1 0.792735
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9485843 0.7692308
2         3     75 1.0000000 0.8215812
3         6     75 1.0000000 0.7425214
4         1    100 0.9594557 0.7702991
5         3    100 1.0000000 0.7991453
6         6    100 1.0000000 0.7553419
7         1    150 0.9747809 0.7692308
8         3    150 1.0000000 0.7991453
9         6    150 1.0000000 0.7799145
10        1    200 0.9853462 0.7724359
11        3    200 1.0000000 0.8023504
12        6    200 1.0000000 0.7991453
13        1   1000 1.0000000 0.7489316
14        3   1000 1.0000000 0.7820513
15        6   1000 1.0000000 0.7927350
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8215812

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 39 controls (perf_imgT1$obs 0) > 24 cases (perf_imgT1$obs 1).
Area under the curve: 0.7714

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 39 controls (perf_imgT1T2$obs 0) > 24 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.7628

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 39 controls (perf_T2wpLMSIR$obs 0) > 24 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8216
```

```
Area under the curve: 0.7714
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5761598 0.5891    0.7179  0.8462 0.6667    0.8333  0.9583
```

```
Area under the curve: 0.7628
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5655913 0.5128    0.6667  0.8205 0.5833      0.75  0.9167
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-7.png)<!-- -->

```
Area under the curve: 0.8216
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5459772 0.7949    0.8974  0.9744 0.5833      0.75  0.9167
   massB    massM nonmassB nonmassM 
     212      143      128       77 
   massB    massM nonmassB nonmassM 
      28       24       14        1 
   massB    massM nonmassB nonmassM 
     212      143      128       77 
   massB    massM nonmassB nonmassM 
      28       24       14        1 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9423486      0.7914286 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9423486 0.7914286
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
          3.00          75.00           1.00           0.76 
  maxdepth mfinal rocTrain rocTest
2        3     75        1    0.76
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7971429 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7971429
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9537067      0.7885714 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9537067 0.7885714
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8114286 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8114286
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7857143 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7857143
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9756834      0.7800000 
  maxdepth mfinal  rocTrain rocTest
7        1    150 0.9756834    0.78
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     150.000000       1.000000       0.807619 
  maxdepth mfinal rocTrain  rocTest
8        3    150        1 0.807619
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7742857 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7742857
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9862197      0.7771429 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9862197 0.7771429
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8009524 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8009524
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7780952 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7780952
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8095238 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8095238
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7971429 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7971429
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7914286 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7914286
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9423486 0.7914286
2         3     75 1.0000000 0.7600000
3         6     75 1.0000000 0.7971429
4         1    100 0.9537067 0.7885714
5         3    100 1.0000000 0.8114286
6         6    100 1.0000000 0.7857143
7         1    150 0.9756834 0.7800000
8         3    150 1.0000000 0.8076190
9         6    150 1.0000000 0.7742857
10        1    200 0.9862197 0.7771429
11        3    200 1.0000000 0.8009524
12        6    200 1.0000000 0.7780952
13        1   1000 1.0000000 0.8095238
14        3   1000 1.0000000 0.7971429
15        6   1000 1.0000000 0.7914286
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8114286

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9494983      0.7647619 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9494983 0.7647619
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7561905 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7561905
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7447619 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7447619
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9594204      0.7704762 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9594204 0.7704762
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7371429 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7371429
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000     100.000000       1.000000       0.747619 
  maxdepth mfinal rocTrain  rocTest
6        6    100        1 0.747619
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9775606      0.7628571 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9775606 0.7628571
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7428571 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7428571
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7561905 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7561905
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9856401      0.7580952 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9856401 0.7580952
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.7409524 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.7409524
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7609524 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7609524
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.7247619 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.7247619
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7771429 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7771429
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7895238 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7895238
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9494983 0.7647619
2         3     75 1.0000000 0.7561905
3         6     75 1.0000000 0.7447619
4         1    100 0.9594204 0.7704762
5         3    100 1.0000000 0.7371429
6         6    100 1.0000000 0.7476190
7         1    150 0.9775606 0.7628571
8         3    150 1.0000000 0.7428571
9         6    150 1.0000000 0.7561905
10        1    200 0.9856401 0.7580952
11        3    200 1.0000000 0.7409524
12        6    200 1.0000000 0.7609524
13        1   1000 1.0000000 0.7247619
14        3   1000 1.0000000 0.7771429
15        6   1000 1.0000000 0.7895238
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7895238

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9579066      0.8342857 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9579066 0.8342857
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000      75.000000       1.000000       0.772381 
  maxdepth mfinal rocTrain  rocTest
2        3     75        1 0.772381
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7742857 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7742857
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9655969      0.8333333 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9655969 0.8333333
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.7819048 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.7819048
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
          6.00         100.00           1.00           0.78 
  maxdepth mfinal rocTrain rocTest
6        6    100        1    0.78
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9775346      0.8352381 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9775346 0.8352381
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7838095 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7838095
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7666667 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7666667
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9852076      0.8380952 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9852076 0.8380952
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     200.000000       1.000000       0.787619 
   maxdepth mfinal rocTrain  rocTest
11        3    200        1 0.787619
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7838095 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7838095
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8295238 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8295238
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8104762 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8104762
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8114286 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8114286
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9579066 0.8342857
2         3     75 1.0000000 0.7723810
3         6     75 1.0000000 0.7742857
4         1    100 0.9655969 0.8333333
5         3    100 1.0000000 0.7819048
6         6    100 1.0000000 0.7800000
7         1    150 0.9775346 0.8352381
8         3    150 1.0000000 0.7838095
9         6    150 1.0000000 0.7666667
10        1    200 0.9852076 0.8380952
11        3    200 1.0000000 0.7876190
12        6    200 1.0000000 0.7838095
13        1   1000 1.0000000 0.8295238
14        3   1000 1.0000000 0.8104762
15        6   1000 1.0000000 0.8114286
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9852076 0.8380952

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 42 controls (perf_imgT1$obs 0) > 25 cases (perf_imgT1$obs 1).
Area under the curve: 0.8114

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 42 controls (perf_imgT1T2$obs 0) > 25 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.7895

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 42 controls (perf_T2wpLMSIR$obs 0) > 25 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8381
```

```
Area under the curve: 0.8114
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5557392  0.619    0.7381  0.8571    0.6      0.76    0.92
```

```
Area under the curve: 0.7895
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5846236 0.7381    0.8571  0.9524   0.44      0.64     0.8
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-8.png)<!-- -->

```
Area under the curve: 0.8381
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5333076 0.7143    0.8333  0.9286   0.64       0.8    0.96
   massB    massM nonmassB nonmassM 
     209      154      125       71 
   massB    massM nonmassB nonmassM 
      31       13       17        7 
   massB    massM nonmassB nonmassM 
     209      154      125       71 
   massB    massM nonmassB nonmassM 
      31       13       17        7 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9249256      0.7322917 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9249256 0.7322917
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.7645833 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.7645833
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7677083 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7677083
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9372871      0.7291667 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9372871 0.7291667
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     100.000000       1.000000       0.734375 
  maxdepth mfinal rocTrain  rocTest
5        3    100        1 0.734375
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.7666667 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.7666667
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9650221      0.7312500 
  maxdepth mfinal  rocTrain rocTest
7        1    150 0.9650221 0.73125
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
        3.0000       150.0000         1.0000         0.7625 
  maxdepth mfinal rocTrain rocTest
8        3    150        1  0.7625
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.7677083 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.7677083
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9787102      0.7375000 
   maxdepth mfinal  rocTrain rocTest
10        1    200 0.9787102  0.7375
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       3.00000      200.00000        1.00000        0.76875 
   maxdepth mfinal rocTrain rocTest
11        3    200        1 0.76875
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
         6.000        200.000          1.000          0.775 
   maxdepth mfinal rocTrain rocTest
12        6    200        1   0.775
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
        1.0000      1000.0000         1.0000         0.7125 
   maxdepth mfinal rocTrain rocTest
13        1   1000        1  0.7125
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000    1000.000000       1.000000       0.753125 
   maxdepth mfinal rocTrain  rocTest
14        3   1000        1 0.753125
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7604167 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7604167
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9249256 0.7322917
2         3     75 1.0000000 0.7645833
3         6     75 1.0000000 0.7677083
4         1    100 0.9372871 0.7291667
5         3    100 1.0000000 0.7343750
6         6    100 1.0000000 0.7666667
7         1    150 0.9650221 0.7312500
8         3    150 1.0000000 0.7625000
9         6    150 1.0000000 0.7677083
10        1    200 0.9787102 0.7375000
11        3    200 1.0000000 0.7687500
12        6    200 1.0000000 0.7750000
13        1   1000 1.0000000 0.7125000
14        3   1000 1.0000000 0.7531250
15        6   1000 1.0000000 0.7604167
   maxdepth mfinal rocTrain rocTest
12        6    200        1   0.775

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9581242      0.6666667 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9581242 0.6666667
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.5989583 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.5989583
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.6864583 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.6864583
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9664384      0.6666667 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9664384 0.6666667
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     100.000000       1.000000       0.621875 
  maxdepth mfinal rocTrain  rocTest
5        3    100        1 0.621875
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.6833333 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.6833333
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9862222      0.6322917 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9862222 0.6322917
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.6427083 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.6427083
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.6885417 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.6885417
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9907849      0.6354167 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9907849 0.6354167
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       3.00000      200.00000        1.00000        0.64375 
   maxdepth mfinal rocTrain rocTest
11        3    200        1 0.64375
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.6895833 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.6895833
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.5979167 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.5979167
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.6291667 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.6291667
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       6.00000     1000.00000        1.00000        0.73125 
   maxdepth mfinal rocTrain rocTest
15        6   1000        1 0.73125
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9581242 0.6666667
2         3     75 1.0000000 0.5989583
3         6     75 1.0000000 0.6864583
4         1    100 0.9664384 0.6666667
5         3    100 1.0000000 0.6218750
6         6    100 1.0000000 0.6833333
7         1    150 0.9862222 0.6322917
8         3    150 1.0000000 0.6427083
9         6    150 1.0000000 0.6885417
10        1    200 0.9907849 0.6354167
11        3    200 1.0000000 0.6437500
12        6    200 1.0000000 0.6895833
13        1   1000 1.0000000 0.5979167
14        3   1000 1.0000000 0.6291667
15        6   1000 1.0000000 0.7312500
   maxdepth mfinal rocTrain rocTest
15        6   1000        1 0.73125

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9420291      0.7197917 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9420291 0.7197917
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000      75.000000       1.000000       0.715625 
  maxdepth mfinal rocTrain  rocTest
2        3     75        1 0.715625
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.7645833 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.7645833
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9588368      0.7343750 
  maxdepth mfinal  rocTrain  rocTest
4        1    100 0.9588368 0.734375
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       3.00000      100.00000        1.00000        0.70625 
  maxdepth mfinal rocTrain rocTest
5        3    100        1 0.70625
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8020833 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8020833
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9779931      0.7312500 
  maxdepth mfinal  rocTrain rocTest
7        1    150 0.9779931 0.73125
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.7333333 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.7333333
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000     150.000000       1.000000       0.771875 
  maxdepth mfinal rocTrain  rocTest
9        6    150        1 0.771875
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9866256      0.7083333 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9866256 0.7083333
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
       3.00000      200.00000        1.00000        0.74375 
   maxdepth mfinal rocTrain rocTest
11        3    200        1 0.74375
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.7697917 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.7697917
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.6770833 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.6770833
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.7114583 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.7114583
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.7697917 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.7697917
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9420291 0.7197917
2         3     75 1.0000000 0.7156250
3         6     75 1.0000000 0.7645833
4         1    100 0.9588368 0.7343750
5         3    100 1.0000000 0.7062500
6         6    100 1.0000000 0.8020833
7         1    150 0.9779931 0.7312500
8         3    150 1.0000000 0.7333333
9         6    150 1.0000000 0.7718750
10        1    200 0.9866256 0.7083333
11        3    200 1.0000000 0.7437500
12        6    200 1.0000000 0.7697917
13        1   1000 1.0000000 0.6770833
14        3   1000 1.0000000 0.7114583
15        6   1000 1.0000000 0.7697917
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8020833

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 48 controls (perf_imgT1$obs 0) > 20 cases (perf_imgT1$obs 1).
Area under the curve: 0.775

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 48 controls (perf_imgT1T2$obs 0) > 20 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.7312

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 48 controls (perf_T2wpLMSIR$obs 0) > 20 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.8021
```

```
Area under the curve: 0.775
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5456363 0.7083    0.8333  0.9375   0.45      0.65    0.85
```

```
Area under the curve: 0.7312
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6148672 0.5208    0.6667  0.7917   0.55      0.75     0.9
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-9.png)<!-- -->

```
Area under the curve: 0.8021
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5295403   0.75    0.8542  0.9375   0.45      0.65    0.85
   massB    massM nonmassB nonmassM 
     218      150      131       65 
   massB    massM nonmassB nonmassM 
      22       17       11       13 
   massB    massM nonmassB nonmassM 
     218      150      131       65 
   massB    massM nonmassB nonmassM 
      22       17       11       13 

============ boosting trees imgT1train 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9507763      0.8727273 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9507763 0.8727273
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8575758 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8575758
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.8757576 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.8757576
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9630011      0.8777778 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9630011 0.8777778
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8717172 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8717172
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.8707071 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.8707071
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9751972      0.8828283 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9751972 0.8828283
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8676768 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8676768
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000     150.000000       1.000000       0.879798 
  maxdepth mfinal rocTrain  rocTest
9        6    150        1 0.879798
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9842202      0.8797980 
   maxdepth mfinal  rocTrain  rocTest
10        1    200 0.9842202 0.879798
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.8787879 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.8787879
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      6.000000     200.000000       1.000000       0.879798 
   maxdepth mfinal rocTrain  rocTest
12        6    200        1 0.879798
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8636364 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8636364
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.8808081 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.8808081
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.8959596 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8959596
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9507763 0.8727273
2         3     75 1.0000000 0.8575758
3         6     75 1.0000000 0.8757576
4         1    100 0.9630011 0.8777778
5         3    100 1.0000000 0.8717172
6         6    100 1.0000000 0.8707071
7         1    150 0.9751972 0.8828283
8         3    150 1.0000000 0.8676768
9         6    150 1.0000000 0.8797980
10        1    200 0.9842202 0.8797980
11        3    200 1.0000000 0.8787879
12        6    200 1.0000000 0.8797980
13        1   1000 1.0000000 0.8636364
14        3   1000 1.0000000 0.8808081
15        6   1000 1.0000000 0.8959596
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.8959596

============ boosting trees imgT1T2 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9431860      0.8909091 
  maxdepth mfinal rocTrain   rocTest
1        1     75 0.943186 0.8909091
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.9010101 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.9010101
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.9222222 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.9222222
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9507065      0.8979798 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9507065 0.8979798
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    100.0000000      1.0000000      0.8989899 
  maxdepth mfinal rocTrain   rocTest
5        3    100        1 0.8989899
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.9080808 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.9080808
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9686292      0.9010101 
  maxdepth mfinal  rocTrain   rocTest
7        1    150 0.9686292 0.9010101
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    150.0000000      1.0000000      0.8969697 
  maxdepth mfinal rocTrain   rocTest
8        3    150        1 0.8969697
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.9040404 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.9040404
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9806159      0.8898990 
   maxdepth mfinal  rocTrain  rocTest
10        1    200 0.9806159 0.889899
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.9111111 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.9111111
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.8969697 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.8969697
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8717172 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8717172
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.9050505 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.9050505
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.9171717 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.9171717
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9431860 0.8909091
2         3     75 1.0000000 0.9010101
3         6     75 1.0000000 0.9222222
4         1    100 0.9507065 0.8979798
5         3    100 1.0000000 0.8989899
6         6    100 1.0000000 0.9080808
7         1    150 0.9686292 0.9010101
8         3    150 1.0000000 0.8969697
9         6    150 1.0000000 0.9040404
10        1    200 0.9806159 0.8898990
11        3    200 1.0000000 0.9111111
12        6    200 1.0000000 0.8969697
13        1   1000 1.0000000 0.8717172
14        3   1000 1.0000000 0.9050505
15        6   1000 1.0000000 0.9171717
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.9222222

============ boosting trees imgT2pLMSIR 
maxdepth:  1 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000     75.0000000      0.9443272      0.9050505 
  maxdepth mfinal  rocTrain   rocTest
1        1     75 0.9443272 0.9050505
maxdepth:  3 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000     75.0000000      1.0000000      0.8828283 
  maxdepth mfinal rocTrain   rocTest
2        3     75        1 0.8828283
maxdepth:  6 mfinal:  75 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000     75.0000000      1.0000000      0.9333333 
  maxdepth mfinal rocTrain   rocTest
3        6     75        1 0.9333333
maxdepth:  1 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    100.0000000      0.9609034      0.8848485 
  maxdepth mfinal  rocTrain   rocTest
4        1    100 0.9609034 0.8848485
maxdepth:  3 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     100.000000       1.000000       0.910101 
  maxdepth mfinal rocTrain  rocTest
5        3    100        1 0.910101
maxdepth:  6 mfinal:  100 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    100.0000000      1.0000000      0.9333333 
  maxdepth mfinal rocTrain   rocTest
6        6    100        1 0.9333333
maxdepth:  1 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    150.0000000      0.9767900      0.8838384 
  maxdepth mfinal rocTrain   rocTest
7        1    150  0.97679 0.8838384
maxdepth:  3 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
      3.000000     150.000000       1.000000       0.910101 
  maxdepth mfinal rocTrain  rocTest
8        3    150        1 0.910101
maxdepth:  6 mfinal:  150 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    150.0000000      1.0000000      0.9383838 
  maxdepth mfinal rocTrain   rocTest
9        6    150        1 0.9383838
maxdepth:  1 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000    200.0000000      0.9864944      0.8808081 
   maxdepth mfinal  rocTrain   rocTest
10        1    200 0.9864944 0.8808081
maxdepth:  3 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000    200.0000000      1.0000000      0.9333333 
   maxdepth mfinal rocTrain   rocTest
11        3    200        1 0.9333333
maxdepth:  6 mfinal:  200 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000    200.0000000      1.0000000      0.9484848 
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.9484848
maxdepth:  1 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     1.0000000   1000.0000000      1.0000000      0.8626263 
   maxdepth mfinal rocTrain   rocTest
13        1   1000        1 0.8626263
maxdepth:  3 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     3.0000000   1000.0000000      1.0000000      0.9161616 
   maxdepth mfinal rocTrain   rocTest
14        3   1000        1 0.9161616
maxdepth:  6 mfinal:  1000 
      maxdepth         mfinal boost.rocTrain  boost.rocTest 
     6.0000000   1000.0000000      1.0000000      0.9363636 
   maxdepth mfinal rocTrain   rocTest
15        6   1000        1 0.9363636
   maxdepth mfinal  rocTrain   rocTest
1         1     75 0.9443272 0.9050505
2         3     75 1.0000000 0.8828283
3         6     75 1.0000000 0.9333333
4         1    100 0.9609034 0.8848485
5         3    100 1.0000000 0.9101010
6         6    100 1.0000000 0.9333333
7         1    150 0.9767900 0.8838384
8         3    150 1.0000000 0.9101010
9         6    150 1.0000000 0.9383838
10        1    200 0.9864944 0.8808081
11        3    200 1.0000000 0.9333333
12        6    200 1.0000000 0.9484848
13        1   1000 1.0000000 0.8626263
14        3   1000 1.0000000 0.9161616
15        6   1000 1.0000000 0.9363636
   maxdepth mfinal rocTrain   rocTest
12        6    200        1 0.9484848

Call:
roc.default(response = perf_imgT1$obs, predictor = perf_imgT1$C)

Data: perf_imgT1$C in 33 controls (perf_imgT1$obs 0) > 30 cases (perf_imgT1$obs 1).
Area under the curve: 0.896

Call:
roc.default(response = perf_imgT1T2$obs, predictor = perf_imgT1T2$C)

Data: perf_imgT1T2$C in 33 controls (perf_imgT1T2$obs 0) > 30 cases (perf_imgT1T2$obs 1).
Area under the curve: 0.9222

Call:
roc.default(response = perf_T2wpLMSIR$obs, predictor = perf_T2wpLMSIR$C)

Data: perf_T2wpLMSIR$C in 33 controls (perf_T2wpLMSIR$obs 0) > 30 cases (perf_T2wpLMSIR$obs 1).
Area under the curve: 0.9485
```

```
Area under the curve: 0.896
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5971968 0.7576    0.8788  0.9697    0.7    0.8333  0.9667
```

```
Area under the curve: 0.9222
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5998114 0.8485    0.9394       1 0.7333    0.8667  0.9667
```

![](Boosting_classification_trees_T2imgvspredLMSIR_wT2wtex_files/figure-html/run-boosting-10.png)<!-- -->

```
Area under the curve: 0.9485
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5795522 0.9091    0.9697       1    0.7    0.8333  0.9667
```


