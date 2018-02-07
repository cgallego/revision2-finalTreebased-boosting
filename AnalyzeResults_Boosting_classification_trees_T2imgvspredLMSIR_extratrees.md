# Analysis of results - T2SI vs predLMSIR
Cristina Gallego  
August 25, 2017  



# Analysis of results by each fold

```r
allcvauc_imgT1_total= c()
allcvauc_T1wT2w_total = c()
allcvauc_T2wpLMSIR_total= c()

## using 3Dtexture first + Boosting  
perfall_imgT1_total = data.frame() 
perfall_imgT1T2_total = data.frame() 
perfall_T2wpLMSIR_total = data.frame() 

varImpall_imgT1w_total = data.frame() 
varImpall_T1wT2w_total = data.frame() 
varImpall_T2wpLMSIR_total = data.frame() 

# perform k-fold-out
kfolds = c(1,3,8:10,1,3,8:10)
for(k in kfolds){  # 1:10f cv
  load(paste0("Outputs/T1T2imgvsT2wtextpredLMSIR_boost_addeddiagvalue_extratrees_cv",k,".RData"))

  # append test results
  perfall_imgT1_total = rbind(perfall_imgT1_total, perf_imgT1) 
  perfall_imgT1T2_total = rbind(perfall_imgT1T2_total, perf_imgT1T2) 
  perfall_T2wpLMSIR_total = rbind(perfall_T2wpLMSIR_total, perf_T2wpLMSIR) 
  
  # append cv results
  allcvauc_imgT1_total = c(allcvauc_imgT1_total, rocperf_imgT1$auc)
  allcvauc_T1wT2w_total = c(allcvauc_T1wT2w_total, rocperf_imgT1T2$auc)
  allcvauc_T2wpLMSIR_total = c(allcvauc_T2wpLMSIR_total, rocperf_T2wpLMSIR$auc)

    
  # plot cv-heldout test cases
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
                    paste0("T1w + T2wText + SI"),
                    paste0("T1w + T2wText + predLMSIR")),
         col = c(colors[2],colors[9],colors[11]), lwd = 2)

  
  ## process varImpts
  # returns the relative importance of each variable in the classification task.   
  # This measure takes into account the gain of the Gini index given by a variable in a tree and the weight of this tree.
  mostImp = treedata_imgT1$forest$importance[treedata_imgT1$forest$importance>0]
  varImp_imgT1w = data.frame(selfeat=names(mostImp),
                             gainGini=mostImp,
                             SelectedFeatureGroup="imgT1")
  row.names(varImp_imgT1w) = NULL
  
  # wimgT2
  mostImp = treedata_imgT1T2$forest$importance[treedata_imgT1T2$forest$importance>0]
  varImp_T2wmLMSIR = data.frame(selfeat=names(mostImp),
                             gainGini=mostImp,
                             SelectedFeatureGroup="imgT2")
  row.names(varImp_T2wmLMSIR) = NULL
  
  # wLMSIR
  mostImp = treedata_T2wpLMSIR$forest$importance[treedata_T2wpLMSIR$forest$importance>0]
  varImp_T2wpLMSIR = data.frame(selfeat=names(mostImp),
                             gainGini=mostImp,
                             SelectedFeatureGroup="T1+T2")
  row.names(varImp_T2wpLMSIR) = NULL
  
  ## group with all of the features spaces combined, most contributing T2w feature
  varImpall_imgT1w_total =  rbind(varImpall_imgT1w_total, cbind(varImp_imgT1w, kfcv=k) )
  varImpall_T1wT2w_total =  rbind(varImpall_T1wT2w_total, cbind(varImp_T2wmLMSIR, kfcv=k) )
  varImpall_T2wpLMSIR_total =  rbind(varImpall_T2wpLMSIR_total, cbind(varImp_T2wpLMSIR, kfcv=k) )
}
```

```
Area under the curve: 0.7406
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4981562   0.44      0.64    0.84 0.6857    0.8286  0.9429
```

```
Area under the curve: 0.776
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4592056   0.64       0.8    0.96 0.5707    0.7143  0.8571
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-1.png)<!-- -->

```
Area under the curve: 0.7691
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5006342   0.52      0.72    0.88 0.6857    0.8286  0.9429
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-2.png)<!-- -->

```
Area under the curve: 0.8246
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.491993 0.3684    0.5789  0.7895 0.8611    0.9444       1
```

```
Area under the curve: 0.841
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3985963  0.759      0.88       1  0.619    0.7619   0.881
```

```
Area under the curve: 0.8038
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4320352    0.8      0.92       1 0.4762     0.619  0.7619
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-3.png)<!-- -->

```
Area under the curve: 0.8562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3802358   0.76      0.88       1 0.5714    0.7143  0.8333
```

```
Area under the curve: 0.7385
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4211708    0.7      0.85       1 0.4792     0.625    0.75
```

```
Area under the curve: 0.7562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4596946   0.55      0.75  0.9012 0.6458    0.7708   0.875
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-4.png)<!-- -->

```
Area under the curve: 0.7802
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4135362   0.75       0.9       1 0.4375    0.5833  0.7292
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-5.png)<!-- -->

```
Area under the curve: 0.9263
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4587403 0.6658       0.8  0.9333 0.8485    0.9394       1
```

```
Area under the curve: 0.7406
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4981562   0.44      0.64    0.84 0.6857    0.8286  0.9429
```

```
Area under the curve: 0.776
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4592056   0.64       0.8    0.96 0.5707    0.7143  0.8571
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-6.png)<!-- -->

```
Area under the curve: 0.7691
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5006342   0.52      0.72    0.88 0.6857    0.8286  0.9429
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-7.png)<!-- -->

```
Area under the curve: 0.8246
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.491993 0.3684    0.5789  0.7895 0.8611    0.9444       1
```

```
Area under the curve: 0.841
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3985963  0.759      0.88       1  0.619    0.7619   0.881
```

```
Area under the curve: 0.8038
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4320352    0.8      0.92       1 0.4762     0.619  0.7619
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-8.png)<!-- -->

```
Area under the curve: 0.8562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.3802358   0.76      0.88       1 0.5714    0.7143  0.8333
```

```
Area under the curve: 0.7385
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4211708    0.7      0.85       1 0.4792     0.625    0.75
```

```
Area under the curve: 0.7562
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4596946   0.55      0.75  0.9012 0.6458    0.7708   0.875
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-9.png)<!-- -->

```
Area under the curve: 0.7802
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4135362   0.75       0.9       1 0.4375    0.5833  0.7292
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/process-cvks-10.png)<!-- -->

```
Area under the curve: 0.9263
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4587403 0.6658       0.8  0.9333 0.8485    0.9394       1
```

# Combine, pool data, and plot final results

```r
print(summary(allcvauc_imgT1_total))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7385  0.7406  0.7471  0.7961  0.8410  0.9131 
```

```r
print(summary(allcvauc_T1wT2w_total))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7562  0.7760  0.7807  0.8046  0.8038  0.9061 
```

```r
print(summary(allcvauc_T2wpLMSIR_total))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7691  0.7802  0.8246  0.8313  0.8562  0.9263 
```

```r
# boxplots of cv-performances
cvperfs = data.frame()
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_imgT1_total, group="T1w"))
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_T1wT2w_total, group="T1w+T2wText+T2w_SI"))
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_T2wpLMSIR_total, group="T1w+T2wText+predictiveLMSIR"))
# find min
minAUC = min(cvperfs$cvAUC)

# plot
p <- ggplot(cvperfs, aes(factor(group), cvAUC))
p + geom_boxplot(aes(fill = factor(group)))
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/Summarize-1.png)<!-- -->

```r
########### 
# plot pooled all cv-heldout test cases
par(mfrow=c(1,1))
n=15
colors = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
# plot 1/4
p1 = calcAUC_plot(perfall_imgT1_total$obs, perfall_imgT1_total$C, 
                           xptext=0.45, yptext=0.75 , 1, colors[2], atitle="")
```

```
Area under the curve: 0.7992
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4156007 0.7185    0.7731  0.8235 0.6546     0.701  0.7474
```

```r
par(new=TRUE)
p2 = calcAUC_plot(perfall_imgT1T2_total$obs, perfall_imgT1T2_total$C, 
                           xptext=0.55, yptext=0.65, 2, colors[9], atitle="")
```

```
Area under the curve: 0.8049
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4592056 0.7185    0.7731  0.8277 0.7088    0.7526  0.7964
```

```r
par(new=TRUE)
p3 = calcAUC_plot(perfall_T2wpLMSIR_total$obs, perfall_T2wpLMSIR_total$C,
                           xptext=0.65, yptext=0.55, 3, colors[11], 
                  atitle="ROCs pooled heldout-patient across ALL folds")
```

```
Area under the curve: 0.8333
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4269815 0.7437    0.7983  0.8487 0.6881     0.732  0.7758
```

```r
legend("bottomright", 
       legend = c(paste0("T1w"),
                    paste0("T1w + T2wText + T2w_SI"),
                    paste0("T1w + T2wText + predictiveLMSIR")),
       col = c(colors[2],colors[9],colors[11]), lty=c(1,2,3), lwd = 2)
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/Summarize-2.png)<!-- -->

```r
# find significants: only imgT1 vs. allT2
p1p2 = roc.test(p1$ROC, p2$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p1p2)
```

```

	Bootstrap test for two correlated ROC curves

data:  p1$ROC and p2$ROC
D = -0.40644, boot.n = 2000, boot.stratified = 1, p-value = 0.6844
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.7991857   0.8048601 
```

```r
# find significants: imgT1 vs wLMSIR
p1p3 = roc.test(p1$ROC, p3$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p1p3)
```

```

	Bootstrap test for two correlated ROC curves

data:  p1$ROC and p3$ROC
D = -2.9071, boot.n = 2000, boot.stratified = 1, p-value = 0.003648
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.7991857   0.8333189 
```

```r
# find significants: imgT1 vs wLMSIR
p2p3 = roc.test(p2$ROC, p3$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p2p3)
```

```

	Bootstrap test for two correlated ROC curves

data:  p2$ROC and p3$ROC
D = -2.4914, boot.n = 2000, boot.stratified = 1, p-value = 0.01273
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8048601   0.8333189 
```
 
 
# Adjust P-values for Multiple Comparisons

The adjustment methods include the Bonferroni correction ("bonferroni") in which the p-values are multiplied by the number of comparisons. Less conservative corrections are also included by Holm (1979) ("holm"), Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"), Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), and Benjamini & Yekutieli (2001) ("BY"), respectively. 
 
The first four methods are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions.


```r
pvals = c(p1p2$p.value, p1p3$p.value, p2p3$p.value)
p.adjust(pvals, method = "bonferroni", n = length(pvals))
```

```
         D          D          D 
1.00000000 0.01094526 0.03817601 
```



# Analysis of the frequency of boosted trees used features

```r
source("functions.R")

# plot features imgT1featsel wimgT2featsel wLMSIRfeatsel
## group with all of the features spaces combined, most contributing T2w feature
# pick frequency of 75% or higher as very common feature
dfimgT1 = data.frame(table(varImpall_imgT1w_total$selfeat))
dfimgT1$high = (dfimgT1$Freq>=0.75*max(varImpall_imgT1w_total$kfcv))*1
dfimgT1 = dfimgT1[order(dfimgT1$Freq, decreasing = TRUE),]

# prune only max frequencies, remove lesser frequencies
# with dictionary assign type of feature
dictimgT1 = feature_dictionary(imgT1train)
dfreqimgT1 = data.frame()
freqf = c()
for(di in 1:nrow(dfimgT1)){
  featname = as.character(dfimgT1[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    dfreqimgT1 = rbind(dfreqimgT1, cbind(dfimgT1[di,],
                       typefeature=dictimgT1$fnnames[dictimgT1$fnnames$f==featname,"type"][1]))}}
print(dfreqimgT1)
```

```
                                      Var1 Freq high     typefeature
1                                A_countor   10    1      T1wdynamic
2                                 A_inside   10    1      T1wdynamic
3                            alpha_countor   10    1      T1wdynamic
4                             alpha_inside   10    1      T1wdynamic
5                              beta_inside   10    1      T1wdynamic
6                              circularity   10    1   T1wmorphology
8                                  dce2SE1   10    1 single-time-Enh
10                                dce2SE11   10    1 single-time-Enh
11                                dce2SE12   10    1 single-time-Enh
12                                dce2SE14   10    1 single-time-Enh
14                                dce2SE16   10    1 single-time-Enh
15                                dce2SE17   10    1 single-time-Enh
18                                 dce2SE2   10    1 single-time-Enh
19                                 dce2SE3   10    1 single-time-Enh
21                                 dce2SE5   10    1 single-time-Enh
22                                 dce2SE8   10    1 single-time-Enh
24                                 dce3SE0   10    1 single-time-Enh
26                                dce3SE10   10    1 single-time-Enh
27                                dce3SE11   10    1 single-time-Enh
28                                dce3SE14   10    1 single-time-Enh
31                                dce3SE19   10    1 single-time-Enh
32                                 dce3SE2   10    1 single-time-Enh
33                                 dce3SE3   10    1 single-time-Enh
34                                 dce3SE4   10    1 single-time-Enh
35                                 dce3SE5   10    1 single-time-Enh
36                                 dce3SE7   10    1 single-time-Enh
37                                 dce3SE8   10    1 single-time-Enh
38                                 dce3SE9   10    1 single-time-Enh
39                                earlySE0   10    1 single-time-Enh
40                                earlySE1   10    1 single-time-Enh
42                               earlySE12   10    1 single-time-Enh
43                               earlySE13   10    1 single-time-Enh
45                               earlySE15   10    1 single-time-Enh
46                               earlySE16   10    1 single-time-Enh
48                               earlySE19   10    1 single-time-Enh
49                                earlySE3   10    1 single-time-Enh
50                                earlySE5   10    1 single-time-Enh
51                                earlySE6   10    1 single-time-Enh
52                                earlySE7   10    1 single-time-Enh
53                                earlySE8   10    1 single-time-Enh
55                         edge_sharp_mean   10    1   T1wmorphology
56                          edge_sharp_std   10    1   T1wmorphology
57                            iAUC1_inside   10    1      T1wdynamic
58                  iiiMax_Margin_Gradient   10    1   T1wmorphology
59            iiMin_change_Variance_uptake   10    1   T1wmorphology
60                    iMax_Variance_uptake   10    1   T1wmorphology
61                            irregularity   10    1   T1wmorphology
62                              ivVariance   10    1   T1wmorphology
64                           Kpeak_countor   10    1      T1wdynamic
65                            Kpeak_inside   10    1      T1wdynamic
66                              kurt_F_r_i   10    1      T1wdynamic
67                                 lateSE0   10    1 single-time-Enh
68                                lateSE10   10    1 single-time-Enh
69                                lateSE11   10    1 single-time-Enh
70                                lateSE12   10    1 single-time-Enh
72                                lateSE16   10    1 single-time-Enh
73                                lateSE17   10    1 single-time-Enh
74                                lateSE19   10    1 single-time-Enh
75                                 lateSE2   10    1 single-time-Enh
76                                 lateSE4   10    1 single-time-Enh
77                                 lateSE5   10    1 single-time-Enh
79                                 lateSE9   10    1 single-time-Enh
80                               max_F_r_i   10    1      T1wdynamic
81                            max_RGH_mean   10    1   T1wmorphology
82                          max_RGH_mean_k   10    1   T1wmorphology
83                             max_RGH_var   10    1   T1wmorphology
85                            maxCr_inside   10    1      T1wdynamic
86                           maxVr_countor   10    1      T1wdynamic
87                            maxVr_inside   10    1      T1wdynamic
88                              mean_F_r_i   10    1      T1wdynamic
89                               min_F_r_i   10    1      T1wdynamic
90                             SER_countor   10    1      T1wdynamic
91                              SER_inside   10    1      T1wdynamic
92                              skew_F_r_i   10    1      T1wdynamic
93                       Slope_ini_countor   10    1      T1wdynamic
94                        Slope_ini_inside   10    1      T1wdynamic
95           texture_contrast_nondir_post2   10    1      T1wtexture
96           texture_contrast_nondir_post3   10    1      T1wtexture
97           texture_contrast_nondir_post4   10    1      T1wtexture
98        texture_correlation_nondir_post1   10    1      T1wtexture
99        texture_correlation_nondir_post2   10    1      T1wtexture
100       texture_correlation_nondir_post3   10    1      T1wtexture
102       texture_diffentropy_nondir_post2   10    1      T1wtexture
103       texture_diffentropy_nondir_post3   10    1      T1wtexture
105      texture_diffvariance_nondir_post1   10    1      T1wtexture
106      texture_diffvariance_nondir_post2   10    1      T1wtexture
107      texture_diffvariance_nondir_post3   10    1      T1wtexture
109            texture_energy_nondir_post1   10    1      T1wtexture
110            texture_energy_nondir_post2   10    1      T1wtexture
111            texture_energy_nondir_post3   10    1      T1wtexture
112           texture_entropy_nondir_post1   10    1      T1wtexture
113           texture_entropy_nondir_post2   10    1      T1wtexture
115 texture_inversediffmoment_nondir_post1   10    1      T1wtexture
116 texture_inversediffmoment_nondir_post2   10    1      T1wtexture
117 texture_inversediffmoment_nondir_post3   10    1      T1wtexture
118 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
119        texture_sumaverage_nondir_post1   10    1      T1wtexture
120        texture_sumaverage_nondir_post2   10    1      T1wtexture
121        texture_sumaverage_nondir_post3   10    1      T1wtexture
122        texture_sumaverage_nondir_post4   10    1      T1wtexture
123        texture_sumentropy_nondir_post1   10    1      T1wtexture
124        texture_sumentropy_nondir_post3   10    1      T1wtexture
125        texture_sumentropy_nondir_post4   10    1      T1wtexture
126       texture_sumvariance_nondir_post1   10    1      T1wtexture
129       texture_sumvariance_nondir_post4   10    1      T1wtexture
130          texture_variance_nondir_post1   10    1      T1wtexture
131          texture_variance_nondir_post2   10    1      T1wtexture
132          texture_variance_nondir_post3   10    1      T1wtexture
134                          Tpeak_countor   10    1      T1wdynamic
135                           Tpeak_inside   10    1      T1wdynamic
136                      UptakeRate_inside   10    1      T1wdynamic
137                                     V0   10    1      dispersion
138                                     V1   10    1      dispersion
139                                    V10   10    1      dispersion
140                                    V11   10    1      dispersion
141                                    V12   10    1      dispersion
142                                    V13   10    1      dispersion
143                                    V14   10    1      dispersion
144                                    V15   10    1      dispersion
145                                    V16   10    1      dispersion
146                                    V17   10    1      dispersion
147                                    V18   10    1      dispersion
148                                    V19   10    1      dispersion
149                                     V2   10    1      dispersion
150                                     V3   10    1      dispersion
151                                     V4   10    1      dispersion
152                                     V5   10    1      dispersion
153                                     V6   10    1      dispersion
154                                     V7   10    1      dispersion
155                                     V8   10    1      dispersion
156                                     V9   10    1      dispersion
157                              var_F_r_i   10    1      T1wdynamic
158              Vr_decreasingRate_countor   10    1      T1wdynamic
159               Vr_decreasingRate_inside   10    1      T1wdynamic
160              Vr_increasingRate_countor   10    1      T1wdynamic
161               Vr_increasingRate_inside   10    1      T1wdynamic
162                      Vr_post_1_countor   10    1      T1wdynamic
163                       Vr_post_1_inside   10    1      T1wdynamic
164                     washoutRate_inside   10    1      T1wdynamic
7                                  dce2SE0    8    1 single-time-Enh
9                                 dce2SE10    8    1 single-time-Enh
13                                dce2SE15    8    1 single-time-Enh
16                                dce2SE18    8    1 single-time-Enh
17                                dce2SE19    8    1 single-time-Enh
20                                 dce2SE4    8    1 single-time-Enh
23                                 dce2SE9    8    1 single-time-Enh
25                                 dce3SE1    8    1 single-time-Enh
29                                dce3SE15    8    1 single-time-Enh
30                                dce3SE18    8    1 single-time-Enh
41                               earlySE11    8    1 single-time-Enh
44                               earlySE14    8    1 single-time-Enh
47                               earlySE17    8    1 single-time-Enh
54                                earlySE9    8    1 single-time-Enh
63                       k_Max_Margin_Grad    8    1   T1wmorphology
71                                lateSE15    8    1 single-time-Enh
78                                 lateSE7    8    1 single-time-Enh
84                           max_RGH_var_k    8    1   T1wmorphology
104       texture_diffentropy_nondir_post4    8    1      T1wtexture
108      texture_diffvariance_nondir_post4    8    1      T1wtexture
114           texture_entropy_nondir_post3    8    1      T1wtexture
128       texture_sumvariance_nondir_post3    8    1      T1wtexture
165                           beta_countor    8    1      T1wdynamic
167                                dce2SE6    8    1 single-time-Enh
168                                dce2SE7    8    1 single-time-Enh
171                               dce3SE16    8    1 single-time-Enh
173                                dce3SE6    8    1 single-time-Enh
174                              earlySE10    8    1 single-time-Enh
175                              earlySE18    8    1 single-time-Enh
176                               earlySE2    8    1 single-time-Enh
177                               earlySE4    8    1 single-time-Enh
178                          iAUC1_countor    8    1      T1wdynamic
179                                lateSE1    8    1 single-time-Enh
180                               lateSE13    8    1 single-time-Enh
182                               lateSE18    8    1 single-time-Enh
183                                lateSE3    8    1 single-time-Enh
184                                lateSE6    8    1 single-time-Enh
185                                lateSE8    8    1 single-time-Enh
188                         peakVr_countor    8    1      T1wdynamic
190          texture_contrast_nondir_post1    8    1      T1wtexture
191       texture_diffentropy_nondir_post1    8    1      T1wtexture
192            texture_energy_nondir_post4    8    1      T1wtexture
194        texture_sumentropy_nondir_post2    8    1      T1wtexture
195                     UptakeRate_countor    8    1      T1wdynamic
196                    washoutRate_countor    8    1      T1wdynamic
101       texture_correlation_nondir_post4    6    0      T1wtexture
127       texture_sumvariance_nondir_post2    6    0      T1wtexture
133          texture_variance_nondir_post4    6    0      T1wtexture
166                               dce2SE13    6    0 single-time-Enh
169                               dce3SE12    6    0 single-time-Enh
172                               dce3SE17    6    0 single-time-Enh
181                               lateSE14    6    0 single-time-Enh
186                          maxCr_countor    6    0      T1wdynamic
187                          peakCr_inside    6    0      T1wdynamic
193           texture_entropy_nondir_post4    6    0      T1wtexture
170                               dce3SE13    4    0 single-time-Enh
189                          peakVr_inside    4    0      T1wdynamic
197                         peakCr_countor    4    0      T1wdynamic
```

```r
#plot
ggplot(dfreqimgT1[dfreqimgT1$high==1,], aes(x=reorder(Var1, order(Freq, decreasing = TRUE)), y=Freq, fill=factor(typefeature))) + 
  geom_bar(stat = "identity") + coord_flip() +
  ggtitle("T1wonly ensembles: feature selection frequency") +
  labs(y="featsel frequency", x=" ") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=18, hjust=0))+
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16))
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/featuressel-1.png)<!-- -->

```r
## varImpall_T1wT2w_total ########### 
# pick frequency of 75% or higher as very common feature
dfallT2 = data.frame(table(varImpall_T1wT2w_total$selfeat))
dfallT2$high = (dfallT2$Freq>=0.75*max(varImpall_T1wT2w_total$kfcv))*1
dfallT2 = dfallT2[order(dfallT2$Freq, decreasing = TRUE),]
print(dfallT2)
```

```
                                      Var1 Freq high
1                                A_countor   10    1
4                             alpha_inside   10    1
5                                  ave_T20   10    1
7                                 ave_T210   10    1
10                                ave_T213   10    1
12                                ave_T215   10    1
18                                 ave_T24   10    1
19                                 ave_T25   10    1
20                                 ave_T26   10    1
22                                 ave_T28   10    1
25                             beta_inside   10    1
26                             circularity   10    1
27                                 dce2SE0   10    1
34                                dce2SE17   10    1
42                                 dce2SE7   10    1
45                                 dce3SE0   10    1
48                                dce3SE11   10    1
53                                dce3SE17   10    1
59                                 dce3SE6   10    1
63                                earlySE0   10    1
69                               earlySE15   10    1
71                               earlySE19   10    1
75                                earlySE6   10    1
77                                earlySE8   10    1
78                                earlySE9   10    1
80                          edge_sharp_std   10    1
81                           iAUC1_countor   10    1
83                  iiiMax_Margin_Gradient   10    1
84            iiMin_change_Variance_uptake   10    1
86                            irregularity   10    1
87                              ivVariance   10    1
89                            Kpeak_inside   10    1
90                              kurt_F_r_i   10    1
91                                 lateSE0   10    1
93                                lateSE10   10    1
94                                lateSE11   10    1
99                                lateSE17   10    1
102                                lateSE2   10    1
103                                lateSE3   10    1
106                                lateSE6   10    1
110                              max_F_r_i   10    1
111                           max_RGH_mean   10    1
113                            max_RGH_var   10    1
117                           maxVr_inside   10    1
118                             mean_F_r_i   10    1
121                            SER_countor   10    1
122                             SER_inside   10    1
123                             skew_F_r_i   10    1
127                         T2_lesionSIstd   10    1
129                      T2grad_margin_var   10    1
130                           T2kurt_F_r_i   10    1
131                            T2max_F_r_i   10    1
133                             T2RGH_mean   10    1
134                              T2RGH_var   10    1
135                           T2skew_F_r_i   10    1
137           T2texture_correlation_nondir   10    1
139          T2texture_diffvariance_nondir   10    1
141               T2texture_entropy_nondir   10    1
142     T2texture_inversediffmoment_nondir   10    1
143            T2texture_sumaverage_nondir   10    1
146              T2texture_variance_nondir   10    1
149          texture_contrast_nondir_post3   10    1
152       texture_correlation_nondir_post2   10    1
153       texture_correlation_nondir_post3   10    1
159      texture_diffvariance_nondir_post1   10    1
160      texture_diffvariance_nondir_post2   10    1
163            texture_energy_nondir_post1   10    1
164            texture_energy_nondir_post3   10    1
167           texture_entropy_nondir_post2   10    1
170 texture_inversediffmoment_nondir_post2   10    1
171 texture_inversediffmoment_nondir_post3   10    1
172 texture_inversediffmoment_nondir_post4   10    1
173        texture_sumaverage_nondir_post1   10    1
174        texture_sumaverage_nondir_post2   10    1
176        texture_sumaverage_nondir_post4   10    1
180       texture_sumvariance_nondir_post1   10    1
184          texture_variance_nondir_post1   10    1
187                          Tpeak_countor   10    1
188                           Tpeak_inside   10    1
189                      UptakeRate_inside   10    1
190                                     V0   10    1
192                                    V10   10    1
193                                    V11   10    1
194                                    V12   10    1
196                                    V14   10    1
198                                    V16   10    1
199                                    V17   10    1
201                                    V19   10    1
202                                     V2   10    1
203                                     V3   10    1
204                                     V4   10    1
205                                     V5   10    1
206                                     V6   10    1
207                                     V7   10    1
208                                     V8   10    1
209                                     V9   10    1
210                              var_F_r_i   10    1
211              Vr_decreasingRate_countor   10    1
212               Vr_decreasingRate_inside   10    1
216                       Vr_post_1_inside   10    1
2                                 A_inside    8    1
3                            alpha_countor    8    1
6                                  ave_T21    8    1
8                                 ave_T211    8    1
9                                 ave_T212    8    1
11                                ave_T214    8    1
15                                ave_T219    8    1
16                                 ave_T22    8    1
17                                 ave_T23    8    1
21                                 ave_T27    8    1
23                                 ave_T29    8    1
24                            beta_countor    8    1
29                                dce2SE11    8    1
30                                dce2SE12    8    1
35                                dce2SE18    8    1
36                                dce2SE19    8    1
37                                 dce2SE2    8    1
38                                 dce2SE3    8    1
40                                 dce2SE5    8    1
41                                 dce2SE6    8    1
43                                 dce2SE8    8    1
46                                 dce3SE1    8    1
50                                dce3SE14    8    1
51                                dce3SE15    8    1
55                                dce3SE19    8    1
57                                 dce3SE3    8    1
58                                 dce3SE4    8    1
60                                 dce3SE7    8    1
61                                 dce3SE8    8    1
62                                 dce3SE9    8    1
64                                earlySE1    8    1
65                               earlySE10    8    1
66                               earlySE12    8    1
67                               earlySE13    8    1
70                               earlySE16    8    1
73                                earlySE4    8    1
76                                earlySE7    8    1
79                         edge_sharp_mean    8    1
82                            iAUC1_inside    8    1
85                    iMax_Variance_uptake    8    1
95                                lateSE12    8    1
97                                lateSE14    8    1
98                                lateSE15    8    1
100                               lateSE18    8    1
104                                lateSE4    8    1
105                                lateSE5    8    1
107                                lateSE7    8    1
108                                lateSE8    8    1
109                                lateSE9    8    1
112                         max_RGH_mean_k    8    1
119                              min_F_r_i    8    1
120                         peakVr_countor    8    1
125                       Slope_ini_inside    8    1
126                            T2_lesionSI    8    1
140                T2texture_energy_nondir    8    1
144            T2texture_sumentropy_nondir    8    1
147          texture_contrast_nondir_post1    8    1
148          texture_contrast_nondir_post2    8    1
151       texture_correlation_nondir_post1    8    1
154       texture_correlation_nondir_post4    8    1
161      texture_diffvariance_nondir_post3    8    1
165            texture_energy_nondir_post4    8    1
166           texture_entropy_nondir_post1    8    1
169 texture_inversediffmoment_nondir_post1    8    1
175        texture_sumaverage_nondir_post3    8    1
177        texture_sumentropy_nondir_post1    8    1
178        texture_sumentropy_nondir_post2    8    1
179        texture_sumentropy_nondir_post3    8    1
182       texture_sumvariance_nondir_post3    8    1
185          texture_variance_nondir_post2    8    1
191                                     V1    8    1
195                                    V13    8    1
197                                    V15    8    1
200                                    V18    8    1
213              Vr_increasingRate_countor    8    1
214               Vr_increasingRate_inside    8    1
215                      Vr_post_1_countor    8    1
223                                dce3SE5    8    1
225                              earlySE18    8    1
228            texture_energy_nondir_post2    8    1
13                                ave_T216    6    0
14                                ave_T217    6    0
28                                 dce2SE1    6    0
31                                dce2SE13    6    0
32                                dce2SE15    6    0
33                                dce2SE16    6    0
39                                 dce2SE4    6    0
47                                dce3SE10    6    0
49                                dce3SE13    6    0
54                                dce3SE18    6    0
56                                 dce3SE2    6    0
72                                earlySE3    6    0
74                                earlySE5    6    0
88                           Kpeak_countor    6    0
92                                 lateSE1    6    0
96                                lateSE13    6    0
101                               lateSE19    6    0
115                           maxCr_inside    6    0
116                          maxVr_countor    6    0
124                      Slope_ini_countor    6    0
128                          T2grad_margin    6    0
136              T2texture_contrast_nondir    6    0
145           T2texture_sumvariance_nondir    6    0
150          texture_contrast_nondir_post4    6    0
155       texture_diffentropy_nondir_post1    6    0
156       texture_diffentropy_nondir_post2    6    0
157       texture_diffentropy_nondir_post3    6    0
162      texture_diffvariance_nondir_post4    6    0
168           texture_entropy_nondir_post3    6    0
181       texture_sumvariance_nondir_post2    6    0
183       texture_sumvariance_nondir_post4    6    0
186          texture_variance_nondir_post4    6    0
217                    washoutRate_countor    6    0
218                     washoutRate_inside    6    0
219                               ave_T218    6    0
220                               dce2SE10    6    0
222                               dce3SE12    6    0
226                               earlySE2    6    0
227                               lateSE16    6    0
229          texture_variance_nondir_post3    6    0
234        texture_sumentropy_nondir_post4    6    0
235                     UptakeRate_countor    6    0
44                                 dce2SE9    4    0
52                                dce3SE16    4    0
68                               earlySE14    4    0
114                          max_RGH_var_k    4    0
132                            T2min_F_r_i    4    0
138           T2texture_diffentropy_nondir    4    0
158       texture_diffentropy_nondir_post4    4    0
221                               dce2SE14    4    0
224                              earlySE17    4    0
230                              earlySE11    4    0
231                      k_Max_Margin_Grad    4    0
232                          maxCr_countor    4    0
233           texture_entropy_nondir_post4    4    0
236                         peakCr_countor    2    0
237                          peakVr_inside    2    0
```

```r
# prune only max frequencies, remove lesser frequencies
# with dictionary assign type of feature
dictimgT2 = feature_dictionary(imgT1T2)
dfreqimgT2 = data.frame()
freqf = c()
for(di in 1:nrow(dfallT2)){
  featname = as.character(dfallT2[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    dfreqimgT2 = rbind(dfreqimgT2, cbind(dfallT2[di,],
                       typefeature=dictimgT2$fnnames[dictimgT2$fnnames$f==featname,"type"][1]))}}
print(dfreqimgT2)
```

```
                                      Var1 Freq high     typefeature
1                                A_countor   10    1      T1wdynamic
4                             alpha_inside   10    1      T1wdynamic
5                                  ave_T20   10    1             T2w
7                                 ave_T210   10    1             T2w
10                                ave_T213   10    1             T2w
12                                ave_T215   10    1             T2w
18                                 ave_T24   10    1             T2w
19                                 ave_T25   10    1             T2w
20                                 ave_T26   10    1             T2w
22                                 ave_T28   10    1             T2w
25                             beta_inside   10    1      T1wdynamic
26                             circularity   10    1   T1wmorphology
27                                 dce2SE0   10    1 single-time-Enh
34                                dce2SE17   10    1 single-time-Enh
42                                 dce2SE7   10    1 single-time-Enh
45                                 dce3SE0   10    1 single-time-Enh
48                                dce3SE11   10    1 single-time-Enh
53                                dce3SE17   10    1 single-time-Enh
59                                 dce3SE6   10    1 single-time-Enh
63                                earlySE0   10    1 single-time-Enh
69                               earlySE15   10    1 single-time-Enh
71                               earlySE19   10    1 single-time-Enh
75                                earlySE6   10    1 single-time-Enh
77                                earlySE8   10    1 single-time-Enh
78                                earlySE9   10    1 single-time-Enh
80                          edge_sharp_std   10    1   T1wmorphology
81                           iAUC1_countor   10    1      T1wdynamic
83                  iiiMax_Margin_Gradient   10    1   T1wmorphology
84            iiMin_change_Variance_uptake   10    1   T1wmorphology
86                            irregularity   10    1   T1wmorphology
87                              ivVariance   10    1   T1wmorphology
89                            Kpeak_inside   10    1      T1wdynamic
90                              kurt_F_r_i   10    1      T1wdynamic
91                                 lateSE0   10    1 single-time-Enh
93                                lateSE10   10    1 single-time-Enh
94                                lateSE11   10    1 single-time-Enh
99                                lateSE17   10    1 single-time-Enh
102                                lateSE2   10    1 single-time-Enh
103                                lateSE3   10    1 single-time-Enh
106                                lateSE6   10    1 single-time-Enh
110                              max_F_r_i   10    1      T1wdynamic
111                           max_RGH_mean   10    1   T1wmorphology
113                            max_RGH_var   10    1   T1wmorphology
117                           maxVr_inside   10    1      T1wdynamic
118                             mean_F_r_i   10    1      T1wdynamic
121                            SER_countor   10    1      T1wdynamic
122                             SER_inside   10    1      T1wdynamic
123                             skew_F_r_i   10    1      T1wdynamic
127                         T2_lesionSIstd   10    1             T2w
129                      T2grad_margin_var   10    1             T2w
130                           T2kurt_F_r_i   10    1             T2w
131                            T2max_F_r_i   10    1             T2w
133                             T2RGH_mean   10    1             T2w
134                              T2RGH_var   10    1             T2w
135                           T2skew_F_r_i   10    1             T2w
137           T2texture_correlation_nondir   10    1             T2w
139          T2texture_diffvariance_nondir   10    1             T2w
141               T2texture_entropy_nondir   10    1             T2w
142     T2texture_inversediffmoment_nondir   10    1             T2w
143            T2texture_sumaverage_nondir   10    1             T2w
146              T2texture_variance_nondir   10    1             T2w
149          texture_contrast_nondir_post3   10    1      T1wtexture
152       texture_correlation_nondir_post2   10    1      T1wtexture
153       texture_correlation_nondir_post3   10    1      T1wtexture
159      texture_diffvariance_nondir_post1   10    1      T1wtexture
160      texture_diffvariance_nondir_post2   10    1      T1wtexture
163            texture_energy_nondir_post1   10    1      T1wtexture
164            texture_energy_nondir_post3   10    1      T1wtexture
167           texture_entropy_nondir_post2   10    1      T1wtexture
170 texture_inversediffmoment_nondir_post2   10    1      T1wtexture
171 texture_inversediffmoment_nondir_post3   10    1      T1wtexture
172 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
173        texture_sumaverage_nondir_post1   10    1      T1wtexture
174        texture_sumaverage_nondir_post2   10    1      T1wtexture
176        texture_sumaverage_nondir_post4   10    1      T1wtexture
180       texture_sumvariance_nondir_post1   10    1      T1wtexture
184          texture_variance_nondir_post1   10    1      T1wtexture
187                          Tpeak_countor   10    1      T1wdynamic
188                           Tpeak_inside   10    1      T1wdynamic
189                      UptakeRate_inside   10    1      T1wdynamic
190                                     V0   10    1      dispersion
192                                    V10   10    1      dispersion
193                                    V11   10    1      dispersion
194                                    V12   10    1      dispersion
196                                    V14   10    1      dispersion
198                                    V16   10    1      dispersion
199                                    V17   10    1      dispersion
201                                    V19   10    1      dispersion
202                                     V2   10    1      dispersion
203                                     V3   10    1      dispersion
204                                     V4   10    1      dispersion
205                                     V5   10    1      dispersion
206                                     V6   10    1      dispersion
207                                     V7   10    1      dispersion
208                                     V8   10    1      dispersion
209                                     V9   10    1      dispersion
210                              var_F_r_i   10    1      T1wdynamic
211              Vr_decreasingRate_countor   10    1      T1wdynamic
212               Vr_decreasingRate_inside   10    1      T1wdynamic
216                       Vr_post_1_inside   10    1      T1wdynamic
2                                 A_inside    8    1      T1wdynamic
3                            alpha_countor    8    1      T1wdynamic
6                                  ave_T21    8    1             T2w
8                                 ave_T211    8    1             T2w
9                                 ave_T212    8    1             T2w
11                                ave_T214    8    1             T2w
15                                ave_T219    8    1             T2w
16                                 ave_T22    8    1             T2w
17                                 ave_T23    8    1             T2w
21                                 ave_T27    8    1             T2w
23                                 ave_T29    8    1             T2w
24                            beta_countor    8    1      T1wdynamic
29                                dce2SE11    8    1 single-time-Enh
30                                dce2SE12    8    1 single-time-Enh
35                                dce2SE18    8    1 single-time-Enh
36                                dce2SE19    8    1 single-time-Enh
37                                 dce2SE2    8    1 single-time-Enh
38                                 dce2SE3    8    1 single-time-Enh
40                                 dce2SE5    8    1 single-time-Enh
41                                 dce2SE6    8    1 single-time-Enh
43                                 dce2SE8    8    1 single-time-Enh
46                                 dce3SE1    8    1 single-time-Enh
50                                dce3SE14    8    1 single-time-Enh
51                                dce3SE15    8    1 single-time-Enh
55                                dce3SE19    8    1 single-time-Enh
57                                 dce3SE3    8    1 single-time-Enh
58                                 dce3SE4    8    1 single-time-Enh
60                                 dce3SE7    8    1 single-time-Enh
61                                 dce3SE8    8    1 single-time-Enh
62                                 dce3SE9    8    1 single-time-Enh
64                                earlySE1    8    1 single-time-Enh
65                               earlySE10    8    1 single-time-Enh
66                               earlySE12    8    1 single-time-Enh
67                               earlySE13    8    1 single-time-Enh
70                               earlySE16    8    1 single-time-Enh
73                                earlySE4    8    1 single-time-Enh
76                                earlySE7    8    1 single-time-Enh
79                         edge_sharp_mean    8    1   T1wmorphology
82                            iAUC1_inside    8    1      T1wdynamic
85                    iMax_Variance_uptake    8    1   T1wmorphology
95                                lateSE12    8    1 single-time-Enh
97                                lateSE14    8    1 single-time-Enh
98                                lateSE15    8    1 single-time-Enh
100                               lateSE18    8    1 single-time-Enh
104                                lateSE4    8    1 single-time-Enh
105                                lateSE5    8    1 single-time-Enh
107                                lateSE7    8    1 single-time-Enh
108                                lateSE8    8    1 single-time-Enh
109                                lateSE9    8    1 single-time-Enh
112                         max_RGH_mean_k    8    1   T1wmorphology
119                              min_F_r_i    8    1      T1wdynamic
120                         peakVr_countor    8    1      T1wdynamic
125                       Slope_ini_inside    8    1      T1wdynamic
126                            T2_lesionSI    8    1             T2w
140                T2texture_energy_nondir    8    1             T2w
144            T2texture_sumentropy_nondir    8    1             T2w
147          texture_contrast_nondir_post1    8    1      T1wtexture
148          texture_contrast_nondir_post2    8    1      T1wtexture
151       texture_correlation_nondir_post1    8    1      T1wtexture
154       texture_correlation_nondir_post4    8    1      T1wtexture
161      texture_diffvariance_nondir_post3    8    1      T1wtexture
165            texture_energy_nondir_post4    8    1      T1wtexture
166           texture_entropy_nondir_post1    8    1      T1wtexture
169 texture_inversediffmoment_nondir_post1    8    1      T1wtexture
175        texture_sumaverage_nondir_post3    8    1      T1wtexture
177        texture_sumentropy_nondir_post1    8    1      T1wtexture
178        texture_sumentropy_nondir_post2    8    1      T1wtexture
179        texture_sumentropy_nondir_post3    8    1      T1wtexture
182       texture_sumvariance_nondir_post3    8    1      T1wtexture
185          texture_variance_nondir_post2    8    1      T1wtexture
191                                     V1    8    1      dispersion
195                                    V13    8    1      dispersion
197                                    V15    8    1      dispersion
200                                    V18    8    1      dispersion
213              Vr_increasingRate_countor    8    1      T1wdynamic
214               Vr_increasingRate_inside    8    1      T1wdynamic
215                      Vr_post_1_countor    8    1      T1wdynamic
223                                dce3SE5    8    1 single-time-Enh
225                              earlySE18    8    1 single-time-Enh
228            texture_energy_nondir_post2    8    1      T1wtexture
13                                ave_T216    6    0             T2w
14                                ave_T217    6    0             T2w
28                                 dce2SE1    6    0 single-time-Enh
31                                dce2SE13    6    0 single-time-Enh
32                                dce2SE15    6    0 single-time-Enh
33                                dce2SE16    6    0 single-time-Enh
39                                 dce2SE4    6    0 single-time-Enh
47                                dce3SE10    6    0 single-time-Enh
49                                dce3SE13    6    0 single-time-Enh
54                                dce3SE18    6    0 single-time-Enh
56                                 dce3SE2    6    0 single-time-Enh
72                                earlySE3    6    0 single-time-Enh
74                                earlySE5    6    0 single-time-Enh
88                           Kpeak_countor    6    0      T1wdynamic
92                                 lateSE1    6    0 single-time-Enh
96                                lateSE13    6    0 single-time-Enh
101                               lateSE19    6    0 single-time-Enh
115                           maxCr_inside    6    0      T1wdynamic
116                          maxVr_countor    6    0      T1wdynamic
124                      Slope_ini_countor    6    0      T1wdynamic
128                          T2grad_margin    6    0             T2w
136              T2texture_contrast_nondir    6    0             T2w
145           T2texture_sumvariance_nondir    6    0             T2w
150          texture_contrast_nondir_post4    6    0      T1wtexture
155       texture_diffentropy_nondir_post1    6    0      T1wtexture
156       texture_diffentropy_nondir_post2    6    0      T1wtexture
157       texture_diffentropy_nondir_post3    6    0      T1wtexture
162      texture_diffvariance_nondir_post4    6    0      T1wtexture
168           texture_entropy_nondir_post3    6    0      T1wtexture
181       texture_sumvariance_nondir_post2    6    0      T1wtexture
183       texture_sumvariance_nondir_post4    6    0      T1wtexture
186          texture_variance_nondir_post4    6    0      T1wtexture
217                    washoutRate_countor    6    0      T1wdynamic
218                     washoutRate_inside    6    0      T1wdynamic
219                               ave_T218    6    0             T2w
220                               dce2SE10    6    0 single-time-Enh
222                               dce3SE12    6    0 single-time-Enh
226                               earlySE2    6    0 single-time-Enh
227                               lateSE16    6    0 single-time-Enh
229          texture_variance_nondir_post3    6    0      T1wtexture
234        texture_sumentropy_nondir_post4    6    0      T1wtexture
235                     UptakeRate_countor    6    0      T1wdynamic
44                                 dce2SE9    4    0 single-time-Enh
52                                dce3SE16    4    0 single-time-Enh
68                               earlySE14    4    0 single-time-Enh
114                          max_RGH_var_k    4    0   T1wmorphology
132                            T2min_F_r_i    4    0             T2w
138           T2texture_diffentropy_nondir    4    0             T2w
158       texture_diffentropy_nondir_post4    4    0      T1wtexture
221                               dce2SE14    4    0 single-time-Enh
224                              earlySE17    4    0 single-time-Enh
230                              earlySE11    4    0 single-time-Enh
231                      k_Max_Margin_Grad    4    0   T1wmorphology
232                          maxCr_countor    4    0      T1wdynamic
233           texture_entropy_nondir_post4    4    0      T1wtexture
236                         peakCr_countor    2    0      T1wdynamic
237                          peakVr_inside    2    0      T1wdynamic
```

```r
#plot
ggplot(dfreqimgT2, aes(x=reorder(Var1, order(Freq, decreasing = TRUE)), y=Freq, fill=factor(typefeature))) + 
  geom_bar(stat = "identity") + coord_flip() +
  ggtitle("T1wonly ensembles: feature selection frequency") +
  labs(y="featsel frequency", x=" ") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=18, hjust=0))+
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16))
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/featuressel-2.png)<!-- -->

```r
## varImpall_T2wpLMSIR_total ########### 
# pick frequency of 75% or higher as very common feature
dfwT1T2 = data.frame(table(varImpall_T2wpLMSIR_total$selfeat))
dfwT1T2$high = (dfwT1T2$Freq>=0.75*max(varImpall_T2wpLMSIR_total$kfcv))*1
dfwT1T2 = dfwT1T2[order(dfwT1T2$Freq, decreasing = TRUE),]
print(dfwT1T2)
```

```
                                      Var1 Freq high
1                                A_countor   10    1
2                                 A_inside   10    1
3                            alpha_countor   10    1
4                             alpha_inside   10    1
5                             beta_countor   10    1
6                              beta_inside   10    1
7                              circularity   10    1
8                                  dce2SE0   10    1
9                                  dce2SE1   10    1
11                                dce2SE11   10    1
12                                dce2SE12   10    1
13                                dce2SE13   10    1
14                                dce2SE14   10    1
15                                dce2SE16   10    1
16                                dce2SE17   10    1
17                                dce2SE18   10    1
19                                 dce2SE2   10    1
20                                 dce2SE4   10    1
21                                 dce2SE5   10    1
22                                 dce2SE7   10    1
23                                 dce2SE8   10    1
24                                 dce2SE9   10    1
25                                 dce3SE0   10    1
26                                 dce3SE1   10    1
27                                dce3SE10   10    1
28                                dce3SE11   10    1
29                                dce3SE12   10    1
30                                dce3SE13   10    1
31                                dce3SE14   10    1
32                                dce3SE15   10    1
33                                dce3SE17   10    1
35                                 dce3SE3   10    1
36                                 dce3SE5   10    1
37                                 dce3SE6   10    1
38                                 dce3SE7   10    1
39                                 dce3SE8   10    1
40                                 dce3SE9   10    1
41                                earlySE0   10    1
42                                earlySE1   10    1
43                               earlySE10   10    1
44                               earlySE11   10    1
45                               earlySE12   10    1
46                               earlySE15   10    1
47                               earlySE16   10    1
48                               earlySE17   10    1
49                               earlySE18   10    1
50                               earlySE19   10    1
51                                earlySE2   10    1
52                                earlySE3   10    1
53                                earlySE5   10    1
54                                earlySE6   10    1
55                                earlySE7   10    1
56                                earlySE8   10    1
57                                earlySE9   10    1
58                         edge_sharp_mean   10    1
59                          edge_sharp_std   10    1
60                           iAUC1_countor   10    1
61                            iAUC1_inside   10    1
62                  iiiMax_Margin_Gradient   10    1
63            iiMin_change_Variance_uptake   10    1
64                    iMax_Variance_uptake   10    1
65                            irregularity   10    1
66                              ivVariance   10    1
67                       k_Max_Margin_Grad   10    1
68                           Kpeak_countor   10    1
69                            Kpeak_inside   10    1
70                              kurt_F_r_i   10    1
71                                 lateSE0   10    1
72                                 lateSE1   10    1
73                                lateSE10   10    1
74                                lateSE11   10    1
75                                lateSE12   10    1
76                                lateSE13   10    1
77                                lateSE14   10    1
79                                lateSE16   10    1
80                                lateSE17   10    1
81                                lateSE18   10    1
83                                 lateSE2   10    1
84                                 lateSE3   10    1
85                                 lateSE4   10    1
86                                 lateSE7   10    1
87                                 lateSE8   10    1
88                                 lateSE9   10    1
89                         LMSIR_predicted   10    1
90                               max_F_r_i   10    1
91                            max_RGH_mean   10    1
92                          max_RGH_mean_k   10    1
93                             max_RGH_var   10    1
94                           max_RGH_var_k   10    1
95                           maxCr_countor   10    1
96                            maxCr_inside   10    1
97                           maxVr_countor   10    1
98                            maxVr_inside   10    1
99                              mean_F_r_i   10    1
100                              min_F_r_i   10    1
101                         peakVr_countor   10    1
102                          peakVr_inside   10    1
103                            SER_countor   10    1
104                             SER_inside   10    1
105                             skew_F_r_i   10    1
106                      Slope_ini_countor   10    1
107                       Slope_ini_inside   10    1
108                          T2grad_margin   10    1
109                      T2grad_margin_var   10    1
110                           T2kurt_F_r_i   10    1
111                            T2max_F_r_i   10    1
112                           T2mean_F_r_i   10    1
113                            T2min_F_r_i   10    1
114                             T2RGH_mean   10    1
115                              T2RGH_var   10    1
116                           T2skew_F_r_i   10    1
117              T2texture_contrast_nondir   10    1
118           T2texture_correlation_nondir   10    1
119           T2texture_diffentropy_nondir   10    1
120          T2texture_diffvariance_nondir   10    1
121                T2texture_energy_nondir   10    1
122               T2texture_entropy_nondir   10    1
123     T2texture_inversediffmoment_nondir   10    1
124            T2texture_sumaverage_nondir   10    1
125            T2texture_sumentropy_nondir   10    1
126           T2texture_sumvariance_nondir   10    1
127              T2texture_variance_nondir   10    1
128                            T2var_F_r_i   10    1
129          texture_contrast_nondir_post1   10    1
130          texture_contrast_nondir_post2   10    1
131          texture_contrast_nondir_post3   10    1
132          texture_contrast_nondir_post4   10    1
133       texture_correlation_nondir_post1   10    1
134       texture_correlation_nondir_post2   10    1
135       texture_correlation_nondir_post3   10    1
136       texture_correlation_nondir_post4   10    1
137       texture_diffentropy_nondir_post1   10    1
138       texture_diffentropy_nondir_post2   10    1
139       texture_diffentropy_nondir_post3   10    1
141      texture_diffvariance_nondir_post1   10    1
142      texture_diffvariance_nondir_post2   10    1
143      texture_diffvariance_nondir_post3   10    1
144      texture_diffvariance_nondir_post4   10    1
145            texture_energy_nondir_post1   10    1
146            texture_energy_nondir_post2   10    1
147            texture_energy_nondir_post3   10    1
148            texture_energy_nondir_post4   10    1
149           texture_entropy_nondir_post1   10    1
150           texture_entropy_nondir_post2   10    1
151           texture_entropy_nondir_post3   10    1
152           texture_entropy_nondir_post4   10    1
153 texture_inversediffmoment_nondir_post1   10    1
154 texture_inversediffmoment_nondir_post2   10    1
155 texture_inversediffmoment_nondir_post3   10    1
156 texture_inversediffmoment_nondir_post4   10    1
157        texture_sumaverage_nondir_post1   10    1
158        texture_sumaverage_nondir_post2   10    1
159        texture_sumaverage_nondir_post3   10    1
160        texture_sumaverage_nondir_post4   10    1
161        texture_sumentropy_nondir_post1   10    1
162        texture_sumentropy_nondir_post2   10    1
163        texture_sumentropy_nondir_post3   10    1
164        texture_sumentropy_nondir_post4   10    1
165       texture_sumvariance_nondir_post1   10    1
166       texture_sumvariance_nondir_post2   10    1
167       texture_sumvariance_nondir_post3   10    1
168       texture_sumvariance_nondir_post4   10    1
169          texture_variance_nondir_post1   10    1
170          texture_variance_nondir_post2   10    1
171          texture_variance_nondir_post3   10    1
172          texture_variance_nondir_post4   10    1
173                          Tpeak_countor   10    1
174                           Tpeak_inside   10    1
175                     UptakeRate_countor   10    1
176                      UptakeRate_inside   10    1
177                                     V0   10    1
178                                     V1   10    1
179                                    V10   10    1
180                                    V11   10    1
181                                    V12   10    1
182                                    V13   10    1
183                                    V14   10    1
184                                    V15   10    1
185                                    V16   10    1
186                                    V17   10    1
187                                    V18   10    1
188                                    V19   10    1
189                                     V2   10    1
190                                     V3   10    1
191                                     V4   10    1
192                                     V5   10    1
193                                     V6   10    1
194                                     V7   10    1
195                                     V8   10    1
196                                     V9   10    1
197                              var_F_r_i   10    1
198              Vr_decreasingRate_countor   10    1
199               Vr_decreasingRate_inside   10    1
200              Vr_increasingRate_countor   10    1
201               Vr_increasingRate_inside   10    1
202                      Vr_post_1_countor   10    1
203                       Vr_post_1_inside   10    1
204                    washoutRate_countor   10    1
10                                dce2SE10    8    1
18                                dce2SE19    8    1
34                                dce3SE18    8    1
78                                lateSE15    8    1
82                                lateSE19    8    1
140       texture_diffentropy_nondir_post4    8    1
205                     washoutRate_inside    8    1
206                               dce2SE15    8    1
207                                dce2SE3    8    1
208                                dce2SE6    8    1
209                               dce3SE16    8    1
210                               dce3SE19    8    1
211                                dce3SE2    8    1
212                                dce3SE4    8    1
213                              earlySE13    8    1
214                              earlySE14    8    1
215                               earlySE4    8    1
216                                lateSE5    8    1
217                                lateSE6    8    1
218                          peakCr_inside    8    1
219                         peakCr_countor    4    0
```

```r
# prune only max frequencies, remove lesser frequencies
# with dictionary assign type of feature
dictimgT1T2 = feature_dictionary(imgT2pLMSIR)
dfreqimgT1T2 = data.frame()
freqf = c()
for(di in 1:nrow(dfwT1T2)){
  featname = as.character(dfwT1T2[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    dfreqimgT1T2 = rbind(dfreqimgT1T2, cbind(dfwT1T2[di,],
                       typefeature=dictimgT1T2$fnnames[dictimgT1T2$fnnames$f==featname,"type"][1]))}
  }
print(dfreqimgT1T2)
```

```
                                      Var1 Freq high     typefeature
1                                A_countor   10    1      T1wdynamic
2                                 A_inside   10    1      T1wdynamic
3                            alpha_countor   10    1      T1wdynamic
4                             alpha_inside   10    1      T1wdynamic
5                             beta_countor   10    1      T1wdynamic
6                              beta_inside   10    1      T1wdynamic
7                              circularity   10    1   T1wmorphology
8                                  dce2SE0   10    1 single-time-Enh
9                                  dce2SE1   10    1 single-time-Enh
11                                dce2SE11   10    1 single-time-Enh
12                                dce2SE12   10    1 single-time-Enh
13                                dce2SE13   10    1 single-time-Enh
14                                dce2SE14   10    1 single-time-Enh
15                                dce2SE16   10    1 single-time-Enh
16                                dce2SE17   10    1 single-time-Enh
17                                dce2SE18   10    1 single-time-Enh
19                                 dce2SE2   10    1 single-time-Enh
20                                 dce2SE4   10    1 single-time-Enh
21                                 dce2SE5   10    1 single-time-Enh
22                                 dce2SE7   10    1 single-time-Enh
23                                 dce2SE8   10    1 single-time-Enh
24                                 dce2SE9   10    1 single-time-Enh
25                                 dce3SE0   10    1 single-time-Enh
26                                 dce3SE1   10    1 single-time-Enh
27                                dce3SE10   10    1 single-time-Enh
28                                dce3SE11   10    1 single-time-Enh
29                                dce3SE12   10    1 single-time-Enh
30                                dce3SE13   10    1 single-time-Enh
31                                dce3SE14   10    1 single-time-Enh
32                                dce3SE15   10    1 single-time-Enh
33                                dce3SE17   10    1 single-time-Enh
35                                 dce3SE3   10    1 single-time-Enh
36                                 dce3SE5   10    1 single-time-Enh
37                                 dce3SE6   10    1 single-time-Enh
38                                 dce3SE7   10    1 single-time-Enh
39                                 dce3SE8   10    1 single-time-Enh
40                                 dce3SE9   10    1 single-time-Enh
41                                earlySE0   10    1 single-time-Enh
42                                earlySE1   10    1 single-time-Enh
43                               earlySE10   10    1 single-time-Enh
44                               earlySE11   10    1 single-time-Enh
45                               earlySE12   10    1 single-time-Enh
46                               earlySE15   10    1 single-time-Enh
47                               earlySE16   10    1 single-time-Enh
48                               earlySE17   10    1 single-time-Enh
49                               earlySE18   10    1 single-time-Enh
50                               earlySE19   10    1 single-time-Enh
51                                earlySE2   10    1 single-time-Enh
52                                earlySE3   10    1 single-time-Enh
53                                earlySE5   10    1 single-time-Enh
54                                earlySE6   10    1 single-time-Enh
55                                earlySE7   10    1 single-time-Enh
56                                earlySE8   10    1 single-time-Enh
57                                earlySE9   10    1 single-time-Enh
58                         edge_sharp_mean   10    1   T1wmorphology
59                          edge_sharp_std   10    1   T1wmorphology
60                           iAUC1_countor   10    1      T1wdynamic
61                            iAUC1_inside   10    1      T1wdynamic
62                  iiiMax_Margin_Gradient   10    1   T1wmorphology
63            iiMin_change_Variance_uptake   10    1   T1wmorphology
64                    iMax_Variance_uptake   10    1   T1wmorphology
65                            irregularity   10    1   T1wmorphology
66                              ivVariance   10    1   T1wmorphology
67                       k_Max_Margin_Grad   10    1   T1wmorphology
68                           Kpeak_countor   10    1      T1wdynamic
69                            Kpeak_inside   10    1      T1wdynamic
70                              kurt_F_r_i   10    1      T1wdynamic
71                                 lateSE0   10    1 single-time-Enh
72                                 lateSE1   10    1 single-time-Enh
73                                lateSE10   10    1 single-time-Enh
74                                lateSE11   10    1 single-time-Enh
75                                lateSE12   10    1 single-time-Enh
76                                lateSE13   10    1 single-time-Enh
77                                lateSE14   10    1 single-time-Enh
79                                lateSE16   10    1 single-time-Enh
80                                lateSE17   10    1 single-time-Enh
81                                lateSE18   10    1 single-time-Enh
83                                 lateSE2   10    1 single-time-Enh
84                                 lateSE3   10    1 single-time-Enh
85                                 lateSE4   10    1 single-time-Enh
86                                 lateSE7   10    1 single-time-Enh
87                                 lateSE8   10    1 single-time-Enh
88                                 lateSE9   10    1 single-time-Enh
89                         LMSIR_predicted   10    1             T2w
90                               max_F_r_i   10    1      T1wdynamic
91                            max_RGH_mean   10    1   T1wmorphology
92                          max_RGH_mean_k   10    1   T1wmorphology
93                             max_RGH_var   10    1   T1wmorphology
94                           max_RGH_var_k   10    1   T1wmorphology
95                           maxCr_countor   10    1      T1wdynamic
96                            maxCr_inside   10    1      T1wdynamic
97                           maxVr_countor   10    1      T1wdynamic
98                            maxVr_inside   10    1      T1wdynamic
99                              mean_F_r_i   10    1      T1wdynamic
100                              min_F_r_i   10    1      T1wdynamic
101                         peakVr_countor   10    1      T1wdynamic
102                          peakVr_inside   10    1      T1wdynamic
103                            SER_countor   10    1      T1wdynamic
104                             SER_inside   10    1      T1wdynamic
105                             skew_F_r_i   10    1      T1wdynamic
106                      Slope_ini_countor   10    1      T1wdynamic
107                       Slope_ini_inside   10    1      T1wdynamic
108                          T2grad_margin   10    1             T2w
109                      T2grad_margin_var   10    1             T2w
110                           T2kurt_F_r_i   10    1             T2w
111                            T2max_F_r_i   10    1             T2w
112                           T2mean_F_r_i   10    1             T2w
113                            T2min_F_r_i   10    1             T2w
114                             T2RGH_mean   10    1             T2w
115                              T2RGH_var   10    1             T2w
116                           T2skew_F_r_i   10    1             T2w
117              T2texture_contrast_nondir   10    1             T2w
118           T2texture_correlation_nondir   10    1             T2w
119           T2texture_diffentropy_nondir   10    1             T2w
120          T2texture_diffvariance_nondir   10    1             T2w
121                T2texture_energy_nondir   10    1             T2w
122               T2texture_entropy_nondir   10    1             T2w
123     T2texture_inversediffmoment_nondir   10    1             T2w
124            T2texture_sumaverage_nondir   10    1             T2w
125            T2texture_sumentropy_nondir   10    1             T2w
126           T2texture_sumvariance_nondir   10    1             T2w
127              T2texture_variance_nondir   10    1             T2w
128                            T2var_F_r_i   10    1             T2w
129          texture_contrast_nondir_post1   10    1      T1wtexture
130          texture_contrast_nondir_post2   10    1      T1wtexture
131          texture_contrast_nondir_post3   10    1      T1wtexture
132          texture_contrast_nondir_post4   10    1      T1wtexture
133       texture_correlation_nondir_post1   10    1      T1wtexture
134       texture_correlation_nondir_post2   10    1      T1wtexture
135       texture_correlation_nondir_post3   10    1      T1wtexture
136       texture_correlation_nondir_post4   10    1      T1wtexture
137       texture_diffentropy_nondir_post1   10    1      T1wtexture
138       texture_diffentropy_nondir_post2   10    1      T1wtexture
139       texture_diffentropy_nondir_post3   10    1      T1wtexture
141      texture_diffvariance_nondir_post1   10    1      T1wtexture
142      texture_diffvariance_nondir_post2   10    1      T1wtexture
143      texture_diffvariance_nondir_post3   10    1      T1wtexture
144      texture_diffvariance_nondir_post4   10    1      T1wtexture
145            texture_energy_nondir_post1   10    1      T1wtexture
146            texture_energy_nondir_post2   10    1      T1wtexture
147            texture_energy_nondir_post3   10    1      T1wtexture
148            texture_energy_nondir_post4   10    1      T1wtexture
149           texture_entropy_nondir_post1   10    1      T1wtexture
150           texture_entropy_nondir_post2   10    1      T1wtexture
151           texture_entropy_nondir_post3   10    1      T1wtexture
152           texture_entropy_nondir_post4   10    1      T1wtexture
153 texture_inversediffmoment_nondir_post1   10    1      T1wtexture
154 texture_inversediffmoment_nondir_post2   10    1      T1wtexture
155 texture_inversediffmoment_nondir_post3   10    1      T1wtexture
156 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
157        texture_sumaverage_nondir_post1   10    1      T1wtexture
158        texture_sumaverage_nondir_post2   10    1      T1wtexture
159        texture_sumaverage_nondir_post3   10    1      T1wtexture
160        texture_sumaverage_nondir_post4   10    1      T1wtexture
161        texture_sumentropy_nondir_post1   10    1      T1wtexture
162        texture_sumentropy_nondir_post2   10    1      T1wtexture
163        texture_sumentropy_nondir_post3   10    1      T1wtexture
164        texture_sumentropy_nondir_post4   10    1      T1wtexture
165       texture_sumvariance_nondir_post1   10    1      T1wtexture
166       texture_sumvariance_nondir_post2   10    1      T1wtexture
167       texture_sumvariance_nondir_post3   10    1      T1wtexture
168       texture_sumvariance_nondir_post4   10    1      T1wtexture
169          texture_variance_nondir_post1   10    1      T1wtexture
170          texture_variance_nondir_post2   10    1      T1wtexture
171          texture_variance_nondir_post3   10    1      T1wtexture
172          texture_variance_nondir_post4   10    1      T1wtexture
173                          Tpeak_countor   10    1      T1wdynamic
174                           Tpeak_inside   10    1      T1wdynamic
175                     UptakeRate_countor   10    1      T1wdynamic
176                      UptakeRate_inside   10    1      T1wdynamic
177                                     V0   10    1      dispersion
178                                     V1   10    1      dispersion
179                                    V10   10    1      dispersion
180                                    V11   10    1      dispersion
181                                    V12   10    1      dispersion
182                                    V13   10    1      dispersion
183                                    V14   10    1      dispersion
184                                    V15   10    1      dispersion
185                                    V16   10    1      dispersion
186                                    V17   10    1      dispersion
187                                    V18   10    1      dispersion
188                                    V19   10    1      dispersion
189                                     V2   10    1      dispersion
190                                     V3   10    1      dispersion
191                                     V4   10    1      dispersion
192                                     V5   10    1      dispersion
193                                     V6   10    1      dispersion
194                                     V7   10    1      dispersion
195                                     V8   10    1      dispersion
196                                     V9   10    1      dispersion
197                              var_F_r_i   10    1      T1wdynamic
198              Vr_decreasingRate_countor   10    1      T1wdynamic
199               Vr_decreasingRate_inside   10    1      T1wdynamic
200              Vr_increasingRate_countor   10    1      T1wdynamic
201               Vr_increasingRate_inside   10    1      T1wdynamic
202                      Vr_post_1_countor   10    1      T1wdynamic
203                       Vr_post_1_inside   10    1      T1wdynamic
204                    washoutRate_countor   10    1      T1wdynamic
10                                dce2SE10    8    1 single-time-Enh
18                                dce2SE19    8    1 single-time-Enh
34                                dce3SE18    8    1 single-time-Enh
78                                lateSE15    8    1 single-time-Enh
82                                lateSE19    8    1 single-time-Enh
140       texture_diffentropy_nondir_post4    8    1      T1wtexture
205                     washoutRate_inside    8    1      T1wdynamic
206                               dce2SE15    8    1 single-time-Enh
207                                dce2SE3    8    1 single-time-Enh
208                                dce2SE6    8    1 single-time-Enh
209                               dce3SE16    8    1 single-time-Enh
210                               dce3SE19    8    1 single-time-Enh
211                                dce3SE2    8    1 single-time-Enh
212                                dce3SE4    8    1 single-time-Enh
213                              earlySE13    8    1 single-time-Enh
214                              earlySE14    8    1 single-time-Enh
215                               earlySE4    8    1 single-time-Enh
216                                lateSE5    8    1 single-time-Enh
217                                lateSE6    8    1 single-time-Enh
218                          peakCr_inside    8    1      T1wdynamic
219                         peakCr_countor    4    0      T1wdynamic
```

```r
#plot
ggplot(dfreqimgT1T2[dfreqimgT1T2$high==1,], aes(x=reorder(Var1, order(Freq, decreasing = TRUE)), y=Freq, fill=factor(typefeature))) + 
  geom_bar(stat = "identity") + coord_flip() +
  ggtitle("T1w+T2wSI_noLMSIR ensembles: feature selection frequency") +
  labs(y="featsel frequency", x=" ") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=18, hjust=0))+
  theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16))
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/featuressel-3.png)<!-- -->


# now based on frequency of selecting a feature

```r
# pick frequency of 75% or higher as very common feature
dfimgT1 = data.frame(table(varImpall_imgT1w_total$selfeat))
dfimgT1$high = (dfimgT1$Freq>=0.75*10)*3 + 
               (dfimgT1$Freq>=0.5*10 & dfimgT1$Freq<0.75*10)*2 + 
               (dfimgT1$Freq<0.5*10)*1
dfimgT1 = dfimgT1[order(dfimgT1$Freq, decreasing = TRUE),]

# with dictionary assign type of feature
dictimgT1 = feature_dictionary(imgT1train)
dfreqimgT1 = data.frame()
freqf = c()
for(di in 1:nrow(dfimgT1)){
  featname = as.character(dfimgT1[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    dfreqimgT1 = rbind(dfreqimgT1, cbind(dfimgT1[di,],
                       typefeature=dictimgT1$fnnames[dictimgT1$fnnames$f==featname,"type"][1]))}}

T1w75orh =  dfreqimgT1[dfreqimgT1$high==3,]
T1w50to75 =  dfreqimgT1[dfreqimgT1$high==2,]
T1w25to50 =  dfreqimgT1[dfreqimgT1$high==1,]

# plot 2
# select freuqutnly selected features
maxFreq_dictT1w = summary(as.factor(dictimgT1$fnnames$type))
pander(maxFreq_dictT1w)
```


------------------------------------------------------------------------
 dispersion   single-time-Enh   T1wdynamic   T1wmorphology   T1wtexture 
------------ ----------------- ------------ --------------- ------------
     20             80              40            14             44     
------------------------------------------------------------------------

```r
T1w75orhF = data.frame(featureGroup = names(summary(T1w75orh$typefeature)),
           freq = summary(T1w75orh$typefeature))
rownames(T1w75orhF)=NULL
T1w75orhF$flagSel = ">75%"
T1w50to75F = data.frame(featureGroup = names(summary(T1w50to75$typefeature)),
           freq = summary(T1w50to75$typefeature))
rownames(T1w50to75F)=NULL
T1w50to75F$flagSel = "50-75%"
T1w25to50F = data.frame(featureGroup = names(summary(T1w25to50$typefeature)),
           freq = summary(T1w25to50$typefeature))
rownames(T1w25to50F)=NULL
T1w25to50F$flagSel = "<50%"

allT1wFfreq = rbind(T1w75orhF, T1w50to75F, T1w25to50F)

maxFf= c()
for(i in 1:nrow(allT1wFfreq)){
  Ff = unname(maxFreq_dictT1w[names(maxFreq_dictT1w) == allT1wFfreq$featureGroup[i]])
  maxFf = c(maxFf, Ff)
}
allT1wFfreq$classifier = "T1wonly"

#plot +scale_fill_manual(values=c(p1,p2))
#p1<- barplot(1, angle = 45, density = 10, col = "black")
#p2<- barplot(1, angle = 145, density = 10, col = "grey")

# Get the positions of the labels
allT1wFfreq = ddply(allT1wFfreq, .(featureGroup), transform, pos = cumsum(freq) - 0.5*freq) 

labelsT1 = allT1wFfreq$freq
labelsT1[labelsT1==0]=""
allT1wFfreq$labels = labelsT1

ggplot(allT1wFfreq, aes(x=featureGroup, y=freq, fill=flagSel)) + 
  geom_bar(stat = "identity") + coord_flip()   +
  geom_text(aes(label = labels, y = pos), size = 3) +
  scale_fill_discrete(guide = guide_legend(title = "cv Selection\nfrequency")) +
  ggtitle("only T1w featsel") +
  labs(y="# features selected", x=" ") 
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/frequency-graphs-1.png)<!-- -->

```r
###########################################
# pick frequency of 75% or higher as very common feature
dfimgT1T2 = data.frame(table(varImpall_T2wpLMSIR_total$selfeat))
dfimgT1T2$high = (dfimgT1T2$Freq>=0.75*10)*3 + 
               (dfimgT1T2$Freq>=0.5*10 & dfimgT1T2$Freq<0.75*10)*2 + 
               (dfimgT1T2$Freq<0.5*10)*1
dfimgT1T2 = dfimgT1T2[order(dfimgT1T2$Freq, decreasing = TRUE),]

# with dictionary assign type of feature
dictimgT1T2 = feature_dictionary(imgT2pLMSIR)
dfreqimgT1T2 = data.frame()
freqf = c()
for(di in 1:nrow(dfimgT1T2)){
  featname = as.character(dfimgT1T2[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    dfreqimgT1T2 = rbind(dfreqimgT1T2, cbind(dfimgT1T2[di,],
                       typefeature=dictimgT1T2$fnnames[dictimgT1T2$fnnames$f==featname,"type"][1]))}
  }

T1T2w75orh =  dfreqimgT1T2[dfreqimgT1T2$high==3,]
T1T2w50to75 =  dfreqimgT1T2[dfreqimgT1T2$high==2,]
T1T2w25to50 =  dfreqimgT1T2[dfreqimgT1T2$high==1,]

# plot 2
# select freuqutnly selected features
maxFreq_dictT1T2 = summary(as.factor(dictimgT1T2$fnnames$type))
pander(maxFreq_dictT1T2)
```


------------------------------------------------------------------------------
 dispersion   single-time-Enh   T1wdynamic   T1wmorphology   T1wtexture   T2w 
------------ ----------------- ------------ --------------- ------------ -----
     20             80              40            14             44       22  
------------------------------------------------------------------------------

```r
T1T2w75orhF = data.frame(featureGroup = names(summary(T1T2w75orh$typefeature)),
           freq = summary(T1T2w75orh$typefeature))
rownames(T1T2w75orhF)=NULL
T1T2w75orhF$flagSel = ">75%"
T1T2w50to75F = data.frame(featureGroup = names(summary(T1T2w50to75$typefeature)),
           freq = summary(T1T2w50to75$typefeature))
rownames(T1T2w50to75F)=NULL
T1T2w50to75F$flagSel = "50-75%"
T1T2w25to50F = data.frame(featureGroup = names(summary(T1T2w25to50$typefeature)),
           freq = summary(T1T2w25to50$typefeature))
rownames(T1T2w25to50F)=NULL
T1T2w25to50F$flagSel = "<50%"

allT1T2wFfreq = rbind(T1T2w75orhF, T1T2w50to75F, T1T2w25to50F)


maxFf= c()
for(i in 1:nrow(allT1T2wFfreq)){
  Ff = unname(maxFreq_dictT1T2[names(maxFreq_dictT1T2) == allT1T2wFfreq$featureGroup[i]])
  maxFf = c(maxFf, Ff)
}
allT1T2wFfreq$classifier = "T1w+T2w"

# finally plot
# Get the positions of the labels
allT1T2wFfreq = ddply(allT1T2wFfreq, .(featureGroup), transform, pos = cumsum(freq) - 0.5*freq)      

labels = allT1T2wFfreq$freq
labels[labels==0]=""
allT1T2wFfreq$labels = labels

ggplot(allT1T2wFfreq, aes(x=featureGroup, y=freq, fill=flagSel)) + 
  geom_bar(stat = "identity") + coord_flip()   +
  geom_text(aes(label = labels, y = pos), size = 3) +
  scale_fill_discrete(guide = guide_legend(title = "cv Selection\nfrequency")) +
  ggtitle("T1w +T2w_predicted LMSIR featsel") +
  labs(y="# features selected", x=" ") 
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/frequency-graphs-2.png)<!-- -->

```r
# plot combingin
allT1T2highF = rbind(allT1wFfreq, allT1T2wFfreq)
# plot
labelsT1T2 = allT1T2highF$freq
labelsT1T2[labelsT1T2==0]=""
allT1T2highF$labels = labelsT1T2

ggplot(allT1T2highF, aes(x=featureGroup, y=freq, fill=flagSel)) + 
  geom_bar(stat = "identity") + coord_flip()  +
  geom_text(aes(label = labels, y = pos), size = 3) +
  scale_fill_discrete(guide = guide_legend(title = "cv Selection\nfrequency")) +
  labs(y="# features selected", x=" ") + facet_grid(. ~ classifier)
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/frequency-graphs-3.png)<!-- -->

```r
#plot +scale_fill_manual(values=c(p1,p2))
#p1<- barplot(1, angle = 45, density = 10, col = "black")
#p2<- barplot(1, angle = 145, density = 10, col = "grey")
```

## Within relevant T2w features

```r
# pick frequency of 75% or higher as very common feature
dfimgT2 = data.frame(table(varImpall_T2wpLMSIR_total$selfeat))
dfimgT2$high = (dfimgT2$Freq>=0.75*10)*3 + 
               (dfimgT2$Freq>=0.5*10 & dfimgT2$Freq<0.75*10)*2 + 
               (dfimgT2$Freq<0.5*10)*1
dfimgT2 = dfimgT2[order(dfimgT2$Freq, decreasing = TRUE),]

# with dictionary assign type of feature
dictimgT2 = feature_dictionary_typesT2(imgT2pLMSIR)
dfreqimgT2 = data.frame()
freqf = c()
for(di in 1:nrow(dfimgT2)){
  featname = as.character(dfimgT2[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    if(any(dictimgT2$fnnames$f==featname)){
        dfreqimgT2 = rbind(dfreqimgT2, cbind(dfimgT2[di,],
                          typefeature=dictimgT2$fnnames[dictimgT2$fnnames$f==featname,"type"][1]))}}}

# only for selected T2w features
T2w75orh =  dfreqimgT2[dfreqimgT2$high==3,]
T2w50to75 =  dfreqimgT2[dfreqimgT2$high==2,]
T2w25to50 =  dfreqimgT2[dfreqimgT2$high==1,]

# plot 2
# select freuqutnly selected features
maxFreq_dictT2 = summary(as.factor(dictimgT2$fnnames$type))
pander(maxFreq_dictT2)
```


---------------------------------------------------------------
 predictedT2wLMSIR   T2wIntensity   T2wmorphology   T2wtexture 
------------------- -------------- --------------- ------------
         1                6               4             11     
---------------------------------------------------------------

```r
T2w75orhF = data.frame(featureGroup = names(summary(T2w75orh$typefeature)),
           freq = summary(T2w75orh$typefeature))
rownames(T2w75orhF)=NULL
T2w75orhF$flagSel = ">75%"
T2w50to75F = data.frame(featureGroup = names(summary(T2w50to75$typefeature)),
           freq = summary(T2w50to75$typefeature))
rownames(T2w50to75F)=NULL
T2w50to75F$flagSel = "50-75%"
T2w25to50F = data.frame(featureGroup = names(summary(T2w25to50$typefeature)),
           freq = summary(T2w25to50$typefeature))
rownames(T2w25to50F)=NULL
T2w25to50F$flagSel = "<50%"

allT2wFfreqpredicted = rbind(T2w75orhF, T2w50to75F, T2w25to50F)


maxFf= c()
for(i in 1:nrow(allT2wFfreqpredicted)){
  Ff = unname(maxFreq_dictT2[names(maxFreq_dictT2) == allT2wFfreqpredicted$featureGroup[i]])
  maxFf = c(maxFf, Ff)
}
allT2wFfreqpredicted$classifier = "selectedT2w"

# finally plot
# Get the positions of the labels
allT2wFfreqpredicted = ddply(allT2wFfreqpredicted, .(featureGroup), transform, pos = cumsum(freq) - 0.5*freq)

labelsT2 = allT2wFfreqpredicted$freq
labelsT2[labelsT2==0]=""
allT2wFfreqpredicted$labels = labelsT2

ggplot(allT2wFfreqpredicted, aes(x=featureGroup, y=freq, fill=flagSel)) + 
  geom_bar(stat = "identity") + coord_flip()  +
  geom_text(aes(label = labelsT2, y = pos), size = 3) +
  scale_fill_discrete(guide = guide_legend(title = "cv Selection\nfrequency")) +
  ggtitle("T1w+T2w_predicted LMSIR ensembles: T2w feature selection frequency") +
  labs(y="# features selected", x=" ") 
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/T2featsel-1.png)<!-- -->

```r
######################################
# pick frequency of 75% or higher as very common feature
dfimgT2 = data.frame(table(varImpall_T1wT2w_total$selfeat))
dfimgT2$high = (dfimgT2$Freq>=0.75*10)*3 + 
               (dfimgT2$Freq>=0.5*10 & dfimgT2$Freq<0.75*10)*2 + 
               (dfimgT2$Freq<0.5*10)*1
dfimgT2 = dfimgT2[order(dfimgT2$Freq, decreasing = TRUE),]

# with dictionary assign type of feature
dictimgT2 = feature_dictionary_typesT2(imgT1T2)
dfreqimgT2 = data.frame()
freqf = c()
for(di in 1:nrow(dfimgT2)){
  featname = as.character(dfimgT2[di,"Var1"])
  if(!(featname %in% freqf)){
    freqf = c(freqf, featname)
    if(any(dictimgT2$fnnames$f==featname)){
        dfreqimgT2 = rbind(dfreqimgT2, cbind(dfimgT2[di,],
                          typefeature=dictimgT2$fnnames[dictimgT2$fnnames$f==featname,"type"][1]))}}}

# only for selected T2w features
T2w75orh =  dfreqimgT2[dfreqimgT2$high==3,]
T2w50to75 =  dfreqimgT2[dfreqimgT2$high==2,]
T2w25to50 =  dfreqimgT2[dfreqimgT2$high==1,]

# plot 2
# select freuqutnly selected features
maxFreq_dictT2 = summary(as.factor(dictimgT2$fnnames$type))
pander(maxFreq_dictT2)
```


-------------------------------------------
 T2wIntensity   T2wmorphology   T2wtexture 
-------------- --------------- ------------
      28              4             11     
-------------------------------------------

```r
T2w75orhF = data.frame(featureGroup = names(summary(T2w75orh$typefeature)),
           freq = summary(T2w75orh$typefeature))
rownames(T2w75orhF)=NULL
T2w75orhF$flagSel = ">75%"
T2w50to75F = data.frame(featureGroup = names(summary(T2w50to75$typefeature)),
           freq = summary(T2w50to75$typefeature))
rownames(T2w50to75F)=NULL
T2w50to75F$flagSel = "50-75%"
T2w25to50F = data.frame(featureGroup = names(summary(T2w25to50$typefeature)),
           freq = summary(T2w25to50$typefeature))
rownames(T2w25to50F)=NULL
T2w25to50F$flagSel = "<50%"

allT2wFfreqmeasured = rbind(T2w75orhF, T2w50to75F, T2w25to50F)


maxFf= c()
for(i in 1:nrow(allT2wFfreqmeasured)){
  Ff = unname(maxFreq_dictT2[names(maxFreq_dictT2) == allT2wFfreqmeasured$featureGroup[i]])
  maxFf = c(maxFf, Ff)
}
allT2wFfreqmeasured$classifier = "selectedT2w"

# finally plot
# Get the positions of the labels
allT2wFfreqmeasured = ddply(allT2wFfreqmeasured, .(featureGroup), transform, pos = cumsum(freq) - 0.5*freq)

labelsT2 = allT2wFfreqmeasured$freq
labelsT2[labelsT2==0]=""
allT2wFfreqmeasured$labels = labelsT2

ggplot(allT2wFfreqmeasured, aes(x=featureGroup, y=freq, fill=flagSel)) + 
  geom_bar(stat = "identity") + coord_flip()  +
  geom_text(aes(label = labelsT2, y = pos), size = 3) +
  scale_fill_discrete(guide = guide_legend(title = "cv Selection\nfrequency")) +
  ggtitle("T1w+T2w_measure LMSIR ensembles: T2w feature selection frequency") +
  labs(y="# features selected", x=" ") 
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_extratrees_files/figure-html/T2featsel-2.png)<!-- -->
               


```r
save.image("Outputs/T2SIvspredLMSIR_extratrees_summaryResults.RData")
```

