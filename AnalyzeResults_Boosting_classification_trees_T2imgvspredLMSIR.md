# Analysis of results - T2SI vs predLMSIR
Cristina Gallego  
August 25, 2017  



# Analysis of results by each fold

```r
allcvauc_imgT1_total= c()
allcvauc_T1wT2w_total = c()
allcvauc_T2wpLMSIR_total= c()

## using 3Dtexture first + Boosting  
perfall_imgT1 = data.frame() 
perfall_imgT1T2 = data.frame() 
perfall_T2wpLMSIR = data.frame() 

varImpall_imgT1w_total = data.frame() 
varImpall_T1wT2w_total = data.frame() 
varImpall_T2wpLMSIR_total = data.frame() 

# perform k-fold-out
for(k in 1:10){  # 1:10f cv
  load(paste0("Outputs/T1T2imgvsT2wtextpredLMSIR_boost_addeddiagvalue_cv",k,".RData"))

  # append test results
  perfall_imgT1 = rbind(perfall_imgT1, perf_imgT1) 
  perfall_imgT1T2 = rbind(perfall_imgT1T2, perf_imgT1T2) 
  perfall_T2wpLMSIR = rbind(perfall_T2wpLMSIR, perf_T2wpLMSIR) 
  
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-1.png)<!-- -->

```
Area under the curve: 0.7463
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5601932 0.6571       0.8  0.9429   0.44      0.64     0.8
```

```
Area under the curve: 0.8177
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5129804   0.65       0.8   0.925 0.5806    0.7419  0.9032
```

```
Area under the curve: 0.8306
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4915374  0.925     0.975       1 0.4194    0.5806  0.7419
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-2.png)<!-- -->

```
Area under the curve: 0.8468
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6506404  0.475     0.625   0.775 0.8387    0.9355       1
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
  0.5121938 0.7222    0.8611  0.9722 0.5263    0.7368  0.9474
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-3.png)<!-- -->

```
Area under the curve: 0.8187
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4627353      1         1       1 0.3158    0.5263  0.7368
```

```
Area under the curve: 0.8153
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5306807    0.6    0.7667     0.9 0.5417      0.75  0.9167
```

```
Area under the curve: 0.8389
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5257692 0.7333    0.8667  0.9667 0.4583    0.6667  0.8333
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-4.png)<!-- -->

```
Area under the curve: 0.8236
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
   0.488011    0.9    0.9667       1 0.4167     0.625  0.7917
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
   0.564311 0.8438    0.9375       1 0.5185    0.7037  0.8889
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-5.png)<!-- -->

```
Area under the curve: 0.8322
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6376378 0.4688    0.6562  0.8125 0.7407    0.8889       1
```

```
Area under the curve: 0.8298
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5750095 0.5745    0.7021  0.8298   0.75       0.9       1
```

```
Area under the curve: 0.85
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5757877 0.5745    0.7021  0.8298   0.75       0.9       1
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-6.png)<!-- -->

```
Area under the curve: 0.8553
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5155473 0.7447    0.8511  0.9362   0.55      0.75     0.9
```

```
Area under the curve: 0.7714
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5761598 0.5641    0.7179  0.8462 0.6667    0.8333  0.9583
```

```
Area under the curve: 0.7628
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5655913 0.5128    0.6667  0.7949 0.5833      0.75  0.9167
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-7.png)<!-- -->

```
Area under the curve: 0.8216
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5459772 0.7949    0.8974  0.9744 0.5417      0.75  0.9167
```

```
Area under the curve: 0.8114
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5557392 0.5952    0.7381  0.8571    0.6      0.76    0.92
```

```
Area under the curve: 0.7895
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5846236 0.7619    0.8571  0.9524   0.44      0.64    0.84
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-8.png)<!-- -->

```
Area under the curve: 0.8381
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5333076 0.7143    0.8333  0.9286   0.64       0.8    0.96
```

```
Area under the curve: 0.775
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5456363 0.7292    0.8333  0.9375   0.45      0.65    0.85
```

```
Area under the curve: 0.7312
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.6148672 0.5417    0.6667  0.7917   0.55      0.75    0.95
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-9.png)<!-- -->

```
Area under the curve: 0.8021
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5295403   0.75    0.8542  0.9375   0.45      0.65    0.85
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/process-cvks-10.png)<!-- -->

```
Area under the curve: 0.9485
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5795522 0.9091    0.9697       1    0.7    0.8333  0.9667
```

# Combine, pool data, and plot final results

```r
print(summary(allcvauc_imgT1))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7383  0.7383  0.7383  0.7383  0.7383  0.7383 
```

```r
print(summary(allcvauc_T1wT2w))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7566  0.7566  0.7566  0.7566  0.7566  0.7566 
```

```r
print(summary(allcvauc_T2wpLMSIR))
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.768   0.768   0.768   0.768   0.768   0.768 
```

```r
# boxplots of cv-performances
cvperfs = data.frame()
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_imgT1, group="T1w"))
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_T1wT2w, group="T1w+T2wText+T2w_SI"))
cvperfs = rbind(cvperfs, data.frame(cvAUC=allcvauc_T2wpLMSIR, group="T1w+T2wText+predictiveLMSIR"))
# find min
minAUC = min(cvperfs$cvAUC)

# plot
p <- ggplot(cvperfs, aes(factor(group), cvAUC))
p + geom_boxplot(aes(fill = factor(group)))
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/Summarize-1.png)<!-- -->

```r
########### 
# plot pooled all cv-heldout test cases
par(mfrow=c(1,1))
n=15
colors = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
# plot 1/4
p1 = calcAUC_plot(perfall_imgT1$obs, perfall_imgT1$C, 
                           xptext=0.45, yptext=0.75 , 1, colors[2], atitle="")
```

```
Area under the curve: 0.8118
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5808038 0.6618    0.7647  0.8529 0.7091    0.8182  0.9091
```

```r
par(new=TRUE)
p2 = calcAUC_plot(perfall_imgT1T2$obs, perfall_imgT1T2$C, 
                           xptext=0.55, yptext=0.65, 2, colors[9], atitle="")
```

```
Area under the curve: 0.8366
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5521249 0.7206    0.8235  0.9118 0.6727    0.7818  0.8909
```

```r
par(new=TRUE)
p3 = calcAUC_plot(perfall_T2wpLMSIR$obs, perfall_T2wpLMSIR$C,
                           xptext=0.65, yptext=0.55, 3, colors[11], 
                  atitle="ROCs pooled heldout-patient across ALL folds")
```

```
Area under the curve: 0.8631
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.5257094 0.7941    0.8824  0.9559 0.6541    0.7636  0.8727
```

```r
legend("bottomright", 
       legend = c(paste0("T1w"),
                    paste0("T1w + T2wText + T2w_SI"),
                    paste0("T1w + T2wText + predictiveLMSIR")),
       col = c(colors[2],colors[9],colors[11]), lty=c(1,2,3), lwd = 2)
```

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/Summarize-2.png)<!-- -->

```r
# find significants: only imgT1 vs. allT2
p1p2 = roc.test(p1$ROC, p2$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p1p2)
```

```

	Bootstrap test for two correlated ROC curves

data:  p1$ROC and p2$ROC
D = -1.7072, boot.n = 2000, boot.stratified = 1, p-value = 0.08778
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8117647   0.8366310 
```

```r
# find significants: imgT1 vs wLMSIR
p1p3 = roc.test(p1$ROC, p3$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p1p3)
```

```

	Bootstrap test for two correlated ROC curves

data:  p1$ROC and p3$ROC
D = -2.7627, boot.n = 2000, boot.stratified = 1, p-value = 0.005732
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8117647   0.8631016 
```

```r
# find significants: imgT1 vs wLMSIR
p2p3 = roc.test(p2$ROC, p3$ROC, method=c("bootstrap"), alternative = c("two.sided"), boot.stratified=TRUE)
print(p2p3)
```

```

	Bootstrap test for two correlated ROC curves

data:  p2$ROC and p3$ROC
D = -1.6445, boot.n = 2000, boot.stratified = 1, p-value = 0.1001
alternative hypothesis: true difference in AUC is not equal to 0
sample estimates:
AUC of roc1 AUC of roc2 
  0.8366310   0.8631016 
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
0.26335435 0.01719611 0.30020780 
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
7                              circularity   10    1   T1wmorphology
37                                dce3SE17   10    1 single-time-Enh
52                               earlySE12   10    1 single-time-Enh
55                               earlySE15   10    1 single-time-Enh
66                                earlySE8   10    1 single-time-Enh
74                    iMax_Variance_uptake   10    1   T1wmorphology
75                            irregularity   10    1   T1wmorphology
90                                lateSE17   10    1 single-time-Enh
97                                 lateSE6   10    1 single-time-Enh
100                                lateSE9   10    1 single-time-Enh
101                              max_F_r_i   10    1      T1wdynamic
116                             SER_inside   10    1      T1wdynamic
139            texture_energy_nondir_post4   10    1      T1wtexture
147 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
148        texture_sumaverage_nondir_post1   10    1      T1wtexture
149        texture_sumaverage_nondir_post2   10    1      T1wtexture
151        texture_sumaverage_nondir_post4   10    1      T1wtexture
156       texture_sumvariance_nondir_post1   10    1      T1wtexture
160          texture_variance_nondir_post1   10    1      T1wtexture
170                                    V10   10    1      dispersion
172                                    V12   10    1      dispersion
183                                     V5   10    1      dispersion
192               Vr_increasingRate_inside   10    1      T1wdynamic
1                                A_countor    9    1      T1wdynamic
2                                 A_inside    9    1      T1wdynamic
3                            alpha_countor    9    1      T1wdynamic
5                             beta_countor    9    1      T1wdynamic
8                                  dce2SE0    9    1 single-time-Enh
9                                  dce2SE1    9    1 single-time-Enh
13                                dce2SE13    9    1 single-time-Enh
23                                 dce2SE5    9    1 single-time-Enh
26                                 dce2SE8    9    1 single-time-Enh
28                                 dce3SE0    9    1 single-time-Enh
32                                dce3SE12    9    1 single-time-Enh
47                                 dce3SE9    9    1 single-time-Enh
53                               earlySE13    9    1 single-time-Enh
54                               earlySE14    9    1 single-time-Enh
56                               earlySE16    9    1 single-time-Enh
59                               earlySE19    9    1 single-time-Enh
61                                earlySE3    9    1 single-time-Enh
64                                earlySE6    9    1 single-time-Enh
68                         edge_sharp_mean    9    1   T1wmorphology
71                            iAUC1_inside    9    1      T1wdynamic
72                  iiiMax_Margin_Gradient    9    1   T1wmorphology
73            iiMin_change_Variance_uptake    9    1   T1wmorphology
79                            Kpeak_inside    9    1      T1wdynamic
80                              kurt_F_r_i    9    1      T1wdynamic
81                                 lateSE0    9    1 single-time-Enh
84                                lateSE11    9    1 single-time-Enh
85                                lateSE12    9    1 single-time-Enh
89                                lateSE16    9    1 single-time-Enh
93                                 lateSE2    9    1 single-time-Enh
94                                 lateSE3    9    1 single-time-Enh
98                                 lateSE7    9    1 single-time-Enh
99                                 lateSE8    9    1 single-time-Enh
102                           max_RGH_mean    9    1   T1wmorphology
104                            max_RGH_var    9    1   T1wmorphology
110                             mean_F_r_i    9    1      T1wdynamic
115                            SER_countor    9    1      T1wdynamic
117                             skew_F_r_i    9    1      T1wdynamic
127       texture_correlation_nondir_post4    9    1      T1wtexture
132      texture_diffvariance_nondir_post1    9    1      T1wtexture
133      texture_diffvariance_nondir_post2    9    1      T1wtexture
134      texture_diffvariance_nondir_post3    9    1      T1wtexture
135      texture_diffvariance_nondir_post4    9    1      T1wtexture
145 texture_inversediffmoment_nondir_post2    9    1      T1wtexture
146 texture_inversediffmoment_nondir_post3    9    1      T1wtexture
152        texture_sumentropy_nondir_post1    9    1      T1wtexture
161          texture_variance_nondir_post2    9    1      T1wtexture
164                          Tpeak_countor    9    1      T1wdynamic
168                                     V0    9    1      dispersion
169                                     V1    9    1      dispersion
171                                    V11    9    1      dispersion
174                                    V14    9    1      dispersion
175                                    V15    9    1      dispersion
176                                    V16    9    1      dispersion
177                                    V17    9    1      dispersion
179                                    V19    9    1      dispersion
180                                     V2    9    1      dispersion
182                                     V4    9    1      dispersion
184                                     V6    9    1      dispersion
186                                     V8    9    1      dispersion
187                                     V9    9    1      dispersion
188                              var_F_r_i    9    1      T1wdynamic
189              Vr_decreasingRate_countor    9    1      T1wdynamic
190               Vr_decreasingRate_inside    9    1      T1wdynamic
191              Vr_increasingRate_countor    9    1      T1wdynamic
193                      Vr_post_1_countor    9    1      T1wdynamic
194                       Vr_post_1_inside    9    1      T1wdynamic
196                     washoutRate_inside    9    1      T1wdynamic
4                             alpha_inside    8    1      T1wdynamic
6                              beta_inside    8    1      T1wdynamic
12                                dce2SE12    8    1 single-time-Enh
14                                dce2SE14    8    1 single-time-Enh
15                                dce2SE15    8    1 single-time-Enh
16                                dce2SE16    8    1 single-time-Enh
17                                dce2SE17    8    1 single-time-Enh
18                                dce2SE18    8    1 single-time-Enh
20                                 dce2SE2    8    1 single-time-Enh
21                                 dce2SE3    8    1 single-time-Enh
22                                 dce2SE4    8    1 single-time-Enh
30                                dce3SE10    8    1 single-time-Enh
34                                dce3SE14    8    1 single-time-Enh
35                                dce3SE15    8    1 single-time-Enh
40                                 dce3SE2    8    1 single-time-Enh
42                                 dce3SE4    8    1 single-time-Enh
43                                 dce3SE5    8    1 single-time-Enh
44                                 dce3SE6    8    1 single-time-Enh
45                                 dce3SE7    8    1 single-time-Enh
49                                earlySE1    8    1 single-time-Enh
50                               earlySE10    8    1 single-time-Enh
57                               earlySE17    8    1 single-time-Enh
60                                earlySE2    8    1 single-time-Enh
62                                earlySE4    8    1 single-time-Enh
65                                earlySE7    8    1 single-time-Enh
69                          edge_sharp_std    8    1   T1wmorphology
70                           iAUC1_countor    8    1      T1wdynamic
76                              ivVariance    8    1   T1wmorphology
78                           Kpeak_countor    8    1      T1wdynamic
82                                 lateSE1    8    1 single-time-Enh
83                                lateSE10    8    1 single-time-Enh
88                                lateSE15    8    1 single-time-Enh
91                                lateSE18    8    1 single-time-Enh
96                                 lateSE5    8    1 single-time-Enh
108                          maxVr_countor    8    1      T1wdynamic
109                           maxVr_inside    8    1      T1wdynamic
111                              min_F_r_i    8    1      T1wdynamic
120          texture_contrast_nondir_post1    8    1      T1wtexture
136            texture_energy_nondir_post1    8    1      T1wtexture
137            texture_energy_nondir_post2    8    1      T1wtexture
138            texture_energy_nondir_post3    8    1      T1wtexture
150        texture_sumaverage_nondir_post3    8    1      T1wtexture
155        texture_sumentropy_nondir_post4    8    1      T1wtexture
157       texture_sumvariance_nondir_post2    8    1      T1wtexture
159       texture_sumvariance_nondir_post4    8    1      T1wtexture
162          texture_variance_nondir_post3    8    1      T1wtexture
163          texture_variance_nondir_post4    8    1      T1wtexture
165                           Tpeak_inside    8    1      T1wdynamic
167                      UptakeRate_inside    8    1      T1wdynamic
173                                    V13    8    1      dispersion
178                                    V18    8    1      dispersion
181                                     V3    8    1      dispersion
195                    washoutRate_countor    8    1      T1wdynamic
10                                dce2SE10    7    0 single-time-Enh
11                                dce2SE11    7    0 single-time-Enh
19                                dce2SE19    7    0 single-time-Enh
24                                 dce2SE6    7    0 single-time-Enh
25                                 dce2SE7    7    0 single-time-Enh
29                                 dce3SE1    7    0 single-time-Enh
31                                dce3SE11    7    0 single-time-Enh
33                                dce3SE13    7    0 single-time-Enh
38                                dce3SE18    7    0 single-time-Enh
39                                dce3SE19    7    0 single-time-Enh
41                                 dce3SE3    7    0 single-time-Enh
46                                 dce3SE8    7    0 single-time-Enh
48                                earlySE0    7    0 single-time-Enh
51                               earlySE11    7    0 single-time-Enh
58                               earlySE18    7    0 single-time-Enh
63                                earlySE5    7    0 single-time-Enh
67                                earlySE9    7    0 single-time-Enh
95                                 lateSE4    7    0 single-time-Enh
103                         max_RGH_mean_k    7    0   T1wmorphology
106                          maxCr_countor    7    0      T1wdynamic
107                           maxCr_inside    7    0      T1wdynamic
118                      Slope_ini_countor    7    0      T1wdynamic
119                       Slope_ini_inside    7    0      T1wdynamic
121          texture_contrast_nondir_post2    7    0      T1wtexture
124       texture_correlation_nondir_post1    7    0      T1wtexture
125       texture_correlation_nondir_post2    7    0      T1wtexture
126       texture_correlation_nondir_post3    7    0      T1wtexture
131       texture_diffentropy_nondir_post4    7    0      T1wtexture
140           texture_entropy_nondir_post1    7    0      T1wtexture
141           texture_entropy_nondir_post2    7    0      T1wtexture
143           texture_entropy_nondir_post4    7    0      T1wtexture
153        texture_sumentropy_nondir_post2    7    0      T1wtexture
154        texture_sumentropy_nondir_post3    7    0      T1wtexture
158       texture_sumvariance_nondir_post3    7    0      T1wtexture
166                     UptakeRate_countor    7    0      T1wdynamic
185                                     V7    7    0      dispersion
77                       k_Max_Margin_Grad    6    0   T1wmorphology
86                                lateSE13    6    0 single-time-Enh
87                                lateSE14    6    0 single-time-Enh
105                          max_RGH_var_k    6    0   T1wmorphology
113                         peakVr_countor    6    0      T1wdynamic
122          texture_contrast_nondir_post3    6    0      T1wtexture
123          texture_contrast_nondir_post4    6    0      T1wtexture
128       texture_diffentropy_nondir_post1    6    0      T1wtexture
130       texture_diffentropy_nondir_post3    6    0      T1wtexture
144 texture_inversediffmoment_nondir_post1    6    0      T1wtexture
27                                 dce2SE9    5    0 single-time-Enh
36                                dce3SE16    5    0 single-time-Enh
92                                lateSE19    5    0 single-time-Enh
114                          peakVr_inside    5    0      T1wdynamic
129       texture_diffentropy_nondir_post2    5    0      T1wtexture
142           texture_entropy_nondir_post3    5    0      T1wtexture
112                          peakCr_inside    4    0      T1wdynamic
197                         peakCr_countor    3    0      T1wdynamic
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/featuressel-1.png)<!-- -->

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
1                                 A_inside   10    1
3                                  ave_T20   10    1
4                                  ave_T21   10    1
5                                 ave_T210   10    1
6                                 ave_T211   10    1
7                                 ave_T212   10    1
8                                 ave_T213   10    1
9                                 ave_T215   10    1
10                                ave_T216   10    1
12                                ave_T219   10    1
13                                 ave_T22   10    1
15                                 ave_T25   10    1
17                                 ave_T27   10    1
18                                 ave_T29   10    1
19                            beta_countor   10    1
20                             beta_inside   10    1
21                             circularity   10    1
22                                 dce2SE0   10    1
24                                dce2SE12   10    1
29                                dce2SE19   10    1
34                                 dce3SE0   10    1
41                                 dce3SE5   10    1
44                                earlySE0   10    1
45                                earlySE1   10    1
48                               earlySE14   10    1
49                               earlySE16   10    1
51                               earlySE19   10    1
53                                earlySE6   10    1
54                                earlySE7   10    1
55                                earlySE8   10    1
56                                earlySE9   10    1
57                         edge_sharp_mean   10    1
58                           iAUC1_countor   10    1
60                  iiiMax_Margin_Gradient   10    1
61            iiMin_change_Variance_uptake   10    1
62                    iMax_Variance_uptake   10    1
63                            irregularity   10    1
64                              ivVariance   10    1
65                            Kpeak_inside   10    1
66                              kurt_F_r_i   10    1
69                                lateSE11   10    1
70                                lateSE12   10    1
71                                lateSE17   10    1
73                                 lateSE3   10    1
75                                 lateSE8   10    1
76                                 lateSE9   10    1
77                               max_F_r_i   10    1
78                            max_RGH_mean   10    1
79                             max_RGH_var   10    1
82                              mean_F_r_i   10    1
84                             SER_countor   10    1
85                              SER_inside   10    1
86                              skew_F_r_i   10    1
91                       T2grad_margin_var   10    1
92                            T2kurt_F_r_i   10    1
93                             T2max_F_r_i   10    1
94                              T2RGH_mean   10    1
95                               T2RGH_var   10    1
97            T2texture_correlation_nondir   10    1
101               T2texture_entropy_nondir   10    1
109       texture_correlation_nondir_post1   10    1
111       texture_correlation_nondir_post4   10    1
114      texture_diffvariance_nondir_post1   10    1
115      texture_diffvariance_nondir_post2   10    1
118            texture_energy_nondir_post3   10    1
119            texture_energy_nondir_post4   10    1
122 texture_inversediffmoment_nondir_post2   10    1
123 texture_inversediffmoment_nondir_post3   10    1
124 texture_inversediffmoment_nondir_post4   10    1
125        texture_sumaverage_nondir_post1   10    1
126        texture_sumaverage_nondir_post2   10    1
127        texture_sumaverage_nondir_post4   10    1
128        texture_sumentropy_nondir_post1   10    1
131       texture_sumvariance_nondir_post1   10    1
133       texture_sumvariance_nondir_post4   10    1
134          texture_variance_nondir_post1   10    1
135          texture_variance_nondir_post2   10    1
136                          Tpeak_countor   10    1
137                           Tpeak_inside   10    1
141                                     V1   10    1
142                                    V10   10    1
143                                    V11   10    1
144                                    V12   10    1
145                                    V13   10    1
147                                    V15   10    1
148                                    V16   10    1
149                                    V17   10    1
150                                    V18   10    1
151                                    V19   10    1
154                                     V5   10    1
156                                     V7   10    1
157                                     V8   10    1
158                                     V9   10    1
159                              var_F_r_i   10    1
160              Vr_decreasingRate_countor   10    1
161               Vr_decreasingRate_inside   10    1
162              Vr_increasingRate_countor   10    1
163               Vr_increasingRate_inside   10    1
164                      Vr_post_1_countor   10    1
165                       Vr_post_1_inside   10    1
2                            alpha_countor    9    1
23                                 dce2SE1    9    1
25                                dce2SE13    9    1
26                                dce2SE14    9    1
27                                dce2SE16    9    1
28                                dce2SE18    9    1
30                                 dce2SE4    9    1
31                                 dce2SE5    9    1
32                                 dce2SE7    9    1
33                                 dce2SE8    9    1
35                                dce3SE11    9    1
36                                dce3SE13    9    1
38                                dce3SE19    9    1
39                                 dce3SE2    9    1
42                                 dce3SE7    9    1
43                                 dce3SE8    9    1
47                               earlySE12    9    1
52                                earlySE5    9    1
59                            iAUC1_inside    9    1
68                                lateSE10    9    1
72                                lateSE19    9    1
74                                 lateSE7    9    1
81                            maxVr_inside    9    1
87                       Slope_ini_countor    9    1
88                        Slope_ini_inside    9    1
89                             T2_lesionSI    9    1
90                           T2grad_margin    9    1
96               T2texture_contrast_nondir    9    1
98            T2texture_diffentropy_nondir    9    1
99           T2texture_diffvariance_nondir    9    1
100                T2texture_energy_nondir    9    1
102     T2texture_inversediffmoment_nondir    9    1
103            T2texture_sumaverage_nondir    9    1
106          texture_contrast_nondir_post1    9    1
108          texture_contrast_nondir_post3    9    1
110       texture_correlation_nondir_post2    9    1
113       texture_diffentropy_nondir_post4    9    1
116      texture_diffvariance_nondir_post3    9    1
120           texture_entropy_nondir_post1    9    1
129        texture_sumentropy_nondir_post2    9    1
139                      UptakeRate_inside    9    1
146                                    V14    9    1
152                                     V3    9    1
153                                     V4    9    1
155                                     V6    9    1
166                    washoutRate_countor    9    1
167                     washoutRate_inside    9    1
168                              A_countor    9    1
170                               ave_T214    9    1
171                               ave_T218    9    1
183                               dce3SE10    9    1
184                               dce3SE12    9    1
185                               dce3SE15    9    1
188                               dce3SE18    9    1
192                              earlySE10    9    1
194                              earlySE15    9    1
199                         edge_sharp_std    9    1
201                          Kpeak_countor    9    1
203                               lateSE13    9    1
209                                lateSE4    9    1
211                                lateSE6    9    1
234        texture_sumaverage_nondir_post3    9    1
236       texture_sumvariance_nondir_post3    9    1
237          texture_variance_nondir_post3    9    1
239                                     V2    9    1
16                                 ave_T26    8    1
37                                dce3SE14    8    1
40                                 dce3SE4    8    1
46                               earlySE11    8    1
50                               earlySE17    8    1
67                                 lateSE0    8    1
80                            maxCr_inside    8    1
104           T2texture_sumvariance_nondir    8    1
107          texture_contrast_nondir_post2    8    1
117            texture_energy_nondir_post1    8    1
121           texture_entropy_nondir_post3    8    1
130        texture_sumentropy_nondir_post4    8    1
132       texture_sumvariance_nondir_post2    8    1
138                     UptakeRate_countor    8    1
140                                     V0    8    1
169                           alpha_inside    8    1
173                                ave_T28    8    1
174                               dce2SE10    8    1
176                               dce2SE15    8    1
177                               dce2SE17    8    1
182                                dce3SE1    8    1
186                               dce3SE16    8    1
190                                dce3SE6    8    1
191                                dce3SE9    8    1
193                              earlySE13    8    1
195                              earlySE18    8    1
196                               earlySE2    8    1
198                               earlySE4    8    1
202                                lateSE1    8    1
205                               lateSE15    8    1
206                               lateSE16    8    1
207                               lateSE18    8    1
208                                lateSE2    8    1
210                                lateSE5    8    1
212                         max_RGH_mean_k    8    1
216                              min_F_r_i    8    1
217                         peakVr_countor    8    1
219                         T2_lesionSIstd    8    1
222                           T2skew_F_r_i    8    1
227       texture_diffentropy_nondir_post1    8    1
229      texture_diffvariance_nondir_post4    8    1
230            texture_energy_nondir_post2    8    1
231           texture_entropy_nondir_post2    8    1
232           texture_entropy_nondir_post4    8    1
233 texture_inversediffmoment_nondir_post1    8    1
238          texture_variance_nondir_post4    8    1
11                                ave_T217    7    0
14                                 ave_T24    7    0
105              T2texture_variance_nondir    7    0
112       texture_diffentropy_nondir_post2    7    0
172                                ave_T23    7    0
175                               dce2SE11    7    0
179                                dce2SE3    7    0
187                               dce3SE17    7    0
189                                dce3SE3    7    0
204                               lateSE14    7    0
215                          maxVr_countor    7    0
221                            T2min_F_r_i    7    0
225          texture_contrast_nondir_post4    7    0
226       texture_correlation_nondir_post3    7    0
228       texture_diffentropy_nondir_post3    7    0
235        texture_sumentropy_nondir_post3    7    0
178                                dce2SE2    6    0
180                                dce2SE6    6    0
181                                dce2SE9    6    0
197                               earlySE3    6    0
200                      k_Max_Margin_Grad    6    0
213                          max_RGH_var_k    6    0
214                          maxCr_countor    6    0
223            T2texture_sumentropy_nondir    6    0
224                            T2var_F_r_i    6    0
83                           peakCr_inside    5    0
218                          peakVr_inside    4    0
220                           T2mean_F_r_i    4    0
240                         peakCr_countor    2    0
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
1                                 A_inside   10    1      T1wdynamic
3                                  ave_T20   10    1             T2w
4                                  ave_T21   10    1             T2w
5                                 ave_T210   10    1             T2w
6                                 ave_T211   10    1             T2w
7                                 ave_T212   10    1             T2w
8                                 ave_T213   10    1             T2w
9                                 ave_T215   10    1             T2w
10                                ave_T216   10    1             T2w
12                                ave_T219   10    1             T2w
13                                 ave_T22   10    1             T2w
15                                 ave_T25   10    1             T2w
17                                 ave_T27   10    1             T2w
18                                 ave_T29   10    1             T2w
19                            beta_countor   10    1      T1wdynamic
20                             beta_inside   10    1      T1wdynamic
21                             circularity   10    1   T1wmorphology
22                                 dce2SE0   10    1 single-time-Enh
24                                dce2SE12   10    1 single-time-Enh
29                                dce2SE19   10    1 single-time-Enh
34                                 dce3SE0   10    1 single-time-Enh
41                                 dce3SE5   10    1 single-time-Enh
44                                earlySE0   10    1 single-time-Enh
45                                earlySE1   10    1 single-time-Enh
48                               earlySE14   10    1 single-time-Enh
49                               earlySE16   10    1 single-time-Enh
51                               earlySE19   10    1 single-time-Enh
53                                earlySE6   10    1 single-time-Enh
54                                earlySE7   10    1 single-time-Enh
55                                earlySE8   10    1 single-time-Enh
56                                earlySE9   10    1 single-time-Enh
57                         edge_sharp_mean   10    1   T1wmorphology
58                           iAUC1_countor   10    1      T1wdynamic
60                  iiiMax_Margin_Gradient   10    1   T1wmorphology
61            iiMin_change_Variance_uptake   10    1   T1wmorphology
62                    iMax_Variance_uptake   10    1   T1wmorphology
63                            irregularity   10    1   T1wmorphology
64                              ivVariance   10    1   T1wmorphology
65                            Kpeak_inside   10    1      T1wdynamic
66                              kurt_F_r_i   10    1      T1wdynamic
69                                lateSE11   10    1 single-time-Enh
70                                lateSE12   10    1 single-time-Enh
71                                lateSE17   10    1 single-time-Enh
73                                 lateSE3   10    1 single-time-Enh
75                                 lateSE8   10    1 single-time-Enh
76                                 lateSE9   10    1 single-time-Enh
77                               max_F_r_i   10    1      T1wdynamic
78                            max_RGH_mean   10    1   T1wmorphology
79                             max_RGH_var   10    1   T1wmorphology
82                              mean_F_r_i   10    1      T1wdynamic
84                             SER_countor   10    1      T1wdynamic
85                              SER_inside   10    1      T1wdynamic
86                              skew_F_r_i   10    1      T1wdynamic
91                       T2grad_margin_var   10    1             T2w
92                            T2kurt_F_r_i   10    1             T2w
93                             T2max_F_r_i   10    1             T2w
94                              T2RGH_mean   10    1             T2w
95                               T2RGH_var   10    1             T2w
97            T2texture_correlation_nondir   10    1             T2w
101               T2texture_entropy_nondir   10    1             T2w
109       texture_correlation_nondir_post1   10    1      T1wtexture
111       texture_correlation_nondir_post4   10    1      T1wtexture
114      texture_diffvariance_nondir_post1   10    1      T1wtexture
115      texture_diffvariance_nondir_post2   10    1      T1wtexture
118            texture_energy_nondir_post3   10    1      T1wtexture
119            texture_energy_nondir_post4   10    1      T1wtexture
122 texture_inversediffmoment_nondir_post2   10    1      T1wtexture
123 texture_inversediffmoment_nondir_post3   10    1      T1wtexture
124 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
125        texture_sumaverage_nondir_post1   10    1      T1wtexture
126        texture_sumaverage_nondir_post2   10    1      T1wtexture
127        texture_sumaverage_nondir_post4   10    1      T1wtexture
128        texture_sumentropy_nondir_post1   10    1      T1wtexture
131       texture_sumvariance_nondir_post1   10    1      T1wtexture
133       texture_sumvariance_nondir_post4   10    1      T1wtexture
134          texture_variance_nondir_post1   10    1      T1wtexture
135          texture_variance_nondir_post2   10    1      T1wtexture
136                          Tpeak_countor   10    1      T1wdynamic
137                           Tpeak_inside   10    1      T1wdynamic
141                                     V1   10    1      dispersion
142                                    V10   10    1      dispersion
143                                    V11   10    1      dispersion
144                                    V12   10    1      dispersion
145                                    V13   10    1      dispersion
147                                    V15   10    1      dispersion
148                                    V16   10    1      dispersion
149                                    V17   10    1      dispersion
150                                    V18   10    1      dispersion
151                                    V19   10    1      dispersion
154                                     V5   10    1      dispersion
156                                     V7   10    1      dispersion
157                                     V8   10    1      dispersion
158                                     V9   10    1      dispersion
159                              var_F_r_i   10    1      T1wdynamic
160              Vr_decreasingRate_countor   10    1      T1wdynamic
161               Vr_decreasingRate_inside   10    1      T1wdynamic
162              Vr_increasingRate_countor   10    1      T1wdynamic
163               Vr_increasingRate_inside   10    1      T1wdynamic
164                      Vr_post_1_countor   10    1      T1wdynamic
165                       Vr_post_1_inside   10    1      T1wdynamic
2                            alpha_countor    9    1      T1wdynamic
23                                 dce2SE1    9    1 single-time-Enh
25                                dce2SE13    9    1 single-time-Enh
26                                dce2SE14    9    1 single-time-Enh
27                                dce2SE16    9    1 single-time-Enh
28                                dce2SE18    9    1 single-time-Enh
30                                 dce2SE4    9    1 single-time-Enh
31                                 dce2SE5    9    1 single-time-Enh
32                                 dce2SE7    9    1 single-time-Enh
33                                 dce2SE8    9    1 single-time-Enh
35                                dce3SE11    9    1 single-time-Enh
36                                dce3SE13    9    1 single-time-Enh
38                                dce3SE19    9    1 single-time-Enh
39                                 dce3SE2    9    1 single-time-Enh
42                                 dce3SE7    9    1 single-time-Enh
43                                 dce3SE8    9    1 single-time-Enh
47                               earlySE12    9    1 single-time-Enh
52                                earlySE5    9    1 single-time-Enh
59                            iAUC1_inside    9    1      T1wdynamic
68                                lateSE10    9    1 single-time-Enh
72                                lateSE19    9    1 single-time-Enh
74                                 lateSE7    9    1 single-time-Enh
81                            maxVr_inside    9    1      T1wdynamic
87                       Slope_ini_countor    9    1      T1wdynamic
88                        Slope_ini_inside    9    1      T1wdynamic
89                             T2_lesionSI    9    1             T2w
90                           T2grad_margin    9    1             T2w
96               T2texture_contrast_nondir    9    1             T2w
98            T2texture_diffentropy_nondir    9    1             T2w
99           T2texture_diffvariance_nondir    9    1             T2w
100                T2texture_energy_nondir    9    1             T2w
102     T2texture_inversediffmoment_nondir    9    1             T2w
103            T2texture_sumaverage_nondir    9    1             T2w
106          texture_contrast_nondir_post1    9    1      T1wtexture
108          texture_contrast_nondir_post3    9    1      T1wtexture
110       texture_correlation_nondir_post2    9    1      T1wtexture
113       texture_diffentropy_nondir_post4    9    1      T1wtexture
116      texture_diffvariance_nondir_post3    9    1      T1wtexture
120           texture_entropy_nondir_post1    9    1      T1wtexture
129        texture_sumentropy_nondir_post2    9    1      T1wtexture
139                      UptakeRate_inside    9    1      T1wdynamic
146                                    V14    9    1      dispersion
152                                     V3    9    1      dispersion
153                                     V4    9    1      dispersion
155                                     V6    9    1      dispersion
166                    washoutRate_countor    9    1      T1wdynamic
167                     washoutRate_inside    9    1      T1wdynamic
168                              A_countor    9    1      T1wdynamic
170                               ave_T214    9    1             T2w
171                               ave_T218    9    1             T2w
183                               dce3SE10    9    1 single-time-Enh
184                               dce3SE12    9    1 single-time-Enh
185                               dce3SE15    9    1 single-time-Enh
188                               dce3SE18    9    1 single-time-Enh
192                              earlySE10    9    1 single-time-Enh
194                              earlySE15    9    1 single-time-Enh
199                         edge_sharp_std    9    1   T1wmorphology
201                          Kpeak_countor    9    1      T1wdynamic
203                               lateSE13    9    1 single-time-Enh
209                                lateSE4    9    1 single-time-Enh
211                                lateSE6    9    1 single-time-Enh
234        texture_sumaverage_nondir_post3    9    1      T1wtexture
236       texture_sumvariance_nondir_post3    9    1      T1wtexture
237          texture_variance_nondir_post3    9    1      T1wtexture
239                                     V2    9    1      dispersion
16                                 ave_T26    8    1             T2w
37                                dce3SE14    8    1 single-time-Enh
40                                 dce3SE4    8    1 single-time-Enh
46                               earlySE11    8    1 single-time-Enh
50                               earlySE17    8    1 single-time-Enh
67                                 lateSE0    8    1 single-time-Enh
80                            maxCr_inside    8    1      T1wdynamic
104           T2texture_sumvariance_nondir    8    1             T2w
107          texture_contrast_nondir_post2    8    1      T1wtexture
117            texture_energy_nondir_post1    8    1      T1wtexture
121           texture_entropy_nondir_post3    8    1      T1wtexture
130        texture_sumentropy_nondir_post4    8    1      T1wtexture
132       texture_sumvariance_nondir_post2    8    1      T1wtexture
138                     UptakeRate_countor    8    1      T1wdynamic
140                                     V0    8    1      dispersion
169                           alpha_inside    8    1      T1wdynamic
173                                ave_T28    8    1             T2w
174                               dce2SE10    8    1 single-time-Enh
176                               dce2SE15    8    1 single-time-Enh
177                               dce2SE17    8    1 single-time-Enh
182                                dce3SE1    8    1 single-time-Enh
186                               dce3SE16    8    1 single-time-Enh
190                                dce3SE6    8    1 single-time-Enh
191                                dce3SE9    8    1 single-time-Enh
193                              earlySE13    8    1 single-time-Enh
195                              earlySE18    8    1 single-time-Enh
196                               earlySE2    8    1 single-time-Enh
198                               earlySE4    8    1 single-time-Enh
202                                lateSE1    8    1 single-time-Enh
205                               lateSE15    8    1 single-time-Enh
206                               lateSE16    8    1 single-time-Enh
207                               lateSE18    8    1 single-time-Enh
208                                lateSE2    8    1 single-time-Enh
210                                lateSE5    8    1 single-time-Enh
212                         max_RGH_mean_k    8    1   T1wmorphology
216                              min_F_r_i    8    1      T1wdynamic
217                         peakVr_countor    8    1      T1wdynamic
219                         T2_lesionSIstd    8    1             T2w
222                           T2skew_F_r_i    8    1             T2w
227       texture_diffentropy_nondir_post1    8    1      T1wtexture
229      texture_diffvariance_nondir_post4    8    1      T1wtexture
230            texture_energy_nondir_post2    8    1      T1wtexture
231           texture_entropy_nondir_post2    8    1      T1wtexture
232           texture_entropy_nondir_post4    8    1      T1wtexture
233 texture_inversediffmoment_nondir_post1    8    1      T1wtexture
238          texture_variance_nondir_post4    8    1      T1wtexture
11                                ave_T217    7    0             T2w
14                                 ave_T24    7    0             T2w
105              T2texture_variance_nondir    7    0             T2w
112       texture_diffentropy_nondir_post2    7    0      T1wtexture
172                                ave_T23    7    0             T2w
175                               dce2SE11    7    0 single-time-Enh
179                                dce2SE3    7    0 single-time-Enh
187                               dce3SE17    7    0 single-time-Enh
189                                dce3SE3    7    0 single-time-Enh
204                               lateSE14    7    0 single-time-Enh
215                          maxVr_countor    7    0      T1wdynamic
221                            T2min_F_r_i    7    0             T2w
225          texture_contrast_nondir_post4    7    0      T1wtexture
226       texture_correlation_nondir_post3    7    0      T1wtexture
228       texture_diffentropy_nondir_post3    7    0      T1wtexture
235        texture_sumentropy_nondir_post3    7    0      T1wtexture
178                                dce2SE2    6    0 single-time-Enh
180                                dce2SE6    6    0 single-time-Enh
181                                dce2SE9    6    0 single-time-Enh
197                               earlySE3    6    0 single-time-Enh
200                      k_Max_Margin_Grad    6    0   T1wmorphology
213                          max_RGH_var_k    6    0   T1wmorphology
214                          maxCr_countor    6    0      T1wdynamic
223            T2texture_sumentropy_nondir    6    0             T2w
224                            T2var_F_r_i    6    0             T2w
83                           peakCr_inside    5    0      T1wdynamic
218                          peakVr_inside    4    0      T1wdynamic
220                           T2mean_F_r_i    4    0             T2w
240                         peakCr_countor    2    0      T1wdynamic
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/featuressel-2.png)<!-- -->

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
6                              beta_inside   10    1
7                              circularity   10    1
17                                dce2SE17   10    1
40                                 dce3SE2   10    1
48                                earlySE0   10    1
50                               earlySE10   10    1
56                               earlySE16   10    1
64                                earlySE6   10    1
66                                earlySE8   10    1
73            iiMin_change_Variance_uptake   10    1
74                    iMax_Variance_uptake   10    1
75                            irregularity   10    1
76                              ivVariance   10    1
79                            Kpeak_inside   10    1
81                                 lateSE0   10    1
84                                lateSE11   10    1
90                                lateSE17   10    1
97                                 lateSE6   10    1
100                                lateSE9   10    1
103                           max_RGH_mean   10    1
105                            max_RGH_var   10    1
108                           maxCr_inside   10    1
110                           maxVr_inside   10    1
111                             mean_F_r_i   10    1
112                              min_F_r_i   10    1
117                             SER_inside   10    1
122                      T2grad_margin_var   10    1
123                           T2kurt_F_r_i   10    1
127                             T2RGH_mean   10    1
128                              T2RGH_var   10    1
129                           T2skew_F_r_i   10    1
131           T2texture_correlation_nondir   10    1
135               T2texture_entropy_nondir   10    1
136     T2texture_inversediffmoment_nondir   10    1
139           T2texture_sumvariance_nondir   10    1
144          texture_contrast_nondir_post3   10    1
146       texture_correlation_nondir_post1   10    1
147       texture_correlation_nondir_post2   10    1
154      texture_diffvariance_nondir_post1   10    1
161            texture_energy_nondir_post4   10    1
164           texture_entropy_nondir_post3   10    1
167 texture_inversediffmoment_nondir_post2   10    1
169 texture_inversediffmoment_nondir_post4   10    1
171        texture_sumaverage_nondir_post2   10    1
178       texture_sumvariance_nondir_post1   10    1
182          texture_variance_nondir_post1   10    1
183          texture_variance_nondir_post2   10    1
186                          Tpeak_countor   10    1
189                      UptakeRate_inside   10    1
192                                    V10   10    1
194                                    V12   10    1
202                                     V2   10    1
204                                     V4   10    1
205                                     V5   10    1
209                                     V9   10    1
210                              var_F_r_i   10    1
212               Vr_decreasingRate_inside   10    1
213              Vr_increasingRate_countor   10    1
214               Vr_increasingRate_inside   10    1
215                      Vr_post_1_countor   10    1
216                       Vr_post_1_inside   10    1
2                                 A_inside    9    1
3                            alpha_countor    9    1
5                             beta_countor    9    1
8                                  dce2SE0    9    1
10                                dce2SE10    9    1
15                                dce2SE15    9    1
18                                dce2SE18    9    1
22                                 dce2SE4    9    1
23                                 dce2SE5    9    1
24                                 dce2SE6    9    1
26                                 dce2SE8    9    1
28                                 dce3SE0    9    1
30                                dce3SE10    9    1
32                                dce3SE12    9    1
38                                dce3SE18    9    1
41                                 dce3SE3    9    1
43                                 dce3SE5    9    1
44                                 dce3SE6    9    1
52                               earlySE12    9    1
55                               earlySE15    9    1
57                               earlySE17    9    1
60                                earlySE2    9    1
63                                earlySE5    9    1
67                                earlySE9    9    1
68                         edge_sharp_mean    9    1
70                           iAUC1_countor    9    1
71                            iAUC1_inside    9    1
72                  iiiMax_Margin_Gradient    9    1
78                           Kpeak_countor    9    1
80                              kurt_F_r_i    9    1
85                                lateSE12    9    1
86                                lateSE13    9    1
87                                lateSE14    9    1
89                                lateSE16    9    1
91                                lateSE18    9    1
92                                lateSE19    9    1
93                                 lateSE2    9    1
94                                 lateSE3    9    1
95                                 lateSE4    9    1
96                                 lateSE5    9    1
98                                 lateSE7    9    1
101                        LMSIR_predicted    9    1
102                              max_F_r_i    9    1
116                            SER_countor    9    1
118                             skew_F_r_i    9    1
119                      Slope_ini_countor    9    1
120                       Slope_ini_inside    9    1
124                            T2max_F_r_i    9    1
125                           T2mean_F_r_i    9    1
138            T2texture_sumentropy_nondir    9    1
141                            T2var_F_r_i    9    1
145          texture_contrast_nondir_post4    9    1
148       texture_correlation_nondir_post3    9    1
152       texture_diffentropy_nondir_post3    9    1
155      texture_diffvariance_nondir_post2    9    1
158            texture_energy_nondir_post1    9    1
160            texture_energy_nondir_post3    9    1
162           texture_entropy_nondir_post1    9    1
168 texture_inversediffmoment_nondir_post3    9    1
170        texture_sumaverage_nondir_post1    9    1
173        texture_sumaverage_nondir_post4    9    1
174        texture_sumentropy_nondir_post1    9    1
180       texture_sumvariance_nondir_post3    9    1
187                           Tpeak_inside    9    1
190                                     V0    9    1
191                                     V1    9    1
193                                    V11    9    1
195                                    V13    9    1
196                                    V14    9    1
198                                    V16    9    1
199                                    V17    9    1
200                                    V18    9    1
201                                    V19    9    1
206                                     V6    9    1
208                                     V8    9    1
211              Vr_decreasingRate_countor    9    1
217                    washoutRate_countor    9    1
218                     washoutRate_inside    9    1
4                             alpha_inside    8    1
9                                  dce2SE1    8    1
11                                dce2SE11    8    1
13                                dce2SE13    8    1
14                                dce2SE14    8    1
16                                dce2SE16    8    1
21                                 dce2SE3    8    1
25                                 dce2SE7    8    1
29                                 dce3SE1    8    1
31                                dce3SE11    8    1
34                                dce3SE14    8    1
35                                dce3SE15    8    1
36                                dce3SE16    8    1
39                                dce3SE19    8    1
42                                 dce3SE4    8    1
45                                 dce3SE7    8    1
53                               earlySE13    8    1
54                               earlySE14    8    1
61                                earlySE3    8    1
62                                earlySE4    8    1
69                          edge_sharp_std    8    1
83                                lateSE10    8    1
88                                lateSE15    8    1
99                                 lateSE8    8    1
106                          max_RGH_var_k    8    1
107                          maxCr_countor    8    1
109                          maxVr_countor    8    1
121                          T2grad_margin    8    1
126                            T2min_F_r_i    8    1
130              T2texture_contrast_nondir    8    1
133          T2texture_diffvariance_nondir    8    1
134                T2texture_energy_nondir    8    1
137            T2texture_sumaverage_nondir    8    1
140              T2texture_variance_nondir    8    1
149       texture_correlation_nondir_post4    8    1
150       texture_diffentropy_nondir_post1    8    1
156      texture_diffvariance_nondir_post3    8    1
157      texture_diffvariance_nondir_post4    8    1
159            texture_energy_nondir_post2    8    1
163           texture_entropy_nondir_post2    8    1
166 texture_inversediffmoment_nondir_post1    8    1
172        texture_sumaverage_nondir_post3    8    1
179       texture_sumvariance_nondir_post2    8    1
181       texture_sumvariance_nondir_post4    8    1
184          texture_variance_nondir_post3    8    1
197                                    V15    8    1
203                                     V3    8    1
207                                     V7    8    1
12                                dce2SE12    7    0
19                                dce2SE19    7    0
27                                 dce2SE9    7    0
37                                dce3SE17    7    0
46                                 dce3SE8    7    0
49                                earlySE1    7    0
51                               earlySE11    7    0
58                               earlySE18    7    0
59                               earlySE19    7    0
65                                earlySE7    7    0
82                                 lateSE1    7    0
104                         max_RGH_mean_k    7    0
114                         peakVr_countor    7    0
142          texture_contrast_nondir_post1    7    0
143          texture_contrast_nondir_post2    7    0
151       texture_diffentropy_nondir_post2    7    0
153       texture_diffentropy_nondir_post4    7    0
165           texture_entropy_nondir_post4    7    0
176        texture_sumentropy_nondir_post3    7    0
177        texture_sumentropy_nondir_post4    7    0
20                                 dce2SE2    6    0
33                                dce3SE13    6    0
47                                 dce3SE9    6    0
132           T2texture_diffentropy_nondir    6    0
175        texture_sumentropy_nondir_post2    6    0
185          texture_variance_nondir_post4    6    0
188                     UptakeRate_countor    6    0
77                       k_Max_Margin_Grad    4    0
113                          peakCr_inside    3    0
115                          peakVr_inside    3    0
219                         peakCr_countor    1    0
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
6                              beta_inside   10    1      T1wdynamic
7                              circularity   10    1   T1wmorphology
17                                dce2SE17   10    1 single-time-Enh
40                                 dce3SE2   10    1 single-time-Enh
48                                earlySE0   10    1 single-time-Enh
50                               earlySE10   10    1 single-time-Enh
56                               earlySE16   10    1 single-time-Enh
64                                earlySE6   10    1 single-time-Enh
66                                earlySE8   10    1 single-time-Enh
73            iiMin_change_Variance_uptake   10    1   T1wmorphology
74                    iMax_Variance_uptake   10    1   T1wmorphology
75                            irregularity   10    1   T1wmorphology
76                              ivVariance   10    1   T1wmorphology
79                            Kpeak_inside   10    1      T1wdynamic
81                                 lateSE0   10    1 single-time-Enh
84                                lateSE11   10    1 single-time-Enh
90                                lateSE17   10    1 single-time-Enh
97                                 lateSE6   10    1 single-time-Enh
100                                lateSE9   10    1 single-time-Enh
103                           max_RGH_mean   10    1   T1wmorphology
105                            max_RGH_var   10    1   T1wmorphology
108                           maxCr_inside   10    1      T1wdynamic
110                           maxVr_inside   10    1      T1wdynamic
111                             mean_F_r_i   10    1      T1wdynamic
112                              min_F_r_i   10    1      T1wdynamic
117                             SER_inside   10    1      T1wdynamic
122                      T2grad_margin_var   10    1             T2w
123                           T2kurt_F_r_i   10    1             T2w
127                             T2RGH_mean   10    1             T2w
128                              T2RGH_var   10    1             T2w
129                           T2skew_F_r_i   10    1             T2w
131           T2texture_correlation_nondir   10    1             T2w
135               T2texture_entropy_nondir   10    1             T2w
136     T2texture_inversediffmoment_nondir   10    1             T2w
139           T2texture_sumvariance_nondir   10    1             T2w
144          texture_contrast_nondir_post3   10    1      T1wtexture
146       texture_correlation_nondir_post1   10    1      T1wtexture
147       texture_correlation_nondir_post2   10    1      T1wtexture
154      texture_diffvariance_nondir_post1   10    1      T1wtexture
161            texture_energy_nondir_post4   10    1      T1wtexture
164           texture_entropy_nondir_post3   10    1      T1wtexture
167 texture_inversediffmoment_nondir_post2   10    1      T1wtexture
169 texture_inversediffmoment_nondir_post4   10    1      T1wtexture
171        texture_sumaverage_nondir_post2   10    1      T1wtexture
178       texture_sumvariance_nondir_post1   10    1      T1wtexture
182          texture_variance_nondir_post1   10    1      T1wtexture
183          texture_variance_nondir_post2   10    1      T1wtexture
186                          Tpeak_countor   10    1      T1wdynamic
189                      UptakeRate_inside   10    1      T1wdynamic
192                                    V10   10    1      dispersion
194                                    V12   10    1      dispersion
202                                     V2   10    1      dispersion
204                                     V4   10    1      dispersion
205                                     V5   10    1      dispersion
209                                     V9   10    1      dispersion
210                              var_F_r_i   10    1      T1wdynamic
212               Vr_decreasingRate_inside   10    1      T1wdynamic
213              Vr_increasingRate_countor   10    1      T1wdynamic
214               Vr_increasingRate_inside   10    1      T1wdynamic
215                      Vr_post_1_countor   10    1      T1wdynamic
216                       Vr_post_1_inside   10    1      T1wdynamic
2                                 A_inside    9    1      T1wdynamic
3                            alpha_countor    9    1      T1wdynamic
5                             beta_countor    9    1      T1wdynamic
8                                  dce2SE0    9    1 single-time-Enh
10                                dce2SE10    9    1 single-time-Enh
15                                dce2SE15    9    1 single-time-Enh
18                                dce2SE18    9    1 single-time-Enh
22                                 dce2SE4    9    1 single-time-Enh
23                                 dce2SE5    9    1 single-time-Enh
24                                 dce2SE6    9    1 single-time-Enh
26                                 dce2SE8    9    1 single-time-Enh
28                                 dce3SE0    9    1 single-time-Enh
30                                dce3SE10    9    1 single-time-Enh
32                                dce3SE12    9    1 single-time-Enh
38                                dce3SE18    9    1 single-time-Enh
41                                 dce3SE3    9    1 single-time-Enh
43                                 dce3SE5    9    1 single-time-Enh
44                                 dce3SE6    9    1 single-time-Enh
52                               earlySE12    9    1 single-time-Enh
55                               earlySE15    9    1 single-time-Enh
57                               earlySE17    9    1 single-time-Enh
60                                earlySE2    9    1 single-time-Enh
63                                earlySE5    9    1 single-time-Enh
67                                earlySE9    9    1 single-time-Enh
68                         edge_sharp_mean    9    1   T1wmorphology
70                           iAUC1_countor    9    1      T1wdynamic
71                            iAUC1_inside    9    1      T1wdynamic
72                  iiiMax_Margin_Gradient    9    1   T1wmorphology
78                           Kpeak_countor    9    1      T1wdynamic
80                              kurt_F_r_i    9    1      T1wdynamic
85                                lateSE12    9    1 single-time-Enh
86                                lateSE13    9    1 single-time-Enh
87                                lateSE14    9    1 single-time-Enh
89                                lateSE16    9    1 single-time-Enh
91                                lateSE18    9    1 single-time-Enh
92                                lateSE19    9    1 single-time-Enh
93                                 lateSE2    9    1 single-time-Enh
94                                 lateSE3    9    1 single-time-Enh
95                                 lateSE4    9    1 single-time-Enh
96                                 lateSE5    9    1 single-time-Enh
98                                 lateSE7    9    1 single-time-Enh
101                        LMSIR_predicted    9    1             T2w
102                              max_F_r_i    9    1      T1wdynamic
116                            SER_countor    9    1      T1wdynamic
118                             skew_F_r_i    9    1      T1wdynamic
119                      Slope_ini_countor    9    1      T1wdynamic
120                       Slope_ini_inside    9    1      T1wdynamic
124                            T2max_F_r_i    9    1             T2w
125                           T2mean_F_r_i    9    1             T2w
138            T2texture_sumentropy_nondir    9    1             T2w
141                            T2var_F_r_i    9    1             T2w
145          texture_contrast_nondir_post4    9    1      T1wtexture
148       texture_correlation_nondir_post3    9    1      T1wtexture
152       texture_diffentropy_nondir_post3    9    1      T1wtexture
155      texture_diffvariance_nondir_post2    9    1      T1wtexture
158            texture_energy_nondir_post1    9    1      T1wtexture
160            texture_energy_nondir_post3    9    1      T1wtexture
162           texture_entropy_nondir_post1    9    1      T1wtexture
168 texture_inversediffmoment_nondir_post3    9    1      T1wtexture
170        texture_sumaverage_nondir_post1    9    1      T1wtexture
173        texture_sumaverage_nondir_post4    9    1      T1wtexture
174        texture_sumentropy_nondir_post1    9    1      T1wtexture
180       texture_sumvariance_nondir_post3    9    1      T1wtexture
187                           Tpeak_inside    9    1      T1wdynamic
190                                     V0    9    1      dispersion
191                                     V1    9    1      dispersion
193                                    V11    9    1      dispersion
195                                    V13    9    1      dispersion
196                                    V14    9    1      dispersion
198                                    V16    9    1      dispersion
199                                    V17    9    1      dispersion
200                                    V18    9    1      dispersion
201                                    V19    9    1      dispersion
206                                     V6    9    1      dispersion
208                                     V8    9    1      dispersion
211              Vr_decreasingRate_countor    9    1      T1wdynamic
217                    washoutRate_countor    9    1      T1wdynamic
218                     washoutRate_inside    9    1      T1wdynamic
4                             alpha_inside    8    1      T1wdynamic
9                                  dce2SE1    8    1 single-time-Enh
11                                dce2SE11    8    1 single-time-Enh
13                                dce2SE13    8    1 single-time-Enh
14                                dce2SE14    8    1 single-time-Enh
16                                dce2SE16    8    1 single-time-Enh
21                                 dce2SE3    8    1 single-time-Enh
25                                 dce2SE7    8    1 single-time-Enh
29                                 dce3SE1    8    1 single-time-Enh
31                                dce3SE11    8    1 single-time-Enh
34                                dce3SE14    8    1 single-time-Enh
35                                dce3SE15    8    1 single-time-Enh
36                                dce3SE16    8    1 single-time-Enh
39                                dce3SE19    8    1 single-time-Enh
42                                 dce3SE4    8    1 single-time-Enh
45                                 dce3SE7    8    1 single-time-Enh
53                               earlySE13    8    1 single-time-Enh
54                               earlySE14    8    1 single-time-Enh
61                                earlySE3    8    1 single-time-Enh
62                                earlySE4    8    1 single-time-Enh
69                          edge_sharp_std    8    1   T1wmorphology
83                                lateSE10    8    1 single-time-Enh
88                                lateSE15    8    1 single-time-Enh
99                                 lateSE8    8    1 single-time-Enh
106                          max_RGH_var_k    8    1   T1wmorphology
107                          maxCr_countor    8    1      T1wdynamic
109                          maxVr_countor    8    1      T1wdynamic
121                          T2grad_margin    8    1             T2w
126                            T2min_F_r_i    8    1             T2w
130              T2texture_contrast_nondir    8    1             T2w
133          T2texture_diffvariance_nondir    8    1             T2w
134                T2texture_energy_nondir    8    1             T2w
137            T2texture_sumaverage_nondir    8    1             T2w
140              T2texture_variance_nondir    8    1             T2w
149       texture_correlation_nondir_post4    8    1      T1wtexture
150       texture_diffentropy_nondir_post1    8    1      T1wtexture
156      texture_diffvariance_nondir_post3    8    1      T1wtexture
157      texture_diffvariance_nondir_post4    8    1      T1wtexture
159            texture_energy_nondir_post2    8    1      T1wtexture
163           texture_entropy_nondir_post2    8    1      T1wtexture
166 texture_inversediffmoment_nondir_post1    8    1      T1wtexture
172        texture_sumaverage_nondir_post3    8    1      T1wtexture
179       texture_sumvariance_nondir_post2    8    1      T1wtexture
181       texture_sumvariance_nondir_post4    8    1      T1wtexture
184          texture_variance_nondir_post3    8    1      T1wtexture
197                                    V15    8    1      dispersion
203                                     V3    8    1      dispersion
207                                     V7    8    1      dispersion
12                                dce2SE12    7    0 single-time-Enh
19                                dce2SE19    7    0 single-time-Enh
27                                 dce2SE9    7    0 single-time-Enh
37                                dce3SE17    7    0 single-time-Enh
46                                 dce3SE8    7    0 single-time-Enh
49                                earlySE1    7    0 single-time-Enh
51                               earlySE11    7    0 single-time-Enh
58                               earlySE18    7    0 single-time-Enh
59                               earlySE19    7    0 single-time-Enh
65                                earlySE7    7    0 single-time-Enh
82                                 lateSE1    7    0 single-time-Enh
104                         max_RGH_mean_k    7    0   T1wmorphology
114                         peakVr_countor    7    0      T1wdynamic
142          texture_contrast_nondir_post1    7    0      T1wtexture
143          texture_contrast_nondir_post2    7    0      T1wtexture
151       texture_diffentropy_nondir_post2    7    0      T1wtexture
153       texture_diffentropy_nondir_post4    7    0      T1wtexture
165           texture_entropy_nondir_post4    7    0      T1wtexture
176        texture_sumentropy_nondir_post3    7    0      T1wtexture
177        texture_sumentropy_nondir_post4    7    0      T1wtexture
20                                 dce2SE2    6    0 single-time-Enh
33                                dce3SE13    6    0 single-time-Enh
47                                 dce3SE9    6    0 single-time-Enh
132           T2texture_diffentropy_nondir    6    0             T2w
175        texture_sumentropy_nondir_post2    6    0      T1wtexture
185          texture_variance_nondir_post4    6    0      T1wtexture
188                     UptakeRate_countor    6    0      T1wdynamic
77                       k_Max_Margin_Grad    4    0   T1wmorphology
113                          peakCr_inside    3    0      T1wdynamic
115                          peakVr_inside    3    0      T1wdynamic
219                         peakCr_countor    1    0      T1wdynamic
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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/featuressel-3.png)<!-- -->


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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/frequency-graphs-1.png)<!-- -->

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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/frequency-graphs-2.png)<!-- -->

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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/frequency-graphs-3.png)<!-- -->

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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/T2featsel-1.png)<!-- -->

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

![](AnalyzeResults_Boosting_classification_trees_T2imgvspredLMSIR_files/figure-html/T2featsel-2.png)<!-- -->
               


```r
save.image("Outputs/T2SIvspredLMSIR_summaryResults.RData")
```

