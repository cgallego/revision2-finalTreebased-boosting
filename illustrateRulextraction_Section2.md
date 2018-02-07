# Rule extraction Illustration with cases - experiments new paper



## Summarize ROC of rules and rulemodel

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
legend("bottomright", legend = c(paste0("probRules"), paste0("probModel")), col = c(colors[15], 
    colors[2]), lty = c(1, 5), lwd = 2)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

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

## find cases using LMSIR rule

```r
LMSIRfeatureflag = c()
for (i in 1:length(eachtestRules)) {
    LMSIRwtop5 = FALSE
    eachtestRules[[i]]$Ruletype
    # recognize whether t2w features are used
    isT2wtop5 = "T2w" %in% unlist(strsplit(eachtestRules[[i]]$Ruletype, split = " & "))
    
    # if isT2wtop5 find if LMSIR
    if (isT2wtop5) {
        LMSIRwtop5 = "LMSIR" %in% unlist(strsplit(unlist(strsplit(eachtestRules[[i]]$condition, 
            split = " & ")), split = " "))
        LMSIRfeatureflag = c(LMSIRfeatureflag, LMSIRwtop5)
    } else {
        LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
    }
}
print(summary(LMSIRfeatureflag))
```

```
   Mode   FALSE    TRUE    NA's 
logical     595      32       0 
```

```r
idxwLMSIR = c(1:length(eachtestRules))[LMSIRfeatureflag]
print(idxwLMSIR)
```

```
 [1]  10  30  54  58  64  66 117 149 166 167 193 205 231 376 386 398 403 416 418 426 451 457
[23] 458 480 489 499 500 532 566 578 610 624
```


## Case by case illustration: Using T2w SI & Texture

```r
# exemplify one # Texture T2w & SI
lesion_id = 41
# or idx= 41
idx = c(1:nrow(allTestinfo))[allTestinfo$lesion_id == lesion_id]
print(allTestinfo[idx, -c(20)])
```

```
   lesion_id          lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
41        41 41_0190_6760690.vtk          0190 1975-04-25 00:00:00.000000              17213
   exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
41           6760690 2011-04-12 00:00:00.000000               Malignant
   cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
41                          BRCA2                     0                        1
   exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
41              Right                 2092 2011-04-28 00:00:00.000000              Right
   proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int proc_lesion_comments_txt
41            Radiology                 US Core Needle Biopsy                     None
   find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
41           None                       N/A                        N/A      4     nonmassM
   lesion_diagnosis
41   InvasiveDuctal
```

```r
# get LMSIR
LMSIR = LMSIR_cv[LMSIR_cv$lesion_id == allTestinfo[idx, "lesion_id"], ]
print(LMSIR)
```

```
   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime BIRADS lesion_label
41        41          0190           6760690 2011-04-12 00:00:00.000000      4     nonmassM
   lesion_diagnosis fold LMSIR_measured LMSIR_predicted mean_sqrtMSEtrain sqrtMSEtest
41   InvasiveDuctal    1       2.277133        2.243625          1.444012  0.03350782
```

```r
selRulesIx = eachtestRules[[idx]]
# average error rate of 5 top scoring rules
print(sum(as.numeric(selRulesIx$err))/nrow(selRulesIx))
```

```
[1] 0.0862
```

```r
# present rules
print(thresh)
```

```
[1] 0.5103875
```

```r
myrules = selRulesIx
myrules$pred = ifelse(myrules$sumtempProb >= thresh, "C", "NC")
print(myrules)
```

```
     len  freq   err
1923   3 0.098 0.074
923    3  0.12 0.076
7803   2 0.084 0.087
593    2 0.248 0.096
6442   3 0.093 0.098
                                                                                               condition
1923        min uptake low to very high  & 2nd post-SE s6 low to very high  & last post-SE s11 very low 
923        min uptake low to very high  & 2nd post-SE s12 low to very high  & last post-SE s11 very low 
7803                                                   min uptake very high  & 3rd post-SE s17 very low 
593                                                 SER(in) low  & max Radial gradient low to very high 
6442 last post-SE s8 low to very high  & last post-SE s11 very low  & T2w Correlation  very low to high 
     pred sumtempProb                     Ruletype
1923   NC     0.45589 T1wdynamic & single-time-Enh
923    NC   0.4745757 T1wdynamic & single-time-Enh
7803    C   0.5477606 T1wdynamic & single-time-Enh
593    NC   0.3858172   T1wdynamic & T1wmorphology
6442    C   0.5173376        single-time-Enh & T2w
```

```r
print(xtable(myrules), include.rownames = FALSE)
```

```
% latex table generated in R 3.3.1 by xtable 1.8-2 package
% Thu Sep 07 11:34:50 2017
\begin{table}[ht]
\centering
\begin{tabular}{lllllrl}
  \hline
len & freq & err & condition & pred & sumtempProb & Ruletype \\ 
  \hline
3 & 0.098 & 0.074 & min uptake low to very high  \& 2nd post-SE s6 low to very high  \& last post-SE s11 very low  & NC & 0.46 & T1wdynamic \& single-time-Enh \\ 
  3 & 0.12 & 0.076 & min uptake low to very high  \& 2nd post-SE s12 low to very high  \& last post-SE s11 very low  & NC & 0.47 & T1wdynamic \& single-time-Enh \\ 
  2 & 0.084 & 0.087 & min uptake very high  \& 3rd post-SE s17 very low  & C & 0.55 & T1wdynamic \& single-time-Enh \\ 
  2 & 0.248 & 0.096 & SER(in) low  \& max Radial gradient low to very high  & NC & 0.39 & T1wdynamic \& T1wmorphology \\ 
  3 & 0.093 & 0.098 & last post-SE s8 low to very high  \& last post-SE s11 very low  \& T2w Correlation  very low to high  & C & 0.52 & single-time-Enh \& T2w \\ 
   \hline
\end{tabular}
\end{table}
```

```r
# averprobTop5[idx] of prediction
print(averprobTop5[idx])
```

```
        probModel 
"0.1924839341693" 
```

```r
# classTop5[idx] based on prob >= theshold found via ROC space
print(classTop5[idx])
```

```
probModel 
     "NC" 
```



![](Z:/Cristina/Section2/papernew_notes/images/41_0190_6760690_nonmassM_None.png)


### Case by case illustration: cases using LMSIR rule

```r
# exemplify one # LMSIR
lesion_id = 10
# or idx= 153 #188 149 or 119 or 331
idx = c(1:nrow(allTestinfo))[allTestinfo$lesion_id == lesion_id]
print(allTestinfo[idx, -c(20)])
```

```
   lesion_id          lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
10        10 10_0102_4755778.vtk          0102 1974-05-19 00:00:00.000000               5142
   exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
10           4755778 2009-03-19 00:00:00.000000     Benign by pathology
   cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
10                      High Risk                     1                        0
   exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
10              Right                 2708 2009-04-16 00:00:00.000000              Right
   proc_proc_source_int proc_proc_guid_int       proc_proc_tp_int proc_lesion_comments_txt
10            Radiology                MRI Vacuum Assisted Biopsy                     None
   find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
10            III                       N/A                    Washout      4        massB
   lesion_diagnosis
10      FIBROCYSTIC
```

```r
selRulesIx = eachtestRules[[idx]]
# average error rate of 5 top scoring rules
print(sum(as.numeric(selRulesIx$err))/nrow(selRulesIx))
```

```
[1] 0.0954
```

```r
# get LMSIR
LMSIR = LMSIR_cv[LMSIR_cv$lesion_id == allTestinfo[idx, "lesion_id"], ]
print(LMSIR)
```

```
   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime BIRADS lesion_label
10        10          0102           4755778 2009-03-19 00:00:00.000000      4        massB
   lesion_diagnosis fold LMSIR_measured LMSIR_predicted mean_sqrtMSEtrain sqrtMSEtest
10      FIBROCYSTIC    5       3.961557        3.547677          1.371497   0.4138799
```

```r
# present rules
print(thresh)
```

```
[1] 0.5103875
```

```r
myrules = selRulesIx
myrules$pred = ifelse(myrules$sumtempProb >= thresh, "C", "NC")
print(myrules)
```

```
     len  freq                err
715    3  0.08 0.0679999999999999
478    3 0.156              0.093
2203   3 0.093              0.098
501    2 0.255              0.107
194    3 0.164              0.111
                                                                                                                    condition
715                          uptake average high  & Inverse difference moment ptime2 very low to high  & 1st post-SE s8 high 
478  irregularity low  & std 3D Sharpness of lesion margin  medium to high  & Inverse difference moment ptime3 low to medium 
2203                         uptake average high  & Inverse difference moment ptime2 very low to high  & 1st post-SE s8 high 
501                                                       irregularity low  & Inverse difference moment ptime2 low to medium 
194               irregularity low  & Sum variance ptime1 low to very high  & Inverse difference moment ptime4 low to medium 
     pred sumtempProb                                  Ruletype
715     C   0.5612487 T1wdynamic & T1wtexture & single-time-Enh
478    NC   0.4570242                T1wmorphology & T1wtexture
2203    C   0.5576801 T1wdynamic & T1wtexture & single-time-Enh
501    NC   0.3823795                T1wmorphology & T1wtexture
194    NC   0.4766391                T1wmorphology & T1wtexture
```

```r
print(xtable(myrules), include.rownames = FALSE)
```

```
% latex table generated in R 3.3.1 by xtable 1.8-2 package
% Thu Sep 07 11:34:51 2017
\begin{table}[ht]
\centering
\begin{tabular}{lllllrl}
  \hline
len & freq & err & condition & pred & sumtempProb & Ruletype \\ 
  \hline
3 & 0.08 & 0.0679999999999999 & uptake average high  \& Inverse difference moment ptime2 very low to high  \& 1st post-SE s8 high  & C & 0.56 & T1wdynamic \& T1wtexture \& single-time-Enh \\ 
  3 & 0.156 & 0.093 & irregularity low  \& std 3D Sharpness of lesion margin  medium to high  \& Inverse difference moment ptime3 low to medium  & NC & 0.46 & T1wmorphology \& T1wtexture \\ 
  3 & 0.093 & 0.098 & uptake average high  \& Inverse difference moment ptime2 very low to high  \& 1st post-SE s8 high  & C & 0.56 & T1wdynamic \& T1wtexture \& single-time-Enh \\ 
  2 & 0.255 & 0.107 & irregularity low  \& Inverse difference moment ptime2 low to medium  & NC & 0.38 & T1wmorphology \& T1wtexture \\ 
  3 & 0.164 & 0.111 & irregularity low  \& Sum variance ptime1 low to very high  \& Inverse difference moment ptime4 low to medium  & NC & 0.48 & T1wmorphology \& T1wtexture \\ 
   \hline
\end{tabular}
\end{table}
```

```r
# averprobTop5[idx] of prediction
print(averprobTop5[idx])
```

```
          probModel 
"0.602796849586872" 
```

```r
# classTop5[idx] based on prob >= theshold found via ROC space
print(classTop5[idx])
```

```
probModel 
      "C" 
```

![](Z:/Cristina/Section2/papernew_notes/images/10_0102_4755778_massB_Hyperintense.png)


## Some Illustrations
### Based on IDC, ISDC, FA, FI ... cannonical examples of each

```r
data = data.frame(cbind(allTraininfo, allTestinfo$lesion_id, allTestinfo$lesion_label, allTestinfo$lesion_diagnosis))

IDC = subset(data, allTestinfo.lesion_diagnosis == "InvasiveDuctal")
summary(IDC$allTestinfo.lesion_label)
```

```
   massB    massM nonmassB nonmassM 
       0      106        0       28 
```

```r
ISDC = subset(data, allTestinfo.lesion_diagnosis == "InsituDuctal")
summary(ISDC$allTestinfo.lesion_label)
```

```
   massB    massM nonmassB nonmassM 
       0       40        0       40 
```

```r
FA = subset(data, allTestinfo.lesion_diagnosis == "FIBROADENOMA")
summary(FA$allTestinfo.lesion_label)
```

```
   massB    massM nonmassB nonmassM 
      59        0       10        0 
```

```r
FI = subset(data, allTestinfo.lesion_diagnosis == "FIBROCYSTIC")
summary(FI$allTestinfo.lesion_label)
```

```
   massB    massM nonmassB nonmassM 
      38        0       32        0 
```

```r
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, "Z:/Cristina/Section2/finalTreebased-boosting/textureUpdatedFeatures.db")

# 2) all T1W features
lesionsQuery <- dbGetQuery(conn, "SELECT *
                           FROM  lesion 
                           INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                           INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                           INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                           INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)
                           INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)
                           INNER JOIN radiologyInfo ON (lesion.lesion_id = radiologyInfo.lesion_id)")

# prune entries and extract feature subsets corresponds to 5 entries lesion info, 34 dynamic, 19
# morpho, 34 texture fueatures
lesioninfo = lesionsQuery[c(1:26, 259)]
# select non foci
lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM")

table(lesioninfo$find_t2_signal_int)
```

```

           Hyperintense Hypointense or not seen                    None 
                     98                     152                     312 
  Slightly hyperintense 
                     65 
```



# 1) Present rules collect top for IDC

```r
library(png)
library(grid)

eachIDCTopRules = list()
allIDCtoprules = c()
for (idx in 1:nrow(IDC)) {
    ############## first analysis
    X = IDC[idx, 2:ncol(imgT2pLMSIRtest)]
    y = IDC[idx, "lesion_label"]
    
    rulesoutputIDC = mynewapplyLearnerxrules(topRulesIDC, X, y, minerr = 0.1, minfrq = 10/627, classes, 
        gbmModel)
    topRulesIDC = rulesoutputIDC[[1]]
    
    if (length(rulesoutputIDC) > 1) {
        eachIDCTopRules[[idx]] = rulesoutputIDC[[2]]
        allIDCtoprules = rbind(allIDCtoprules, rulesoutputIDC[[3]])
    } else {
        eachIDCTopRules[[idx]] = list()
        allIDCtoprules = rbind(allIDCtoprules, 1:nrow(topRulesIDC) * 0)
    }
}
```

```
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 25 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 30 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 29 out of 50"
[1] "test complies with # rules: 28 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 23 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 29 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 28 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 21 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 26 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 28 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 30 out of 50"
[1] "test complies with # rules: 28 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 31 out of 50"
[1] "test complies with # rules: 32 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 29 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 22 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 25 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 1 out of 50"
```

```r
# take the occurance of each explaining rule (columns) per case(row) + lesion identifies
df = data.frame(cbind(allIDCtoprules, as.character(IDC$allTestinfo.lesion_id), as.character(IDC$allTestinfo.lesion_label), 
    as.character(IDC$allTestinfo.lesion_diagnosis)), stringsAsFactors = FALSE)
colnames(df) <- c(as.character(1:nrow(topRulesIDC)), "lesion_id", "lesion_label", "lesion_diagnosis")
cols = 1:nrow(topRulesIDC)
df[, cols] = apply(df[, cols], 2, function(x) as.numeric(x))


cat(paste0("\n============== Top five explaining rules for InvasiveDuctal (IDC): n = ", nrow(IDC) + 
    2, "\n"))
```

```

============== Top five explaining rules for InvasiveDuctal (IDC): n = 136
```

```r
topIDC = sort(colSums(df[, 1:nrow(topRulesIDC)]), decreasing = TRUE)
print(topIDC[1:5])
```

```
43 38 45 21 44 
53 52 50 49 49 
```

```r
freqtopIDC = topIDC[1:5]
totalIDC = nrow(IDC) + 2
print(freqtopIDC/totalIDC)
```

```
       43        38        45        21        44 
0.3897059 0.3823529 0.3676471 0.3602941 0.3602941 
```

```r
rulesTopIDC = as.numeric(names(topIDC[1:5]))
preserulesTopIDC = mypresentRules(topRulesIDC[rulesTopIDC, ], colnames(X), fnames)
rownames(preserulesTopIDC) <- NULL
print(preserulesTopIDC)
```

```
  len  freq   err
1   4 0.159 0.093
2   3  0.16  0.09
3   2 0.156 0.094
4   2 0.144 0.076
5   2 0.155 0.094
                                                                                                                                             condition
1 irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
2                     max relative signal enhancement(rim) medium to high  & Variance ptime1 very high  & Difference variance ptime2 very low to high 
3                                                                                                         SER(in) medium to high  & irregularity high 
4                                                                                                         SER(in) medium to high  & irregularity high 
5                                                                                          Variance ptime1 very high  & dispersion s16 medium to high 
  pred                                       prob sumtempProb                         Ruletype
1    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589 T1wmorphology & T1wtexture & T2w
2    C            0.8237662, 0.4610768, 0.7663881   0.6837437          T1wdynamic & T1wtexture
3    C                       0.7579719, 0.4714506   0.6147112       T1wdynamic & T1wmorphology
4    C                       0.5006988, 0.6626664   0.5816826       T1wdynamic & T1wmorphology
5    C                       0.4291609, 0.7692691    0.599215          T1wtexture & dispersion
```

```r
# display a case that meets them all
casesIDC = df[, c(rulesTopIDC, nrow(topRulesIDC) + 1, nrow(topRulesIDC) + 2, nrow(topRulesIDC) + 
    3)]
print("Sorted cases meeting top five explaining rules for InvasiveDuctal (IDC):")
```

```
[1] "Sorted cases meeting top five explaining rules for InvasiveDuctal (IDC):"
```

```r
topcasesIDC = casesIDC[sort(rowSums(casesIDC[, 1:5]), index = TRUE, decreasing = TRUE)$ix, ]
print(topcasesIDC)
```

```
    43 38 45 21 44 lesion_id lesion_label lesion_diagnosis
10   1  1  1  1  1       515     nonmassM   InvasiveDuctal
11   1  1  1  1  1       538        massM   InvasiveDuctal
19   1  1  1  1  1       235        massM   InvasiveDuctal
21   1  1  1  1  1       299        massM   InvasiveDuctal
44   1  1  1  1  1       168        massM   InvasiveDuctal
49   1  1  1  1  1       320        massM   InvasiveDuctal
51   1  1  1  1  1       397        massM   InvasiveDuctal
52   1  1  1  1  1       430        massM   InvasiveDuctal
58   1  1  1  1  1       199        massM   InvasiveDuctal
59   1  1  1  1  1       200        massM   InvasiveDuctal
63   1  1  1  1  1       544        massM   InvasiveDuctal
75   1  1  1  1  1       341        massM   InvasiveDuctal
77   1  1  1  1  1       417        massM   InvasiveDuctal
80   1  1  1  1  1       154        massM   InvasiveDuctal
83   1  1  1  1  1       169        massM   InvasiveDuctal
92   1  1  1  1  1       560        massM   InvasiveDuctal
93   1  1  1  1  1       561        massM   InvasiveDuctal
94   1  1  1  1  1       567        massM   InvasiveDuctal
103  1  1  1  1  1       241        massM   InvasiveDuctal
104  1  1  1  1  1       280        massM   InvasiveDuctal
109  1  1  1  1  1       554        massM   InvasiveDuctal
119  1  1  1  1  1       160        massM   InvasiveDuctal
124  1  1  1  1  1       392     nonmassM   InvasiveDuctal
128  1  1  1  1  1       520        massM   InvasiveDuctal
130  1  1  1  1  1       526        massM   InvasiveDuctal
131  1  1  1  1  1       527     nonmassM   InvasiveDuctal
133  1  1  1  1  1       529        massM   InvasiveDuctal
5    1  1  1  1  0       217        massM   InvasiveDuctal
24   1  0  1  1  1       398     nonmassM   InvasiveDuctal
37   1  0  1  1  1       404        massM   InvasiveDuctal
45   0  1  1  1  1       182        massM   InvasiveDuctal
53   1  0  1  1  1       559        massM   InvasiveDuctal
88   1  1  1  0  1       394        massM   InvasiveDuctal
117  1  1  0  1  1       504        massM   InvasiveDuctal
120  1  0  1  1  1       161     nonmassM   InvasiveDuctal
125  1  1  1  1  0       393        massM   InvasiveDuctal
3    1  1  0  0  1       153        massM   InvasiveDuctal
9    1  1  0  0  1       514        massM   InvasiveDuctal
12   1  1  0  0  1       539        massM   InvasiveDuctal
15   0  1  1  1  0       139        massM   InvasiveDuctal
18   1  1  0  0  1       209     nonmassM   InvasiveDuctal
31   1  1  0  0  1       121        massM   InvasiveDuctal
34   1  1  0  0  1       171        massM   InvasiveDuctal
41   1  1  0  0  1        12     nonmassM   InvasiveDuctal
42   1  1  0  0  1       149        massM   InvasiveDuctal
60   1  1  0  0  1       210        massM   InvasiveDuctal
105  1  1  0  0  1       281     nonmassM   InvasiveDuctal
110  1  1  0  0  1       569        massM   InvasiveDuctal
8    0  0  1  1  0       424        massM   InvasiveDuctal
16   0  0  1  1  0       172        massM   InvasiveDuctal
20   0  0  1  0  1       282     nonmassM   InvasiveDuctal
29   1  1  0  0  0       513        massM   InvasiveDuctal
30   0  0  1  1  0       570     nonmassM   InvasiveDuctal
33   0  1  0  0  1       124        massM   InvasiveDuctal
36   0  0  1  1  0       303        massM   InvasiveDuctal
38   0  0  1  1  0       498     nonmassM   InvasiveDuctal
43   1  1  0  0  0       150        massM   InvasiveDuctal
47   0  0  1  1  0       252        massM   InvasiveDuctal
70   0  0  1  1  0       144        massM   InvasiveDuctal
71   0  0  1  1  0       177     nonmassM   InvasiveDuctal
74   0  0  1  1  0       257     nonmassM   InvasiveDuctal
81   1  1  0  0  0       162        massM   InvasiveDuctal
84   0  0  1  1  0       251        massM   InvasiveDuctal
102  1  1  0  0  0       240        massM   InvasiveDuctal
113  1  1  0  0  0         8     nonmassM   InvasiveDuctal
114  0  0  1  1  0         9     nonmassM   InvasiveDuctal
118  1  0  0  0  1       546        massM   InvasiveDuctal
129  0  0  1  1  0       523        massM   InvasiveDuctal
4    0  1  0  0  0       179        massM   InvasiveDuctal
17   0  0  1  0  0       188        massM   InvasiveDuctal
23   0  1  0  0  0       384        massM   InvasiveDuctal
48   0  0  0  1  0       276        massM   InvasiveDuctal
72   1  0  0  0  0       190        massM   InvasiveDuctal
96   0  0  0  0  1       628        massM   InvasiveDuctal
1    0  0  0  0  0        40        massM   InvasiveDuctal
2    0  0  0  0  0        41     nonmassM   InvasiveDuctal
6    0  0  0  0  0       229        massM   InvasiveDuctal
7    0  0  0  0  0       380        massM   InvasiveDuctal
13   0  0  0  0  0       540     nonmassM   InvasiveDuctal
14   0  0  0  0  0        89     nonmassM   InvasiveDuctal
22   0  0  0  0  0       360        massM   InvasiveDuctal
25   0  0  0  0  0       457        massM   InvasiveDuctal
26   0  0  0  0  0       487        massM   InvasiveDuctal
27   0  0  0  0  0       490        massM   InvasiveDuctal
28   0  0  0  0  0       506     nonmassM   InvasiveDuctal
32   0  0  0  0  0       122     nonmassM   InvasiveDuctal
35   0  0  0  0  0       230        massM   InvasiveDuctal
39   0  0  0  0  0       552        massM   InvasiveDuctal
40   0  0  0  0  0       578        massM   InvasiveDuctal
46   0  0  0  0  0       189        massM   InvasiveDuctal
50   0  0  0  0  0       321     nonmassM   InvasiveDuctal
54   0  0  0  0  0       621        massM   InvasiveDuctal
55   0  0  0  0  0        81        massM   InvasiveDuctal
56   0  0  0  0  0        95        massM   InvasiveDuctal
57   0  0  0  0  0       108     nonmassM   InvasiveDuctal
61   0  0  0  0  0       247        massM   InvasiveDuctal
62   0  0  0  0  0       401        massM   InvasiveDuctal
64   0  0  0  0  0       545        massM   InvasiveDuctal
65   0  0  0  0  0       555        massM   InvasiveDuctal
66   0  0  0  0  0       556        massM   InvasiveDuctal
67   0  0  0  0  0       620        massM   InvasiveDuctal
68   0  0  0  0  0       142        massM   InvasiveDuctal
69   0  0  0  0  0       143     nonmassM   InvasiveDuctal
73   0  0  0  0  0       215        massM   InvasiveDuctal
76   0  0  0  0  0       403        massM   InvasiveDuctal
78   0  0  0  0  0       494     nonmassM   InvasiveDuctal
79   0  0  0  0  0        88        massM   InvasiveDuctal
82   0  0  0  0  0       163        massM   InvasiveDuctal
85   0  0  0  0  0       260        massM   InvasiveDuctal
86   0  0  0  0  0       329     nonmassM   InvasiveDuctal
87   0  0  0  0  0       350        massM   InvasiveDuctal
89   0  0  0  0  0       541        massM   InvasiveDuctal
90   0  0  0  0  0       542        massM   InvasiveDuctal
91   0  0  0  0  0       543        massM   InvasiveDuctal
95   0  0  0  0  0       580        massM   InvasiveDuctal
97   0  0  0  0  0        67        massM   InvasiveDuctal
98   0  0  0  0  0        68        massM   InvasiveDuctal
99   0  0  0  0  0       111        massM   InvasiveDuctal
100  0  0  0  0  0       126        massM   InvasiveDuctal
101  0  0  0  0  0       136        massM   InvasiveDuctal
106  0  0  0  0  0       444        massM   InvasiveDuctal
107  0  0  0  0  0       447        massM   InvasiveDuctal
108  0  0  0  0  0       553        massM   InvasiveDuctal
111  0  0  0  0  0       589        massM   InvasiveDuctal
112  0  0  0  0  0       626        massM   InvasiveDuctal
115  0  0  0  0  0       106        massM   InvasiveDuctal
116  0  0  0  0  0       495        massM   InvasiveDuctal
121  0  0  0  0  0       238        massM   InvasiveDuctal
122  0  0  0  0  0       265        massM   InvasiveDuctal
123  0  0  0  0  0       337     nonmassM   InvasiveDuctal
126  0  0  0  0  0       441     nonmassM   InvasiveDuctal
127  0  0  0  0  0       512        massM   InvasiveDuctal
132  0  0  0  0  0       528        massM   InvasiveDuctal
134  0  0  0  0  0       530     nonmassM   InvasiveDuctal
```

```r
print("Fine explaining rules for InvasiveDuctal (IDC) with LMSIR:")
```

```
[1] "Fine explaining rules for InvasiveDuctal (IDC) with LMSIR:"
```

```r
LMSIRfeatureflag = c()
idLMSIR = c()
for (i in 1:length(eachIDCTopRules)) {
    if (length(eachIDCTopRules[[i]]) > 0) {
        rules = mypresentRules(eachIDCTopRules[[i]], colnames(X), fnames)
        LMSIRwtop5 = "LMSIR" %in% unlist(strsplit(unlist(strsplit(rules$condition, split = " & ")), 
            split = " "))
        if (LMSIRwtop5) {
            LMSIRfeatureflag = c(LMSIRfeatureflag, LMSIRwtop5)
            print(i)
            print(rules)
            idLMSIR = c(idLMSIR, i)
        } else {
            LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
        }
    } else {
        LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
    }
}
```

```
[1] 46
     len  freq                err
3868   3 0.031 0.0590000000000001
300    3 0.115              0.081
                                                                               condition pred
3868 dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low     C
300        SER(in) medium to high  & uptake skewness low  & irregularity medium to high     C
                                prob sumtempProb                           Ruletype
3868 0.4853810, 0.6962119, 0.4631950   0.5482627 dispersion & single-time-Enh & T2w
300  0.5949279, 0.7518615, 0.4569649   0.6012514         T1wdynamic & T1wmorphology
[1] 56
     len  freq                err
3868   3 0.031 0.0590000000000001
4999   3 0.059 0.0620000000000001
442    3 0.237               0.07
858    2 0.062              0.088
                                                                                condition pred
3868  dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low     C
4999                  uptake average high  & 1st post-SE s8 high  & T2w Correlation  low     C
442  irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high    NC
858                      SER(in) medium to high  & Max Radial gradient variance very low    NC
                                prob sumtempProb                           Ruletype
3868 0.4853810, 0.6962119, 0.4631950   0.5482627 dispersion & single-time-Enh & T2w
4999 0.3727659, 0.5435277, 0.4635614   0.4599517 T1wdynamic & single-time-Enh & T2w
442  0.4378788, 0.2352466, 0.2721108   0.3150787         T1wmorphology & T1wtexture
858             0.4812519, 0.5830327   0.5321423         T1wdynamic & T1wmorphology
[1] 87
     len  freq                err
68     3 0.111              0.049
3868   3 0.031 0.0590000000000001
4999   3 0.059 0.0620000000000001
300    3 0.115              0.081
                                                                                condition pred
68   Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high     C
3868  dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low     C
4999                  uptake average high  & 1st post-SE s8 high  & T2w Correlation  low     C
300         SER(in) medium to high  & uptake skewness low  & irregularity medium to high     C
                                prob sumtempProb                           Ruletype
68   0.4479334, 0.7561858, 0.6379134   0.6140108         T1wdynamic & T1wmorphology
3868 0.4853810, 0.6962119, 0.4631950   0.5482627 dispersion & single-time-Enh & T2w
4999 0.3727659, 0.5435277, 0.4635614   0.4599517 T1wdynamic & single-time-Enh & T2w
300  0.5949279, 0.7518615, 0.4569649   0.6012514         T1wdynamic & T1wmorphology
[1] 89
      len  freq                err
10957   5 0.051                  0
68      3 0.111              0.049
3868    3 0.031 0.0590000000000001
4999    3 0.059 0.0620000000000001
2791    3 0.117              0.078
2963    3 0.091               0.08
300     3 0.115              0.081
124     3 0.133              0.097
                                                                                                                                                               condition
10957 Time-to-peak(in) low  & Variance of spatial Margin Gradient very low to high  & circularity low & dispersion s12 very low to high  & dispersion s17 low to medium 
68                                                                                  Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high 
3868                                                                                 dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
4999                                                                                                 uptake average high  & 1st post-SE s8 high  & T2w Correlation  low 
2791                                                                            Time-to-peak(in) low  & circularity low & Max Radial gradient variance low to very high 
2963                                                                                      uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high 
300                                                                                        SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
124                                                                  Time-to-peak(in) low  & Rate of signal increase(rim) medium to high  & irregularity medium to high 
      pred                                                  prob sumtempProb
10957    C 0.5682319, 0.6838746, 0.4879725, 0.4530955, 0.4300886   0.5246526
68       C                       0.4479334, 0.7561858, 0.6379134   0.6140108
3868     C                       0.4853810, 0.6962119, 0.4631950   0.5482627
4999     C                       0.3727659, 0.5435277, 0.4635614   0.4599517
2791     C                       0.6633012, 0.5110174, 0.6863620   0.6202269
2963     C                       0.4627302, 0.4869596, 0.5895322    0.513074
300      C                       0.5949279, 0.7518615, 0.4569649   0.6012514
124      C                       0.7280304, 0.4252003, 0.5877582   0.5803296
                                       Ruletype
10957   T1wdynamic & T1wmorphology & dispersion
68                   T1wdynamic & T1wmorphology
3868         dispersion & single-time-Enh & T2w
4999         T1wdynamic & single-time-Enh & T2w
2791                 T1wdynamic & T1wmorphology
2963  T1wdynamic & T1wtexture & single-time-Enh
300                  T1wdynamic & T1wmorphology
124                  T1wdynamic & T1wmorphology
[1] 95
     len  freq                err
3868   3 0.031 0.0590000000000001
                                                                               condition pred
3868 dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low     C
                                prob sumtempProb                           Ruletype
3868 0.4853810, 0.6962119, 0.4631950   0.5482627 dispersion & single-time-Enh & T2w
[1] 111
      len  freq                err
10957   5 0.051                  0
3868    3 0.031 0.0590000000000001
2791    3 0.117              0.078
454     3 0.132              0.097
124     3 0.133              0.097
                                                                                                                                                               condition
10957 Time-to-peak(in) low  & Variance of spatial Margin Gradient very low to high  & circularity low & dispersion s12 very low to high  & dispersion s17 low to medium 
3868                                                                                 dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
2791                                                                            Time-to-peak(in) low  & circularity low & Max Radial gradient variance low to very high 
454                                                                               Time-to-peak(in) low  & irregularity medium to high  & 1st post-SE s12 medium to high 
124                                                                  Time-to-peak(in) low  & Rate of signal increase(rim) medium to high  & irregularity medium to high 
      pred                                                  prob sumtempProb
10957    C 0.5682319, 0.6838746, 0.4879725, 0.4530955, 0.4300886   0.5246526
3868     C                       0.4853810, 0.6962119, 0.4631950   0.5482627
2791     C                       0.6633012, 0.5110174, 0.6863620   0.6202269
454      C                       0.7569990, 0.5009352, 0.6581838    0.638706
124      C                       0.7280304, 0.4252003, 0.5877582   0.5803296
                                          Ruletype
10957      T1wdynamic & T1wmorphology & dispersion
3868            dispersion & single-time-Enh & T2w
2791                    T1wdynamic & T1wmorphology
454   T1wdynamic & T1wmorphology & single-time-Enh
124                     T1wdynamic & T1wmorphology
[1] 113
      len  freq                err
420     3 0.148              0.049
3868    3 0.031 0.0590000000000001
2963    3 0.091               0.08
3524    2 0.043              0.087
182     2 0.123              0.088
62      3  0.16               0.09
10681   2  0.08              0.091
7286    4 0.159              0.093
2589    2 0.095              0.096
                                                                                                                                                 condition
420                                                         Variance ptime1 very high  & dispersion s16 low to very high  & 3rd post-SE s6 medium to high 
3868                                                                   dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
2963                                                                        uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high 
3524                                                                        change in Variance of spatial Margin Gradient low  & dispersion s17 very high 
182                                                                                             Variance ptime1 very high  & dispersion s9 medium to high 
62                        max relative signal enhancement(rim) medium to high  & Variance ptime1 very high  & Difference variance ptime2 very low to high 
10681                                                                 Difference variance ptime1 very high  & Difference variance ptime2 very low to high 
7286  irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
2589                                                                             Variance ptime2 very high  & Difference variance ptime2 very low to high 
      pred                                       prob sumtempProb
420      C            0.4705374, 0.8140437, 0.7544592   0.6796801
3868     C            0.4853810, 0.6962119, 0.4631950   0.5482627
2963     C            0.4627302, 0.4869596, 0.5895322    0.513074
3524     C                       0.3885357, 0.5447962    0.466666
182      C                       0.4279656, 0.7763369   0.6021512
62       C            0.8237662, 0.4610768, 0.7663881   0.6837437
10681    C                       0.4832205, 0.4638078   0.4735142
7286     C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
2589     C                       0.4973411, 0.4768212   0.4870812
                                       Ruletype
420   T1wtexture & dispersion & single-time-Enh
3868         dispersion & single-time-Enh & T2w
2963  T1wdynamic & T1wtexture & single-time-Enh
3524                 T1wmorphology & dispersion
182                     T1wtexture & dispersion
62                      T1wdynamic & T1wtexture
10681                                T1wtexture
7286           T1wmorphology & T1wtexture & T2w
2589                                 T1wtexture
[1] 114
     len  freq                err
1193   3  0.03                  0
68     3 0.111              0.049
3868   3 0.031 0.0590000000000001
310    2 0.144              0.076
300    3 0.115              0.081
1063   3 0.067              0.081
3524   2 0.043              0.087
588    2 0.121               0.09
615    2  0.14              0.091
184    2 0.156              0.094
                                                                                            condition
1193 Time-to-peak(rim) medium to high  & max Radial gradient low  & T2w radial gradient variance low 
68               Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high 
3868              dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
310                                                      SER(in) medium to high  & irregularity high 
300                     SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
1063      irregularity high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
3524                   change in Variance of spatial Margin Gradient low  & dispersion s17 very high 
588                                                                SER(in) high  & irregularity high 
615                                              irregularity high  & 1st post-SE s15 medium to high 
184                                                      SER(in) medium to high  & irregularity high 
     pred                            prob sumtempProb
1193   NC 0.4662758, 0.3631702, 0.5977509   0.4757323
68      C 0.4479334, 0.7561858, 0.6379134   0.6140108
3868    C 0.4853810, 0.6962119, 0.4631950   0.5482627
310     C            0.5006988, 0.6626664   0.5816826
300     C 0.5949279, 0.7518615, 0.4569649   0.6012514
1063    C 0.5128461, 0.6452170, 0.4944342   0.5508324
3524    C            0.3885357, 0.5447962    0.466666
588     C            0.4565346, 0.6397499   0.5481423
615     C            0.4819968, 0.7471799   0.6145883
184     C            0.7579719, 0.4714506   0.6147112
                                         Ruletype
1193             T1wdynamic & T1wmorphology & T2w
68                     T1wdynamic & T1wmorphology
3868           dispersion & single-time-Enh & T2w
310                    T1wdynamic & T1wmorphology
300                    T1wdynamic & T1wmorphology
1063 T1wmorphology & T1wtexture & single-time-Enh
3524                   T1wmorphology & dispersion
588                    T1wdynamic & T1wmorphology
615               T1wmorphology & single-time-Enh
184                    T1wdynamic & T1wmorphology
[1] 126
      len  freq                err
3499    4 0.031                  0
4984    2 0.042              0.043
11396   3 0.063 0.0590000000000001
3868    3 0.031 0.0590000000000001
4999    3 0.059 0.0620000000000001
2963    3 0.091               0.08
7460    3 0.062              0.088
                                                                                                                                  condition
3499  max uptake very high  & Sum average ptime2 very low to high  & last post-SE s8 low to very high  & T2w Sum average  low to very high 
4984                                                                              max uptake very high  & last post-SE s8 low to very high 
11396                                               Washout rate(in) very high  & Sum Entropy ptime1 very high  & Entropy ptime1 very high 
3868                                                    dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
4999                                                                    uptake average high  & 1st post-SE s8 high  & T2w Correlation  low 
2963                                                         uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high 
7460                                     uptake average high  & std 3D Sharpness of lesion margin  high  & last post-SE s18 medium to high 
      pred                                       prob sumtempProb
3499     C 0.4648681, 0.4791564, 0.4145251, 0.4431299   0.4504199
4984     C                       0.5415453, 0.4910949   0.5163201
11396    C            0.5191539, 0.4959777, 0.3986381   0.4712566
3868     C            0.4853810, 0.6962119, 0.4631950   0.5482627
4999     C            0.3727659, 0.5435277, 0.4635614   0.4599517
2963     C            0.4627302, 0.4869596, 0.5895322    0.513074
7460     C            0.2956963, 0.4250158, 0.5423104   0.4210075
                                             Ruletype
3499  T1wdynamic & T1wtexture & single-time-Enh & T2w
4984                     T1wdynamic & single-time-Enh
11396                         T1wdynamic & T1wtexture
3868               dispersion & single-time-Enh & T2w
4999               T1wdynamic & single-time-Enh & T2w
2963        T1wdynamic & T1wtexture & single-time-Enh
7460     T1wdynamic & T1wmorphology & single-time-Enh
```

```r
print(summary(LMSIRfeatureflag))
```

```
   Mode   FALSE    TRUE    NA's 
logical     125       9       0 
```

```r
### display
lesion_id = 241

idx = c(1:nrow(IDC))[c(IDC$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
241       241 241_0817_5363917.vtk          0817 1970-07-04 00:00:00.000000              13171
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
241           5363917 2010-08-24 00:00:00.000000                 Unknown
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
241                          Other                     1                        0
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
241              Right                 2499 2010-09-08 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int       proc_proc_tp_int
241            Radiology                MRI Vacuum Assisted Biopsy
                                                                                                                                           proc_lesion_comments_txt
241 Investigation in Bangladesh of palpable abnormality 1\no'clock right breast with FNAB positive for malignancy. FNAB of\nprominent node negative for malignancy.
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
241             Ia                     Rapid                 Persistent      6        massM
    lesion_diagnosis
241   InvasiveDuctal
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
?<Report title=""Surgical Pathology Report""><Sessions>
Date: Sept 08, 2010

AP Report: 
 (NOTE)
 
 SURGICAL PATHOLOGY REPORT
 

 HFN: 2693587
 Encounter #: 7323410R
 Specimen #: S10-23904
 
 CLINICAL INFORMATION
 A) 14 G 4 cores. Probable FA. B) 14 G 4 cores. Known Ca by FNA and
 receptors. Receptors please; ER/PR/Her2Neu.
 
 SPECIMENS SUBMITTED
 A: Right breast, core biopsy 9 o'clock
 B: Right breast, core biopsy 1 o'clock
 
 
 
 DIAGNOSIS
 
 A. Breast, right, 9 o'clock, core biopsy:
 - FIBROADENOMA
 - PROLIFERATIVE FIBROCYSTIC CHANGES
 - NO EVIDENCE OF MALIGNANCY
 
 B. Breast, right, 1 o'clock, core biopsy:
 - INVASIVE DUCTAL CARCINOMA, MICROPAPILLARY TYPE, INTERMEDIATE
 DIFFERENTIATION
 
 
 
 
 MACROSCOPIC DESCRIPTION
 A. The specimen container is labeled with the patient name and 'right
 breast biopsies site A 9H'. The accompanying requisition matches the
 container's label.The specimen consists of 4 cores of tan hemorrhagic
 and fatty tissue, each with a diameter of 0.1 cm, ranging from 0.8 to
 1.5 cm in length. Submitted in toto in one block.
 
 B. The specimen container is labeled with the patient name and 'right
 breast biopsies site B 1H'. The accompanying requisition matches the
 container's label.The specimen consists of 4 cores of tan hemorrhagic
 and fatty tissue, each with a diameter of 0.1 cm, ranging from 0.6 to
 1.5 cm in length. Submitted in toto in one block.
 LS
 Dictated 09/08/2010
 
 
 
 
 
 Reda Saad, MD, FRCPC
 Report Electronically Signed
 2010/09/09 16:30


</Sessions></Report>
```

```r
print(mypresentRules(eachIDCTopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
10957   5 0.051                  0
322     2 0.113              0.049
420     3 0.148              0.049
68      3 0.111              0.049
9295    3 0.073               0.05
11396   3 0.063 0.0590000000000001
112     2 0.089 0.0610000000000001
4999    3 0.059 0.0620000000000001
40      2 0.132              0.069
670     3 0.079               0.07
310     2 0.144              0.076
12068   3 0.094              0.077
2791    3 0.117              0.078
2963    3 0.091               0.08
300     3 0.115              0.081
1063    3 0.067              0.081
455     2 0.128              0.085
7460    3 0.062              0.088
62      3  0.16               0.09
588     2 0.121               0.09
615     2  0.14              0.091
10681   2  0.08              0.091
85      2 0.138              0.092
7286    4 0.159              0.093
155     2 0.155              0.094
184     2 0.156              0.094
135     2 0.115              0.095
2589    2 0.095              0.096
511     2  0.13              0.097
454     3 0.132              0.097
124     3 0.133              0.097
                                                                                                                                                               condition
10957 Time-to-peak(in) low  & Variance of spatial Margin Gradient very low to high  & circularity low & dispersion s12 very low to high  & dispersion s17 low to medium 
322                                                                                                                           irregularity high  & Variance ptime1 high 
420                                                                       Variance ptime1 very high  & dispersion s16 low to very high  & 3rd post-SE s6 medium to high 
68                                                                                  Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high 
9295                                                                                                  Time-to-peak(in) low  & irregularity high  & 1st post-SE s10 high 
11396                                                                            Washout rate(in) very high  & Sum Entropy ptime1 very high  & Entropy ptime1 very high 
112                                                                                                                  irregularity very high  & Sum variance ptime1 high 
4999                                                                                                 uptake average high  & 1st post-SE s8 high  & T2w Correlation  low 
40                                                                                                                            irregularity high  & Variance ptime1 high 
670                                                                                           irregularity high  & Variance ptime2 high  & Sum Entropy ptime2 very high 
310                                                                                                                         SER(in) medium to high  & irregularity high 
12068                                                                   irregularity low to very high  & Sum variance ptime1 very high  & dispersion s18 medium to high 
2791                                                                            Time-to-peak(in) low  & circularity low & Max Radial gradient variance low to very high 
2963                                                                                      uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high 
300                                                                                        SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
1063                                                                         irregularity high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
455                                                                                                                         SER(in) medium to high  & irregularity high 
7460                                                                  uptake average high  & std 3D Sharpness of lesion margin  high  & last post-SE s18 medium to high 
62                                      max relative signal enhancement(rim) medium to high  & Variance ptime1 very high  & Difference variance ptime2 very low to high 
588                                                                                                                                   SER(in) high  & irregularity high 
615                                                                                                                 irregularity high  & 1st post-SE s15 medium to high 
10681                                                                               Difference variance ptime1 very high  & Difference variance ptime2 very low to high 
85                                                                                                           Variance ptime1 very high  & dispersion s16 medium to high 
7286                irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
155                                                                                                          Variance ptime1 very high  & dispersion s16 medium to high 
184                                                                                                                         SER(in) medium to high  & irregularity high 
135                                                                                                                         SER(in) medium to high  & irregularity high 
2589                                                                                           Variance ptime2 very high  & Difference variance ptime2 very low to high 
511                                                                                                                         SER(in) medium to high  & irregularity high 
454                                                                               Time-to-peak(in) low  & irregularity medium to high  & 1st post-SE s12 medium to high 
124                                                                  Time-to-peak(in) low  & Rate of signal increase(rim) medium to high  & irregularity medium to high 
      pred                                                  prob sumtempProb
10957    C 0.5682319, 0.6838746, 0.4879725, 0.4530955, 0.4300886   0.5246526
322      C                                  0.4675507, 0.7419192    0.604735
420      C                       0.4705374, 0.8140437, 0.7544592   0.6796801
68       C                       0.4479334, 0.7561858, 0.6379134   0.6140108
9295     C                       0.7741849, 0.4788727, 0.7091933   0.6540836
11396    C                       0.5191539, 0.4959777, 0.3986381   0.4712566
112      C                                  0.7075908, 0.4421720   0.5748814
4999     C                       0.3727659, 0.5435277, 0.4635614   0.4599517
40       C                                  0.5017755, 0.7615645     0.63167
670      C                       0.5627317, 0.5135349, 0.6720011   0.5827559
310      C                                  0.5006988, 0.6626664   0.5816826
12068    C                       0.4815758, 0.4739428, 0.5346694   0.4967293
2791     C                       0.6633012, 0.5110174, 0.6863620   0.6202269
2963     C                       0.4627302, 0.4869596, 0.5895322    0.513074
300      C                       0.5949279, 0.7518615, 0.4569649   0.6012514
1063     C                       0.5128461, 0.6452170, 0.4944342   0.5508324
455      C                                  0.7662063, 0.4502787   0.6082425
7460     C                       0.2956963, 0.4250158, 0.5423104   0.4210075
62       C                       0.8237662, 0.4610768, 0.7663881   0.6837437
588      C                                  0.4565346, 0.6397499   0.5481423
615      C                                  0.4819968, 0.7471799   0.6145883
10681    C                                  0.4832205, 0.4638078   0.4735142
85       C                                  0.4344588, 0.7756773   0.6050681
7286     C            0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
155      C                                  0.4291609, 0.7692691    0.599215
184      C                                  0.7579719, 0.4714506   0.6147112
135      C                                  0.7467186, 0.4575242   0.6021214
2589     C                                  0.4973411, 0.4768212   0.4870812
511      C                                  0.7969179, 0.5002271   0.6485725
454      C                       0.7569990, 0.5009352, 0.6581838    0.638706
124      C                       0.7280304, 0.4252003, 0.5877582   0.5803296
                                          Ruletype
10957      T1wdynamic & T1wmorphology & dispersion
322                     T1wmorphology & T1wtexture
420      T1wtexture & dispersion & single-time-Enh
68                      T1wdynamic & T1wmorphology
9295  T1wdynamic & T1wmorphology & single-time-Enh
11396                      T1wdynamic & T1wtexture
112                     T1wmorphology & T1wtexture
4999            T1wdynamic & single-time-Enh & T2w
40                      T1wmorphology & T1wtexture
670                     T1wmorphology & T1wtexture
310                     T1wdynamic & T1wmorphology
12068      T1wmorphology & T1wtexture & dispersion
2791                    T1wdynamic & T1wmorphology
2963     T1wdynamic & T1wtexture & single-time-Enh
300                     T1wdynamic & T1wmorphology
1063  T1wmorphology & T1wtexture & single-time-Enh
455                     T1wdynamic & T1wmorphology
7460  T1wdynamic & T1wmorphology & single-time-Enh
62                         T1wdynamic & T1wtexture
588                     T1wdynamic & T1wmorphology
615                T1wmorphology & single-time-Enh
10681                                   T1wtexture
85                         T1wtexture & dispersion
7286              T1wmorphology & T1wtexture & T2w
155                        T1wtexture & dispersion
184                     T1wdynamic & T1wmorphology
135                     T1wdynamic & T1wmorphology
2589                                    T1wtexture
511                     T1wdynamic & T1wmorphology
454   T1wdynamic & T1wmorphology & single-time-Enh
124                     T1wdynamic & T1wmorphology
```

```r
# print(IDC[IDC$allTestinfo.lesion_id==241,])
imgIDC <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgIDC)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
## with LMSIR
if (length(idLMSIR) > 0) {
    for (i in 1:length(idLMSIR)) {
        idx = idLMSIR[i]
        lesion_id = IDC[idx, ]$allTestinfo.lesion_id
        print(lesion_id)
    }
}
```

```
[1] 189
[1] 95
[1] 350
[1] 541
[1] 580
[1] 589
[1] 8
[1] 9
[1] 441
```

```r
idx = c(1:nrow(IDC))[c(IDC$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
441       441 441_3053_7449310.vtk          3053 1960-08-21 00:00:00.000000              31654
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
441           7449310 2013-04-07 00:00:00.000000                 Unknown
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
441                           None                     0                        1
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
441              Right                 4075 2013-04-23 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int proc_lesion_comments_txt
441            Radiology             Stereo Core Needle Biopsy                     None
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
441           None                     Rapid                    Washout      5     nonmassM
    lesion_diagnosis
441   InvasiveDuctal
```

```r
print(mypresentRules(eachIDCTopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
3499    4 0.031                  0
4984    2 0.042              0.043
11396   3 0.063 0.0590000000000001
3868    3 0.031 0.0590000000000001
4999    3 0.059 0.0620000000000001
2963    3 0.091               0.08
7460    3 0.062              0.088
                                                                                                                                  condition
3499  max uptake very high  & Sum average ptime2 very low to high  & last post-SE s8 low to very high  & T2w Sum average  low to very high 
4984                                                                              max uptake very high  & last post-SE s8 low to very high 
11396                                               Washout rate(in) very high  & Sum Entropy ptime1 very high  & Entropy ptime1 very high 
3868                                                    dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
4999                                                                    uptake average high  & 1st post-SE s8 high  & T2w Correlation  low 
2963                                                         uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high 
7460                                     uptake average high  & std 3D Sharpness of lesion margin  high  & last post-SE s18 medium to high 
      pred                                       prob sumtempProb
3499     C 0.4648681, 0.4791564, 0.4145251, 0.4431299   0.4504199
4984     C                       0.5415453, 0.4910949   0.5163201
11396    C            0.5191539, 0.4959777, 0.3986381   0.4712566
3868     C            0.4853810, 0.6962119, 0.4631950   0.5482627
4999     C            0.3727659, 0.5435277, 0.4635614   0.4599517
2963     C            0.4627302, 0.4869596, 0.5895322    0.513074
7460     C            0.2956963, 0.4250158, 0.5423104   0.4210075
                                             Ruletype
3499  T1wdynamic & T1wtexture & single-time-Enh & T2w
4984                     T1wdynamic & single-time-Enh
11396                         T1wdynamic & T1wtexture
3868               dispersion & single-time-Enh & T2w
4999               T1wdynamic & single-time-Enh & T2w
2963        T1wdynamic & T1wtexture & single-time-Enh
7460     T1wdynamic & T1wmorphology & single-time-Enh
```

![](Z:/Cristina/Section2/papernew_notes/images/241_0817_5363917_massM_Hyperintense.png)


# 2) Present rules collect top for InsituDuctal ISDC

```r
eachISDCTopRules = list()
allISDCtoprules = c()
for (idx in 1:nrow(ISDC)) {
    ############## first analysis
    X = ISDC[idx, 2:ncol(imgT2pLMSIRtest)]
    y = ISDC[idx, "lesion_label"]
    
    rulesoutputISDC = mynewapplyLearnerxrules(topRulesISDC, X, y, minerr = 0.1, minfrq = 10/627, 
        classes, gbmModel)
    topRulesISDC = rulesoutputISDC[[1]]
    
    if (length(rulesoutputISDC) > 1) {
        eachISDCTopRules[[idx]] = rulesoutputISDC[[2]]
        allISDCtoprules = rbind(allISDCtoprules, rulesoutputISDC[[3]])
    } else {
        eachISDCTopRules[[idx]] = list()
        allISDCtoprules = rbind(allISDCtoprules, 1:nrow(topRulesISDC) * 0)
    }
}
```

```
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 10 out of 50"
```

```r
# take the occurance of each explaining rule (columns) per case(row) + lesion identifies
df = data.frame(cbind(allISDCtoprules, as.character(ISDC$allTestinfo.lesion_id), as.character(ISDC$allTestinfo.lesion_label), 
    as.character(ISDC$allTestinfo.lesion_diagnosis)), stringsAsFactors = FALSE)
colnames(df) <- c(as.character(1:nrow(topRulesISDC)), "lesion_id", "lesion_label", "lesion_diagnosis")
cols = 1:nrow(topRulesISDC)
df[, cols] = apply(df[, cols], 2, function(x) as.numeric(x))


cat(paste0("\n============== Top five explaining rules for InsituDuctal (ISDC): n = ", nrow(ISDC), 
    "\n"))
```

```

============== Top five explaining rules for InsituDuctal (ISDC): n = 80
```

```r
topISDC = sort(colSums(df[, 1:nrow(topRulesISDC)]), decreasing = TRUE)
print(topISDC[1:5])
```

```
42 50 21 48 41 
27 21 20 19 18 
```

```r
freqtopISDC = topISDC[1:5]
print(freqtopISDC/nrow(ISDC))
```

```
    42     50     21     48     41 
0.3375 0.2625 0.2500 0.2375 0.2250 
```

```r
rulesTopISDC = as.numeric(names(topISDC[1:5]))
preserulesTopISDC = mypresentRules(topRulesISDC[rulesTopISDC, ], colnames(X), fnames)
rownames(preserulesTopISDC) <- NULL
print(preserulesTopISDC)
```

```
  len  freq                err
1   4 0.159              0.093
2   2 0.152              0.098
3   2 0.135 0.0679999999999999
4   2 0.095              0.096
5   3 0.101              0.091
                                                                                                                                             condition
1 irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
2                                                                                                           irregularity high  & Variance ptime1 high 
3                                                                                                       irregularity high  & Sum variance ptime1 high 
4                                                                            Variance ptime2 very high  & Difference variance ptime2 very low to high 
5                                                                       SER(in) medium to high  & irregularity high  & Difference entropy ptime3 high 
  pred                                       prob sumtempProb
1    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
2    C                       0.4852856, 0.7650575   0.6251715
3    C                       0.5087730, 0.8044078   0.6565904
4    C                       0.4973411, 0.4768212   0.4870812
5    C            0.4542278, 0.3886270, 0.5466488   0.4631679
                                 Ruletype
1        T1wmorphology & T1wtexture & T2w
2              T1wmorphology & T1wtexture
3              T1wmorphology & T1wtexture
4                              T1wtexture
5 T1wdynamic & T1wmorphology & T1wtexture
```

```r
# display a case that meets them all
casesISDC = df[, c(rulesTopISDC, nrow(topRulesISDC) + 1, nrow(topRulesISDC) + 2, nrow(topRulesISDC) + 
    3)]
print("Sorted cases meeting top five explaining rules for InsituDuctal (ISDC):")
```

```
[1] "Sorted cases meeting top five explaining rules for InsituDuctal (ISDC):"
```

```r
topcasesISDC = casesISDC[sort(rowSums(casesISDC[, 1:5]), index = TRUE, decreasing = TRUE)$ix, ]
print(topcasesISDC)
```

```
   42 50 21 48 41 lesion_id lesion_label lesion_diagnosis
9   1  1  1  1  1       294        massM     InsituDuctal
10  1  1  1  1  1       295     nonmassM     InsituDuctal
31  1  1  1  1  1       395        massM     InsituDuctal
32  1  1  1  1  1       612        massM     InsituDuctal
45  1  1  1  1  1       246        massM     InsituDuctal
59  1  1  1  1  1       564        massM     InsituDuctal
7   1  1  1  1  0       466     nonmassM     InsituDuctal
16  1  1  1  0  1       571        massM     InsituDuctal
41  1  1  1  0  1       551        massM     InsituDuctal
57  1  1  1  0  1       562        massM     InsituDuctal
76  1  1  1  1  0       595        massM     InsituDuctal
78  1  1  1  1  0       608        massM     InsituDuctal
79  1  1  1  0  1       609        massM     InsituDuctal
80  1  1  1  0  1       610     nonmassM     InsituDuctal
13  0  1  1  0  1       359     nonmassM     InsituDuctal
38  0  1  1  0  1       454     nonmassM     InsituDuctal
40  0  1  1  0  1       550     nonmassM     InsituDuctal
46  1  0  0  1  1       340        massM     InsituDuctal
69  0  1  1  0  1       402        massM     InsituDuctal
72  0  1  1  1  0       581     nonmassM     InsituDuctal
3   1  0  0  1  0       232     nonmassM     InsituDuctal
6   0  1  1  0  0       416     nonmassM     InsituDuctal
18  1  0  0  1  0       573        massM     InsituDuctal
19  1  0  0  1  0       574        massM     InsituDuctal
20  1  0  0  1  0       575        massM     InsituDuctal
37  1  0  0  1  0       382        massM     InsituDuctal
58  1  0  0  1  0       563        massM     InsituDuctal
60  1  0  0  1  0       565        massM     InsituDuctal
70  1  0  0  1  0       459        massM     InsituDuctal
5   0  1  0  0  0       414     nonmassM     InsituDuctal
11  1  0  0  0  0       296     nonmassM     InsituDuctal
17  1  0  0  0  0       572        massM     InsituDuctal
25  0  0  0  0  1       198     nonmassM     InsituDuctal
64  1  0  0  0  0       170        massM     InsituDuctal
68  0  0  0  0  1       376     nonmassM     InsituDuctal
71  1  0  0  0  0       460        massM     InsituDuctal
1   0  0  0  0  0         1     nonmassM     InsituDuctal
2   0  0  0  0  0       134     nonmassM     InsituDuctal
4   0  0  0  0  0       285        massM     InsituDuctal
8   0  0  0  0  0       250     nonmassM     InsituDuctal
12  0  0  0  0  0       300     nonmassM     InsituDuctal
14  0  0  0  0  0       374     nonmassM     InsituDuctal
15  0  0  0  0  0       458     nonmassM     InsituDuctal
21  0  0  0  0  0        33        massM     InsituDuctal
22  0  0  0  0  0        34     nonmassM     InsituDuctal
23  0  0  0  0  0        83     nonmassM     InsituDuctal
24  0  0  0  0  0       133     nonmassM     InsituDuctal
26  0  0  0  0  0       318        massM     InsituDuctal
27  0  0  0  0  0       118     nonmassM     InsituDuctal
28  0  0  0  0  0       138     nonmassM     InsituDuctal
29  0  0  0  0  0       195     nonmassM     InsituDuctal
30  0  0  0  0  0       316     nonmassM     InsituDuctal
33  0  0  0  0  0        69     nonmassM     InsituDuctal
34  0  0  0  0  0       184        massM     InsituDuctal
35  0  0  0  0  0       275        massM     InsituDuctal
36  0  0  0  0  0       381     nonmassM     InsituDuctal
39  0  0  0  0  0       549        massM     InsituDuctal
42  0  0  0  0  0       131        massM     InsituDuctal
43  0  0  0  0  0       132     nonmassM     InsituDuctal
44  0  0  0  0  0       173        massM     InsituDuctal
47  0  0  0  0  0       365        massM     InsituDuctal
48  0  0  0  0  0       531     nonmassM     InsituDuctal
49  0  0  0  0  0       113     nonmassM     InsituDuctal
50  0  0  0  0  0       164     nonmassM     InsituDuctal
51  0  0  0  0  0       483     nonmassM     InsituDuctal
52  0  0  0  0  0       606     nonmassM     InsituDuctal
53  0  0  0  0  0       129        massM     InsituDuctal
54  0  0  0  0  0       218        massM     InsituDuctal
55  0  0  0  0  0       234        massM     InsituDuctal
56  0  0  0  0  0       277        massM     InsituDuctal
61  0  0  0  0  0         5     nonmassM     InsituDuctal
62  0  0  0  0  0        30        massM     InsituDuctal
63  0  0  0  0  0       120        massM     InsituDuctal
65  0  0  0  0  0       268     nonmassM     InsituDuctal
66  0  0  0  0  0       298     nonmassM     InsituDuctal
67  0  0  0  0  0       375        massM     InsituDuctal
73  0  0  0  0  0        37     nonmassM     InsituDuctal
74  0  0  0  0  0       140     nonmassM     InsituDuctal
75  0  0  0  0  0       464     nonmassM     InsituDuctal
77  0  0  0  0  0       596        massM     InsituDuctal
```

```r
print("Fine explaining rules for InsituDuctal (ISDC) with LMSIR:")
```

```
[1] "Fine explaining rules for InsituDuctal (ISDC) with LMSIR:"
```

```r
LMSIRfeatureflag = c()
idLMSIR = c()
for (i in 1:length(eachISDCTopRules)) {
    if (length(eachISDCTopRules[[i]]) > 0) {
        rules = mypresentRules(eachISDCTopRules[[i]], colnames(X), fnames)
        LMSIRwtop5 = "LMSIR" %in% unlist(strsplit(unlist(strsplit(rules$condition, split = " & ")), 
            split = " "))
        if (LMSIRwtop5) {
            LMSIRfeatureflag = c(LMSIRfeatureflag, LMSIRwtop5)
            print(i)
            print(rules)
            idLMSIR = c(idLMSIR, i)
        } else {
            LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
        }
    } else {
        LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
    }
}
```

```
[1] 2
      len  freq                err
4715    4 0.027                  0
6488    3 0.102              0.036
13591   4  0.09 0.0610000000000001
                                                                                                                                           condition
4715                        2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
6488                                                                Time-to-peak(in) medium to high  & uptake skewness low  & irregularity very low 
13591 enhancement-variance decreasing rate(rim) very low to high  & irregularity low  & Sum average ptime1 high  & 2nd post-SE s14 very low to high 
      pred                                       prob sumtempProb
4715     C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264
6488    NC            0.4742468, 0.3732442, 0.4938563   0.4471158
13591   NC 0.4442097, 0.5323714, 0.4614969, 0.5905932   0.5071678
                                                       Ruletype
4715                                      single-time-Enh & T2w
6488                                 T1wdynamic & T1wmorphology
13591 T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
[1] 27
     len  freq   err
4715   4 0.027     0
593    2 0.248 0.096
                                                                                                                    condition
4715 2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
593                                                                        SER(in) low  & max Radial gradient medium to high 
     pred                                       prob sumtempProb                   Ruletype
4715    C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264      single-time-Enh & T2w
593    NC                       0.5088292, 0.2628051   0.3858172 T1wdynamic & T1wmorphology
[1] 33
      len  freq                err
4715    4 0.027                  0
13591   4  0.09 0.0610000000000001
442     3 0.237               0.07
478     3 0.156              0.093
                                                                                                                                           condition
4715                        2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
13591 enhancement-variance decreasing rate(rim) very low to high  & irregularity low  & Sum average ptime1 high  & 2nd post-SE s14 very low to high 
442                                                             irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
478                         irregularity low  & std 3D Sharpness of lesion margin  medium to high  & Inverse difference moment ptime3 low to medium 
      pred                                       prob sumtempProb
4715     C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264
13591   NC 0.4442097, 0.5323714, 0.4614969, 0.5905932   0.5071678
442     NC            0.4378788, 0.2352466, 0.2721108   0.3150787
478     NC            0.5702187, 0.4669343, 0.3339196   0.4570242
                                                       Ruletype
4715                                      single-time-Enh & T2w
13591 T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
442                                  T1wmorphology & T1wtexture
478                                  T1wmorphology & T1wtexture
[1] 52
     len  freq                err
4715   4 0.027                  0
9188   5 0.038              0.048
715    3  0.08 0.0679999999999999
442    3 0.237               0.07
2980   3 0.076              0.071
4592   2 0.091              0.082
2203   3 0.093              0.098
                                                                                                                                                                   condition
4715                                                2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
9188 max uptake medium to high  & Entropy ptime2 low to very high  & 1st post-SE s16 low to medium  & 3rd post-SE s10 low to very high  & T2w Difference variance  very low 
715                                                                            uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
442                                                                                     irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
2980                                                                                                    uptake average high  & 1st post-SE s8 high  & last post-SE s11 high 
4592                                                                                                                             uptake average high  & 2nd post-SE s2 high 
2203                                                                           uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
     pred                                                  prob sumtempProb
4715    C            0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264
9188   NC 0.4555243, 0.4955115, 0.5378557, 0.4726280, 0.4208891   0.4764817
715     C                       0.5436611, 0.6397459, 0.5003391   0.5612487
442    NC                       0.4378788, 0.2352466, 0.2721108   0.3150787
2980    C                       0.5757935, 0.5294741, 0.6741425   0.5931367
4592    C                                  0.4996731, 0.5976669     0.54867
2203    C                       0.5260448, 0.6336959, 0.5132996   0.5576801
                                            Ruletype
4715                           single-time-Enh & T2w
9188 T1wdynamic & T1wtexture & single-time-Enh & T2w
715        T1wdynamic & T1wtexture & single-time-Enh
442                       T1wmorphology & T1wtexture
2980                    T1wdynamic & single-time-Enh
4592                    T1wdynamic & single-time-Enh
2203       T1wdynamic & T1wtexture & single-time-Enh
[1] 73
     len  freq err
4715   4 0.027   0
                                                                                                                    condition
4715 2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
     pred                                       prob sumtempProb              Ruletype
4715    C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264 single-time-Enh & T2w
[1] 74
     len  freq err
4715   4 0.027   0
                                                                                                                    condition
4715 2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
     pred                                       prob sumtempProb              Ruletype
4715    C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264 single-time-Enh & T2w
```

```r
print(summary(LMSIRfeatureflag))
```

```
   Mode   FALSE    TRUE    NA's 
logical      74       6       0 
```

```r
### display
lesion_id = 595
idx = c(1:nrow(ISDC))[c(ISDC$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
595       595 595_7018_6803089.vtk          7018 1956-11-07 00:00:00.000000              17969
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
595           6803089 2011-05-24 00:00:00.000000               Malignant
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
595                           None                     1                        0
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
595              Right                 3315 2012-04-26 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int
595            Radiology                 US Core Needle Biopsy
                                                           proc_lesion_comments_txt
595 CLINICAL INDICATION:  Suspicious segmental microcalcifications on the mammogram
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
595           None                       N/A                    Plateau      4        massM
    lesion_diagnosis
595     InsituDuctal
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
ULTRASOUND GUIDED CORE BIOPSY RIGHT BREAST

CLINICAL HISTORY: suspicious for recurrence medial right breast 3 o'clock 2 cm from nipple

PROCEDURE:

Targeted ultrasound was obtained and the previously described abnormality was identified.

Following informed consent and using sterile technique and local anaesthesia (1% xylocaine), core biopsies were obtained with ultrasound guidance using a ML approach. A skin incision was made and cores were obtained with a 14 gauge tru-cut needle in
an automated gun. Through the same incision several passes were made in order to sample different areas. 4 cores were obtained.

There were no immediate complications. The patient was given a sheet of post biopsy instructions.

The patient was requested to contact her referring physician for the results.

CORE PATHOLOGY: highly suspicious for ductal carcinoma in situ. This result is concordant.

MRI is recommended for further assessment.
_____________

This report was electronically signed by WRIGHT, BARBARA, Staff Radiologist on 2012/05/03 at 17:26






--- Addendum ---
ADDENDUM:

This case was discussed at breast multidisciplinary rounds May 31, 2012.  Repeat stereotactic guided biopsy with vacuum assistance has been recommended to try to obtain a larger sample of tissue and confirm a diagnosis of DCIS definitively prior to
planning surgical treatment.


This addendum was electronically signed by HACK, KALESHA, Staff Radiologist on 2012/06/01 at 11:28

--- Addendum end ---
```

```r
print(mypresentRules(eachISDCTopRules[[idx]], colnames(X), fnames))
```

```
     len  freq                err
5728   3  0.05              0.037
8280   4 0.035              0.053
90     2 0.135 0.0679999999999999
670    3 0.079               0.07
6675   3 0.087              0.085
4150   5 0.104              0.089
7286   4 0.159              0.093
135    2 0.115              0.095
2589   2 0.095              0.096
315    2 0.152              0.098
                                                                                                                                                                      condition
5728                                                                                          uptake variance very high  & last post-SE s6 low to very high  & T2w Energy  low 
8280                                                          Rate of signal decrease(in) high & uptake variance very high  & circularity low & last post-SE s8 medium to high 
90                                                                                                                               irregularity high  & Sum variance ptime1 high 
670                                                                                                  irregularity high  & Variance ptime2 high  & Sum Entropy ptime2 very high 
6675                                                               Sum variance ptime1 high  & Inverse difference moment ptime4 medium to high  & dispersion s2 medium to high 
4150 uptake average high  & Max Radial gradient variance low to very high  & Energy ptime4 very low to high  & dispersion s19 medium to high  & 1st post-SE s12 medium to high 
7286                       irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
135                                                                                                                                SER(in) medium to high  & irregularity high 
2589                                                                                                  Variance ptime2 very high  & Difference variance ptime2 very low to high 
315                                                                                                                                  irregularity high  & Variance ptime1 high 
     pred                                                  prob sumtempProb
5728    C                       0.7332721, 0.4739669, 0.5252209   0.5774866
8280    C            0.4768186, 0.4505621, 0.7410536, 0.6456018    0.578509
90      C                                  0.5087730, 0.8044078   0.6565904
670     C                       0.5627317, 0.5135349, 0.6720011   0.5827559
6675    C                       0.3978654, 0.7406423, 0.5150453   0.5511844
4150    C 0.5303059, 0.5125338, 0.5635975, 0.6702745, 0.6316332    0.581669
7286    C            0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
135     C                                  0.7467186, 0.4575242   0.6021214
2589    C                                  0.4973411, 0.4768212   0.4870812
315     C                                  0.4852856, 0.7650575   0.6251715
                                                                   Ruletype
5728                                     T1wdynamic & single-time-Enh & T2w
8280                           T1wdynamic & T1wmorphology & single-time-Enh
90                                               T1wmorphology & T1wtexture
670                                              T1wmorphology & T1wtexture
6675                                                T1wtexture & dispersion
4150 T1wdynamic & T1wmorphology & T1wtexture & dispersion & single-time-Enh
7286                                       T1wmorphology & T1wtexture & T2w
135                                              T1wdynamic & T1wmorphology
2589                                                             T1wtexture
315                                              T1wmorphology & T1wtexture
```

```r
# print(ISDC[ISDC$allTestinfo.lesion_id==241,])
imgISDC <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgISDC)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
## with LMSIR
if (length(idLMSIR) > 0) {
    for (i in 1:length(idLMSIR)) {
        idx = idLMSIR[i]
        lesion_id = ISDC[idx, ]$allTestinfo.lesion_id
        print(lesion_id)
    }
}
```

```
[1] 134
[1] 118
[1] 69
[1] 606
[1] 37
[1] 140
```

```r
idx = c(1:nrow(ISDC))[c(ISDC$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
140       140 140_0684_5266209.vtk          0684 1956-10-14 00:00:00.000000              12766
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
140           5266209 2010-07-24 00:00:00.000000               Malignant
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
140                      High Risk                     0                        1
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
140              Right                  173 2010-08-05 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int
140            Radiology                 US Core Needle Biopsy
                                                                                            proc_lesion_comments_txt
140 Life time risk > 25%. Left lumpectomy for ADH in\n2004. Family history of breast cancer. LMP: Early June 2010.\n
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
140           None                       N/A                        N/A      4     nonmassM
    lesion_diagnosis
140     InsituDuctal
```

```r
print(mypresentRules(eachISDCTopRules[[idx]], colnames(X), fnames))
```

```
     len  freq err
4715   4 0.027   0
                                                                                                                    condition
4715 2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
     pred                                       prob sumtempProb              Ruletype
4715    C 0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264 single-time-Enh & T2w
```

![](Z:/Cristina/Section2/papernew_notes/images/595_7018_6803089_massM_None.png)

# 3) Present rules collect top for FIBROADENOMA FA

```r
eachFATopRules = list()
allFAtoprules = c()
for (idx in 1:nrow(FA)) {
    ############## first analysis
    X = FA[idx, 2:ncol(imgT2pLMSIRtest)]
    y = FA[idx, "lesion_label"]
    
    rulesoutputFA = mynewapplyLearnerxrules(topRulesFA, X, y, minerr = 0.1, minfrq = 10/627, classes, 
        gbmModel)
    topRulesFA = rulesoutputFA[[1]]
    
    if (length(rulesoutputFA) > 1) {
        eachFATopRules[[idx]] = rulesoutputFA[[2]]
        allFAtoprules = rbind(allFAtoprules, rulesoutputFA[[3]])
    } else {
        eachFATopRules[[idx]] = list()
        allFAtoprules = rbind(allFAtoprules, 1:nrow(topRulesFA) * 0)
    }
}
```

```
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 21 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 23 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 18 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 24 out of 50"
[1] "test complies with # rules: 27 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 18 out of 50"
[1] "test complies with # rules: 26 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 6 out of 50"
```

```r
# take the occurance of each explaining rule (columns) per case(row) + lesion identifies
df = data.frame(cbind(allFAtoprules, as.character(FA$allTestinfo.lesion_id), as.character(FA$allTestinfo.lesion_label), 
    as.character(FA$allTestinfo.lesion_diagnosis)), stringsAsFactors = FALSE)
colnames(df) <- c(as.character(1:nrow(topRulesFA)), "lesion_id", "lesion_label", "lesion_diagnosis")
cols = 1:nrow(topRulesFA)
df[, cols] = apply(df[, cols], 2, function(x) as.numeric(x))


cat(paste0("\n============== Top five explaining rules for FIBROADENOMA (FA): n = ", nrow(FA) + 
    1, "\n"))
```

```

============== Top five explaining rules for FIBROADENOMA (FA): n = 70
```

```r
topFA = sort(colSums(df[, 1:nrow(topRulesFA)]), decreasing = TRUE)
print(topFA[1:5])
```

```
28 35 26 30 32 
32 32 28 28 28 
```

```r
freqtopFA = topFA[1:5]
totalFA = nrow(FA) + 1
print(freqtopFA/totalFA)
```

```
       28        35        26        30        32 
0.4571429 0.4571429 0.4000000 0.4000000 0.4000000 
```

```r
rulesTopFA = as.numeric(names(topFA[1:5]))
preserulesTopFA = mypresentRules(topRulesFA[rulesTopFA, ], colnames(X), fnames)
rownames(preserulesTopFA) <- NULL
print(preserulesTopFA)
```

```
  len  freq   err
1   2 0.228 0.072
2   2 0.246 0.081
3   3 0.237  0.07
4   2   0.2 0.073
5   2 0.193 0.077
                                                                             condition pred
1                                           SER(in) low to medium  & irregularity low    NC
2                                           SER(in) low to medium  & irregularity low    NC
3 irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high    NC
4                                     irregularity low  & Energy ptime4 low to medium    NC
5                                     irregularity low  & Energy ptime4 low to medium    NC
                             prob sumtempProb                   Ruletype
1            0.4405761, 0.3043518    0.372464 T1wdynamic & T1wmorphology
2            0.5053493, 0.3706932   0.4380213 T1wdynamic & T1wmorphology
3 0.4378788, 0.2352466, 0.2721108   0.3150787 T1wmorphology & T1wtexture
4            0.4867494, 0.2974405    0.392095 T1wmorphology & T1wtexture
5            0.4991863, 0.2996197    0.399403 T1wmorphology & T1wtexture
```

```r
# display a case that meets them all
casesFA = df[, c(rulesTopFA, nrow(topRulesFA) + 1, nrow(topRulesFA) + 2, nrow(topRulesFA) + 3)]
print("Sorted cases meeting top five explaining rules for FIBROADENOMA (FA):")
```

```
[1] "Sorted cases meeting top five explaining rules for FIBROADENOMA (FA):"
```

```r
topcasesFA = casesFA[sort(rowSums(casesFA[, 1:5]), index = TRUE, decreasing = TRUE)$ix, ]
print(topcasesFA)
```

```
   28 35 26 30 32 lesion_id lesion_label lesion_diagnosis
1   1  1  1  1  1       104        massB     FIBROADENOMA
9   1  1  1  1  1       112        massB     FIBROADENOMA
14  1  1  1  1  1       436        massB     FIBROADENOMA
16  1  1  1  1  1       592        massB     FIBROADENOMA
17  1  1  1  1  1       604        massB     FIBROADENOMA
23  1  1  1  1  1       326        massB     FIBROADENOMA
24  1  1  1  1  1       408        massB     FIBROADENOMA
28  1  1  1  1  1       377     nonmassB     FIBROADENOMA
35  1  1  1  1  1       425        massB     FIBROADENOMA
37  1  1  1  1  1        73        massB     FIBROADENOMA
45  1  1  1  1  1       366        massB     FIBROADENOMA
46  1  1  1  1  1        78        massB     FIBROADENOMA
52  1  1  1  1  1       272        massB     FIBROADENOMA
53  1  1  1  1  1       618        massB     FIBROADENOMA
54  1  1  1  1  1        49        massB     FIBROADENOMA
63  1  1  1  1  1       151        massB     FIBROADENOMA
64  1  1  1  1  1       310        massB     FIBROADENOMA
30  1  1  0  1  1       502        massB     FIBROADENOMA
33  1  1  0  1  1       584        massB     FIBROADENOMA
38  1  1  0  1  1       100        massB     FIBROADENOMA
47  1  1  0  1  1       167     nonmassB     FIBROADENOMA
5   1  1  1  0  0       535        massB     FIBROADENOMA
7   0  0  1  1  1        55        massB     FIBROADENOMA
19  0  0  1  1  1       616        massB     FIBROADENOMA
20  1  1  1  0  0       127        massB     FIBROADENOMA
22  1  1  1  0  0       325        massB     FIBROADENOMA
26  1  1  1  0  0       601        massB     FIBROADENOMA
32  1  1  1  0  0       583        massB     FIBROADENOMA
55  0  0  1  1  1        60        massB     FIBROADENOMA
56  0  0  1  1  1        61        massB     FIBROADENOMA
60  1  1  1  0  0       439        massB     FIBROADENOMA
3   0  0  0  1  1       378        massB     FIBROADENOMA
31  0  0  0  1  1       517        massB     FIBROADENOMA
36  1  1  0  0  0        36     nonmassB     FIBROADENOMA
40  1  1  0  0  0       130        massB     FIBROADENOMA
48  0  0  0  1  1       208        massB     FIBROADENOMA
62  1  1  0  0  0        53        massB     FIBROADENOMA
65  1  1  0  0  0       314     nonmassB     FIBROADENOMA
67  1  1  0  0  0       467        massB     FIBROADENOMA
49  0  0  1  0  0        18     nonmassB     FIBROADENOMA
2   0  0  0  0  0       309        massB     FIBROADENOMA
4   0  0  0  0  0       534        massB     FIBROADENOMA
6   0  0  0  0  0       536        massB     FIBROADENOMA
8   0  0  0  0  0        90        massB     FIBROADENOMA
10  0  0  0  0  0       323     nonmassB     FIBROADENOMA
11  0  0  0  0  0       327        massB     FIBROADENOMA
12  0  0  0  0  0       373        massB     FIBROADENOMA
13  0  0  0  0  0       435        massB     FIBROADENOMA
15  0  0  0  0  0       437        massB     FIBROADENOMA
18  0  0  0  0  0       615        massB     FIBROADENOMA
21  0  0  0  0  0       304        massB     FIBROADENOMA
25  0  0  0  0  0       597        massB     FIBROADENOMA
27  0  0  0  0  0       635     nonmassB     FIBROADENOMA
29  0  0  0  0  0       501        massB     FIBROADENOMA
34  0  0  0  0  0       345        massB     FIBROADENOMA
39  0  0  0  0  0       125        massB     FIBROADENOMA
41  0  0  0  0  0       174        massB     FIBROADENOMA
42  0  0  0  0  0       176        massB     FIBROADENOMA
43  0  0  0  0  0       342        massB     FIBROADENOMA
44  0  0  0  0  0       343     nonmassB     FIBROADENOMA
50  0  0  0  0  0       224        massB     FIBROADENOMA
51  0  0  0  0  0       271        massB     FIBROADENOMA
57  0  0  0  0  0       157        massB     FIBROADENOMA
58  0  0  0  0  0       349     nonmassB     FIBROADENOMA
59  0  0  0  0  0       431        massB     FIBROADENOMA
61  0  0  0  0  0       496        massB     FIBROADENOMA
66  0  0  0  0  0       336        massB     FIBROADENOMA
68  0  0  0  0  0       522     nonmassB     FIBROADENOMA
69  0  0  0  0  0       629        massB     FIBROADENOMA
```

```r
print("Fine explaining rules for FIBROADENOMA (FA) with LMSIR:")
```

```
[1] "Fine explaining rules for FIBROADENOMA (FA) with LMSIR:"
```

```r
LMSIRfeatureflag = c()
idLMSIR = c()
for (i in 1:length(eachFATopRules)) {
    if (length(eachFATopRules[[i]]) > 0) {
        rules = mypresentRules(eachFATopRules[[i]], colnames(X), fnames)
        LMSIRwtop5 = "LMSIR" %in% unlist(strsplit(unlist(strsplit(rules$condition, split = " & ")), 
            split = " "))
        if (LMSIRwtop5) {
            LMSIRfeatureflag = c(LMSIRfeatureflag, LMSIRwtop5)
            print(i)
            print(rules)
            idLMSIR = c(idLMSIR, i)
        } else {
            LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
        }
    } else {
        LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
    }
}
print(summary(LMSIRfeatureflag))
```

```
   Mode   FALSE    NA's 
logical      69       0 
```

```r
### display
lesion_id = 436
idx = c(1:nrow(FA))[c(FA$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
436       436 436_3046_7289130.vtk          3046 1970-11-16 00:00:00.000000              29427
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
436           7289130 2012-12-22 00:00:00.000000     Benign by pathology
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
436                           None                     1                        0
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
436              Right                 2691 2013-01-22 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int proc_lesion_comments_txt
436            Radiology                 US Core Needle Biopsy                     None
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
436           None                       N/A                        N/A      4        massB
    lesion_diagnosis
436     FIBROADENOMA
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
AP Report: 
     (NOTE)
     
     SURGICAL  PATHOLOGY  REPORT
     

     Encounter #: 2695007R
     Submitting Physician: MESSNER,SANDRA
     Specimen #: S13-1891
     
     CLINICAL INFORMATION
     MRI lesion with slight enhancement, hypoechoic on ultrasound - 1cm.
     ?Fibroadenoma - family history of breast cancer.
     
     SPECIMENS SUBMITTED
     Right breast 1 oclock x4 cores
     
     
     
     DIAGNOSIS
     Right breast 1 oclock x4 cores: FIBROADENOMA, SEE NOTE.
     
     Note: Sections show benign fibroepithelial lesion with mildly
     increased stromal cellularity, no cytologic atypia and absent mitoses,
     most consistent with a fibroadenoma in this material.
     
     
     
     
     MACROSCOPIC DESCRIPTION
     The specimen container is labeled with the patient name and 'right
     breast core biopsy'. The accompanying requisition matches the
     container's label. The specimen consists of 4 cores of tan and fatty
     tissue, each with a diameter of 0.1 cm, ranging from 0.6 to 1.3 cm in
     length. Submitted in toto in one block
     JP
     Dictated 01/24/2013
     
     
     
     
     
     Elzbieta Slodkowska, MD
     Report Electronically Signed
     2013/01/25 17:00
```

```r
print(mypresentRules(eachFATopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
6093    6 0.031                  0
12940   4 0.049                  0
519     2 0.176              0.053
1581    3 0.096 0.0570000000000001
459     2 0.064 0.0570000000000001
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
331     3 0.106               0.07
442     3 0.237               0.07
3987    5 0.078              0.071
411     2 0.228              0.072
766     2   0.2              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
183     2 0.232              0.079
336     2 0.246              0.081
11043   3 0.111              0.082
1128    2 0.215              0.085
1538    4 0.085              0.087
2330    3 0.216              0.093
722     4 0.136              0.095
                                                                                                                                                                                                                                               condition
6093  SER(in) low to medium  & Variance of spatial Margin Gradient low to very high  & std 3D Sharpness of lesion margin  low to medium  & Inverse difference moment ptime4 medium to high  & 2nd post-SE s12 very low  & 2nd post-SE s17 low to medium 
12940                                                                                                                    Uptake(in) curve amplitude medium to high  & SER(in) low to medium  & Variance ptime3 low to medium  & dispersion s11 very low 
519                                                                                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1581                                                                                                                                                          Difference variance ptime1 low  & dispersion s16 low to medium  & T2w margin gradient low 
459                                                                                                                                                                                               irregularity very low  & dispersion s9 medium to high 
1541                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
331                                                                                                                                                                           SER(in) medium to high  & irregularity low  & Energy ptime4 low to medium 
442                                                                                                                                                                 irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
3987                                                              Inverse difference moment ptime1 very low to high  & dispersion s5 low to very high  & dispersion s7 low to very high  & 1st post-SE s19 low to very high  & 3rd post-SE s12 very low 
411                                                                                                                                                                                                           SER(in) low to medium  & irregularity low 
766                                                                                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
515                                                                                                                       SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
183                                                                                                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
336                                                                                                                                                                                                           SER(in) low to medium  & irregularity low 
11043                                                                                                                                                         Sum Entropy ptime1 low to medium  & 1st post-SE s8 low to medium  & last post-SE s17 high 
1128                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1538                                                                                                             Uptake(in) curve amplitude medium to high  & Sum Entropy ptime4 low to medium  & dispersion s17 low to very high  & 3rd post-SE s4 low 
2330                                                                                                                                                                  irregularity low  & Energy ptime3 low to medium  & dispersion s5 low to very high 
722                                                                                                                               uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
      pred                                                             prob sumtempProb
6093    NC 0.6855391, 0.5934397, 0.5476637, 0.4856720, 0.6810844, 0.5714698   0.5941448
12940   NC                       0.4874018, 0.5424656, 0.3437655, 0.4102728   0.4459764
519     NC                                             0.4703946, 0.2461935    0.358294
1581    NC                                  0.5724939, 0.3448565, 0.4659352   0.4610952
459     NC                                             0.6755135, 0.4776196   0.5765665
1541    NC                                             0.4792871, 0.2735542   0.3764207
1996    NC                                             0.5007048, 0.3009807   0.4008427
331     NC                                  0.4648444, 0.5897234, 0.3071130   0.4538936
442     NC                                  0.4378788, 0.2352466, 0.2721108   0.3150787
3987    NC            0.5232709, 0.5661026, 0.5040197, 0.5410826, 0.4760760   0.5221104
411     NC                                             0.4405761, 0.3043518    0.372464
766     NC                                             0.4867494, 0.2974405    0.392095
515     NC                       0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                                             0.4991863, 0.2996197    0.399403
183     NC                                             0.4863329, 0.2652150   0.3757739
336     NC                                             0.5053493, 0.3706932   0.4380213
11043   NC                                  0.4335317, 0.4984163, 0.3944782   0.4421421
1128    NC                                             0.4553747, 0.2413562   0.3483655
1538    NC                       0.6057217, 0.5218727, 0.4647596, 0.4578509   0.5125512
2330    NC                                  0.4678680, 0.2644552, 0.4902931   0.4075388
722     NC                       0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                                       Ruletype
6093  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
12940                      T1wdynamic & T1wtexture & dispersion
519                                  T1wmorphology & T1wtexture
1581                              T1wtexture & dispersion & T2w
459                                  T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
1996                                 T1wmorphology & T1wtexture
331                     T1wdynamic & T1wmorphology & T1wtexture
442                                  T1wmorphology & T1wtexture
3987                  T1wtexture & dispersion & single-time-Enh
411                                  T1wdynamic & T1wmorphology
766                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
183                                  T1wmorphology & T1wtexture
336                                  T1wdynamic & T1wmorphology
11043                              T1wtexture & single-time-Enh
1128                                 T1wmorphology & T1wtexture
1538     T1wdynamic & T1wtexture & dispersion & single-time-Enh
2330                    T1wmorphology & T1wtexture & dispersion
722          T1wdynamic & T1wmorphology & single-time-Enh & T2w
```

```r
# print(FA[FA$allTestinfo.lesion_id==241,])
imgFA <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgFA)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
## with LMSIR
if (length(idLMSIR) > 0) {
    for (i in 1:length(idLMSIR)) {
        idx = idLMSIR[i]
        lesion_id = FA[idx, ]$allTestinfo.lesion_id
        print(lesion_id)
    }
}
idx = c(1:nrow(FA))[c(FA$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
    lesion_id           lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
436       436 436_3046_7289130.vtk          3046 1970-11-16 00:00:00.000000              29427
    exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
436           7289130 2012-12-22 00:00:00.000000     Benign by pathology
    cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
436                           None                     1                        0
    exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
436              Right                 2691 2013-01-22 00:00:00.000000              Right
    proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int proc_lesion_comments_txt
436            Radiology                 US Core Needle Biopsy                     None
    find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
436           None                       N/A                        N/A      4        massB
    lesion_diagnosis
436     FIBROADENOMA
```

```r
print(mypresentRules(eachFATopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
6093    6 0.031                  0
12940   4 0.049                  0
519     2 0.176              0.053
1581    3 0.096 0.0570000000000001
459     2 0.064 0.0570000000000001
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
331     3 0.106               0.07
442     3 0.237               0.07
3987    5 0.078              0.071
411     2 0.228              0.072
766     2   0.2              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
183     2 0.232              0.079
336     2 0.246              0.081
11043   3 0.111              0.082
1128    2 0.215              0.085
1538    4 0.085              0.087
2330    3 0.216              0.093
722     4 0.136              0.095
                                                                                                                                                                                                                                               condition
6093  SER(in) low to medium  & Variance of spatial Margin Gradient low to very high  & std 3D Sharpness of lesion margin  low to medium  & Inverse difference moment ptime4 medium to high  & 2nd post-SE s12 very low  & 2nd post-SE s17 low to medium 
12940                                                                                                                    Uptake(in) curve amplitude medium to high  & SER(in) low to medium  & Variance ptime3 low to medium  & dispersion s11 very low 
519                                                                                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1581                                                                                                                                                          Difference variance ptime1 low  & dispersion s16 low to medium  & T2w margin gradient low 
459                                                                                                                                                                                               irregularity very low  & dispersion s9 medium to high 
1541                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
331                                                                                                                                                                           SER(in) medium to high  & irregularity low  & Energy ptime4 low to medium 
442                                                                                                                                                                 irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
3987                                                              Inverse difference moment ptime1 very low to high  & dispersion s5 low to very high  & dispersion s7 low to very high  & 1st post-SE s19 low to very high  & 3rd post-SE s12 very low 
411                                                                                                                                                                                                           SER(in) low to medium  & irregularity low 
766                                                                                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
515                                                                                                                       SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
183                                                                                                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
336                                                                                                                                                                                                           SER(in) low to medium  & irregularity low 
11043                                                                                                                                                         Sum Entropy ptime1 low to medium  & 1st post-SE s8 low to medium  & last post-SE s17 high 
1128                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1538                                                                                                             Uptake(in) curve amplitude medium to high  & Sum Entropy ptime4 low to medium  & dispersion s17 low to very high  & 3rd post-SE s4 low 
2330                                                                                                                                                                  irregularity low  & Energy ptime3 low to medium  & dispersion s5 low to very high 
722                                                                                                                               uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
      pred                                                             prob sumtempProb
6093    NC 0.6855391, 0.5934397, 0.5476637, 0.4856720, 0.6810844, 0.5714698   0.5941448
12940   NC                       0.4874018, 0.5424656, 0.3437655, 0.4102728   0.4459764
519     NC                                             0.4703946, 0.2461935    0.358294
1581    NC                                  0.5724939, 0.3448565, 0.4659352   0.4610952
459     NC                                             0.6755135, 0.4776196   0.5765665
1541    NC                                             0.4792871, 0.2735542   0.3764207
1996    NC                                             0.5007048, 0.3009807   0.4008427
331     NC                                  0.4648444, 0.5897234, 0.3071130   0.4538936
442     NC                                  0.4378788, 0.2352466, 0.2721108   0.3150787
3987    NC            0.5232709, 0.5661026, 0.5040197, 0.5410826, 0.4760760   0.5221104
411     NC                                             0.4405761, 0.3043518    0.372464
766     NC                                             0.4867494, 0.2974405    0.392095
515     NC                       0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                                             0.4991863, 0.2996197    0.399403
183     NC                                             0.4863329, 0.2652150   0.3757739
336     NC                                             0.5053493, 0.3706932   0.4380213
11043   NC                                  0.4335317, 0.4984163, 0.3944782   0.4421421
1128    NC                                             0.4553747, 0.2413562   0.3483655
1538    NC                       0.6057217, 0.5218727, 0.4647596, 0.4578509   0.5125512
2330    NC                                  0.4678680, 0.2644552, 0.4902931   0.4075388
722     NC                       0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                                       Ruletype
6093  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
12940                      T1wdynamic & T1wtexture & dispersion
519                                  T1wmorphology & T1wtexture
1581                              T1wtexture & dispersion & T2w
459                                  T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
1996                                 T1wmorphology & T1wtexture
331                     T1wdynamic & T1wmorphology & T1wtexture
442                                  T1wmorphology & T1wtexture
3987                  T1wtexture & dispersion & single-time-Enh
411                                  T1wdynamic & T1wmorphology
766                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
183                                  T1wmorphology & T1wtexture
336                                  T1wdynamic & T1wmorphology
11043                              T1wtexture & single-time-Enh
1128                                 T1wmorphology & T1wtexture
1538     T1wdynamic & T1wtexture & dispersion & single-time-Enh
2330                    T1wmorphology & T1wtexture & dispersion
722          T1wdynamic & T1wmorphology & single-time-Enh & T2w
```


![](Z:/Cristina/Section2/papernew_notes/images/436_3046_7289130_massB_Hyperintense.png)

# 4) Present rules collect top for FIBROCYSTIC FI

```r
eachFITopRules = list()
allFItoprules = c()
for (idx in 1:nrow(FI)) {
    ############## first analysis
    X = FI[idx, 2:ncol(imgT2pLMSIRtest)]
    y = FI[idx, "lesion_label"]
    
    rulesoutputFI = mynewapplyLearnerxrules(topRulesFI, X, y, minerr = 0.1, minfrq = 10/627, classes, 
        gbmModel)
    topRulesFI = rulesoutputFI[[1]]
    
    if (length(rulesoutputFI) > 1) {
        eachFITopRules[[idx]] = rulesoutputFI[[2]]
        allFItoprules = rbind(allFItoprules, rulesoutputFI[[3]])
    } else {
        eachFITopRules[[idx]] = list()
        allFItoprules = rbind(allFItoprules, 1:nrow(topRulesFI) * 0)
    }
}
```

```
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 16 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 19 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 20 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 2 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 17 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 15 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 10 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 8 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 4 out of 50"
[1] "test complies with # rules: 13 out of 50"
[1] "test complies with # rules: 14 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 9 out of 50"
[1] "test complies with # rules: 12 out of 50"
[1] "test complies with # rules: 6 out of 50"
[1] "test complies with # rules: 11 out of 50"
[1] "test complies with # rules: 1 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 3 out of 50"
[1] "test complies with # rules: 5 out of 50"
[1] "test complies with # rules: 7 out of 50"
[1] "test complies with # rules: 7 out of 50"
```

```r
# take the occurance of each explaining rule (columns) per case(row) + lesion identifies
df = data.frame(cbind(allFItoprules, as.character(FI$allTestinfo.lesion_id), as.character(FI$allTestinfo.lesion_label), 
    as.character(FI$allTestinfo.lesion_diagnosis)), stringsAsFactors = FALSE)
colnames(df) <- c(as.character(1:nrow(topRulesFI)), "lesion_id", "lesion_label", "lesion_diagnosis")
cols = 1:nrow(topRulesFI)
df[, cols] = apply(df[, cols], 2, function(x) as.numeric(x))


cat(paste0("\n============== Top five explaining rules for FIBROCYSTIC (FI): n = ", nrow(FA), "\n"))
```

```

============== Top five explaining rules for FIBROCYSTIC (FI): n = 69
```

```r
topFI = sort(colSums(df[, 1:nrow(topRulesFI)]), decreasing = TRUE)
print(topFI[1:5])
```

```
42 35 39 38 40 
27 25 25 24 24 
```

```r
freqtopFI = topFI[1:5]
print(freqtopFI/nrow(FI))
```

```
       42        35        39        38        40 
0.3857143 0.3571429 0.3571429 0.3428571 0.3428571 
```

```r
rulesTopFI = as.numeric(names(topFI))
preserulesTopFI = mypresentRules(topRulesFI[rulesTopFI, ], colnames(X), fnames)
rownames(preserulesTopFI) <- NULL
print(preserulesTopFI)
```

```
   len  freq                err
1    2  0.25              0.094
2    5  0.17              0.087
3    4 0.202              0.092
4    3 0.185               0.09
5    2 0.215              0.093
6    4 0.195              0.074
7    2 0.193              0.077
8    2 0.197              0.073
9    2 0.184              0.069
10   2 0.162 0.0570000000000001
11   5 0.168              0.075
12   3 0.127              0.043
13   3 0.177              0.093
14   4 0.094              0.077
15   2 0.124              0.088
16   2 0.091              0.082
17   4 0.136              0.095
18   3 0.053              0.034
19   4 0.103              0.088
20   4 0.066              0.028
21   6 0.096 0.0580000000000001
22   3  0.06 0.0610000000000001
23   2 0.052              0.071
24   3 0.095              0.096
25   3 0.116              0.095
26   4 0.041                  0
27   6 0.047              0.038
28   5  0.07              0.077
29   3 0.108              0.083
30   2 0.057              0.097
31   4 0.057              0.097
32   2 0.071              0.051
33   3 0.076              0.095
34   4 0.057              0.032
35   4 0.055              0.033
36   3 0.048              0.038
37   2 0.048              0.038
38   5 0.038              0.048
39   3 0.057 0.0649999999999999
40   3 0.057              0.097
41   4 0.031 0.0590000000000001
42   2 0.045               0.08
43   4 0.066              0.083
44   2 0.044              0.042
45   3 0.038              0.095
46   3 0.031                  0
47   2 0.135 0.0679999999999999
48   3 0.107 0.0679999999999999
49   3  0.08 0.0679999999999999
50   2  0.14              0.078
                                                                                                                                                                                                               condition
1                                                                                                                                                    irregularity low  & Inverse difference moment ptime4 low to medium 
2                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
4                                                                                                                   Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
5                                                                                                                                                                             SER(in) low  & T2w Energy  medium to high 
6                                                                                         SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
7                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
8                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
9                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
10                                                                                                                                                                      irregularity low  & Energy ptime4 low to medium 
11                                     Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
12                                                                                                                                       uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
13                                                                                                                                  SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
14                                                                  Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
15                                                                                                                                                       circularity high & T2w radial gradient variance medium to high 
16                                                                                                                                                                            circularity high & Sum average ptime1 low 
17                                                                                                uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
18                                                                                                        enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
19                                                                                                                     SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
20                                                                      enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
21 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
22                                                                                                      circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
23                                                                                                                                                     SER(in) low  & enhancement-variance increasing rate(in) very low 
24                                                                                               enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
25                                                                          SER(in) very low  & enhancement-variance increasing rate(in) very low to high  & enhancement-variance decreasing rate(rim) very low to high 
26                                                                                          Initial Uptake slope(in) low  & uptake average high  & 1st post-SE s15 low to very high  & 1st post-SE s16 low to very high 
27                                    circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
28                                                                       SER(rim) very low to high  & dispersion s14 very low  & dispersion s15 very low to high  & 1st post-SE s4 low to very high  & T2w Entropy  low 
29                                                                                                                                irregularity low to medium  & 3rd post-SE s10 low to very high  & 3rd post-SE s16 low 
30                                                                                                       max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
31                                                             enhancement-variance at first time-point(rim) low to medium  & uptake variance medium to high  & max Radial gradient low  & dispersion s2 medium to high 
32                                                                                                                                                        3rd post-SE s10 low to very high  & last post-SE s11 very low 
33                                                                                                           enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
34                                                    Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
35                                                                                     min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
36                                                                                                                                               SER(in) low  & irregularity medium to high  & last post-SE s2 very low 
37                                                                                                                                                            Energy ptime1 very low to high  & Entropy ptime2 very low 
38                                              Max Radial gradient variance very low  & Sum average ptime3 medium to high  & Entropy ptime4 very low to high  & dispersion s18 low  & 2nd post-SE s15 very low to high 
39                                                                                                  irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
40                                                                                                                                    SER(rim) low to medium  & irregularity medium to high  & last post-SE s2 very low 
41                                                                Rate of signal increase(in) low  & Rate of signal increase(in) medium to high  & 1st post-SE s15 low to very high  & last post-SE s8 low to very high 
42                                                                                                                                                                        dispersion s5 very low  & 1st post-SE s11 low 
43                                                                            uptake skewness medium to high  & irregularity medium to high  & kurtosis T2w SI low to very high  & T2w average radial gradient very low 
44                                                                                                                                                 irregularity very low  & T2w radial gradient variance medium to high 
45                                                                                                             Rate of signal increase(rim) very low  & dispersion s17 low to very high  & T2w mean SI low to very high 
46                                                                                                                                             SER(rim) very low  & 1st post-SE s6 very low  & last post-SE s7 very low 
47                                                                                                                                                                        irregularity high  & Sum variance ptime1 high 
48                                                                                                              irregularity medium to high  & Max Radial gradient variance medium to high  & Variance ptime1 very high 
49                                                                                                                         uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
50                                                                                                                                                                        irregularity high  & Sum variance ptime1 high 
   pred                                                             prob sumtempProb
1    NC                                             0.4236079, 0.2408969   0.3322524
2    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
4    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
5    NC                                             0.4881946, 0.2745578   0.3813762
6    NC                       0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
7    NC                                             0.4991863, 0.2996197    0.399403
8    NC                                             0.4694149, 0.2405692   0.3549921
9    NC                                             0.5007048, 0.3009807   0.4008427
10   NC                                             0.4792871, 0.2735542   0.3764207
11   NC            0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
12   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
13   NC                                  0.5505793, 0.3818076, 0.5017071   0.4780313
14   NC                       0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
15   NC                                             0.5514954, 0.4847561   0.5181258
16   NC                                             0.5779883, 0.4539037    0.515946
17   NC                       0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
18   NC                                  0.4848342, 0.5546800, 0.6219773   0.5538305
19   NC                       0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
20   NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
21   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
22   NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
23   NC                                             0.6114913, 0.4662233   0.5388573
24   NC                                  0.4806604, 0.4285987, 0.5949174   0.5013922
25   NC                                  0.5192031, 0.5443866, 0.4975629   0.5203842
26   NC                       0.6920559, 0.4833696, 0.7240647, 0.5912338    0.622681
27   NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
28   NC            0.4499065, 0.4352475, 0.4759258, 0.4943733, 0.5792785   0.4869463
29   NC                                  0.5052503, 0.4047316, 0.3668232   0.4256017
30   NC                                             0.6236805, 0.4747261   0.5492033
31   NC                       0.5904771, 0.6088937, 0.5063315, 0.6877946   0.5983742
32   NC                                             0.5144900, 0.4830707   0.4987804
33   NC                                  0.3132262, 0.5008952, 0.5531483   0.4557566
34   NC                       0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
35   NC                       0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
36   NC                                  0.6131108, 0.4701726, 0.3805070   0.4879301
37   NC                                             0.5127610, 0.4835754   0.4981682
38   NC            0.4289456, 0.6226296, 0.5110246, 0.5385531, 0.4767274   0.5155761
39   NC                                  0.5914519, 0.4717202, 0.4894149    0.517529
40   NC                                  0.5743398, 0.4464800, 0.3780448   0.4662882
41   NC                       0.5271721, 0.3903370, 0.5439983, 0.4908403    0.488087
42   NC                                             0.4735013, 0.6948949   0.5841981
43   NC                       0.6209921, 0.4870514, 0.5075204, 0.4635613   0.5197813
44   NC                                             0.5748740, 0.4661645   0.5205192
45   NC                                  0.5582364, 0.5165599, 0.5742535   0.5496832
46   NC                                  0.6420296, 0.5124139, 0.4394809   0.5313081
47    C                                             0.5087730, 0.8044078   0.6565904
48    C                                  0.4969203, 0.6008289, 0.5469925   0.5482472
49    C                                  0.5436611, 0.6397459, 0.5003391   0.5612487
50    C                                             0.4957177, 0.7767643    0.636241
                                                    Ruletype
1                                 T1wmorphology & T1wtexture
2              T1wmorphology & T1wtexture & dispersion & T2w
3                           T1wmorphology & dispersion & T2w
4                                 T1wdynamic & T1wmorphology
5                                           T1wdynamic & T2w
6  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
7                                 T1wmorphology & T1wtexture
8                                 T1wmorphology & T1wtexture
9                                 T1wmorphology & T1wtexture
10                                T1wmorphology & T1wtexture
11                        T1wdynamic & single-time-Enh & T2w
12                   T1wdynamic & T1wmorphology & dispersion
13                             T1wdynamic & T1wtexture & T2w
14                 T1wdynamic & dispersion & single-time-Enh
15                                       T1wmorphology & T2w
16                                T1wmorphology & T1wtexture
17        T1wdynamic & T1wmorphology & single-time-Enh & T2w
18                          T1wdynamic & T1wmorphology & T2w
19                 T1wdynamic & dispersion & single-time-Enh
20                   T1wdynamic & T1wmorphology & T1wtexture
21    T1wdynamic & T1wtexture & dispersion & single-time-Enh
22                                       T1wmorphology & T2w
23                                                T1wdynamic
24                 T1wdynamic & T1wtexture & single-time-Enh
25                                                T1wdynamic
26                              T1wdynamic & single-time-Enh
27              T1wmorphology & dispersion & single-time-Enh
28           T1wdynamic & dispersion & single-time-Enh & T2w
29                           T1wmorphology & single-time-Enh
30                                T1wdynamic & T1wmorphology
31                   T1wdynamic & T1wmorphology & dispersion
32                                           single-time-Enh
33                 T1wdynamic & dispersion & single-time-Enh
34        T1wmorphology & T1wtexture & single-time-Enh & T2w
35      T1wdynamic & T1wmorphology & T1wtexture & dispersion
36              T1wdynamic & T1wmorphology & single-time-Enh
37                                                T1wtexture
38 T1wmorphology & T1wtexture & dispersion & single-time-Enh
39                          T1wmorphology & T1wtexture & T2w
40              T1wdynamic & T1wmorphology & single-time-Enh
41                              T1wdynamic & single-time-Enh
42                              dispersion & single-time-Enh
43                          T1wdynamic & T1wmorphology & T2w
44                                       T1wmorphology & T2w
45                             T1wdynamic & dispersion & T2w
46                              T1wdynamic & single-time-Enh
47                                T1wmorphology & T1wtexture
48                                T1wmorphology & T1wtexture
49                 T1wdynamic & T1wtexture & single-time-Enh
50                                T1wmorphology & T1wtexture
```

```r
# display a case that meets them all
casesFI = df[, c(rulesTopFI, nrow(topRulesFI) + 1, nrow(topRulesFI) + 2, nrow(topRulesFI) + 3)]
print("Sorted cases meeting top five explaining rules for FIBROCYSTIC (FI):")
```

```
[1] "Sorted cases meeting top five explaining rules for FIBROCYSTIC (FI):"
```

```r
topcasesFI = casesFI[sort(rowSums(casesFI[, 1:5]), index = TRUE, decreasing = TRUE)$ix, ]
print(topcasesFI)
```

```
   42 35 39 38 40 25 27 24 22 14 26 11 41 29 37 32 46 6 36 3 15 17 23 47 45 1 8 28 33 49 50 13
9   1  1  1  1  1  0  1  1  1  1  0  1  0  0  1  0  1 1  0 0  0  0  1  0  1 0 0  0  0  0  0  0
14  1  1  1  1  1  1  1  1  1  1  0  1  1  0  1  1  0 0  0 1  1  1  0  1  0 0 0  0  0  1  0  0
15  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  0 0  0 1  1  1  0  0  0 0 0  0  0  1  0  0
21  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  1  1 1  0 1  0  0  1  0  0 1 0  0  1  0  0  0
3   1  1  1  1  0  0  1  1  1  1  0  0  0  0  1  1  0 0  0 0  1  1  0  0  0 0 0  0  0  1  0  0
8   1  1  1  1  0  0  0  0  0  0  0  0  0  1  0  1  0 0  0 1  0  0  0  0  0 0 0  1  0  0  0  0
13  1  1  1  1  0  0  1  1  1  1  0  1  1  0  1  1  0 1  0 1  1  0  0  0  0 0 0  0  0  0  0  0
18  1  1  0  1  1  0  1  1  1  1  0  0  0  0  0  0  0 0  1 0  0  0  1  0  1 0 0  1  0  1  0  0
26  1  1  1  0  1  1  1  1  1  1  1  0  0  1  0  0  1 0  1 0  0  0  1  0  0 1 0  0  0  0  0  0
50  1  0  1  1  1  1  1  1  1  1  0  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
54  1  0  1  1  1  1  1  1  1  1  1  0  0  1  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
55  1  1  1  1  0  0  1  1  1  1  0  1  0  1  0  1  0 0  1 0  0  0  0  0  0 0 1  0  0  1  0  0
5   0  1  1  1  0  0  0  0  0  0  0  0  0  0  1  1  0 0  0 0  0  1  0  0  0 0 0  0  0  0  0  0
16  1  1  1  0  0  0  1  1  1  1  0  0  0  0  0  0  0 1  1 0  0  0  0  0  0 0 1  0  1  0  0  0
19  1  0  1  0  1  1  1  1  1  1  0  1  0  0  0  0  0 1  0 0  0  0  1  0  1 0 0  0  1  0  0  0
28  1  1  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0 1  1 0  0  0  0  0  0 0 0  1  0  0  0  0
30  1  1  1  0  0  0  1  1  1  1  0  0  0  0  1  1  0 0  0 0  1  0  0  1  0 0 0  0  0  0  0  0
32  1  0  1  1  0  0  0  0  0  0  1  1  1  1  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
33  1  0  1  1  0  0  1  1  1  1  0  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
44  1  1  1  0  0  0  0  0  0  0  0  0  1  1  1  1  0 0  0 1  0  1  0  0  1 0 0  0  0  1  0  0
45  1  0  0  1  1  0  0  0  0  0  0  0  1  0  0  0  1 1  1 1  0  0  1  0  0 1 0  0  0  0  0  0
57  1  1  1  0  0  0  1  1  1  1  0  1  0  0  1  0  0 0  0 0  0  1  0  0  0 0 1  0  0  0  0  1
64  1  1  1  0  0  0  1  1  1  1  0  0  0  1  1  0  0 0  0 1  1  0  0  0  0 0 0  0  0  0  0  0
2   0  1  1  0  0  0  1  1  1  0  0  1  0  0  1  0  0 0  0 0  0  1  0  1  0 0 0  0  0  0  0  0
12  0  1  1  0  0  0  0  0  0  0  1  0  1  0  0  1  1 0  0 1  0  0  0  0  0 0 0  0  0  0  0  0
35  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  1  0 1  0 0  0  1  0  0  0 0 0  0  0  0  0  1
36  0  0  0  1  1  1  0  0  0  0  1  1  1  0  1  0  0 1  0 0  0  1  1  0  0 1 0  0  0  0  0  0
42  1  1  0  0  0  0  0  0  0  0  1  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 1  1  0  0  0  0
43  1  0  0  1  0  0  1  1  1  1  1  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  1  0  0
48  0  1  0  1  0  0  0  0  0  0  1  1  1  1  0  0  0 1  0 1  1  0  0  0  0 0 1  0  0  0  0  0
51  0  0  1  1  0  0  0  0  0  0  0  1  0  0  1  0  0 1  0 0  0  0  0  0  0 0 1  0  0  0  0  0
56  1  1  0  0  0  0  1  1  1  1  0  0  0  0  1  0  0 0  0 0  0  1  0  0  0 0 0  0  0  0  0  0
62  1  0  0  0  1  1  0  0  0  0  1  0  1  1  0  0  0 0  0 0  1  0  1  0  1 1 1  0  0  0  0  0
63  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 0  1 0  0  0  0  1  0 0 0  1  0  0  0  0
68  0  1  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  1  0  0  1
4   0  0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  1
11  0  0  0  0  1  1  0  0  0  0  1  0  0  0  0  0  0 0  1 0  0  0  0  0  0 1 0  0  0  0  0  0
17  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
20  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
23  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  1  0 0  0 0  1  1  0  0  0 0 0  0  0  0  0  0
27  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  1  0 1  1 1  0  0  0  0  0 0 0  1  1  0  0  1
37  0  0  0  0  1  1  0  0  0  0  1  0  0  0  0  0  0 0  0 0  0  0  0  0  0 1 0  0  0  0  0  0
39  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  1 0  0  0  0  0  0 0 0  1  0  0  0  0
46  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  1  0  0
47  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  1  0 0  0 0  0  0  0  0  0 0 1  0  0  0  0  0
49  0  0  0  0  1  1  0  0  0  0  1  0  1  0  0  0  0 0  1 0  1  0  0  0  1 0 0  1  1  0  0  0
52  0  0  0  0  1  1  0  0  0  0  0  0  0  1  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
61  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  1 0 0  0  1  0  1  1
66  0  0  0  0  1  1  0  0  0  0  1  0  0  0  0  0  1 0  0 0  0  0  0  0  0 1 0  0  0  0  1  0
69  0  0  0  0  1  1  0  0  0  0  1  0  1  0  0  0  1 0  0 0  0  0  1  0  0 0 0  0  0  0  1  0
70  0  0  0  0  1  1  0  0  0  0  0  0  1  0  0  0  1 0  0 0  0  0  0  0  1 0 0  0  0  0  0  1
1   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
6   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
7   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
10  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  1 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
22  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1 0  0 0  0  0  1  0  1 0 0  0  0  0  1  0
24  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  1  0
25  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 0  1 0  0  0  0  0  0 0 0  0  0  0  0  0
29  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  1  0
31  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0  0  1 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
34  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
38  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
40  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
41  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1 0  0 0  0  0  0  0  0 0 0  0  0  0  1  0
53  0  0  0  0  0  1  0  0  0  0  1  0  0  0  0  0  0 0  0 0  0  0  1  0  0 1 0  0  0  0  0  0
58  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
59  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0  0  0 0  0 0  1  0  0  0  0 0 1  0  1  0  0  0
60  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0 0  0 0  0  0  0  0  1 0 0  0  0  0  0  0
65  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 0  0 0  0  0  0  1  0 0 0  0  0  0  0  0
67  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 0  0 0  0  0  0  0  0 0 0  0  0  0  1  0
   43 4 5 7 9 12 18 48 16 31 34 10 44 2 19 20 21 30 lesion_id lesion_label lesion_diagnosis
9   0 0 0 0 0  1  1  0  0  0  0  0  0 0  0  0  0  0        74        massB      FIBROCYSTIC
14  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       411     nonmassB      FIBROCYSTIC
15  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       412     nonmassB      FIBROCYSTIC
21  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       311        massB      FIBROCYSTIC
3   0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       301     nonmassB      FIBROCYSTIC
8   1 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        66     nonmassB      FIBROCYSTIC
13  0 0 0 0 0  0  1  0  0  0  0  0  0 0  0  0  0  0       383     nonmassB      FIBROCYSTIC
18  1 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       212        massB      FIBROCYSTIC
26  1 0 0 0 0  1  0  0  0  0  0  0  0 0  0  0  0  0       135        massB      FIBROCYSTIC
50  0 0 0 0 0  0  0  0  1  0  0  0  0 0  0  0  0  0       627        massB      FIBROCYSTIC
54  0 0 1 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        48        massB      FIBROCYSTIC
55  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       263        massB      FIBROCYSTIC
5   0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       503        massB      FIBROCYSTIC
16  1 1 0 0 0  1  0  0  0  0  0  0  0 0  0  0  0  0        93        massB      FIBROCYSTIC
19  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       213        massB      FIBROCYSTIC
28  1 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       317     nonmassB      FIBROCYSTIC
30  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       396     nonmassB      FIBROCYSTIC
32  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       613     nonmassB      FIBROCYSTIC
33  0 1 1 0 0  0  0  0  0  0  0  0  0 0  0  0  1  0        10        massB      FIBROCYSTIC
44  0 0 0 0 1  0  1  0  0  0  0  1  0 0  0  0  0  0       166        massB      FIBROCYSTIC
45  0 1 0 0 1  1  0  0  0  1  0  0  0 0  0  0  0  0       193     nonmassB      FIBROCYSTIC
57  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       267     nonmassB      FIBROCYSTIC
64  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        29        massB      FIBROCYSTIC
2   0 0 1 0 0  0  1  0  0  0  0  1  0 0  0  0  0  0       228        massB      FIBROCYSTIC
12  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       288        massB      FIBROCYSTIC
35  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       256        massB      FIBROCYSTIC
36  0 0 0 0 1  0  1  0  0  1  0  1  0 0  0  0  0  0       274        massB      FIBROCYSTIC
42  0 0 0 0 0  0  0  0  0  0  0  0  1 0  0  0  0  0        79     nonmassB      FIBROCYSTIC
43  0 1 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       137        massB      FIBROCYSTIC
48  0 0 0 0 0  0  0  0  1  0  0  0  0 0  0  0  0  0       579     nonmassB      FIBROCYSTIC
51  0 0 0 0 1  0  0  0  1  0  0  0  0 0  0  0  0  0        26        massB      FIBROCYSTIC
56  0 0 0 0 0  1  0  0  0  0  0  0  0 0  0  0  0  0       266     nonmassB      FIBROCYSTIC
62  0 0 0 0 0  0  0  0  0  0  0  0  1 0  0  0  0  0       590     nonmassB      FIBROCYSTIC
63  0 1 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        25        massB      FIBROCYSTIC
68  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       330     nonmassB      FIBROCYSTIC
4   0 0 0 1 0  0  0  1  0  0  0  0  0 0  1  1  0  1       415     nonmassB      FIBROCYSTIC
11  1 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       287        massB      FIBROCYSTIC
17  0 0 0 1 0  0  0  1  0  0  0  0  0 0  0  0  0  0       146     nonmassB      FIBROCYSTIC
20  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       233        massB      FIBROCYSTIC
23  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       419     nonmassB      FIBROCYSTIC
27  1 0 0 0 0  1  0  0  0  1  0  0  0 0  0  0  0  0       302        massB      FIBROCYSTIC
37  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       636     nonmassB      FIBROCYSTIC
39  0 1 0 0 0  0  1  0  0  1  0  1  1 0  0  0  0  0        76        massB      FIBROCYSTIC
46  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       348     nonmassB      FIBROCYSTIC
47  0 0 1 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       356     nonmassB      FIBROCYSTIC
49  0 0 0 0 1  0  0  0  0  1  0  0  0 0  0  0  0  0       211        massB      FIBROCYSTIC
52  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        27     nonmassB      FIBROCYSTIC
61  0 0 0 1 0  0  0  1  0  0  1  0  0 1  0  0  0  0       537     nonmassB      FIBROCYSTIC
66  0 0 0 0 0  0  0  0  0  0  1  0  0 0  0  0  0  0       141        massB      FIBROCYSTIC
69  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       468        massB      FIBROCYSTIC
70  0 0 0 0 0  0  0  0  0  0  1  0  0 0  0  0  0  0       469     nonmassB      FIBROCYSTIC
1   0 0 0 0 0  0  0  0  1  0  0  0  0 0  1  1  1  1       165        massB      FIBROCYSTIC
6   0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        57     nonmassB      FIBROCYSTIC
7   0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        58     nonmassB      FIBROCYSTIC
10  0 0 1 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       123        massB      FIBROCYSTIC
22  0 0 0 1 0  0  0  1  0  0  0  0  0 1  0  0  0  0       339     nonmassB      FIBROCYSTIC
24  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       461     nonmassB      FIBROCYSTIC
25  0 0 1 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       462        massB      FIBROCYSTIC
29  0 0 0 1 0  0  0  1  0  0  0  0  0 1  0  0  0  0       338     nonmassB      FIBROCYSTIC
31  0 0 0 0 0  0  0  0  0  0  1  0  0 0  0  0  0  0       558        massB      FIBROCYSTIC
34  0 0 0 0 0  0  0  0  0  0  0  0  1 0  0  0  0  0       255        massB      FIBROCYSTIC
38  0 0 0 1 0  0  0  1  0  0  0  0  0 0  0  0  0  0        75        massB      FIBROCYSTIC
40  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       236     nonmassB      FIBROCYSTIC
41  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       524        massB      FIBROCYSTIC
53  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        31        massB      FIBROCYSTIC
58  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       307        massB      FIBROCYSTIC
59  0 0 0 0 1  0  0  0  1  0  0  0  0 0  0  0  0  0       426        massB      FIBROCYSTIC
60  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0       432        massB      FIBROCYSTIC
65  0 0 0 0 0  0  0  0  0  0  0  0  0 0  0  0  0  0        56     nonmassB      FIBROCYSTIC
67  0 0 0 0 0  0  0  0  0  0  1  0  0 0  0  0  0  0       231     nonmassB      FIBROCYSTIC
```

```r
print("Fine explaining rules for FIBROCYSTIC (FI) with LMSIR:")
```

```
[1] "Fine explaining rules for FIBROCYSTIC (FI) with LMSIR:"
```

```r
LMSIRfeatureflag = c()
idLMSIR = c()
for (i in 1:length(eachFITopRules)) {
    if (length(eachFITopRules[[i]]) > 0) {
        rules = mypresentRules(eachFITopRules[[i]], colnames(X), fnames)
        LMSIRwtop5 = "LMSIR" %in% unlist(strsplit(unlist(strsplit(rules$condition, split = " & ")), 
            split = " "))
        if (LMSIRwtop5) {
            LMSIRfeatureflag = c(LMSIRfeatureflag, LMSIRwtop5)
            print(i)
            print(rules)
            idLMSIR = c(idLMSIR, i)
        } else {
            LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
        }
    } else {
        LMSIRfeatureflag = c(LMSIRfeatureflag, FALSE)
    }
}
```

```
[1] 2
      len  freq                err
4971    4 0.055              0.033
7103    2 0.044              0.042
12721   3 0.127              0.043
1073    3  0.06 0.0610000000000001
1578    3 0.057 0.0649999999999999
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
4809    3 0.095              0.096
                                                                                                                                                       condition
4971                           min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
7103                                                                                       irregularity very low  & T2w radial gradient variance medium to high 
12721                                                                            uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1073                                            circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1578                                        irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
1996                                                                                                            irregularity low  & Energy ptime4 low to medium 
883                                                                                                             irregularity low  & Energy ptime4 low to medium 
1338                                                                                                            irregularity low  & Energy ptime4 low to medium 
5547  circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                             circularity high & T2w radial gradient variance medium to high 
6375                                  irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
4809                                     enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
      pred                                                  prob sumtempProb
4971    NC            0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
7103    NC                                  0.5748740, 0.4661645   0.5205192
12721   NC                       0.4677355, 0.5728223, 0.2321752   0.4242443
1073    NC                       0.6781025, 0.4757408, 0.5921319   0.5819917
1578    NC                       0.5914519, 0.4717202, 0.4894149    0.517529
1996    NC                                  0.5007048, 0.3009807   0.4008427
883     NC                                  0.4694149, 0.2405692   0.3549921
1338    NC                                  0.4991863, 0.2996197    0.399403
5547    NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                  0.5514954, 0.4847561   0.5181258
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
4809    NC                       0.4806604, 0.4285987, 0.5949174   0.5013922
                                                  Ruletype
4971  T1wdynamic & T1wmorphology & T1wtexture & dispersion
7103                                   T1wmorphology & T2w
12721              T1wdynamic & T1wmorphology & dispersion
1073                                   T1wmorphology & T2w
1578                      T1wmorphology & T1wtexture & T2w
1996                            T1wmorphology & T1wtexture
883                             T1wmorphology & T1wtexture
1338                            T1wmorphology & T1wtexture
5547         T1wmorphology & T1wtexture & dispersion & T2w
3902                                   T1wmorphology & T2w
6375                      T1wmorphology & dispersion & T2w
4809             T1wdynamic & T1wtexture & single-time-Enh
[1] 3
      len  freq                err
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1073    3  0.06 0.0610000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
750     2  0.25              0.094
9527    2 0.057              0.097
                                                                                                                                                                                                                  condition
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1073                                                                                                       circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
8424                                                                                                                   Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
9527                                                                                                        max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
      pred                                                             prob sumtempProb
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1073    NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
9527    NC                                             0.6236805, 0.4747261   0.5492033
                                                    Ruletype
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1073                                     T1wmorphology & T2w
1996                              T1wmorphology & T1wtexture
883                               T1wmorphology & T1wtexture
1338                              T1wmorphology & T1wtexture
7750                              T1wmorphology & T1wtexture
5547           T1wmorphology & T1wtexture & dispersion & T2w
3902                                     T1wmorphology & T2w
8424                              T1wdynamic & T1wmorphology
6375                        T1wmorphology & dispersion & T2w
750                               T1wmorphology & T1wtexture
9527                              T1wdynamic & T1wmorphology
[1] 5
     len  freq                err
1073   3  0.06 0.0610000000000001
7750   2 0.091              0.082
5547   5  0.17              0.087
3902   2 0.124              0.088
8424   3 0.185               0.09
6375   4 0.202              0.092
                                                                                                                                                      condition
1073                                           circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
7750                                                                                                                 circularity high & Sum average ptime1 low 
5547 circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                            circularity high & T2w radial gradient variance medium to high 
8424                                                       Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                 irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
     pred                                                  prob sumtempProb
1073   NC                       0.6781025, 0.4757408, 0.5921319   0.5819917
7750   NC                                  0.5779883, 0.4539037    0.515946
5547   NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902   NC                                  0.5514954, 0.4847561   0.5181258
8424   NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
                                          Ruletype
1073                           T1wmorphology & T2w
7750                    T1wmorphology & T1wtexture
5547 T1wmorphology & T1wtexture & dispersion & T2w
3902                           T1wmorphology & T2w
8424                    T1wdynamic & T1wmorphology
6375              T1wmorphology & dispersion & T2w
[1] 8
     len  freq   err
1233   4 0.066 0.028
3314   5  0.07 0.077
7111   4 0.094 0.077
7750   2 0.091 0.082
5547   5  0.17 0.087
8424   3 0.185  0.09
6375   4 0.202 0.092
750    2  0.25 0.094
2209   3 0.076 0.095
                                                                                                                                                      condition
1233           enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
3314            SER(rim) very low to high  & dispersion s14 very low  & dispersion s15 very low to high  & 1st post-SE s4 low to very high  & T2w Entropy  low 
7111       Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
7750                                                                                                                 circularity high & Sum average ptime1 low 
5547 circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
8424                                                       Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                 irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                         irregularity low  & Inverse difference moment ptime4 low to medium 
2209                                                enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
     pred                                                  prob sumtempProb
1233   NC            0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
3314   NC 0.4499065, 0.4352475, 0.4759258, 0.4943733, 0.5792785   0.4869463
7111   NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
7750   NC                                  0.5779883, 0.4539037    0.515946
5547   NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
8424   NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750    NC                                  0.4236079, 0.2408969   0.3322524
2209   NC                       0.3132262, 0.5008952, 0.5531483   0.4557566
                                            Ruletype
1233         T1wdynamic & T1wmorphology & T1wtexture
3314 T1wdynamic & dispersion & single-time-Enh & T2w
7111       T1wdynamic & dispersion & single-time-Enh
7750                      T1wmorphology & T1wtexture
5547   T1wmorphology & T1wtexture & dispersion & T2w
8424                      T1wdynamic & T1wmorphology
6375                T1wmorphology & dispersion & T2w
750                       T1wmorphology & T1wtexture
2209       T1wdynamic & dispersion & single-time-Enh
[1] 9
      len  freq                err
1342    3 0.053              0.034
12721   3 0.127              0.043
4833    5 0.038              0.048
1541    2 0.162 0.0570000000000001
1578    3 0.057 0.0649999999999999
1996    2 0.184              0.069
2076    2 0.052              0.071
883     2 0.197              0.073
1338    2 0.193              0.077
5547    5  0.17              0.087
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
726     2 0.215              0.093
750     2  0.25              0.094
5386    3 0.116              0.095
722     4 0.136              0.095
                                                                                                                                                                     condition
1342                                                            enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
12721                                                                                          uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
4833  Max Radial gradient variance very low  & Sum average ptime3 medium to high  & Entropy ptime4 very low to high  & dispersion s18 low  & 2nd post-SE s15 very low to high 
1541                                                                                                                          irregularity low  & Energy ptime4 low to medium 
1578                                                      irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
1996                                                                                                                          irregularity low  & Energy ptime4 low to medium 
2076                                                                                                         SER(in) low  & enhancement-variance increasing rate(in) very low 
883                                                                                                                           irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                          irregularity low  & Energy ptime4 low to medium 
5547                circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                           circularity high & T2w radial gradient variance medium to high 
8424                                                                      Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                 SER(in) low  & T2w Energy  medium to high 
750                                                                                                        irregularity low  & Inverse difference moment ptime4 low to medium 
5386                              SER(in) very low  & enhancement-variance increasing rate(in) very low to high  & enhancement-variance decreasing rate(rim) very low to high 
722                                                     uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
      pred                                                  prob sumtempProb
1342    NC                       0.4848342, 0.5546800, 0.6219773   0.5538305
12721   NC                       0.4677355, 0.5728223, 0.2321752   0.4242443
4833    NC 0.4289456, 0.6226296, 0.5110246, 0.5385531, 0.4767274   0.5155761
1541    NC                                  0.4792871, 0.2735542   0.3764207
1578    NC                       0.5914519, 0.4717202, 0.4894149    0.517529
1996    NC                                  0.5007048, 0.3009807   0.4008427
2076    NC                                  0.6114913, 0.4662233   0.5388573
883     NC                                  0.4694149, 0.2405692   0.3549921
1338    NC                                  0.4991863, 0.2996197    0.399403
5547    NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                  0.5514954, 0.4847561   0.5181258
8424    NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                                  0.4881946, 0.2745578   0.3813762
750     NC                                  0.4236079, 0.2408969   0.3322524
5386    NC                       0.5192031, 0.5443866, 0.4975629   0.5203842
722     NC            0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                                       Ruletype
1342                           T1wdynamic & T1wmorphology & T2w
12721                   T1wdynamic & T1wmorphology & dispersion
4833  T1wmorphology & T1wtexture & dispersion & single-time-Enh
1541                                 T1wmorphology & T1wtexture
1578                           T1wmorphology & T1wtexture & T2w
1996                                 T1wmorphology & T1wtexture
2076                                                 T1wdynamic
883                                  T1wmorphology & T1wtexture
1338                                 T1wmorphology & T1wtexture
5547              T1wmorphology & T1wtexture & dispersion & T2w
3902                                        T1wmorphology & T2w
8424                                 T1wdynamic & T1wmorphology
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
750                                  T1wmorphology & T1wtexture
5386                                                 T1wdynamic
722          T1wdynamic & T1wmorphology & single-time-Enh & T2w
[1] 12
     len  freq   err
1233   4 0.066 0.028
5316   5 0.168 0.075
7750   2 0.091 0.082
5547   5  0.17 0.087
6375   4 0.202 0.092
5658   3 0.177 0.093
722    4 0.136 0.095
                                                                                                                                                                             condition
1233                                  enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
5316 Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
7750                                                                                                                                        circularity high & Sum average ptime1 low 
5547                        circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
6375                                                        irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
5658                                                                                              SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
722                                                             uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
     pred                                                  prob sumtempProb
1233   NC            0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
7750   NC                                  0.5779883, 0.4539037    0.515946
5547   NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
5658   NC                       0.5505793, 0.3818076, 0.5017071   0.4780313
722    NC            0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                               Ruletype
1233            T1wdynamic & T1wmorphology & T1wtexture
5316                 T1wdynamic & single-time-Enh & T2w
7750                         T1wmorphology & T1wtexture
5547      T1wmorphology & T1wtexture & dispersion & T2w
6375                   T1wmorphology & dispersion & T2w
5658                      T1wdynamic & T1wtexture & T2w
722  T1wdynamic & T1wmorphology & single-time-Enh & T2w
[1] 13
      len  freq                err
1233    4 0.066              0.028
1342    3 0.053              0.034
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1578    3 0.057 0.0649999999999999
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
5658    3 0.177              0.093
750     2  0.25              0.094
                                                                                                                                                                                                                  condition
1233                                                                       enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
1342                                                                                                         enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
12721                                                                                                                                       uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1578                                                                                                   irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
8424                                                                                                                   Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
5658                                                                                                                                   SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                                             prob sumtempProb
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
1342    NC                                  0.4848342, 0.5546800, 0.6219773   0.5538305
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1578    NC                                  0.5914519, 0.4717202, 0.4894149    0.517529
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
5658    NC                                  0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                             0.4236079, 0.2408969   0.3322524
                                                    Ruletype
1233                 T1wdynamic & T1wmorphology & T1wtexture
1342                        T1wdynamic & T1wmorphology & T2w
12721                T1wdynamic & T1wmorphology & dispersion
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1578                        T1wmorphology & T1wtexture & T2w
1996                              T1wmorphology & T1wtexture
883                               T1wmorphology & T1wtexture
1338                              T1wmorphology & T1wtexture
7750                              T1wmorphology & T1wtexture
5547           T1wmorphology & T1wtexture & dispersion & T2w
3902                                     T1wmorphology & T2w
8424                              T1wdynamic & T1wmorphology
6375                        T1wmorphology & dispersion & T2w
5658                           T1wdynamic & T1wtexture & T2w
750                               T1wmorphology & T1wtexture
[1] 14
      len  freq                err
1233    4 0.066              0.028
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1073    3  0.06 0.0610000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
726     2 0.215              0.093
5658    3 0.177              0.093
750     2  0.25              0.094
4809    3 0.095              0.096
9527    2 0.057              0.097
                                                                                                                                                                                                                  condition
1233                                                                       enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
12721                                                                                                                                       uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1073                                                                                                       circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
515                                                                                          SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
8424                                                                                                                   Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                                                              SER(in) low  & T2w Energy  medium to high 
5658                                                                                                                                   SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
4809                                                                                                enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
9527                                                                                                        max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
      pred                                                             prob sumtempProb
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1073    NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
515     NC                       0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                                             0.4991863, 0.2996197    0.399403
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                                             0.4881946, 0.2745578   0.3813762
5658    NC                                  0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                             0.4236079, 0.2408969   0.3322524
4809    NC                                  0.4806604, 0.4285987, 0.5949174   0.5013922
9527    NC                                             0.6236805, 0.4747261   0.5492033
                                                       Ruletype
1233                    T1wdynamic & T1wmorphology & T1wtexture
12721                   T1wdynamic & T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
12632    T1wdynamic & T1wtexture & dispersion & single-time-Enh
1073                                        T1wmorphology & T2w
1996                                 T1wmorphology & T1wtexture
883                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
7750                                 T1wmorphology & T1wtexture
5547              T1wmorphology & T1wtexture & dispersion & T2w
3902                                        T1wmorphology & T2w
8424                                 T1wdynamic & T1wmorphology
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
5658                              T1wdynamic & T1wtexture & T2w
750                                  T1wmorphology & T1wtexture
4809                  T1wdynamic & T1wtexture & single-time-Enh
9527                                 T1wdynamic & T1wmorphology
[1] 15
      len  freq                err
1233    4 0.066              0.028
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1073    3  0.06 0.0610000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
7111    4 0.094              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
726     2 0.215              0.093
5658    3 0.177              0.093
750     2  0.25              0.094
9527    2 0.057              0.097
                                                                                                                                                                                                                  condition
1233                                                                       enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
12721                                                                                                                                       uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1073                                                                                                       circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
515                                                                                          SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7111                                                                   Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
8424                                                                                                                   Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                                                              SER(in) low  & T2w Energy  medium to high 
5658                                                                                                                                   SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
9527                                                                                                        max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
      pred                                                             prob sumtempProb
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1073    NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
515     NC                       0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                                             0.4991863, 0.2996197    0.399403
7111    NC                       0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                                             0.4881946, 0.2745578   0.3813762
5658    NC                                  0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                             0.4236079, 0.2408969   0.3322524
9527    NC                                             0.6236805, 0.4747261   0.5492033
                                                       Ruletype
1233                    T1wdynamic & T1wmorphology & T1wtexture
12721                   T1wdynamic & T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
12632    T1wdynamic & T1wtexture & dispersion & single-time-Enh
1073                                        T1wmorphology & T2w
1996                                 T1wmorphology & T1wtexture
883                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
7111                  T1wdynamic & dispersion & single-time-Enh
7750                                 T1wmorphology & T1wtexture
5547              T1wmorphology & T1wtexture & dispersion & T2w
3902                                        T1wmorphology & T2w
8424                                 T1wdynamic & T1wmorphology
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
5658                              T1wdynamic & T1wtexture & T2w
750                                  T1wmorphology & T1wtexture
9527                                 T1wdynamic & T1wmorphology
[1] 16
      len  freq                err
10203   4 0.057              0.032
1342    3 0.053              0.034
1593    6 0.047              0.038
4833    5 0.038              0.048
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
723     3 0.108              0.083
5547    5  0.17              0.087
8534    4 0.103              0.088
6375    4 0.202              0.092
750     2  0.25              0.094
2209    3 0.076              0.095
                                                                                                                                                                               condition
10203                 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
1342                                                                      enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
1593  circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
4833            Max Radial gradient variance very low  & Sum average ptime3 medium to high  & Entropy ptime4 very low to high  & dispersion s18 low  & 2nd post-SE s15 very low to high 
1541                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
723                                                                                               irregularity low to medium  & 3rd post-SE s10 low to very high  & 3rd post-SE s16 low 
5547                          circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
8534                                                                                   SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
6375                                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
2209                                                                         enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
      pred                                                             prob sumtempProb
10203   NC                       0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
1342    NC                                  0.4848342, 0.5546800, 0.6219773   0.5538305
1593    NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
4833    NC            0.4289456, 0.6226296, 0.5110246, 0.5385531, 0.4767274   0.5155761
1541    NC                                             0.4792871, 0.2735542   0.3764207
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
723     NC                                  0.5052503, 0.4047316, 0.3668232   0.4256017
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
8534    NC                       0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
2209    NC                                  0.3132262, 0.5008952, 0.5531483   0.4557566
                                                       Ruletype
10203        T1wmorphology & T1wtexture & single-time-Enh & T2w
1342                           T1wdynamic & T1wmorphology & T2w
1593               T1wmorphology & dispersion & single-time-Enh
4833  T1wmorphology & T1wtexture & dispersion & single-time-Enh
1541                                 T1wmorphology & T1wtexture
1996                                 T1wmorphology & T1wtexture
883                                  T1wmorphology & T1wtexture
1338                                 T1wmorphology & T1wtexture
723                             T1wmorphology & single-time-Enh
5547              T1wmorphology & T1wtexture & dispersion & T2w
8534                  T1wdynamic & dispersion & single-time-Enh
6375                           T1wmorphology & dispersion & T2w
750                                  T1wmorphology & T1wtexture
2209                  T1wdynamic & dispersion & single-time-Enh
[1] 19
      len  freq                err
1342    3 0.053              0.034
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
2076    2 0.052              0.071
883     2 0.197              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
723     3 0.108              0.083
6375    4 0.202              0.092
726     2 0.215              0.093
750     2  0.25              0.094
5386    3 0.116              0.095
                                                                                                                                         condition
1342                                enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
12721                                                              uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                              irregularity low  & Energy ptime4 low to medium 
1996                                                                                              irregularity low  & Energy ptime4 low to medium 
2076                                                                             SER(in) low  & enhancement-variance increasing rate(in) very low 
883                                                                                               irregularity low  & Energy ptime4 low to medium 
515                 SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                              irregularity low  & Energy ptime4 low to medium 
723                                                         irregularity low to medium  & 3rd post-SE s10 low to very high  & 3rd post-SE s16 low 
6375                    irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                     SER(in) low  & T2w Energy  medium to high 
750                                                                            irregularity low  & Inverse difference moment ptime4 low to medium 
5386  SER(in) very low  & enhancement-variance increasing rate(in) very low to high  & enhancement-variance decreasing rate(rim) very low to high 
      pred                                       prob sumtempProb
1342    NC            0.4848342, 0.5546800, 0.6219773   0.5538305
12721   NC            0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                       0.4792871, 0.2735542   0.3764207
1996    NC                       0.5007048, 0.3009807   0.4008427
2076    NC                       0.6114913, 0.4662233   0.5388573
883     NC                       0.4694149, 0.2405692   0.3549921
515     NC 0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                       0.4991863, 0.2996197    0.399403
723     NC            0.5052503, 0.4047316, 0.3668232   0.4256017
6375    NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                       0.4881946, 0.2745578   0.3813762
750     NC                       0.4236079, 0.2408969   0.3322524
5386    NC            0.5192031, 0.5443866, 0.4975629   0.5203842
                                                       Ruletype
1342                           T1wdynamic & T1wmorphology & T2w
12721                   T1wdynamic & T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
1996                                 T1wmorphology & T1wtexture
2076                                                 T1wdynamic
883                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
723                             T1wmorphology & single-time-Enh
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
750                                  T1wmorphology & T1wtexture
5386                                                 T1wdynamic
[1] 21
      len  freq                err
9109    4 0.041                  0
1233    4 0.066              0.028
1342    3 0.053              0.034
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
2076    2 0.052              0.071
883     2 0.197              0.073
515     4 0.195              0.074
5316    5 0.168              0.075
1338    2 0.193              0.077
7750    2 0.091              0.082
723     3 0.108              0.083
5547    5  0.17              0.087
8424    3 0.185               0.09
6375    4 0.202              0.092
726     2 0.215              0.093
5658    3 0.177              0.093
750     2  0.25              0.094
722     4 0.136              0.095
                                                                                                                                                                              condition
9109                                                       Initial Uptake slope(in) low  & uptake average high  & 1st post-SE s15 low to very high  & 1st post-SE s16 low to very high 
1233                                   enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
1342                                                                     enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
12721                                                                                                   uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
2076                                                                                                                  SER(in) low  & enhancement-variance increasing rate(in) very low 
883                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
515                                                      SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
5316  Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
1338                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                                         circularity high & Sum average ptime1 low 
723                                                                                              irregularity low to medium  & 3rd post-SE s10 low to very high  & 3rd post-SE s16 low 
5547                         circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
8424                                                                               Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                         irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                          SER(in) low  & T2w Energy  medium to high 
5658                                                                                               SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                                                 irregularity low  & Inverse difference moment ptime4 low to medium 
722                                                              uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
      pred                                                  prob sumtempProb
9109    NC            0.6920559, 0.4833696, 0.7240647, 0.5912338    0.622681
1233    NC            0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
1342    NC                       0.4848342, 0.5546800, 0.6219773   0.5538305
12721   NC                       0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                  0.4792871, 0.2735542   0.3764207
1996    NC                                  0.5007048, 0.3009807   0.4008427
2076    NC                                  0.6114913, 0.4662233   0.5388573
883     NC                                  0.4694149, 0.2405692   0.3549921
515     NC            0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
5316    NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
1338    NC                                  0.4991863, 0.2996197    0.399403
7750    NC                                  0.5779883, 0.4539037    0.515946
723     NC                       0.5052503, 0.4047316, 0.3668232   0.4256017
5547    NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
8424    NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                                  0.4881946, 0.2745578   0.3813762
5658    NC                       0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                  0.4236079, 0.2408969   0.3322524
722     NC            0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                                       Ruletype
9109                               T1wdynamic & single-time-Enh
1233                    T1wdynamic & T1wmorphology & T1wtexture
1342                           T1wdynamic & T1wmorphology & T2w
12721                   T1wdynamic & T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
1996                                 T1wmorphology & T1wtexture
2076                                                 T1wdynamic
883                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
5316                         T1wdynamic & single-time-Enh & T2w
1338                                 T1wmorphology & T1wtexture
7750                                 T1wmorphology & T1wtexture
723                             T1wmorphology & single-time-Enh
5547              T1wmorphology & T1wtexture & dispersion & T2w
8424                                 T1wdynamic & T1wmorphology
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
5658                              T1wdynamic & T1wtexture & T2w
750                                  T1wmorphology & T1wtexture
722          T1wdynamic & T1wmorphology & single-time-Enh & T2w
[1] 23
      len  freq                err
12632   6 0.096 0.0580000000000001
1073    3  0.06 0.0610000000000001
7750    2 0.091              0.082
3902    2 0.124              0.088
6375    4 0.202              0.092
                                                                                                                                                                                                                  condition
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1073                                                                                                       circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
      pred                                                             prob sumtempProb
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1073    NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
7750    NC                                             0.5779883, 0.4539037    0.515946
3902    NC                                             0.5514954, 0.4847561   0.5181258
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
                                                    Ruletype
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1073                                     T1wmorphology & T2w
7750                              T1wmorphology & T1wtexture
3902                                     T1wmorphology & T2w
6375                        T1wmorphology & dispersion & T2w
[1] 26
     len  freq                err
9109   4 0.041                  0
4833   5 0.038              0.048
1541   2 0.162 0.0570000000000001
1996   2 0.184              0.069
2076   2 0.052              0.071
883    2 0.197              0.073
515    4 0.195              0.074
5316   5 0.168              0.075
1338   2 0.193              0.077
7111   4 0.094              0.077
5547   5  0.17              0.087
8534   4 0.103              0.088
6375   4 0.202              0.092
726    2 0.215              0.093
750    2  0.25              0.094
2209   3 0.076              0.095
722    4 0.136              0.095
                                                                                                                                                                             condition
9109                                                      Initial Uptake slope(in) low  & uptake average high  & 1st post-SE s15 low to very high  & 1st post-SE s16 low to very high 
4833          Max Radial gradient variance very low  & Sum average ptime3 medium to high  & Entropy ptime4 very low to high  & dispersion s18 low  & 2nd post-SE s15 very low to high 
1541                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
2076                                                                                                                 SER(in) low  & enhancement-variance increasing rate(in) very low 
883                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
515                                                     SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
5316 Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
1338                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
7111                              Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
5547                        circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
8534                                                                                 SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
6375                                                        irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                         SER(in) low  & T2w Energy  medium to high 
750                                                                                                                irregularity low  & Inverse difference moment ptime4 low to medium 
2209                                                                       enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
722                                                             uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
     pred                                                  prob sumtempProb
9109   NC            0.6920559, 0.4833696, 0.7240647, 0.5912338    0.622681
4833   NC 0.4289456, 0.6226296, 0.5110246, 0.5385531, 0.4767274   0.5155761
1541   NC                                  0.4792871, 0.2735542   0.3764207
1996   NC                                  0.5007048, 0.3009807   0.4008427
2076   NC                                  0.6114913, 0.4662233   0.5388573
883    NC                                  0.4694149, 0.2405692   0.3549921
515    NC            0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
1338   NC                                  0.4991863, 0.2996197    0.399403
7111   NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
5547   NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
8534   NC            0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726    NC                                  0.4881946, 0.2745578   0.3813762
750    NC                                  0.4236079, 0.2408969   0.3322524
2209   NC                       0.3132262, 0.5008952, 0.5531483   0.4557566
722    NC            0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
                                                      Ruletype
9109                              T1wdynamic & single-time-Enh
4833 T1wmorphology & T1wtexture & dispersion & single-time-Enh
1541                                T1wmorphology & T1wtexture
1996                                T1wmorphology & T1wtexture
2076                                                T1wdynamic
883                                 T1wmorphology & T1wtexture
515  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
5316                        T1wdynamic & single-time-Enh & T2w
1338                                T1wmorphology & T1wtexture
7111                 T1wdynamic & dispersion & single-time-Enh
5547             T1wmorphology & T1wtexture & dispersion & T2w
8534                 T1wdynamic & dispersion & single-time-Enh
6375                          T1wmorphology & dispersion & T2w
726                                           T1wdynamic & T2w
750                                 T1wmorphology & T1wtexture
2209                 T1wdynamic & dispersion & single-time-Enh
722         T1wdynamic & T1wmorphology & single-time-Enh & T2w
[1] 30
      len  freq                err
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
750     2  0.25              0.094
4809    3 0.095              0.096
                                                                                                                                                                                                                  condition
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                                                                             circularity high & Sum average ptime1 low 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
4809                                                                                                enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
      pred                                                             prob sumtempProb
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
4809    NC                                  0.4806604, 0.4285987, 0.5949174   0.5013922
                                                    Ruletype
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1996                              T1wmorphology & T1wtexture
883                               T1wmorphology & T1wtexture
1338                              T1wmorphology & T1wtexture
7750                              T1wmorphology & T1wtexture
5547           T1wmorphology & T1wtexture & dispersion & T2w
3902                                     T1wmorphology & T2w
6375                        T1wmorphology & dispersion & T2w
750                               T1wmorphology & T1wtexture
4809               T1wdynamic & T1wtexture & single-time-Enh
[1] 32
      len  freq   err
12721   3 0.127 0.043
5316    5 0.168 0.075
7111    4 0.094 0.077
8424    3 0.185  0.09
6375    4 0.202 0.092
5658    3 0.177 0.093
750     2  0.25 0.094
4809    3 0.095 0.096
                                                                                                                                                                              condition
12721                                                                                                   uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
5316  Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
7111                               Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
8424                                                                               Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                         irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
5658                                                                                               SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                                                 irregularity low  & Inverse difference moment ptime4 low to medium 
4809                                                            enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
      pred                                                  prob sumtempProb
12721   NC                       0.4677355, 0.5728223, 0.2321752   0.4242443
5316    NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
7111    NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
8424    NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
5658    NC                       0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                  0.4236079, 0.2408969   0.3322524
4809    NC                       0.4806604, 0.4285987, 0.5949174   0.5013922
                                       Ruletype
12721   T1wdynamic & T1wmorphology & dispersion
5316         T1wdynamic & single-time-Enh & T2w
7111  T1wdynamic & dispersion & single-time-Enh
8424                 T1wdynamic & T1wmorphology
6375           T1wmorphology & dispersion & T2w
5658              T1wdynamic & T1wtexture & T2w
750                  T1wmorphology & T1wtexture
4809  T1wdynamic & T1wtexture & single-time-Enh
[1] 33
      len  freq                err
10203   4 0.057              0.032
4971    4 0.055              0.033
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
715     3  0.08 0.0679999999999999
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
8424    3 0.185               0.09
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                               condition
10203 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
4971                                   min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
12721                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                    irregularity low  & Energy ptime4 low to medium 
715                                                                        uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
1996                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                    irregularity low  & Energy ptime4 low to medium 
8424                                                                Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                       prob sumtempProb
10203   NC 0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
4971    NC 0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
12721   NC            0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                       0.4792871, 0.2735542   0.3764207
715      C            0.5436611, 0.6397459, 0.5003391   0.5612487
1996    NC                       0.5007048, 0.3009807   0.4008427
883     NC                       0.4694149, 0.2405692   0.3549921
1338    NC                       0.4991863, 0.2996197    0.399403
8424    NC            0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                       0.4236079, 0.2408969   0.3322524
                                                  Ruletype
10203   T1wmorphology & T1wtexture & single-time-Enh & T2w
4971  T1wdynamic & T1wmorphology & T1wtexture & dispersion
12721              T1wdynamic & T1wmorphology & dispersion
1541                            T1wmorphology & T1wtexture
715              T1wdynamic & T1wtexture & single-time-Enh
1996                            T1wmorphology & T1wtexture
883                             T1wmorphology & T1wtexture
1338                            T1wmorphology & T1wtexture
8424                            T1wdynamic & T1wmorphology
6375                      T1wmorphology & dispersion & T2w
750                             T1wmorphology & T1wtexture
[1] 39
      len  freq                err
10203   4 0.057              0.032
7103    2 0.044              0.042
1578    3 0.057 0.0649999999999999
3314    5  0.07              0.077
1977    2 0.045               0.08
8534    4 0.103              0.088
6375    4 0.202              0.092
2931    3 0.038              0.095
                                                                                                                                                               condition
10203 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
7103                                                                                               irregularity very low  & T2w radial gradient variance medium to high 
1578                                                irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
3314                     SER(rim) very low to high  & dispersion s14 very low  & dispersion s15 very low to high  & 1st post-SE s4 low to very high  & T2w Entropy  low 
1977                                                                                                                      dispersion s5 very low  & 1st post-SE s11 low 
8534                                                                   SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
6375                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
2931                                                           Rate of signal increase(rim) very low  & dispersion s17 low to very high  & T2w mean SI low to very high 
      pred                                                  prob sumtempProb
10203   NC            0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
7103    NC                                  0.5748740, 0.4661645   0.5205192
1578    NC                       0.5914519, 0.4717202, 0.4894149    0.517529
3314    NC 0.4499065, 0.4352475, 0.4759258, 0.4943733, 0.5792785   0.4869463
1977    NC                                  0.4735013, 0.6948949   0.5841981
8534    NC            0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
2931    NC                       0.5582364, 0.5165599, 0.5742535   0.5496832
                                                Ruletype
10203 T1wmorphology & T1wtexture & single-time-Enh & T2w
7103                                 T1wmorphology & T2w
1578                    T1wmorphology & T1wtexture & T2w
3314     T1wdynamic & dispersion & single-time-Enh & T2w
1977                        dispersion & single-time-Enh
8534           T1wdynamic & dispersion & single-time-Enh
6375                    T1wmorphology & dispersion & T2w
2931                       T1wdynamic & dispersion & T2w
[1] 44
      len  freq                err
1233    4 0.066              0.028
10158   2 0.048              0.038
7103    2 0.044              0.042
1073    3  0.06 0.0610000000000001
1578    3 0.057 0.0649999999999999
7111    4 0.094              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
5658    3 0.177              0.093
750     2  0.25              0.094
5386    3 0.116              0.095
9527    2 0.057              0.097
                                                                                                                                                       condition
1233            enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
10158                                                                                                 Energy ptime1 very low to high  & Entropy ptime2 very low 
7103                                                                                       irregularity very low  & T2w radial gradient variance medium to high 
1073                                            circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1578                                        irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
7111        Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
7750                                                                                                                  circularity high & Sum average ptime1 low 
5547  circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                             circularity high & T2w radial gradient variance medium to high 
6375                                  irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
5658                                                                        SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
750                                                                                          irregularity low  & Inverse difference moment ptime4 low to medium 
5386                SER(in) very low  & enhancement-variance increasing rate(in) very low to high  & enhancement-variance decreasing rate(rim) very low to high 
9527                                             max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
      pred                                                  prob sumtempProb
1233    NC            0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
10158   NC                                  0.5127610, 0.4835754   0.4981682
7103    NC                                  0.5748740, 0.4661645   0.5205192
1073    NC                       0.6781025, 0.4757408, 0.5921319   0.5819917
1578    NC                       0.5914519, 0.4717202, 0.4894149    0.517529
7111    NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
7750    NC                                  0.5779883, 0.4539037    0.515946
5547    NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                  0.5514954, 0.4847561   0.5181258
6375    NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
5658    NC                       0.5505793, 0.3818076, 0.5017071   0.4780313
750     NC                                  0.4236079, 0.2408969   0.3322524
5386    NC                       0.5192031, 0.5443866, 0.4975629   0.5203842
9527    NC                                  0.6236805, 0.4747261   0.5492033
                                           Ruletype
1233        T1wdynamic & T1wmorphology & T1wtexture
10158                                    T1wtexture
7103                            T1wmorphology & T2w
1073                            T1wmorphology & T2w
1578               T1wmorphology & T1wtexture & T2w
7111      T1wdynamic & dispersion & single-time-Enh
7750                     T1wmorphology & T1wtexture
5547  T1wmorphology & T1wtexture & dispersion & T2w
3902                            T1wmorphology & T2w
6375               T1wmorphology & dispersion & T2w
5658                  T1wdynamic & T1wtexture & T2w
750                      T1wmorphology & T1wtexture
5386                                     T1wdynamic
9527                     T1wdynamic & T1wmorphology
[1] 50
      len  freq                err
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
6591    4 0.031 0.0590000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
515     4 0.195              0.074
1338    2 0.193              0.077
8424    3 0.185               0.09
6375    4 0.202              0.092
726     2 0.215              0.093
750     2  0.25              0.094
                                                                                                                                                   condition
12721                                                                        uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                        irregularity low  & Energy ptime4 low to medium 
6591  Rate of signal increase(in) low  & Rate of signal increase(in) medium to high  & 1st post-SE s15 low to very high  & last post-SE s8 low to very high 
1996                                                                                                        irregularity low  & Energy ptime4 low to medium 
883                                                                                                         irregularity low  & Energy ptime4 low to medium 
515                           SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
1338                                                                                                        irregularity low  & Energy ptime4 low to medium 
8424                                                    Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                              irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                               SER(in) low  & T2w Energy  medium to high 
750                                                                                      irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                       prob sumtempProb
12721   NC            0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                       0.4792871, 0.2735542   0.3764207
6591    NC 0.5271721, 0.3903370, 0.5439983, 0.4908403    0.488087
1996    NC                       0.5007048, 0.3009807   0.4008427
883     NC                       0.4694149, 0.2405692   0.3549921
515     NC 0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
1338    NC                       0.4991863, 0.2996197    0.399403
8424    NC            0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726     NC                       0.4881946, 0.2745578   0.3813762
750     NC                       0.4236079, 0.2408969   0.3322524
                                                       Ruletype
12721                   T1wdynamic & T1wmorphology & dispersion
1541                                 T1wmorphology & T1wtexture
6591                               T1wdynamic & single-time-Enh
1996                                 T1wmorphology & T1wtexture
883                                  T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
1338                                 T1wmorphology & T1wtexture
8424                                 T1wdynamic & T1wmorphology
6375                           T1wmorphology & dispersion & T2w
726                                            T1wdynamic & T2w
750                                  T1wmorphology & T1wtexture
[1] 51
      len  freq                err
1342    3 0.053              0.034
1593    6 0.047              0.038
10158   2 0.048              0.038
12721   3 0.127              0.043
6591    4 0.031 0.0590000000000001
3902    2 0.124              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
                                                                                                                                                                               condition
1342                                                                      enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
1593  circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
10158                                                                                                                         Energy ptime1 very low to high  & Entropy ptime2 very low 
12721                                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
6591                              Rate of signal increase(in) low  & Rate of signal increase(in) medium to high  & 1st post-SE s15 low to very high  & last post-SE s8 low to very high 
3902                                                                                                                     circularity high & T2w radial gradient variance medium to high 
8424                                                                                Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
      pred                                                             prob sumtempProb
1342    NC                                  0.4848342, 0.5546800, 0.6219773   0.5538305
1593    NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
10158   NC                                             0.5127610, 0.4835754   0.4981682
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
6591    NC                       0.5271721, 0.3903370, 0.5439983, 0.4908403    0.488087
3902    NC                                             0.5514954, 0.4847561   0.5181258
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
                                          Ruletype
1342              T1wdynamic & T1wmorphology & T2w
1593  T1wmorphology & dispersion & single-time-Enh
10158                                   T1wtexture
12721      T1wdynamic & T1wmorphology & dispersion
6591                  T1wdynamic & single-time-Enh
3902                           T1wmorphology & T2w
8424                    T1wdynamic & T1wmorphology
6375              T1wmorphology & dispersion & T2w
[1] 54
     len  freq                err
4971   4 0.055              0.033
1541   2 0.162 0.0570000000000001
1996   2 0.184              0.069
883    2 0.197              0.073
515    4 0.195              0.074
5316   5 0.168              0.075
1338   2 0.193              0.077
7111   4 0.094              0.077
8424   3 0.185               0.09
6375   4 0.202              0.092
726    2 0.215              0.093
750    2  0.25              0.094
4809   3 0.095              0.096
                                                                                                                                                                             condition
4971                                                 min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
1541                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
515                                                     SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
5316 Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
1338                                                                                                                                  irregularity low  & Energy ptime4 low to medium 
7111                              Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
8424                                                                              Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                        irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
726                                                                                                                                         SER(in) low  & T2w Energy  medium to high 
750                                                                                                                irregularity low  & Inverse difference moment ptime4 low to medium 
4809                                                           enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
     pred                                                  prob sumtempProb
4971   NC            0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
1541   NC                                  0.4792871, 0.2735542   0.3764207
1996   NC                                  0.5007048, 0.3009807   0.4008427
883    NC                                  0.4694149, 0.2405692   0.3549921
515    NC            0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
1338   NC                                  0.4991863, 0.2996197    0.399403
7111   NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
8424   NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
726    NC                                  0.4881946, 0.2745578   0.3813762
750    NC                                  0.4236079, 0.2408969   0.3322524
4809   NC                       0.4806604, 0.4285987, 0.5949174   0.5013922
                                                      Ruletype
4971      T1wdynamic & T1wmorphology & T1wtexture & dispersion
1541                                T1wmorphology & T1wtexture
1996                                T1wmorphology & T1wtexture
883                                 T1wmorphology & T1wtexture
515  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
5316                        T1wdynamic & single-time-Enh & T2w
1338                                T1wmorphology & T1wtexture
7111                 T1wdynamic & dispersion & single-time-Enh
8424                                T1wdynamic & T1wmorphology
6375                          T1wmorphology & dispersion & T2w
726                                           T1wdynamic & T2w
750                                 T1wmorphology & T1wtexture
4809                 T1wdynamic & T1wtexture & single-time-Enh
[1] 55
      len  freq                err
1593    6 0.047              0.038
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7111    4 0.094              0.077
7750    2 0.091              0.082
5547    5  0.17              0.087
8534    4 0.103              0.088
8424    3 0.185               0.09
6375    4 0.202              0.092
750     2  0.25              0.094
9527    2 0.057              0.097
                                                                                                                                                                               condition
1593  circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
12721                                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1996                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
7111                                Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
7750                                                                                                                                          circularity high & Sum average ptime1 low 
5547                          circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
8534                                                                                   SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
8424                                                                                Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
9527                                                                     max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
      pred                                                             prob sumtempProb
1593    NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                             0.4792871, 0.2735542   0.3764207
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7111    NC                       0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
7750    NC                                             0.5779883, 0.4539037    0.515946
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
8534    NC                       0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
8424    NC                                  0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
9527    NC                                             0.6236805, 0.4747261   0.5492033
                                           Ruletype
1593   T1wmorphology & dispersion & single-time-Enh
12721       T1wdynamic & T1wmorphology & dispersion
1541                     T1wmorphology & T1wtexture
1996                     T1wmorphology & T1wtexture
883                      T1wmorphology & T1wtexture
1338                     T1wmorphology & T1wtexture
7111      T1wdynamic & dispersion & single-time-Enh
7750                     T1wmorphology & T1wtexture
5547  T1wmorphology & T1wtexture & dispersion & T2w
8534      T1wdynamic & dispersion & single-time-Enh
8424                     T1wdynamic & T1wmorphology
6375               T1wmorphology & dispersion & T2w
750                      T1wmorphology & T1wtexture
9527                     T1wdynamic & T1wmorphology
[1] 57
      len  freq                err
1593    6 0.047              0.038
12721   3 0.127              0.043
2822    2 0.071              0.051
1541    2 0.162 0.0570000000000001
1073    3  0.06 0.0610000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                                               condition
1593  circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
12721                                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
2822                                                                                                                      3rd post-SE s10 low to very high  & last post-SE s11 very low 
1541                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
1073                                                                    circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
1996                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
5547                          circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                     circularity high & T2w radial gradient variance medium to high 
6375                                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                                             prob sumtempProb
1593    NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
2822    NC                                             0.5144900, 0.4830707   0.4987804
1541    NC                                             0.4792871, 0.2735542   0.3764207
1073    NC                                  0.6781025, 0.4757408, 0.5921319   0.5819917
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
                                           Ruletype
1593   T1wmorphology & dispersion & single-time-Enh
12721       T1wdynamic & T1wmorphology & dispersion
2822                                single-time-Enh
1541                     T1wmorphology & T1wtexture
1073                            T1wmorphology & T2w
1996                     T1wmorphology & T1wtexture
883                      T1wmorphology & T1wtexture
1338                     T1wmorphology & T1wtexture
5547  T1wmorphology & T1wtexture & dispersion & T2w
3902                            T1wmorphology & T2w
6375               T1wmorphology & dispersion & T2w
750                      T1wmorphology & T1wtexture
[1] 64
      len  freq                err
1233    4 0.066              0.028
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7111    4 0.094              0.077
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                                                                                  condition
1233                                                                       enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7111                                                                   Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                                             prob sumtempProb
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7111    NC                       0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
                                                    Ruletype
1233                 T1wdynamic & T1wmorphology & T1wtexture
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1996                              T1wmorphology & T1wtexture
883                               T1wmorphology & T1wtexture
1338                              T1wmorphology & T1wtexture
7111               T1wdynamic & dispersion & single-time-Enh
5547           T1wmorphology & T1wtexture & dispersion & T2w
3902                                     T1wmorphology & T2w
6375                        T1wmorphology & dispersion & T2w
750                               T1wmorphology & T1wtexture
```

```r
print(summary(LMSIRfeatureflag))
```

```
   Mode   FALSE    TRUE    NA's 
logical      45      25       0 
```

```r
### display
lesion_id = 10
idx = c(1:nrow(FI))[c(FI$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
   lesion_id          lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
10        10 10_0102_4755778.vtk          0102 1974-05-19 00:00:00.000000               5142
   exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
10           4755778 2009-03-19 00:00:00.000000     Benign by pathology
   cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
10                      High Risk                     1                        0
   exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
10              Right                 2708 2009-04-16 00:00:00.000000              Right
   proc_proc_source_int proc_proc_guid_int       proc_proc_tp_int proc_lesion_comments_txt
10            Radiology                MRI Vacuum Assisted Biopsy                     None
   find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
10            III                       N/A                    Washout      4        massB
   lesion_diagnosis
10      FIBROCYSTIC
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
AP Report: 
     (NOTE)
     
     SURGICAL  PATHOLOGY  REPORT

     Encounter #: 4050499R
     Specimen #: S09-10059
     
     CLINICAL INFORMATION
     MRI enhancement. Rule out cancer.
     
     SPECIMENS SUBMITTED
     Right breast biopsies
     
     
     
     DIAGNOSIS
     Right breast, core biopsy:
     - FIBROCYSTIC AND COLUMNAR CELL CHANGES
     - FIBROEPITHELIAL LESION
     
     COMMENT:
     Sections show a benign fibroepithelial lesion with fibrocystic changes
     as well as columnar cell metaplasia and hyperplasia.  There are
     fragments of a fibroadenoma.  There is no evidence of atypia or
     malignancy.
     
     
     
     
     MACROSCOPIC DESCRIPTION
     The specimen container is labeled with the patient name and 'right
     breast samples'. The accompanying requisition matches the container' s
     label. The specimen consists of approx..2 cc's of soft tan and fatty
     tissue pieces which are submitted in toto in blocks 1-4.
     AMM  Dictated 4/17/2009
     
     
     
     
     
     Judit Zubovits, MD, FRCPC
     Report Electronically Signed
     4/20/2009 15:43
```

```r
print(mypresentRules(eachFITopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
10203   4 0.057              0.032
4971    4 0.055              0.033
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
715     3  0.08 0.0679999999999999
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
8424    3 0.185               0.09
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                               condition
10203 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
4971                                   min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
12721                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                    irregularity low  & Energy ptime4 low to medium 
715                                                                        uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
1996                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                    irregularity low  & Energy ptime4 low to medium 
8424                                                                Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                       prob sumtempProb
10203   NC 0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
4971    NC 0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
12721   NC            0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                       0.4792871, 0.2735542   0.3764207
715      C            0.5436611, 0.6397459, 0.5003391   0.5612487
1996    NC                       0.5007048, 0.3009807   0.4008427
883     NC                       0.4694149, 0.2405692   0.3549921
1338    NC                       0.4991863, 0.2996197    0.399403
8424    NC            0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                       0.4236079, 0.2408969   0.3322524
                                                  Ruletype
10203   T1wmorphology & T1wtexture & single-time-Enh & T2w
4971  T1wdynamic & T1wmorphology & T1wtexture & dispersion
12721              T1wdynamic & T1wmorphology & dispersion
1541                            T1wmorphology & T1wtexture
715              T1wdynamic & T1wtexture & single-time-Enh
1996                            T1wmorphology & T1wtexture
883                             T1wmorphology & T1wtexture
1338                            T1wmorphology & T1wtexture
8424                            T1wdynamic & T1wmorphology
6375                      T1wmorphology & dispersion & T2w
750                             T1wmorphology & T1wtexture
```

```r
# print(FI[FI$allTestinfo.lesion_id==241,])
imgFI <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgFI)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
## with LMSIR
if (length(idLMSIR) > 0) {
    for (i in 1:length(idLMSIR)) {
        idx = idLMSIR[i]
        lesion_id = FI[idx, ]$allTestinfo.lesion_id
        print(lesion_id)
    }
}
```

```
[1] 228
[1] 301
[1] 503
[1] 66
[1] 74
[1] 288
[1] 383
[1] 411
[1] 412
[1] 93
[1] 213
[1] 311
[1] 419
[1] 135
[1] 396
[1] 613
[1] 10
[1] 76
[1] 166
[1] 627
[1] 26
[1] 48
[1] 263
[1] 267
[1] 29
```

```r
idx = c(1:nrow(FI))[c(FI$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
   lesion_id          lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
29        29 29_0168_5240535.vtk          0168 1958-03-27 00:00:00.000000              11164
   exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
29           5240535 2010-04-12 00:00:00.000000     Benign by pathology
   cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
29                          Other                     1                        0
   exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
29               Left                 2798 2010-04-27 00:00:00.000000               Left
   proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int
29            Radiology                 US Core Needle Biopsy
                                                                                                                                                                                                                       proc_lesion_comments_txt
29 CLINICAL INDICATION: Intraductal papillary lesion on recent\nbreast ultrasound 11 o'clock subareolar position left breast,\nrecommended for surgical excision biopsy. Recommended for MRI\nevaluation. Bilateral reduction mammoplasty 1990.
   find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
29             II                     Rapid                    Plateau      4        massB
   lesion_diagnosis
29      FIBROCYSTIC
```

```r
print(mypresentRules(eachFITopRules[[idx]], colnames(X), fnames))
```

```
      len  freq                err
1233    4 0.066              0.028
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
7111    4 0.094              0.077
5547    5  0.17              0.087
3902    2 0.124              0.088
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                                                                                  condition
1233                                                                       enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
1996                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
883                                                                                                                                                                        irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
7111                                                                   Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
5547                                                             circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
3902                                                                                                                                                        circularity high & T2w radial gradient variance medium to high 
6375                                                                                             irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                                                                     irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                                             prob sumtempProb
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
1996    NC                                             0.5007048, 0.3009807   0.4008427
883     NC                                             0.4694149, 0.2405692   0.3549921
1338    NC                                             0.4991863, 0.2996197    0.399403
7111    NC                       0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
5547    NC            0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
3902    NC                                             0.5514954, 0.4847561   0.5181258
6375    NC                       0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                                             0.4236079, 0.2408969   0.3322524
                                                    Ruletype
1233                 T1wdynamic & T1wmorphology & T1wtexture
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
1996                              T1wmorphology & T1wtexture
883                               T1wmorphology & T1wtexture
1338                              T1wmorphology & T1wtexture
7111               T1wdynamic & dispersion & single-time-Enh
5547           T1wmorphology & T1wtexture & dispersion & T2w
3902                                     T1wmorphology & T2w
6375                        T1wmorphology & dispersion & T2w
750                               T1wmorphology & T1wtexture
```


![](Z:/Cristina/Section2/papernew_notes/images/10_0102_4755778_massB_Hyperintense.png)

# Hand picked cases for illustration:for IDC

```r
print(topcasesIDC[topcasesIDC$lesion_id == 9, ])
```

```
    43 38 45 21 44 lesion_id lesion_label lesion_diagnosis
114  0  0  1  1  0         9     nonmassM   InvasiveDuctal
```

```r
lesion_id = 9
idx = c(1:nrow(IDC))[c(IDC$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
  lesion_id         lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
9         9 9_0093_7156466.vtk          0093 1976-07-24 00:00:00.000000              28100
  exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
9           7156466 2012-10-16 00:00:00.000000                 Unknown
  cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
9                          BRCA2                     0                        1
  exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
9               Left                 2692 2012-11-28 00:00:00.000000               Left
  proc_proc_source_int proc_proc_guid_int   proc_proc_tp_int proc_lesion_comments_txt
9            Radiology             Stereo Core Needle Biopsy                     None
  find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
9           None                       N/A                 Persistent      4     nonmassM
  lesion_diagnosis
9   InvasiveDuctal
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
Nov 28, 2012

AP Report: 
     (NOTE)
     
     SURGICAL  PATHOLOGY  REPORT

     Encounter #: 14925508R
     Submitting Physician: MESSNER,SANDRA
     Specimen #: S12-30804
     
     CLINICAL INFORMATION
     Calcifications.
     
     SPECIMENS SUBMITTED
     Left breast stereo core biopsy x4
     
     
     
     DIAGNOSIS
     Left breast, core biopsy: INVASIVE AND IN SITU CARCINOMA, SEE NOTE.
     
     Note: Sections show predominantly carcinoma in situ with intermediate
     nuclear grade, solid and cribriform type, with moderate necrosis and
     calcifications, as well as focal stromal invasion. The invasive
     carcinoma grows as single cell files and is negative for e-cadherin,
     which is in keeping with invasive lobular carcinoma, classical type in
     this material. The in situ carcinoma shows aberrant staining with
     e-cadherin ranging from absent to circumferential membranous,
     therefore would be best classified as mammary carcinoma in situ with
     mixed ductal and lobular features.
     
     
     
     
     MACROSCOPIC DESCRIPTION
     The specimen container is labeled with the patient name and 'left
     breast core biopsy'. The accompanying requisition matches the
     container's label. The specimen consists of 4 cores of tan and fatty
     tissue, each with a diameter of 0.1 cm, ranging from 1.0 to 1.2 cm in
     length. Submitted in toto in one block
     JP
     Dictated 12/03/2012
     
     
     
     
     
     Elzbieta Slodkowska, MD
     Report Electronically Signed
     2012/12/06 14:53
```

```r
rulesaIDC = mypresentRules(eachIDCTopRules[[idx]], colnames(X), fnames)
print(rulesaIDC)
```

```
     len  freq                err
1193   3  0.03                  0
68     3 0.111              0.049
3868   3 0.031 0.0590000000000001
310    2 0.144              0.076
300    3 0.115              0.081
1063   3 0.067              0.081
3524   2 0.043              0.087
588    2 0.121               0.09
615    2  0.14              0.091
184    2 0.156              0.094
                                                                                            condition
1193 Time-to-peak(rim) medium to high  & max Radial gradient low  & T2w radial gradient variance low 
68               Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high 
3868              dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
310                                                      SER(in) medium to high  & irregularity high 
300                     SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
1063      irregularity high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
3524                   change in Variance of spatial Margin Gradient low  & dispersion s17 very high 
588                                                                SER(in) high  & irregularity high 
615                                              irregularity high  & 1st post-SE s15 medium to high 
184                                                      SER(in) medium to high  & irregularity high 
     pred                            prob sumtempProb
1193   NC 0.4662758, 0.3631702, 0.5977509   0.4757323
68      C 0.4479334, 0.7561858, 0.6379134   0.6140108
3868    C 0.4853810, 0.6962119, 0.4631950   0.5482627
310     C            0.5006988, 0.6626664   0.5816826
300     C 0.5949279, 0.7518615, 0.4569649   0.6012514
1063    C 0.5128461, 0.6452170, 0.4944342   0.5508324
3524    C            0.3885357, 0.5447962    0.466666
588     C            0.4565346, 0.6397499   0.5481423
615     C            0.4819968, 0.7471799   0.6145883
184     C            0.7579719, 0.4714506   0.6147112
                                         Ruletype
1193             T1wdynamic & T1wmorphology & T2w
68                     T1wdynamic & T1wmorphology
3868           dispersion & single-time-Enh & T2w
310                    T1wdynamic & T1wmorphology
300                    T1wdynamic & T1wmorphology
1063 T1wmorphology & T1wtexture & single-time-Enh
3524                   T1wmorphology & dispersion
588                    T1wdynamic & T1wmorphology
615               T1wmorphology & single-time-Enh
184                    T1wdynamic & T1wmorphology
```

```r
# print(FI[FI$allTestinfo.lesion_id==241,])
imgaIDC <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgaIDC)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


# Hand picked cases for illustration:for FI

```r
print(topcasesFI[topcasesFI$lesion_id == 10, ])  ## 301, 66
```

```
   42 35 39 38 40 25 27 24 22 14 26 11 41 29 37 32 46 6 36 3 15 17 23 47 45 1 8 28 33 49 50 13
33  1  0  1  1  0  0  1  1  1  1  0  1  0  0  0  0  0 0  0 0  0  0  0  0  0 0 0  0  0  0  0  0
   43 4 5 7 9 12 18 48 16 31 34 10 44 2 19 20 21 30 lesion_id lesion_label lesion_diagnosis
33  0 1 1 0 0  0  0  0  0  0  0  0  0 0  0  0  1  0        10        massB      FIBROCYSTIC
```

```r
lesion_id = 10
idx = c(1:nrow(FI))[c(FI$allTestinfo.lesion_id == lesion_id)]
infoidx = allTestinfo[allTestinfo$lesion_id == lesion_id, -c(20)]
print(infoidx)
```

```
   lesion_id          lesionfile cad_pt_no_txt         anony_dob_datetime exam_img_dicom_txt
10        10 10_0102_4755778.vtk          0102 1974-05-19 00:00:00.000000               5142
   exam_a_number_txt           exam_dt_datetime exam_mri_cad_status_txt
10           4755778 2009-03-19 00:00:00.000000     Benign by pathology
   cad_latest_mutation_status_int exam_find_mri_mass_yn exam_find_mri_nonmass_yn
10                      High Risk                     1                        0
   exam_find_side_int proc_pt_procedure_id      proc_proc_dt_datetime proc_proc_side_int
10              Right                 2708 2009-04-16 00:00:00.000000              Right
   proc_proc_source_int proc_proc_guid_int       proc_proc_tp_int proc_lesion_comments_txt
10            Radiology                MRI Vacuum Assisted Biopsy                     None
   find_curve_int find_mri_dce_init_enh_int find_mri_dce_delay_enh_int BIRADS lesion_label
10            III                       N/A                    Washout      4        massB
   lesion_diagnosis
10      FIBROCYSTIC
```

```r
cat(allTestinfo[allTestinfo$lesion_id == lesion_id, c(20)])
```

```
AP Report: 
     (NOTE)
     
     SURGICAL  PATHOLOGY  REPORT

     Encounter #: 4050499R
     Specimen #: S09-10059
     
     CLINICAL INFORMATION
     MRI enhancement. Rule out cancer.
     
     SPECIMENS SUBMITTED
     Right breast biopsies
     
     
     
     DIAGNOSIS
     Right breast, core biopsy:
     - FIBROCYSTIC AND COLUMNAR CELL CHANGES
     - FIBROEPITHELIAL LESION
     
     COMMENT:
     Sections show a benign fibroepithelial lesion with fibrocystic changes
     as well as columnar cell metaplasia and hyperplasia.  There are
     fragments of a fibroadenoma.  There is no evidence of atypia or
     malignancy.
     
     
     
     
     MACROSCOPIC DESCRIPTION
     The specimen container is labeled with the patient name and 'right
     breast samples'. The accompanying requisition matches the container' s
     label. The specimen consists of approx..2 cc's of soft tan and fatty
     tissue pieces which are submitted in toto in blocks 1-4.
     AMM  Dictated 4/17/2009
     
     
     
     
     
     Judit Zubovits, MD, FRCPC
     Report Electronically Signed
     4/20/2009 15:43
```

```r
rulesaFI = mypresentRules(eachFITopRules[[idx]], colnames(X), fnames)
print(rulesaFI)
```

```
      len  freq                err
10203   4 0.057              0.032
4971    4 0.055              0.033
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
715     3  0.08 0.0679999999999999
1996    2 0.184              0.069
883     2 0.197              0.073
1338    2 0.193              0.077
8424    3 0.185               0.09
6375    4 0.202              0.092
750     2  0.25              0.094
                                                                                                                                                               condition
10203 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
4971                                   min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
12721                                                                                    uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                    irregularity low  & Energy ptime4 low to medium 
715                                                                        uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
1996                                                                                                                    irregularity low  & Energy ptime4 low to medium 
883                                                                                                                     irregularity low  & Energy ptime4 low to medium 
1338                                                                                                                    irregularity low  & Energy ptime4 low to medium 
8424                                                                Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
6375                                          irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
750                                                                                                  irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                       prob sumtempProb
10203   NC 0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
4971    NC 0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
12721   NC            0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                       0.4792871, 0.2735542   0.3764207
715      C            0.5436611, 0.6397459, 0.5003391   0.5612487
1996    NC                       0.5007048, 0.3009807   0.4008427
883     NC                       0.4694149, 0.2405692   0.3549921
1338    NC                       0.4991863, 0.2996197    0.399403
8424    NC            0.2901904, 0.4699588, 0.5212986   0.4271493
6375    NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
750     NC                       0.4236079, 0.2408969   0.3322524
                                                  Ruletype
10203   T1wmorphology & T1wtexture & single-time-Enh & T2w
4971  T1wdynamic & T1wmorphology & T1wtexture & dispersion
12721              T1wdynamic & T1wmorphology & dispersion
1541                            T1wmorphology & T1wtexture
715              T1wdynamic & T1wtexture & single-time-Enh
1996                            T1wmorphology & T1wtexture
883                             T1wmorphology & T1wtexture
1338                            T1wmorphology & T1wtexture
8424                            T1wdynamic & T1wmorphology
6375                      T1wmorphology & dispersion & T2w
750                             T1wmorphology & T1wtexture
```

```r
# print(FI[FI$allTestinfo.lesion_id==241,])
imgaIF <- readPNG(paste0("Z:/Cristina/Section2/papernew_notes/images/", as.character(lesion_id), 
    "_", infoidx$cad_pt_no_txt, "_", infoidx$exam_a_number_txt, "_", infoidx$lesion_label, "_", 
    lesioninfo[lesioninfo$lesion_id == lesion_id, "find_t2_signal_int"], ".png"))
grid.raster(imgaIF)
```

![](illustrateRulextraction_Section2_files/figure-html/unnamed-chunk-11-1.png)<!-- -->
