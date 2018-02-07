# OBSP - Sensitivity - Analysis revision 2




# OBSP High risk screening
Screening outcome summary, from July 2011 to June 2012: 

* 6,863 women aged 29-69 referred and registered in the OBSP High Risk Screening Program

* stratified in 2 risks groups: Category A: Known risk 964 (14.0%) and Category B: Referred to genetic assessment 5899 (86.0%) Completed

* From Category B: 3349 (56.8%) Completed genetic counselling only (1189 were elegible for OBSP high risk screening), and 1852 (31.4%) Completed genetic counselling and testing (440 were elegible for OBSP high risk screening)

* For a total of 2359 elegible woman. 

* 2207 Women had at least an MRI/ultrasound screen

* 611 had an Abnormal screens (abnormal call rate = 27.7%)

* 554 (90.7% of 611) Women had a final result for diagnosis

* 35 women had breast cancer (positive predictive value = 6.3%)

* 27 of 2207 were Invasive (cancer detection rate = 12.6 per 1,000) 

* 8 of 2207 Ductal carcinoma in situ (cancer detection rate = 3.7 per 1,000)


After running 10 folds of cross-validation (cv), below is the AUC distributions achieved on held-out test sets:

## Performance of combined T1w+T2w vs. only T1w classifiers 
(all lesions -pooled data across cv-folds and plot pooled results)

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7385  0.7406  0.7471  0.7961  0.8410  0.9131 
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7562  0.7760  0.7807  0.8046  0.8038  0.9061 
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7691  0.7802  0.8246  0.8313  0.8562  0.9263 
```

![](OBSP_Sensitivity_vsCAD_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```
[1] "Results for T1w-only features classifier:"
```

```
Area under the curve: 0.7992
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4156007 0.7185    0.7731  0.8277 0.6546     0.701  0.7474
```

```
[1] "Results for T2w + T2w_SI features classifier:"
```

```
Area under the curve: 0.8049
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4592056 0.7185    0.7731  0.8277 0.7088    0.7526  0.7938
```

```
[1] "Results for T2w T2wText + predictiveLMSIR features classifier:"
```

```
Area under the curve: 0.8333
95% CI (2000 stratified bootstrap replicates):
 thresholds sp.low sp.median sp.high se.low se.median se.high
  0.4269815 0.7437    0.7983  0.8487 0.6881    0.7345  0.7784
```

![](OBSP_Sensitivity_vsCAD_files/figure-html/unnamed-chunk-1-2.png)<!-- -->
 

## Confusion matrices of performance at optimal operating point

```r
library(caret)

print("Results for T1w-only features classifier (group 1):")
```

```
[1] "Results for T1w-only features classifier (group 1):"
```

```r
print(p1$best_thr$sensitivity)
```

```
                      2.5%       50%     97.5%
0.41560071718722 0.6546392 0.7010309 0.7474227
```

```r
th1 = as.numeric(row.names(p1$best_thr$sensitivity))
confusionMatrix(ifelse(perfall_imgT1_total$C >= th1, "C", "NC"), perfall_imgT1_total$obs, positive = "C")
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  184 116
        NC  54 272
                                          
               Accuracy : 0.7284          
                 95% CI : (0.6918, 0.7629)
    No Information Rate : 0.6198          
    P-Value [Acc > NIR] : 6.495e-09       
                                          
                  Kappa : 0.4514          
 Mcnemar's Test P-Value : 2.890e-06       
                                          
            Sensitivity : 0.7731          
            Specificity : 0.7010          
         Pos Pred Value : 0.6133          
         Neg Pred Value : 0.8344          
             Prevalence : 0.3802          
         Detection Rate : 0.2939          
   Detection Prevalence : 0.4792          
      Balanced Accuracy : 0.7371          
                                          
       'Positive' Class : C               
                                          
```

```r
print("Results for T2w text + SI features classifier (group 2):")
```

```
[1] "Results for T2w text + SI features classifier (group 2):"
```

```r
print(p2$best_thr$sensitivity)
```

```
                       2.5%       50%     97.5%
0.459205558580196 0.7087629 0.7525773 0.7938144
```

```r
th2 = as.numeric(row.names(p2$best_thr$sensitivity))
confusionMatrix(ifelse(perfall_imgT1T2_total$C >= th2, "C", "NC"), perfall_imgT1T2_total$obs, positive = "C")
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  184  96
        NC  54 292
                                         
               Accuracy : 0.7604         
                 95% CI : (0.725, 0.7933)
    No Information Rate : 0.6198         
    P-Value [Acc > NIR] : 4.967e-14      
                                         
                  Kappa : 0.5083         
 Mcnemar's Test P-Value : 0.000815       
                                         
            Sensitivity : 0.7731         
            Specificity : 0.7526         
         Pos Pred Value : 0.6571         
         Neg Pred Value : 0.8439         
             Prevalence : 0.3802         
         Detection Rate : 0.2939         
   Detection Prevalence : 0.4473         
      Balanced Accuracy : 0.7628         
                                         
       'Positive' Class : C              
                                         
```

```r
print("Results for T2w predicted LMSIR features classifier (group 3):")
```

```
[1] "Results for T2w predicted LMSIR features classifier (group 3):"
```

```r
print(p3$best_thr$sensitivity)
```

```
                       2.5%       50%     97.5%
0.426981526555143 0.6881443 0.7345361 0.7783505
```

```r
th3 = as.numeric(row.names(p3$best_thr$sensitivity))
confusionMatrix(ifelse(perfall_T2wpLMSIR_total$C >= th3, "C", "NC"), perfall_T2wpLMSIR_total$obs, positive = "C")
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  190 104
        NC  48 284
                                          
               Accuracy : 0.7572          
                 95% CI : (0.7216, 0.7903)
    No Information Rate : 0.6198          
    P-Value [Acc > NIR] : 1.875e-13       
                                          
                  Kappa : 0.5072          
 Mcnemar's Test P-Value : 8.154e-06       
                                          
            Sensitivity : 0.7983          
            Specificity : 0.7320          
         Pos Pred Value : 0.6463          
         Neg Pred Value : 0.8554          
             Prevalence : 0.3802          
         Detection Rate : 0.3035          
   Detection Prevalence : 0.4696          
      Balanced Accuracy : 0.7651          
                                          
       'Positive' Class : C               
                                          
```


# comparing Sensitivity/Specificity at a reference value th1 for only T1w features

```r
T1wcoords <- coords(roc = p1$ROC, x = "all")
rownames(T1wcoords) <- c("threshold", "sensitivity", "specificity")
senT1w = T1wcoords[, T1wcoords["sensitivity", ] >= 0.7010309]
print(senT1w)
```

```
                  all       all       all       all       all       all       all       all
threshold   0.4354237 0.4336886 0.4323613 0.4317271 0.4305428 0.4293794 0.4271350 0.4249132
sensitivity 0.7058824 0.7142857 0.7226891 0.7226891 0.7226891 0.7226891 0.7226891 0.7310924
specificity 0.7474227 0.7474227 0.7474227 0.7422680 0.7371134 0.7319588 0.7268041 0.7268041
                  all       all       all       all       all       all       all       all
threshold   0.4243530 0.4237932 0.4232289 0.4223987 0.4216606 0.4211708 0.4202907 0.4187521
sensitivity 0.7394958 0.7394958 0.7394958 0.7394958 0.7394958 0.7478992 0.7478992 0.7563025
specificity 0.7268041 0.7216495 0.7164948 0.7113402 0.7061856 0.7061856 0.7010309 0.7010309
                  all       all       all       all       all       all       all       all
threshold   0.4168958 0.4156007 0.4135140 0.4119091 0.4117375 0.4112678 0.4108991 0.4106538
sensitivity 0.7647059 0.7731092 0.7731092 0.7731092 0.7731092 0.7815126 0.7815126 0.7899160
specificity 0.7010309 0.7010309 0.6958763 0.6907216 0.6855670 0.6855670 0.6804124 0.6804124
                  all       all       all       all       all       all       all       all
threshold   0.4101642 0.4096650 0.4094110 0.4090691 0.4082209 0.4072931 0.4062231 0.4047941
sensitivity 0.7899160 0.7899160 0.7899160 0.7899160 0.7899160 0.7899160 0.7899160 0.7899160
specificity 0.6752577 0.6701031 0.6649485 0.6597938 0.6546392 0.6494845 0.6443299 0.6391753
                  all       all       all       all       all       all       all       all
threshold   0.4038113 0.4033512 0.4029966 0.4028697 0.4027178 0.4018121 0.4008059 0.4004872
sensitivity 0.7899160 0.7983193 0.7983193 0.7983193 0.8067227 0.8067227 0.8151261 0.8235294
specificity 0.6340206 0.6340206 0.6288660 0.6237113 0.6237113 0.6185567 0.6185567 0.6185567
                  all       all       all       all       all       all       all       all
threshold   0.4002302 0.4000184 0.3996863 0.3993021 0.3991765 0.3978816 0.3959194 0.3948864
sensitivity 0.8319328 0.8319328 0.8403361 0.8403361 0.8403361 0.8403361 0.8403361 0.8403361
specificity 0.6185567 0.6134021 0.6134021 0.6082474 0.6030928 0.5979381 0.5927835 0.5876289
                  all       all       all       all       all       all       all       all
threshold   0.3941368 0.3934826 0.3921405 0.3900182 0.3878066 0.3865541 0.3860894 0.3852591
sensitivity 0.8403361 0.8403361 0.8487395 0.8487395 0.8487395 0.8571429 0.8571429 0.8571429
specificity 0.5824742 0.5773196 0.5773196 0.5721649 0.5670103 0.5670103 0.5618557 0.5567010
                  all       all       all       all       all       all       all       all
threshold   0.3844868 0.3834628 0.3814569 0.3797943 0.3792933 0.3786353 0.3775267 0.3768561
sensitivity 0.8571429 0.8571429 0.8571429 0.8571429 0.8571429 0.8655462 0.8655462 0.8655462
specificity 0.5515464 0.5463918 0.5412371 0.5360825 0.5309278 0.5309278 0.5257732 0.5206186
                  all       all       all       all       all       all       all       all
threshold   0.3767632 0.3762328 0.3751886 0.3743944 0.3740064 0.3732198 0.3723501 0.3715567
sensitivity 0.8739496 0.8739496 0.8739496 0.8823529 0.8823529 0.8823529 0.8823529 0.8823529
specificity 0.5206186 0.5154639 0.5103093 0.5103093 0.5051546 0.5000000 0.4948454 0.4896907
                  all       all       all       all       all       all       all       all
threshold   0.3700382 0.3684096 0.3675405 0.3667413 0.3647826 0.3632885 0.3625054 0.3616759
sensitivity 0.8823529 0.8823529 0.8823529 0.8823529 0.8823529 0.8823529 0.8907563 0.8907563
specificity 0.4845361 0.4793814 0.4742268 0.4690722 0.4639175 0.4587629 0.4587629 0.4536082
                  all       all       all       all       all       all       all       all
threshold   0.3613660 0.3609245 0.3601013 0.3595306 0.3593603 0.3592213 0.3591110 0.3590031
sensitivity 0.8907563 0.8907563 0.8907563 0.8907563 0.8991597 0.8991597 0.8991597 0.8991597
specificity 0.4484536 0.4432990 0.4381443 0.4329897 0.4329897 0.4278351 0.4226804 0.4175258
                  all       all       all       all       all       all       all       all
threshold   0.3585619 0.3575238 0.3565964 0.3557384 0.3550403 0.3547626 0.3539514 0.3532116
sensitivity 0.8991597 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630
specificity 0.4123711 0.4123711 0.4072165 0.4020619 0.3969072 0.3917526 0.3865979 0.3814433
                  all       all       all       all       all       all       all       all
threshold   0.3531240 0.3524369 0.3493043 0.3452080 0.3428371 0.3417941 0.3413219 0.3407544
sensitivity 0.9075630 0.9159664 0.9243697 0.9243697 0.9243697 0.9243697 0.9243697 0.9243697
specificity 0.3762887 0.3762887 0.3762887 0.3711340 0.3659794 0.3608247 0.3556701 0.3505155
                  all       all       all       all       all       all       all       all
threshold   0.3395499 0.3385460 0.3382907 0.3377480 0.3362302 0.3351043 0.3349220 0.3344205
sensitivity 0.9243697 0.9243697 0.9243697 0.9243697 0.9327731 0.9411765 0.9411765 0.9411765
specificity 0.3453608 0.3402062 0.3350515 0.3298969 0.3298969 0.3298969 0.3247423 0.3195876
                  all       all       all       all       all       all       all       all
threshold   0.3335780 0.3307822 0.3271563 0.3244627 0.3223957 0.3210336 0.3190966 0.3177645
sensitivity 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765
specificity 0.3144330 0.3092784 0.3041237 0.2989691 0.2938144 0.2886598 0.2835052 0.2783505
                  all       all       all       all       all       all       all       all
threshold   0.3167001 0.3153688 0.3141513 0.3115508 0.3094783 0.3082008 0.3068723 0.3064394
sensitivity 0.9411765 0.9411765 0.9411765 0.9495798 0.9495798 0.9579832 0.9579832 0.9579832
specificity 0.2731959 0.2680412 0.2628866 0.2628866 0.2577320 0.2577320 0.2525773 0.2474227
                  all       all       all       all       all       all       all       all
threshold   0.3060430 0.3046440 0.3032570 0.3028027 0.3025668 0.3022130 0.3018139 0.3012570
sensitivity 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9663866
specificity 0.2422680 0.2371134 0.2319588 0.2268041 0.2216495 0.2164948 0.2113402 0.2113402
                  all       all       all       all       all       all       all       all
threshold   0.3005337 0.2995512 0.2980326 0.2966159 0.2954221 0.2924609 0.2870837 0.2835092
sensitivity 0.9663866 0.9663866 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899
specificity 0.2061856 0.2010309 0.2010309 0.1958763 0.1907216 0.1855670 0.1804124 0.1752577
                  all       all       all       all       all       all       all       all
threshold   0.2809228 0.2785866 0.2765220 0.2734538 0.2708614 0.2675064 0.2650719 0.2628688
sensitivity 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899
specificity 0.1701031 0.1649485 0.1597938 0.1546392 0.1494845 0.1443299 0.1391753 0.1340206
                  all       all       all       all       all       all       all        all
threshold   0.2609266 0.2604769 0.2586921 0.2556324 0.2536687 0.2514290 0.2486503 0.24616731
sensitivity 0.9747899 0.9747899 0.9747899 0.9747899 0.9831933 0.9831933 0.9831933 0.98319328
specificity 0.1288660 0.1237113 0.1185567 0.1134021 0.1134021 0.1082474 0.1030928 0.09793814
                   all        all        all        all        all        all        all        all
threshold   0.24411307 0.24320323 0.24029535 0.23730432 0.23593792 0.23469474 0.23193763 0.22888450
sensitivity 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664
specificity 0.09793814 0.09278351 0.08762887 0.08247423 0.07731959 0.07216495 0.06701031 0.06185567
                   all        all        all        all        all        all        all       all
threshold   0.22550994 0.22041898 0.21688884 0.21111968 0.20401743 0.20127896 0.19910105 0.1965687
sensitivity 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664 0.99159664 1.00000000 1.0000000
specificity 0.05670103 0.05154639 0.04639175 0.04123711 0.03608247 0.03092784 0.03092784 0.0257732
                   all        all        all         all  all
threshold   0.17857193 0.16008928 0.15225815 0.144576461 -Inf
sensitivity 1.00000000 1.00000000 1.00000000 1.000000000    1
specificity 0.02061856 0.01546392 0.01030928 0.005154639    0
```

```r
bestsenT1w = senT1w[, 1]
print("Results for T1w-only features classifier:")
```

```
[1] "Results for T1w-only features classifier:"
```

```r
print(senT1w[, 1])
```

```
  threshold sensitivity specificity 
  0.4354237   0.7058824   0.7474227 
```

```r
T1T1coords <- coords(roc = p2$ROC, x = "all")
rownames(T1T1coords) <- c("threshold", "sensitivity", "specificity")
senT1T2 = T1T1coords[, T1T1coords["sensitivity", ] >= 0.7010309]
print(senT1T2)
```

```
                  all       all       all       all       all       all       all       all
threshold   0.4680429 0.4671396 0.4664100 0.4659551 0.4650569 0.4638576 0.4625323 0.4617103
sensitivity 0.7058824 0.7058824 0.7142857 0.7226891 0.7310924 0.7310924 0.7394958 0.7478992
specificity 0.7628866 0.7577320 0.7577320 0.7577320 0.7577320 0.7525773 0.7525773 0.7525773
                  all       all       all       all       all       all       all       all
threshold   0.4613745 0.4605833 0.4592056 0.4579810 0.4567407 0.4557261 0.4550164 0.4544421
sensitivity 0.7563025 0.7647059 0.7731092 0.7731092 0.7731092 0.7731092 0.7731092 0.7731092
specificity 0.7525773 0.7525773 0.7525773 0.7474227 0.7422680 0.7371134 0.7319588 0.7268041
                  all       all       all       all       all       all       all       all
threshold   0.4540575 0.4536351 0.4533452 0.4530538 0.4526710 0.4520977 0.4516081 0.4506814
sensitivity 0.7731092 0.7731092 0.7815126 0.7815126 0.7815126 0.7899160 0.7899160 0.7899160
specificity 0.7216495 0.7164948 0.7164948 0.7113402 0.7061856 0.7061856 0.7010309 0.6958763
                  all       all       all       all       all       all       all       all
threshold   0.4486560 0.4465095 0.4451021 0.4443868 0.4440302 0.4437237 0.4434316 0.4430048
sensitivity 0.7899160 0.7983193 0.7983193 0.7983193 0.8067227 0.8067227 0.8067227 0.8067227
specificity 0.6907216 0.6907216 0.6855670 0.6804124 0.6804124 0.6752577 0.6701031 0.6649485
                  all       all       all       all       all       all       all       all
threshold   0.4426972 0.4426178 0.4418901 0.4407250 0.4401092 0.4392686 0.4383980 0.4377659
sensitivity 0.8151261 0.8235294 0.8235294 0.8235294 0.8235294 0.8235294 0.8235294 0.8235294
specificity 0.6649485 0.6649485 0.6597938 0.6546392 0.6494845 0.6443299 0.6391753 0.6340206
                  all       all       all       all       all       all       all       all
threshold   0.4371875 0.4370112 0.4364265 0.4355750 0.4347590 0.4338066 0.4332669 0.4332149
sensitivity 0.8235294 0.8319328 0.8319328 0.8319328 0.8319328 0.8319328 0.8403361 0.8487395
specificity 0.6288660 0.6288660 0.6237113 0.6185567 0.6134021 0.6082474 0.6082474 0.6082474
                  all       all       all       all       all       all       all       all
threshold   0.4325802 0.4317236 0.4314413 0.4311666 0.4308937 0.4308363 0.4301045 0.4290302
sensitivity 0.8487395 0.8487395 0.8487395 0.8487395 0.8487395 0.8487395 0.8487395 0.8487395
specificity 0.6030928 0.5979381 0.5927835 0.5876289 0.5824742 0.5773196 0.5721649 0.5670103
                  all       all       all       all       all       all       all       all
threshold   0.4283785 0.4276473 0.4263920 0.4243262 0.4228951 0.4219582 0.4203754 0.4191938
sensitivity 0.8487395 0.8571429 0.8571429 0.8571429 0.8571429 0.8571429 0.8655462 0.8739496
specificity 0.5618557 0.5618557 0.5567010 0.5515464 0.5463918 0.5412371 0.5412371 0.5412371
                  all       all       all       all       all       all       all       all
threshold   0.4184976 0.4176518 0.4169672 0.4156569 0.4138057 0.4128315 0.4125358 0.4124069
sensitivity 0.8739496 0.8739496 0.8739496 0.8739496 0.8823529 0.8823529 0.8823529 0.8907563
specificity 0.5360825 0.5309278 0.5257732 0.5206186 0.5206186 0.5154639 0.5103093 0.5103093
                  all       all       all       all       all       all       all       all
threshold   0.4105524 0.4070602 0.4052190 0.4050903 0.4044779 0.4029774 0.4019088 0.4016707
sensitivity 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563 0.8991597 0.8991597
specificity 0.5051546 0.5000000 0.4948454 0.4896907 0.4845361 0.4793814 0.4793814 0.4742268
                  all       all       all       all       all       all       all       all
threshold   0.4014484 0.4012543 0.4009321 0.4003171 0.3999422 0.3984215 0.3963126 0.3955732
sensitivity 0.8991597 0.8991597 0.8991597 0.8991597 0.8991597 0.8991597 0.8991597 0.8991597
specificity 0.4690722 0.4639175 0.4587629 0.4536082 0.4484536 0.4432990 0.4381443 0.4329897
                  all       all       all       all       all       all       all       all
threshold   0.3949186 0.3941832 0.3939927 0.3929500 0.3918375 0.3917507 0.3916036 0.3912513
sensitivity 0.8991597 0.8991597 0.8991597 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630
specificity 0.4278351 0.4226804 0.4175258 0.4175258 0.4123711 0.4072165 0.4020619 0.3969072
                  all       all       all       all       all       all       all       all
threshold   0.3910328 0.3906962 0.3898941 0.3887966 0.3874807 0.3865708 0.3860149 0.3856311
sensitivity 0.9075630 0.9075630 0.9159664 0.9159664 0.9243697 0.9243697 0.9243697 0.9327731
specificity 0.3917526 0.3865979 0.3865979 0.3814433 0.3814433 0.3762887 0.3711340 0.3711340
                  all       all       all       all       all       all       all       all
threshold   0.3844666 0.3828267 0.3822036 0.3811799 0.3797476 0.3790705 0.3781507 0.3773098
sensitivity 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731
specificity 0.3659794 0.3608247 0.3556701 0.3505155 0.3453608 0.3402062 0.3350515 0.3298969
                  all       all       all       all       all       all       all       all
threshold   0.3767636 0.3751509 0.3731155 0.3708915 0.3684639 0.3673642 0.3669724 0.3663531
sensitivity 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731 0.9327731
specificity 0.3247423 0.3195876 0.3144330 0.3092784 0.3041237 0.2989691 0.2938144 0.2886598
                  all       all       all       all       all       all       all       all
threshold   0.3658147 0.3654898 0.3651708 0.3648775 0.3645278 0.3632396 0.3618279 0.3607711
sensitivity 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765 0.9411765 0.9495798 0.9495798
specificity 0.2886598 0.2835052 0.2783505 0.2731959 0.2680412 0.2628866 0.2628866 0.2577320
                  all       all       all       all       all       all       all       all
threshold   0.3595353 0.3586961 0.3565204 0.3535417 0.3510565 0.3486411 0.3473542 0.3462929
sensitivity 0.9495798 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832
specificity 0.2525773 0.2525773 0.2474227 0.2422680 0.2371134 0.2319588 0.2268041 0.2216495
                  all       all       all       all       all       all       all       all
threshold   0.3451090 0.3445561 0.3444819 0.3440616 0.3433769 0.3404071 0.3369545 0.3360761
sensitivity 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832
specificity 0.2164948 0.2113402 0.2061856 0.2010309 0.1958763 0.1907216 0.1855670 0.1804124
                  all       all       all       all       all       all       all       all
threshold   0.3349824 0.3337408 0.3326706 0.3310333 0.3266412 0.3223884 0.3212934 0.3174735
sensitivity 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9663866 0.9663866
specificity 0.1752577 0.1701031 0.1649485 0.1597938 0.1546392 0.1494845 0.1494845 0.1443299
                  all       all       all       all       all       all       all       all
threshold   0.3113076 0.3080133 0.3071604 0.3063832 0.3055621 0.3046437 0.3039177 0.3037168
sensitivity 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866
specificity 0.1391753 0.1340206 0.1288660 0.1237113 0.1185567 0.1134021 0.1082474 0.1030928
                   all        all        all        all        all        all        all        all
threshold   0.30237987 0.29953915 0.29682860 0.29379125 0.29119219 0.28958586 0.28800053 0.28666099
sensitivity 0.96638655 0.96638655 0.96638655 0.96638655 0.96638655 0.96638655 0.96638655 0.96638655
specificity 0.09793814 0.09278351 0.08762887 0.08247423 0.07731959 0.07216495 0.06701031 0.06185567
                   all        all        all        all        all        all        all        all
threshold   0.28549963 0.28171975 0.27607623 0.27140909 0.26846785 0.26770932 0.26721418 0.26677108
sensitivity 0.96638655 0.97478992 0.98319328 0.98319328 0.99159664 0.99159664 1.00000000 1.00000000
specificity 0.05670103 0.05670103 0.05670103 0.05154639 0.05154639 0.04639175 0.04639175 0.04123711
                   all        all       all        all        all        all         all  all
threshold   0.26643184 0.26560539 0.2640810 0.25980344 0.24683435 0.23713757 0.231707545 -Inf
sensitivity 1.00000000 1.00000000 1.0000000 1.00000000 1.00000000 1.00000000 1.000000000    1
specificity 0.03608247 0.03092784 0.0257732 0.02061856 0.01546392 0.01030928 0.005154639    0
```

```r
bestsensenT1T2 = senT1T2[, 1]
print("Results for T1w + T2wText + T2w_SI features classifier:")
```

```
[1] "Results for T1w + T2wText + T2w_SI features classifier:"
```

```r
print(senT1T2[, 1])
```

```
  threshold sensitivity specificity 
  0.4680429   0.7058824   0.7628866 
```

```r
T1T2pLMSIRcoords <- coords(roc = p3$ROC, x = "all")
rownames(T1T2pLMSIRcoords) <- c("threshold", "sensitivity", "specificity")
senT1T2pLMSIR = T1T2pLMSIRcoords[, T1T2pLMSIRcoords["sensitivity", ] >= 0.7010309]
print(senT1T2pLMSIR)
```

```
                  all       all       all       all       all       all       all       all
threshold   0.4522711 0.4499933 0.4496323 0.4494524 0.4490592 0.4476646 0.4461738 0.4456368
sensitivity 0.7058824 0.7058824 0.7058824 0.7058824 0.7058824 0.7142857 0.7142857 0.7226891
specificity 0.7989691 0.7938144 0.7886598 0.7835052 0.7783505 0.7783505 0.7731959 0.7731959
                  all       all       all       all       all       all       all       all
threshold   0.4449380 0.4444320 0.4441055 0.4435558 0.4431745 0.4423245 0.4406265 0.4385546
sensitivity 0.7226891 0.7310924 0.7310924 0.7310924 0.7394958 0.7394958 0.7394958 0.7478992
specificity 0.7680412 0.7680412 0.7628866 0.7577320 0.7577320 0.7525773 0.7474227 0.7474227
                  all       all       all       all       all       all       all       all
threshold   0.4373798 0.4369902 0.4366322 0.4353561 0.4330413 0.4309674 0.4288755 0.4277062
sensitivity 0.7563025 0.7647059 0.7731092 0.7731092 0.7815126 0.7815126 0.7815126 0.7899160
specificity 0.7474227 0.7474227 0.7474227 0.7422680 0.7422680 0.7371134 0.7319588 0.7319588
                  all       all       all       all       all       all       all       all
threshold   0.4269815 0.4259642 0.4246981 0.4233049 0.4220641 0.4212520 0.4209776 0.4204775
sensitivity 0.7983193 0.7983193 0.7983193 0.7983193 0.7983193 0.7983193 0.8067227 0.8067227
specificity 0.7319588 0.7268041 0.7216495 0.7164948 0.7113402 0.7061856 0.7061856 0.7010309
                  all       all       all       all       all       all       all       all
threshold   0.4200247 0.4190347 0.4176183 0.4165110 0.4144848 0.4121809 0.4112323 0.4109159
sensitivity 0.8067227 0.8067227 0.8151261 0.8151261 0.8235294 0.8235294 0.8319328 0.8319328
specificity 0.6958763 0.6907216 0.6907216 0.6855670 0.6855670 0.6804124 0.6804124 0.6752577
                  all       all       all       all       all       all       all       all
threshold   0.4106368 0.4104695 0.4103616 0.4099565 0.4091690 0.4077596 0.4066037 0.4057166
sensitivity 0.8319328 0.8319328 0.8319328 0.8319328 0.8403361 0.8403361 0.8403361 0.8487395
specificity 0.6701031 0.6649485 0.6597938 0.6546392 0.6546392 0.6494845 0.6443299 0.6443299
                  all       all       all       all       all       all       all       all
threshold   0.4047971 0.4039606 0.4029872 0.4012004 0.3995722 0.3992930 0.3992325 0.3991964
sensitivity 0.8487395 0.8487395 0.8487395 0.8571429 0.8571429 0.8571429 0.8571429 0.8655462
specificity 0.6391753 0.6340206 0.6288660 0.6288660 0.6237113 0.6185567 0.6134021 0.6134021
                  all       all       all       all       all       all       all       all
threshold   0.3982827 0.3970846 0.3959890 0.3942234 0.3925624 0.3918165 0.3915566 0.3893486
sensitivity 0.8739496 0.8739496 0.8739496 0.8739496 0.8823529 0.8823529 0.8823529 0.8907563
specificity 0.6134021 0.6082474 0.6030928 0.5979381 0.5979381 0.5927835 0.5876289 0.5876289
                  all       all       all       all       all       all       all       all
threshold   0.3868690 0.3862030 0.3860006 0.3853562 0.3844261 0.3838931 0.3832914 0.3824939
sensitivity 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563 0.8907563
specificity 0.5824742 0.5773196 0.5721649 0.5670103 0.5618557 0.5567010 0.5515464 0.5463918
                  all       all       all       all       all       all       all       all
threshold   0.3820204 0.3818201 0.3815723 0.3810836 0.3798032 0.3788754 0.3785469 0.3782345
sensitivity 0.8991597 0.8991597 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630 0.9075630
specificity 0.5463918 0.5412371 0.5412371 0.5360825 0.5309278 0.5257732 0.5206186 0.5154639
                  all       all       all       all       all       all       all       all
threshold   0.3781313 0.3778973 0.3771675 0.3762058 0.3755103 0.3751577 0.3748709 0.3741364
sensitivity 0.9075630 0.9159664 0.9159664 0.9159664 0.9159664 0.9159664 0.9159664 0.9159664
specificity 0.5103093 0.5103093 0.5051546 0.5000000 0.4948454 0.4896907 0.4845361 0.4793814
                  all       all       all       all       all       all       all       all
threshold   0.3735284 0.3732910 0.3731015 0.3729151 0.3721505 0.3712314 0.3704068 0.3698198
sensitivity 0.9243697 0.9243697 0.9243697 0.9243697 0.9327731 0.9327731 0.9327731 0.9327731
specificity 0.4793814 0.4742268 0.4690722 0.4639175 0.4639175 0.4587629 0.4536082 0.4484536
                  all       all       all       all       all       all       all       all
threshold   0.3693749 0.3687424 0.3679345 0.3666613 0.3655365 0.3643474 0.3629612 0.3621820
sensitivity 0.9327731 0.9411765 0.9411765 0.9411765 0.9411765 0.9495798 0.9495798 0.9495798
specificity 0.4432990 0.4432990 0.4381443 0.4329897 0.4278351 0.4278351 0.4226804 0.4175258
                  all       all       all       all       all       all       all       all
threshold   0.3619208 0.3617374 0.3614007 0.3609236 0.3601367 0.3590841 0.3576431 0.3562847
sensitivity 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798
specificity 0.4123711 0.4072165 0.4020619 0.3969072 0.3917526 0.3865979 0.3814433 0.3762887
                  all       all       all       all       all       all       all       all
threshold   0.3547172 0.3533771 0.3527050 0.3519540 0.3510940 0.3491384 0.3474755 0.3468742
sensitivity 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798
specificity 0.3711340 0.3659794 0.3608247 0.3556701 0.3505155 0.3453608 0.3402062 0.3350515
                  all       all       all       all       all       all       all       all
threshold   0.3464090 0.3457552 0.3447976 0.3435113 0.3411039 0.3393053 0.3382369 0.3369591
sensitivity 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798
specificity 0.3298969 0.3247423 0.3195876 0.3144330 0.3092784 0.3041237 0.2989691 0.2938144
                  all       all       all       all       all       all       all       all
threshold   0.3357047 0.3347952 0.3341216 0.3333999 0.3327074 0.3312866 0.3303840 0.3292923
sensitivity 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9495798 0.9579832
specificity 0.2886598 0.2835052 0.2783505 0.2731959 0.2680412 0.2628866 0.2577320 0.2577320
                  all       all       all       all       all       all       all       all
threshold   0.3280827 0.3273232 0.3267269 0.3263706 0.3258835 0.3245648 0.3226320 0.3216622
sensitivity 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832
specificity 0.2525773 0.2474227 0.2422680 0.2371134 0.2319588 0.2268041 0.2216495 0.2164948
                  all       all       all       all       all       all       all       all
threshold   0.3212763 0.3202135 0.3176872 0.3153009 0.3138197 0.3128829 0.3107336 0.3080830
sensitivity 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832 0.9579832
specificity 0.2113402 0.2061856 0.2010309 0.1958763 0.1907216 0.1855670 0.1804124 0.1752577
                  all       all       all       all       all       all       all       all
threshold   0.3069983 0.3063810 0.3038318 0.2978599 0.2931119 0.2907524 0.2886414 0.2868972
sensitivity 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866 0.9663866
specificity 0.1752577 0.1701031 0.1649485 0.1597938 0.1546392 0.1494845 0.1443299 0.1391753
                  all       all       all       all       all       all       all       all
threshold   0.2840562 0.2816127 0.2803735 0.2788575 0.2775695 0.2763861 0.2750835 0.2737451
sensitivity 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899 0.9747899
specificity 0.1391753 0.1340206 0.1288660 0.1237113 0.1185567 0.1134021 0.1082474 0.1030928
                   all        all        all        all        all        all        all        all
threshold   0.27243653 0.27092422 0.26870071 0.26652472 0.26476519 0.26332772 0.25869652 0.25236152
sensitivity 0.97478992 0.97478992 0.98319328 0.98319328 0.98319328 0.98319328 0.98319328 0.98319328
specificity 0.09793814 0.09278351 0.09278351 0.08762887 0.08247423 0.07731959 0.07216495 0.06701031
                   all        all        all        all        all        all        all        all
threshold   0.24716631 0.24421238 0.24417933 0.24124964 0.23318498 0.22737254 0.21974179 0.21141915
sensitivity 0.98319328 0.98319328 0.98319328 0.98319328 0.98319328 0.99159664 0.99159664 0.99159664
specificity 0.06185567 0.05670103 0.05154639 0.04639175 0.04123711 0.04123711 0.03608247 0.03092784
                  all        all        all        all         all       all  all
threshold   0.2046557 0.19814339 0.19266265 0.18784291 0.181462010 0.1738346 -Inf
sensitivity 0.9915966 0.99159664 0.99159664 0.99159664 0.991596639 0.9915966    1
specificity 0.0257732 0.02061856 0.01546392 0.01030928 0.005154639 0.0000000    0
```

```r
bestsensenT1T2pLMSIR = senT1T2pLMSIR[, 1]
print("Results for T1w + T2wText + pLMSIR features classifier:")
```

```
[1] "Results for T1w + T2wText + pLMSIR features classifier:"
```

```r
print(senT1T2pLMSIR[, 1])
```

```
  threshold sensitivity specificity 
  0.4522711   0.7058824   0.7989691 
```


# Sensitivity and Specificity Analysis

source: Cancer Care Ontario, 2011. Ontario Breast Screening Program 2011 Report, 2011.

OBSP program sensitivity (percentage of women diagnosed with breast cancer (DCIS or invasive) within a year of the mammogram date who had an abnormal OBSP screening mammogram result followed by a final diagnosis of breast cancer after completion of diagnostic assessment) has remained relatively high over time and was 86.1% for 2009. Therefore, 13.9% of women with breast cancer diagnosed within a year after the OBSP screen date did not have their breast cancer detected by the program.

OBSP program specificity (percentage of women without a breast cancer diagnosis (DCIS and/or invasive) who had a normal screening mammogram result) has also remained relatively high over time and was 93.1% for 2009. Therefore, 6.9% of women without breast cancer had a false-positive result.

Sensitivity and specificity are affected by a number of factors, including the radiologist's level of experience, the number of previous screens, and the woman's age, breast density and hormone replacement therapy use.

In 2009, sensitivity was 83.0% in women aged 50 to 54, compared with 86.6% to 88.7% in women aged 60 and older. Sensitivity was greater in older women because their breasts are less dense and cancer detection rates are higher for this age group.
Specificity was 90.6% in women aged 50 to 54, compared with 94.6% in women aged 70 to 74 in 2009. The specificity of older women's current screens is improved because these women have more previous screens for comparison.
 
## Confusion matrices of performance compared to OBSP performing sensitivity:

```r
T1wcoords <- coords(roc = p1$ROC, x = "all")
rownames(T1wcoords) <- c("threshold", "sensitivity", "specificity")
senT1w = T1wcoords[, T1wcoords["sensitivity", ] >= 0.861]
bestsenT1w = senT1w[, 1]
print("Results for T1w-only features classifier:")
```

[1] "Results for T1w-only features classifier:"

```r
print(bestsenT1w)
```

  threshold sensitivity specificity 
  0.3786353   0.8655462   0.5309278 

```r
mT1w = cbind(as.numeric(c(round(100 * bestsenT1w[2]), round(100 * (1 - bestsenT1w[2])))), as.numeric(c(round(100 * 
    (1 - bestsenT1w[3])), round(100 * bestsenT1w[3]))))
colnames(mT1w) = c("C", "NC")
rownames(mT1w) = c("predC", "predNC")
pandoc.table(mT1w, keep.line.breaks = TRUE)
```


---------------------
   &nbsp;     C   NC 
------------ --- ----
 **predC**   87   47 

 **predNC**  13   53 
---------------------

```r
T2T1wcoords <- coords(roc = p3$ROC, x = "all")
rownames(T2T1wcoords) <- c("threshold", "sensitivity", "specificity")
senT2T1w = T2T1wcoords[, T2T1wcoords["sensitivity", ] >= 0.861]
bestsenT2T1w = senT2T1w[, 1]
print("Results for T2w predicted LMSIR features classifier:")
```

[1] "Results for T2w predicted LMSIR features classifier:"

```r
print(bestsenT2T1w)
```

  threshold sensitivity specificity 
  0.3991964   0.8655462   0.6134021 

```r
mT2T1w = cbind(as.numeric(c(round(100 * bestsenT2T1w[2]), round(100 * (1 - bestsenT2T1w[2])))), as.numeric(c(round(100 * 
    (1 - bestsenT2T1w[3])), round(100 * bestsenT2T1w[3]))))
colnames(mT2T1w) = c("C", "NC")
rownames(mT2T1w) = c("predC", "predNC")
pandoc.table(mT2T1w, keep.line.breaks = TRUE)
```


---------------------
   &nbsp;     C   NC 
------------ --- ----
 **predC**   87   39 

 **predNC**  13   61 
---------------------


## Confusion matrices of performance compared to OBSP performing specificity:

```r
T1wcoords <- coords(roc = p1$ROC, x = "all")
rownames(T1wcoords) <- c("threshold", "sensitivity", "specificity")
senT1w = T1wcoords[, T1wcoords["specificity", ] >= 0.931]
bestsenT1w = senT1w[, ncol(senT1w)]
print("Results for T1w-only features classifier:")
```

```
[1] "Results for T1w-only features classifier:"
```

```r
print(bestsenT1w)
```

```
  threshold sensitivity specificity 
  0.5235310   0.4201681   0.9329897 
```

```r
mT1w = cbind(as.numeric(c(round(100 * bestsenT1w[2]), round(100 * (1 - bestsenT1w[2])))), as.numeric(c(round(100 * 
    (1 - bestsenT1w[3])), round(100 * bestsenT1w[3]))))
colnames(mT1w) = c("C", "NC")
rownames(mT1w) = c("predC", "predNC")
pandoc.table(mT1w, keep.line.breaks = TRUE)
```

```

---------------------
   &nbsp;     C   NC 
------------ --- ----
 **predC**   42   7  

 **predNC**  58   93 
---------------------
```

```r
T2T1wcoords <- coords(roc = p3$ROC, x = "all")
rownames(T2T1wcoords) <- c("threshold", "sensitivity", "specificity")
senT2T1w = T2T1wcoords[, T2T1wcoords["specificity", ] >= 0.931]
bestsenT2T1w = senT2T1w[, ncol(senT2T1w)]
print("Results for T2w predicted LMSIR features classifier:")
```

```
[1] "Results for T2w predicted LMSIR features classifier:"
```

```r
print(bestsenT2T1w)
```

```
  threshold sensitivity specificity 
  0.5076700   0.5210084   0.9329897 
```

```r
mT2T1w = cbind(as.numeric(c(round(100 * bestsenT2T1w[2]), round(100 * (1 - bestsenT2T1w[2])))), as.numeric(c(round(100 * 
    (1 - bestsenT2T1w[3])), round(100 * bestsenT2T1w[3]))))
colnames(mT2T1w) = c("C", "NC")
rownames(mT2T1w) = c("predC", "predNC")
pandoc.table(mT2T1w, keep.line.breaks = TRUE)
```

```

---------------------
   &nbsp;     C   NC 
------------ --- ----
 **predC**   52   7  

 **predNC**  48   93 
---------------------
```

# Conclusion:
Since high predictive performance is one of the pre-requisites for the applicability of a CAD system to human reading studies, results of this study are encouraging. The Ontario Screening Program (OBSP) defines program sensitivity as the percentage of women diagnosed with breast cancer (DCIS or invasive) within a year of the screening date who had an abnormal OBSP screening result followed by a final diagnosis of breast cancer after completion of diagnostic assessment. Similarly, specificity of OBSP screening is defined as the percentage of women without a breast cancer diagnosis (DCIS and/or invasive) who had a normal screening result. Currently, the reported program sensitivity is 86.1\% and specificity is 93.1\% \cite{CancerCareOntario2011}. Therefore, 13.9\% of women with breast cancer diagnosed within a year after screening did not have their breast cancer detected by the program and 6.9\% of women without breast cancer had a false-positive result. Increasing sensitivity while maintaining high specificity will result in a reduction of breast cancers missed by screening breast MRI, and increasing specificity while maintaining high sensitivity will result in a reduction of false-positive results with screening breast MRI. 

In conclusion, our results showed that CAD predictive performance improved significantly with the introduction of T2w MRI lesion characterization. This improvement is important because CAD can potentially increase both sensitivity and specificity of screening breast MRI by aiding medical profesionals. The assistance of CAD in screening breast MRI is still a matter of further research, but we can argue in favor of the implications of a significant increase in CAD predictive performance. CAD outcomes can be positive (classifying the lesion as malignant) or negative (classifying the lesion as benign). Our CAD sensitivity is estimated as the percentage of diagnosed breast cancer lesions who had an positive CAD diagnosis and CAD specificity as the percentage of diagnosed benign lesions who had an negative CAD diagnosis. If we use the OBSP sensitivity as the minimum operating sensitivity for CAD, the resulting CAD specificity is 53.1\% for CAD lesion characterization based on T1w DCE-MRI only, and specificity increases to 61.3\% with the inclusion of T2w MRI in CAD lesion characterization. This increase in specificity is equivalent to a reduction of 8 false-positive results per 100 women screened with breast MRI who are found to not have cancer. 

If instead we use the OBSP specificity of 93.1\%  as the minimum operating specificity for CAD, the resulting CAD sensitivities are 42.1\% for CAD lesion characterization with T1w features only and increases to 52.1\% with additional T2w features. This increase in sensitivity is equivalent to finding 10 additional cancers per 100 women screened who are found to have cancer. In conclusion, our results showed that CAD predictive performance improved significantly with the introduction of T2w MRI lesion characterization. This improvement is important because CAD can potentially increase both sensitivity and specificity of screening breast MRI by aiding medical professionals.
