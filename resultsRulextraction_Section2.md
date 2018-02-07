# Results Rule extraction - experiments revision 2




## Rules: Extract and summarize rules per k cross

```r
# read k-fold-out
avererrorTop5 = c()
averprobTop5 = c()
classTop5 = c()
T2wfeatureflag = c()
eachtestRules = list()

# allinfo
allTestinfo = c()
allTraininfo = c()
allidx = 0

## for 2nd analysis
alltopRules = c()

for (k in 1:10) {
    # 1:10f cv
    load(paste0("Outputs/gmbRule_noT2SIpredLMSIR_boost_addeddiagvalue_cv", k, ".RData"))
    source("Z:/Cristina/Section2/revision2-finalTreebased-boosting/FunctionsRules.R")
    
    # form feature dictionary
    fdict = feature_dictionary(imgT2pLMSIR)
    fnames = fdict$fnnames
    thresh = as.numeric(row.names(p4$best_thr$sensitivity))
    
    allTestinfo = rbind(allTestinfo, T1T2testinfo)
    allTraininfo = rbind(allTraininfo, imgT2pLMSIRtest)
    print(paste0("Reading cv-fold: ", k, " left out cases"))
    
    # for 2ns analys8is
    eachtestTopRules = list()
    allseltoprules = c()
    
    # present rules collect top scoring rules per k cross
    for (idx in 1:nrow(imgT2pLMSIRtest)) {
        ############## first analysis
        X = imgT2pLMSIRtest[idx, 2:ncol(imgT2pLMSIRtest)]
        y = imgT2pLMSIRtest[idx, "lesion_label"]
        rulesoutput = myapplyLearner(alearner, X, y, minerr = 0.1, minfrq = 10/627, classes, gbmModel)
        
        selRulesIx = rulesoutput[[1]]
        resRulesIx = rulesoutput[[2]]
        
        # present rules
        myrules = mypresentRules(selRulesIx[, c(1:5, 7)], colnames(X), fnames)
        eachtestRules[[allidx + idx]] = myrules
        
        if (length(rulesoutput) == 3) {
            selRuleswT2Ix = rulesoutput[[3]]
            # present rules
            myruleswT2 = mypresentRules(selRuleswT2Ix[, c(1:5, 7)], colnames(X), fnames)
        }
        
        # average error rate of 5 top scoring rules
        avererrorTop5 = c(avererrorTop5, sum(as.numeric(selRulesIx$err))/nrow(selRulesIx))
        
        # recognize whether t2w features are used
        isT2wtop5 = "T2w" %in% unlist(strsplit(myrules$Ruletype, split = " & "))
        T2wfeatureflag = c(T2wfeatureflag, isT2wtop5)
        
        # consider error by rules if less than row.names(p5$best_thr$sensitivity)
        thresh = as.numeric(row.names(p5$best_thr$sensitivity))
        probtop5rules = resRulesIx["probModel"]  #sum(unlist(selRulesIx$sumtempProb))/nrow(selRulesIx)
        averprobTop5 = c(averprobTop5, probtop5rules)
        classTop5 = c(classTop5, ifelse(probtop5rules >= thresh, "C", "NC"))
        
        ############## second analysis
        rulesoutput2 = mynewapplyLearnerxrules(alearner, X, y, minerr = 0.1, minfrq = 10/627, classes, 
            gbmModel)
        topRules = rulesoutput2[[1]]
        
        if (length(rulesoutput2) > 1) {
            eachtestTopRules[[allidx + idx]] = rulesoutput2[[2]]
            allseltoprules = rbind(allseltoprules, rulesoutput2[[3]])
        } else {
            eachtestTopRules[[allidx + idx]] = list()
            allseltoprules = rbind(allseltoprules, 1:nrow(topRules) * 0)
        }
    }
    
    # increment accourdingly
    allidx = allidx + nrow(imgT2pLMSIRtest)
    # append per fold
    alltopRules[[k]] = list(info = T1T2testinfo, data = imgT2pLMSIRtest, topRules = topRules, eachtestTopRules = eachtestTopRules, 
        allseltoprules = allseltoprules)
}
```

```
[1] "Reading cv-fold: 1 left out cases"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
9679   4 0.069 0.158
                                                                                                                             condition
9679 V12<=4.15872171672308 & dce3SE15>0.464838174729024 & lateSE17>0.858933269740971 & T2texture_correlation_nondir>0.0841295312637462
     pred                                       prob sumtempProb
9679   NC 0.6143961, 0.5155817, 0.5006583, 0.5675936   0.5495574
[1] "test complies with # rules: 2 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4959   3 0.086 0.106
                                                                               condition pred
4959 irregularity<=0.909859616438443 & V13>5.24330958981904 & T2RGH_var>345.831042045511   NC
                                prob sumtempProb
4959 0.5741232, 0.6082200, 0.5286033   0.5703155
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 27 out of 112"
     len  freq   err
4815   3 0.056 0.097
                                                                                                condition
4815 texture_variance_nondir_post1>192.620387925557 & V5>3.51255119080561 & T2skew_F_r_i>1.22986842858167
     pred                            prob sumtempProb
4815    C 0.4011533, 0.4281161, 0.6413092   0.4901929
[1] "test complies with # rules: 27 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
2014   3 0.264 0.152
                                                                                                    condition
2014 texture_variance_nondir_post1<=153.171430752191 & lateSE9>1.2281372126798 & T2_lesionSI>53.3142075787255
     pred                            prob sumtempProb
2014   NC 0.4799925, 0.3713933, 0.5005795   0.4506551
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 13 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 3 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 20 out of 112"
[1] "test complies with # rules: 20 out of 112"
[1] "test complies with # rules: 15 out of 112"
     len  freq   err
8608   3 0.062 0.088
                                                                                   condition
8608 dce2SE14>0.432624096558318 & dce3SE13<=0.660588526384995 & T2_lesionSI>43.1284281668948
     pred                            prob sumtempProb
8608   NC 0.4585595, 0.4432739, 0.4748801   0.4589045
[1] "test complies with # rules: 15 out of 112"
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4241   3 0.155 0.118
                                                                                                             condition
4241 earlySE8<=0.956225206128341 & T2grad_margin_var<=1523.4375605976 & T2texture_sumvariance_nondir<=972.655393130357
     pred                            prob sumtempProb
4241   NC 0.3919588, 0.5853889, 0.5176618   0.4983365
[1] "test complies with # rules: 2 out of 112"
[1] "test complies with # rules: 18 out of 112"
[1] "test complies with # rules: 18 out of 112"
[1] "test complies with # rules: 19 out of 112"
[1] "test complies with # rules: 19 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 3 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 2 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 18 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 18 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 14 out of 112"
     len  freq   err
5329   4 0.044 0.083
                                                                                                       condition
5329 V3>4.17885241787875 & earlySE19>0.421479508719241 & lateSE11<=0.810454267166442 & T2RGH_var>345.89795954425
     pred                                       prob sumtempProb
5329   NC 0.5252401, 0.5621318, 0.4823440, 0.4617543   0.5078676
[1] "test complies with # rules: 14 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
1713   3 0.217 0.118
                                                                                           condition
1713 SER_inside<=0.6680364195439 & T2kurt_F_r_i<=4.46086165558844 & LMSIR_predicted>1.97039398616693
     pred                            prob sumtempProb
1713   NC 0.5304901, 0.4858707, 0.5115495   0.5093034
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 3 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
6859   4 0.066 0.139
                                                                                                                                                condition
6859 Kpeak_inside>-0.00173418493839903 & Vr_increasingRate_inside<=0.158266757432716 & V3>3.59029918667748 & T2texture_sumaverage_nondir>41.1799811180647
     pred                                       prob sumtempProb
6859   NC 0.4450344, 0.4982416, 0.4757835, 0.5864314   0.5013727
[1] "test complies with # rules: 13 out of 112"
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 16 out of 112"
[1] "test complies with # rules: 16 out of 112"
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 7 out of 112"
     len  freq   err
4815   3 0.056 0.097
                                                                                                condition
4815 texture_variance_nondir_post1>192.620387925557 & V5>3.51255119080561 & T2skew_F_r_i>1.22986842858167
     pred                            prob sumtempProb
4815    C 0.4011533, 0.4281161, 0.6413092   0.4901929
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 4 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 10 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 3 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
8259   3 0.097 0.132
                                                                                                         condition
8259 circularity<=0.855726514911507 & texture_sumaverage_nondir_post1>125.817629568624 & T2RGH_var>345.89795954425
     pred                            prob sumtempProb
8259   NC 0.6441597, 0.5996009, 0.4813968   0.5750525
[1] "test complies with # rules: 19 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq  err
1042   3 0.078 0.14
                                                                                                                               condition
1042 texture_diffvariance_nondir_post1<=124.098515348777 & texture_energy_nondir_post4<=0.00567725672571473 & T2RGH_var>706.065771712568
     pred                            prob sumtempProb
1042   NC 0.4661382, 0.4903879, 0.3816503   0.4460588
[1] "test complies with # rules: 28 out of 112"
[1] "test complies with # rules: 7 out of 112"
     len  freq   err
5637   3 0.259 0.092
                                                                                        condition
5637 SER_inside<=0.684074475459813 & skew_F_r_i<=0.670444046675986 & T2_lesionSI>47.8757379930724
     pred                            prob sumtempProb
5637   NC 0.3961950, 0.4296275, 0.4587135   0.4281787
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 6 out of 112"
    len  freq   err
884   4 0.168 0.098
                                                                                                                  condition
884 irregularity<=0.920362096693739 & V2>5.72871946106871 & V5>4.48442653681277 & T2texture_entropy_nondir>2.63686154430888
    pred                                       prob sumtempProb
884   NC 0.4926052, 0.4893595, 0.5431961, 0.5185153   0.5109191
[1] "test complies with # rules: 6 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len  freq   err
4241   3 0.155 0.118
                                                                                                             condition
4241 earlySE8<=0.956225206128341 & T2grad_margin_var<=1523.4375605976 & T2texture_sumvariance_nondir<=972.655393130357
     pred                            prob sumtempProb
4241   NC 0.3919588, 0.5853889, 0.5176618   0.4983365
[1] "test complies with # rules: 7 out of 112"
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 11 out of 112"
[1] "test complies with # rules: 2000 out of 112"
[1] "test complies with # rules: 3 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len freq   err
7014   3 0.12 0.167
                                                                                                   condition
7014 mean_F_r_i>760.47219110698 & dce3SE4>1.22621084356217 & T2texture_diffvariance_nondir<=519.719345448646
     pred                            prob sumtempProb
7014    C 0.5166803, 0.6039917, 0.6390194   0.5865638
[1] "test complies with # rules: 12 out of 112"
[1] "test complies with # rules: 2000 out of 112"
     len freq   err
7014   3 0.12 0.167
                                                                                                   condition
7014 mean_F_r_i>760.47219110698 & dce3SE4>1.22621084356217 & T2texture_diffvariance_nondir<=519.719345448646
     pred                            prob sumtempProb
7014    C 0.5166803, 0.6039917, 0.6390194   0.5865638
[1] "test complies with # rules: 5 out of 112"
[1] "test complies with # rules: 8 out of 112"
     len  freq   err
5637   3 0.259 0.092
                                                                                        condition
5637 SER_inside<=0.684074475459813 & skew_F_r_i<=0.670444046675986 & T2_lesionSI>47.8757379930724
     pred                            prob sumtempProb
5637   NC 0.3961950, 0.4296275, 0.4587135   0.4281787
[1] "test complies with # rules: 8 out of 112"
[1] "test complies with # rules: 15 out of 112"
[1] "test complies with # rules: 15 out of 112"
[1] "test complies with # rules: 9 out of 112"
     len freq   err
7086   3 0.04 0.091
                                                                                               condition
7086 max_RGH_var<=0.0830499307748255 & V9<=3.5187690383568 & T2texture_variance_nondir<=325.568685326362
     pred                            prob sumtempProb
7086   NC 0.5486425, 0.6723696, 0.4721203   0.5643775
[1] "test complies with # rules: 9 out of 112"
[1] "test complies with # rules: 11 out of 112"
     len  freq   err
8608   3 0.062 0.088
                                                                                   condition
8608 dce2SE14>0.432624096558318 & dce3SE13<=0.660588526384995 & T2_lesionSI>43.1284281668948
     pred                            prob sumtempProb
8608   NC 0.4585595, 0.4432739, 0.4748801   0.4589045
[1] "test complies with # rules: 11 out of 112"
[1] "Reading cv-fold: 2 left out cases"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 1 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 4 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq   err
2266   3 0.174 0.138
                                                                                           condition
2266 SER_countor>0.626636929139987 & irregularity>0.968408057848163 & T2kurt_F_r_i<=4.66181507751613
     pred                            prob sumtempProb
2266    C 0.7980711, 0.5773883, 0.5266421   0.6340339
[1] "test complies with # rules: 1 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 1 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 45 out of 194"
[1] "test complies with # rules: 45 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 2000 out of 194"
    len  freq   err
834   3 0.131 0.113
                                                                                                      condition
834 iiMin_change_Variance_uptake<=0.232882294446878 & earlySE8<=0.887942203214057 & T2RGH_var<=348.454579378676
    pred                            prob sumtempProb
834   NC 0.4963651, 0.3699555, 0.5985809   0.4883005
[1] "test complies with # rules: 3 out of 194"
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 10 out of 194"
[1] "test complies with # rules: 10 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 27 out of 194"
[1] "test complies with # rules: 27 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 4 out of 194"
[1] "test complies with # rules: 8 out of 194"
     len  freq   err
7286   4 0.159 0.093
                                                                                                                                                                                 condition
7286 irregularity>0.925762714811578 & texture_variance_nondir_post1>188.263791397188 & texture_diffvariance_nondir_post4<=352.460877999769 & T2texture_diffentropy_nondir>1.24666771645525
     pred                                       prob sumtempProb
7286    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq err
6656   5 0.074 0.1
                                                                                                                                                        condition
6656 max_RGH_mean_k<=3.5 & texture_sumvariance_nondir_post1>415.025384577296 & dce2SE2>1.69971466829934 & dce2SE4>1.34623918866759 & T2_lesionSI>52.7200517318938
     pred                                                  prob sumtempProb
6656    C 0.4868222, 0.5391032, 0.6519079, 0.7495536, 0.5621248   0.5979023
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 21 out of 194"
[1] "test complies with # rules: 21 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 2000 out of 194"
    len  freq   err
834   3 0.131 0.113
                                                                                                      condition
834 iiMin_change_Variance_uptake<=0.232882294446878 & earlySE8<=0.887942203214057 & T2RGH_var<=348.454579378676
    pred                            prob sumtempProb
834   NC 0.4963651, 0.3699555, 0.5985809   0.4883005
[1] "test complies with # rules: 4 out of 194"
[1] "test complies with # rules: 21 out of 194"
     len  freq                err
1560   4 0.059 0.0620000000000001
                                                                                                                                  condition
1560 washoutRate_inside>0.00337769286510059 & max_RGH_var>0.0706830683129879 & T2min_F_r_i<=2.5 & T2texture_entropy_nondir>3.44135312373849
     pred                                       prob sumtempProb
1560    C 0.5070210, 0.6912999, 0.7587984, 0.6857119   0.6607078
[1] "test complies with # rules: 21 out of 194"
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 22 out of 194"
[1] "test complies with # rules: 23 out of 194"
[1] "test complies with # rules: 23 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 18 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len freq   err
3455   4 0.03 0.125
                                                                                                                                                                condition
3455 irregularity>0.976454144711226 & texture_sumvariance_nondir_post3<=332.960432019186 & T2kurt_F_r_i<=4.66181507751613 & T2texture_diffentropy_nondir>1.24354568156965
     pred                                       prob sumtempProb
3455    C 0.5264760, 0.7506956, 0.5101984, 0.5449949   0.5830912
[1] "test complies with # rules: 5 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 30 out of 194"
[1] "test complies with # rules: 30 out of 194"
[1] "test complies with # rules: 26 out of 194"
[1] "test complies with # rules: 26 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 10 out of 194"
      len  freq   err
13991   4 0.106 0.088
                                                                                                             condition
13991 mean_F_r_i>781.648597277964 & V19>4.23672233467332 & dce3SE5<=0.913188263359108 & T2skew_F_r_i<=1.38531917646111
      pred                                       prob sumtempProb
13991   NC 0.4503873, 0.5643563, 0.5306730, 0.2458475   0.4478161
[1] "test complies with # rules: 10 out of 194"
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 14 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 20 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 16 out of 194"
[1] "test complies with # rules: 16 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 17 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 4 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 1 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 3 out of 194"
[1] "test complies with # rules: 2000 out of 194"
[1] "test complies with # rules: 2 out of 194"
[1] "test complies with # rules: 32 out of 194"
[1] "test complies with # rules: 32 out of 194"
[1] "test complies with # rules: 9 out of 194"
[1] "test complies with # rules: 9 out of 194"
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 19 out of 194"
      len  freq                err
10069   3 0.063 0.0590000000000001
                                                                                                                               condition
10069 earlySE0<=0.384462022449602 & T2texture_inversediffmoment_nondir<=0.200298613838423 & T2texture_sumaverage_nondir>23.6389015999873
      pred                            prob sumtempProb
10069   NC 0.5055803, 0.4846939, 0.5071630   0.4991458
[1] "test complies with # rules: 19 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 13 out of 194"
[1] "test complies with # rules: 15 out of 194"
     len  freq   err
7286   4 0.159 0.093
                                                                                                                                                                                 condition
7286 irregularity>0.925762714811578 & texture_variance_nondir_post1>188.263791397188 & texture_diffvariance_nondir_post4<=352.460877999769 & T2texture_diffentropy_nondir>1.24666771645525
     pred                                       prob sumtempProb
7286    C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
[1] "test complies with # rules: 15 out of 194"
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 24 out of 194"
[1] "test complies with # rules: 2000 out of 194"
     len  freq   err
8026   3 0.061 0.121
                                                                                              condition
8026 V12<=4.26731764525766 & dce3SE8>0.889522814281942 & T2texture_correlation_nondir>0.124732725643069
     pred                            prob sumtempProb
8026   NC 0.5140166, 0.4364616, 0.4869693   0.4791492
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 6 out of 194"
     len  freq   err
6375   4 0.202 0.092
                                                                                                           condition
6375 irregularity<=0.925751210761359 & V5>3.55748428286012 & V10>5.29508712685388 & LMSIR_predicted>1.98623476838631
     pred                                       prob sumtempProb
6375   NC 0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
[1] "test complies with # rules: 6 out of 194"
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 8 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 11 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "test complies with # rules: 7 out of 194"
[1] "Reading cv-fold: 3 left out cases"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 17 out of 170"
     len  freq                err
2355   5 0.166 0.0649999999999999
                                                                                                                                                                                             condition
2355 maxVr_inside>0.132933970351532 & irregularity<=0.926301266037192 & texture_inversediffmoment_nondir_post4<=0.181096333105491 & T2max_F_r_i>237.5 & T2texture_sumvariance_nondir<=1892.48143275042
     pred                                                  prob sumtempProb
2355   NC 0.5671989, 0.5243535, 0.2820157, 0.4817899, 0.5310127   0.4772741
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 18 out of 170"
[1] "test complies with # rules: 18 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
1023   4 0.146 0.136
                                                                                                                                            condition
1023 washoutRate_countor<=0.00237263537187332 & var_F_r_i<=60724.8179962257 & dce2SE4<=0.86655964324336 & T2texture_variance_nondir<=380.057119029906
     pred                                       prob sumtempProb
1023   NC 0.4613493, 0.5235084, 0.4794116, 0.3752399   0.4598773
[1] "test complies with # rules: 3 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 10 out of 170"
      len  freq   err
10557   3 0.103 0.088
                                                                                     condition
10557 SER_countor<=0.713066433977772 & lateSE3<=1.00154434150194 & T2RGH_var<=348.899261925448
      pred                            prob sumtempProb
10557   NC 0.2830341, 0.4749445, 0.5770729   0.4450172
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 20 out of 170"
[1] "test complies with # rules: 20 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 7 out of 170"
     len  freq   err
7487   5 0.128 0.099
                                                                                                                                                                                  condition
7487 iAUC1_inside<=1867.20963919127 & max_F_r_i<=1997 & texture_variance_nondir_post1>210.870898937361 & texture_diffvariance_nondir_post2<=324.174200101235 & T2_lesionSI>51.1411574636996
     pred                                                  prob sumtempProb
7487    C 0.4257563, 0.4854569, 0.4236016, 0.4099355, 0.4474628   0.4384426
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 19 out of 170"
[1] "test complies with # rules: 19 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 2 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 2 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 3 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 9 out of 170"
[1] "test complies with # rules: 9 out of 170"
[1] "test complies with # rules: 11 out of 170"
[1] "test complies with # rules: 11 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 3 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
9620   4 0.069 0.105
                                                                                                                                               condition
9620 texture_sumaverage_nondir_post4<=168.539244885183 & V11<=51.2306200640228 & earlySE8<=0.974715284142456 & T2texture_entropy_nondir>2.47386447490301
     pred                                       prob sumtempProb
9620   NC 0.5130526, 0.4818617, 0.5684597, 0.5257756   0.5222874
[1] "test complies with # rules: 4 out of 170"
[1] "test complies with # rules: 21 out of 170"
[1] "test complies with # rules: 21 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 8 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 16 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 6 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
1812   4 0.146 0.111
                                                                                                                               condition
1812 texture_variance_nondir_post1>162.468635587987 & lateSE4>0.90412425359233 & lateSE19>1.04549581201775 & T2RGH_var<=605.680930822889
     pred                                       prob sumtempProb
1812    C 0.5779666, 0.4707040, 0.6166417, 0.5418299   0.5517855
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 15 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq   err
5968   5 0.354 0.168
                                                                                                                                                                                             condition
5968 var_F_r_i<=73622.3679819383 & texture_sumvariance_nondir_post1<=262.013409085089 & texture_energy_nondir_post4<=0.00584870387417049 & dce3SE17>0.522817659910163 & T2kurt_F_r_i>-0.55567242234486
     pred                                                  prob sumtempProb
5968   NC 0.4572754, 0.3762719, 0.4043077, 0.4322990, 0.4773361    0.429498
[1] "test complies with # rules: 5 out of 170"
[1] "test complies with # rules: 26 out of 170"
[1] "test complies with # rules: 26 out of 170"
[1] "test complies with # rules: 24 out of 170"
      len  freq   err
13111   5 0.074 0.073
                                                                                                                                                                            condition
13111 texture_sumentropy_nondir_post1>1.94006995925894 & dce2SE7>0.524160803062216 & lateSE5>0.945952590906062 & lateSE8>1.13819699152219 & T2texture_entropy_nondir>2.33506839869461
      pred                                                  prob sumtempProb
13111    C 0.4301397, 0.4668065, 0.4695489, 0.5370841, 0.4446505    0.469646
[1] "test complies with # rules: 24 out of 170"
[1] "test complies with # rules: 17 out of 170"
     len  freq   err
7487   5 0.128 0.099
                                                                                                                                                                                  condition
7487 iAUC1_inside<=1867.20963919127 & max_F_r_i<=1997 & texture_variance_nondir_post1>210.870898937361 & texture_diffvariance_nondir_post2<=324.174200101235 & T2_lesionSI>51.1411574636996
     pred                                                  prob sumtempProb
7487    C 0.4257563, 0.4854569, 0.4236016, 0.4099355, 0.4474628   0.4384426
[1] "test complies with # rules: 17 out of 170"
[1] "test complies with # rules: 12 out of 170"
[1] "test complies with # rules: 12 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 2000 out of 170"
     len  freq err
7061   3 0.036 0.1
                                                                         condition pred
7061 A_countor>0.878898973669725 & max_F_r_i>2132 & T2kurt_F_r_i<=5.44474459334273    C
                                prob sumtempProb
7061 0.5297514, 0.5577010, 0.5004046   0.5292857
[1] "test complies with # rules: 7 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 10 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 13 out of 170"
[1] "test complies with # rules: 2000 out of 170"
[1] "test complies with # rules: 1 out of 170"
[1] "test complies with # rules: 24 out of 170"
[1] "test complies with # rules: 24 out of 170"
[1] "test complies with # rules: 9 out of 170"
[1] "test complies with # rules: 9 out of 170"
[1] "Reading cv-fold: 4 left out cases"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
3745   4 0.031 0.118
                                                                                                                                                                                       condition
3745 Vr_increasingRate_inside>0.119387287517364 & texture_variance_nondir_post1>158.282915050042 & texture_inversediffmoment_nondir_post2>0.144212576138197 & T2grad_margin_var>3296.95677657013
     pred                                       prob sumtempProb
3745    C 0.4794050, 0.3560235, 0.3728296, 0.5111716   0.4298574
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 16 out of 187"
[1] "test complies with # rules: 16 out of 187"
[1] "test complies with # rules: 17 out of 187"
     len  freq   err
3603   5 0.124 0.087
                                                                                                                                condition
3603 irregularity<=0.92495165286538 & V0>7.25458872783576 & V6>5.44452206416059 & lateSE6>0.535319521686134 & T2RGH_mean>21.8641883788852
     pred                                                  prob sumtempProb
3603   NC 0.5119074, 0.4878885, 0.5521666, 0.4338532, 0.4158351   0.4803301
[1] "test complies with # rules: 17 out of 187"
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 2 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 15 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err
12196   3 0.124 0.101
                                                                                                                 condition
12196 SER_inside<=0.656166577515658 & texture_sumaverage_nondir_post4<=207.067310504744 & LMSIR_predicted>1.72917461715247
      pred                            prob sumtempProb
12196   NC 0.4866332, 0.4729611, 0.5208541   0.4934828
[1] "test complies with # rules: 7 out of 187"
[1] "test complies with # rules: 10 out of 187"
     len  freq   err
3603   5 0.124 0.087
                                                                                                                                condition
3603 irregularity<=0.92495165286538 & V0>7.25458872783576 & V6>5.44452206416059 & lateSE6>0.535319521686134 & T2RGH_mean>21.8641883788852
     pred                                                  prob sumtempProb
3603   NC 0.5119074, 0.4878885, 0.5521666, 0.4338532, 0.4158351   0.4803301
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 21 out of 187"
     len  freq                err
3190   4 0.029 0.0620000000000001
                                                                                                                                                          condition
3190 skew_F_r_i>-0.0687768048959314 & irregularity>0.925751210761359 & texture_inversediffmoment_nondir_post4<=0.0923255421666012 & T2kurt_F_r_i>-0.349415185608738
     pred                                       prob sumtempProb
3190   NC 0.5586323, 0.4173910, 0.4476405, 0.4538664   0.4693826
[1] "test complies with # rules: 21 out of 187"
[1] "test complies with # rules: 24 out of 187"
[1] "test complies with # rules: 24 out of 187"
[1] "test complies with # rules: 25 out of 187"
[1] "test complies with # rules: 25 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 33 out of 187"
[1] "test complies with # rules: 33 out of 187"
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 2 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len freq   err
6273   5 0.05 0.107
                                                                                                                                   condition
6273 A_inside<=1.61072218817226 & circularity<=0.855734481103919 & dce2SE11<=0.904584163031592 & lateSE0<=1.15150229077305 & T2min_F_r_i<=23
     pred                                                  prob sumtempProb
6273   NC 0.4461920, 0.6252745, 0.7221661, 0.6488880, 0.5635662   0.6012174
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 4 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 7 out of 187"
     len  freq   err
7128   5 0.137 0.092
                                                                                                                                                                                                      condition
7128 UptakeRate_inside<=1.00418756143313 & texture_energy_nondir_post1<=0.0130423566941265 & texture_variance_nondir_post1<=53.5902830040147 & earlySE10>0.34987914202482 & T2grad_margin_var<=6402.45036359788
     pred                                                  prob sumtempProb
7128   NC 0.4445401, 0.4908341, 0.4556805, 0.5115279, 0.4697661   0.4744697
[1] "test complies with # rules: 7 out of 187"
[1] "test complies with # rules: 10 out of 187"
      len  freq   err
12569   5 0.061 0.088
                                                                                                                                                                                                                         condition
12569 var_F_r_i<=71220.2954156667 & texture_inversediffmoment_nondir_post2>0.200974618005213 & earlySE8<=0.86237226360737 & T2texture_inversediffmoment_nondir>0.0996454758962476 & T2texture_diffentropy_nondir<=1.52849412755374
      pred                                                  prob sumtempProb
12569   NC 0.4320504, 0.3620002, 0.5822682, 0.4883600, 0.4268523   0.4583062
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 26 out of 187"
[1] "test complies with # rules: 26 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 13 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
7574   4 0.079 0.114
                                                                                                                   condition
7574 SER_inside>0.742111684078729 & Vr_post_1_inside>0.171768754142706 & V13>25.2275266554215 & T2_lesionSI>69.3437696866711
     pred                                       prob sumtempProb
7574    C 0.5765578, 0.4723026, 0.4088882, 0.6560729   0.5284554
[1] "test complies with # rules: 16 out of 187"
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 8 out of 187"
[1] "test complies with # rules: 12 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 12 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 21 out of 187"
     len  freq   err
5316   5 0.168 0.075
                                                                                                                                                                   condition
5316 Slope_ini_inside<=0.604112962427516 & mean_F_r_i>533.030523255814 & earlySE10>0.390044237891578 & T2_lesionSIstd<=91.561372458236 & T2grad_margin_var<=6396.57136507791
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
[1] "test complies with # rules: 21 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 18 out of 187"
[1] "test complies with # rules: 7 out of 187"
[1] "test complies with # rules: 7 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq  err
2610   4 0.045 0.12
                                                                                                                                                                               condition
2610 circularity>0.855735583523288 & texture_sumaverage_nondir_post1>128.049051372998 & texture_energy_nondir_post4<=0.00765036336922699 & T2texture_diffentropy_nondir>1.75785784833312
     pred                                       prob sumtempProb
2610    C 0.5030484, 0.4860747, 0.5198417, 0.6825833    0.547887
[1] "test complies with # rules: 2 out of 187"
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 12 out of 187"
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 15 out of 187"
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err
12138   5 0.047 0.115
                                                                                                                                               condition
12138 A_inside>2.41863937468322 & SER_inside>0.717590940320822 & irregularity>0.879980421830873 & lateSE19>1.10405688473339 & T2RGH_var>354.155127809895
      pred                                                  prob sumtempProb
12138    C 0.5931262, 0.5791689, 0.5281811, 0.7094207, 0.4091477   0.5638089
[1] "test complies with # rules: 3 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 4 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 6 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 4 out of 187"
[1] "test complies with # rules: 2000 out of 187"
     len  freq   err
1427   4 0.032 0.111
                                                                                                                                                                              condition
1427 irregularity<=0.984668990886858 & texture_variance_nondir_post1>191.786000082219 & T2texture_correlation_nondir>0.263972379118135 & T2texture_correlation_nondir>0.361597155528048
     pred                                       prob sumtempProb
1427    C 0.4937985, 0.4742500, 0.4468187, 0.5846211   0.4998721
[1] "test complies with # rules: 10 out of 187"
[1] "test complies with # rules: 19 out of 187"
[1] "test complies with # rules: 19 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 14 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 11 out of 187"
[1] "test complies with # rules: 2000 out of 187"
      len  freq   err                                                      condition pred
13176   2 0.032 0.111 max_RGH_mean<=0.495991412376266 & T2_lesionSI>113.683860342556   NC
                      prob sumtempProb
13176 0.5332697, 0.3199769   0.4266233
[1] "test complies with # rules: 5 out of 187"
[1] "test complies with # rules: 9 out of 187"
[1] "test complies with # rules: 9 out of 187"
[1] "test complies with # rules: 2000 out of 187"
[1] "test complies with # rules: 3 out of 187"
[1] "Reading cv-fold: 5 left out cases"
[1] "test complies with # rules: 1000 out of 25"
     len  freq  err
3815   3 0.195 0.14
                                                                                                                           condition
3815 irregularity<=0.929432771288307 & texture_sumaverage_nondir_post4>170.690447525835 & T2texture_variance_nondir>134.392530089492
     pred                            prob sumtempProb
3815   NC 0.5956767, 0.5393015, 0.4608842   0.5319541
[1] "test complies with # rules: 3 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq  err
3815   3 0.195 0.14
                                                                                                                           condition
3815 irregularity<=0.929432771288307 & texture_sumaverage_nondir_post4>170.690447525835 & T2texture_variance_nondir>134.392530089492
     pred                            prob sumtempProb
3815   NC 0.5956767, 0.5393015, 0.4608842   0.5319541
[1] "test complies with # rules: 1 out of 25"
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
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                        condition pred
1560   2 0.093 0.157 irregularity>0.984788961473386 & T2kurt_F_r_i>-0.785487357074651    C
                     prob sumtempProb
1560 0.4785715, 0.5020898   0.4903307
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
2283   3 0.058 0.156
                                                                                            condition
2283 Vr_post_1_inside<=1.9595435001623 & earlySE8>1.84525127427783 & LMSIR_predicted>1.95958273598223
     pred                            prob sumtempProb
2283    C 0.4964885, 0.5170830, 0.5275262   0.5136992
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
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
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len freq   err
1210   3 0.06 0.152
                                                                                                          condition
1210 texture_sumaverage_nondir_post1>178.317115217582 & V19<=4.69302425134638 & T2grad_margin_var<=5974.84492572428
     pred                            prob sumtempProb
1210   NC 0.6718742, 0.5561653, 0.5057248   0.5779214
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 6 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 978 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 3 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
2734   3 0.175 0.115
                                                                                                            condition
2734 iiMin_change_Variance_uptake>0.318089063969771 & T2kurt_F_r_i<=4.55628566003179 & T2kurt_F_r_i<=2.86325409387345
     pred                            prob sumtempProb
2734   NC 0.4691266, 0.4615858, 0.5104885   0.4804003
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1289   3 0.424 0.193
                                                                                                                         condition
1289 washoutRate_inside<=0.00395395574465152 & var_F_r_i<=41313.4935531794 & T2texture_inversediffmoment_nondir<=0.189300676674374
     pred                            prob sumtempProb
1289   NC 0.4808645, 0.3469245, 0.3847784   0.4041892
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 3 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1 out of 25"
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
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
    len  freq   err                                                  condition pred
941   2 0.233 0.141 earlySE10<=0.887901346641671 & T2RGH_var<=348.454579378676   NC
                    prob sumtempProb
941 0.3420054, 0.4877525   0.4148789
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1995   3 0.031 0.176
                                                                                                          condition
1995 Vr_increasingRate_countor>0.0473468061227342 & max_RGH_mean<=0.535957127734903 & T2kurt_F_r_i>4.79798690590228
     pred                            prob sumtempProb
1995   NC 0.5666348, 0.4564168, 0.7178750   0.5803089
[1] "test complies with # rules: 5 out of 25"
[1] "test complies with # rules: 5 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 1 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 4 out of 25"
[1] "test complies with # rules: 7 out of 25"
[1] "test complies with # rules: 7 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err                                                       condition pred
1827   2 0.169 0.183 irregularity>0.977785790197412 & T2kurt_F_r_i<=4.65160033948678    C
                     prob sumtempProb
1827 0.5495094, 0.5000662   0.5247878
[1] "test complies with # rules: 5 out of 25"
[1] "test complies with # rules: 1000 out of 25"
     len  freq   err
1605   2 0.462 0.232
                                                                        condition pred
1605 SER_countor<=0.705775175050863 & T2texture_variance_nondir<=648.653091721668   NC
                     prob sumtempProb
1605 0.4425944, 0.4801895    0.461392
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 1000 out of 25"
[1] "test complies with # rules: 2 out of 25"
[1] "test complies with # rules: 6 out of 25"
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
[1] "test complies with # rules: 1 out of 25"
[1] "Reading cv-fold: 6 left out cases"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                                 condition
4163   3 0.084 0.109 V5>3.51255119080561 & V11<=5.03787177277192 & T2_lesionSI>52.876161736488
     pred                            prob sumtempProb
4163   NC 0.4877236, 0.4626612, 0.5053751   0.4852533
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
4779   3 0.092 0.12
                                                                                                      condition
4779 Tpeak_inside>11.755177337987 & mean_F_r_i<=776.259473090249 & T2texture_sumaverage_nondir>39.9589002202509
     pred                            prob sumtempProb
4779   NC 0.4713799, 0.4513984, 0.5841410   0.5023064
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
1959   3 0.073 0.15
                                                                                                                                                       condition
1959 iiMin_change_Variance_uptake>0.331545153209991 & texture_inversediffmoment_nondir_post2>0.179893974878372 & T2texture_correlation_nondir<=0.370601717909582
     pred                            prob sumtempProb
1959   NC 0.6521837, 0.4705059, 0.5648111   0.5625002
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
5631   3 0.125 0.162
                                                                                                condition
5631 min_F_r_i>386.5 & texture_sumvariance_nondir_post2<=1048.75233234583 & T2kurt_F_r_i>1.14883296665627
     pred                            prob sumtempProb
5631   NC 0.5786830, 0.4887894, 0.4385449   0.5020058
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
    len  freq   err
999   3 0.057 0.129
                                                                                 condition
999 earlySE10>0.398275630549396 & dce3SE11<=0.736907120014966 & T2RGH_var>289.255846492511
    pred                            prob sumtempProb
999   NC 0.5753032, 0.5538037, 0.4956976   0.5416015
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
4215   3 0.205 0.143
                                                                                   condition
4215 irregularity<=0.907737082350515 & V4>2.45154174506599 & T2kurt_F_r_i>-0.521339102394693
     pred                            prob sumtempProb
4215   NC 0.4551003, 0.5220241, 0.5059258   0.4943501
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
3746   3 0.178 0.175
                                                                                    condition
3746 mean_F_r_i<=1117.32137367602 & lateSE5<=0.952831756279501 & T2_lesionSI>52.9125674273991
     pred                            prob sumtempProb
3746   NC 0.4856048, 0.4487932, 0.5044204   0.4796061
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
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1060   3 0.114 0.177
                                                                             condition pred
1060 max_F_r_i>1454.125 & dce2SE9>0.644619457119457 & LMSIR_predicted>2.91344151959955    C
                                prob sumtempProb
1060 0.4692342, 0.6110102, 0.7437256     0.60799
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 6 out of 43"
     len  freq   err
6322   3 0.038 0.095
                                                                                          condition
6322 Kpeak_countor>-0.00631091014572267 & dce2SE6<=0.621422517899093 & T2RGH_mean<=28.5072354996839
     pred                            prob sumtempProb
6322   NC 0.4934857, 0.3293870, 0.5828066   0.4685598
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
6488   3 0.081 0.205
                                                                                                               condition
6488 SER_inside>0.824019225011374 & Vr_post_1_countor>0.225922745750899 & T2texture_correlation_nondir>0.220883358585512
     pred                            prob sumtempProb
6488    C 0.5431932, 0.7044468, 0.5513300   0.5996567
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 8 out of 43"
[1] "test complies with # rules: 1857 out of 43"
    len  freq   err
270   3 0.033 0.167
                                                                                                            condition
270 max_RGH_mean<=0.531406259436042 & texture_variance_nondir_post1<=189.286036058763 & T2kurt_F_r_i>4.73267665300798
    pred                            prob sumtempProb
270   NC 0.4224090, 0.5101070, 0.6310541     0.52119
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1075   3 0.115 0.175
                                                                                               condition
1075 max_F_r_i<=1472.625 & irregularity>0.976213715508256 & T2texture_sumentropy_nondir>1.61393912493361
     pred                            prob sumtempProb
1075    C 0.4853057, 0.4139845, 0.5153649   0.4715517
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len freq   err
5725   2 0.19 0.135
                                                                           condition pred
5725 circularity<=0.855726888231947 & T2texture_correlation_nondir<=0.32564507505023   NC
                     prob sumtempProb
5725 0.4705375, 0.2400161   0.3552768
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len freq   err
5725   2 0.19 0.135
                                                                           condition pred
5725 circularity<=0.855726888231947 & T2texture_correlation_nondir<=0.32564507505023   NC
                     prob sumtempProb
5725 0.4705375, 0.2400161   0.3552768
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
1060   3 0.114 0.177
                                                                             condition pred
1060 max_F_r_i>1454.125 & dce2SE9>0.644619457119457 & LMSIR_predicted>2.91344151959955    C
                                prob sumtempProb
1060 0.4692342, 0.6110102, 0.7437256     0.60799
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                          condition pred
5422   2 0.062 0.147 beta_countor>0.0714608053643474 & T2_lesionSIstd<=54.2682625055716   NC
                     prob sumtempProb
5422 0.5066059, 0.6093130   0.5579594
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 6 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err                                                 condition pred
5545   2 0.081 0.159 dce2SE10<=0.735088866532529 & T2RGH_var<=346.429594497573   NC
                     prob sumtempProb
5545 0.4999955, 0.6814744    0.590735
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 2 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
2730   2 0.121 0.121
                                                                                          condition
2730 texture_variance_nondir_post1<=48.1430868084606 & T2texture_sumaverage_nondir>44.8353718717041
     pred                 prob sumtempProb
2730   NC 0.4030437, 0.4857393   0.4443915
[1] "test complies with # rules: 3 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq  err
2499   3 0.042 0.13
                                                                                                   condition
2499 texture_sumaverage_nondir_post1>131.395713357342 & V10<=3.29645166883552 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
2499   NC 0.4583520, 0.4879422, 0.4728997   0.4730646
[1] "test complies with # rules: 7 out of 43"
[1] "test complies with # rules: 1798 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 1664 out of 43"
[1] "test complies with # rules: 4 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 5 out of 43"
[1] "test complies with # rules: 2000 out of 43"
[1] "test complies with # rules: 1 out of 43"
[1] "test complies with # rules: 2000 out of 43"
     len  freq   err
5631   3 0.125 0.162
                                                                                                condition
5631 min_F_r_i>386.5 & texture_sumvariance_nondir_post2<=1048.75233234583 & T2kurt_F_r_i>1.14883296665627
     pred                            prob sumtempProb
5631   NC 0.5786830, 0.4887894, 0.4385449   0.5020058
[1] "test complies with # rules: 2 out of 43"
[1] "Reading cv-fold: 7 left out cases"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
2358   3 0.104 0.105
                                                                                                         condition
2358 iiMin_change_Variance_uptake>0.283856875469795 & irregularity<=0.970416148881159 & T2RGH_var>353.975029351905
     pred                            prob sumtempProb
2358   NC 0.5610004, 0.2744690, 0.4693122   0.4349272
[1] "test complies with # rules: 3 out of 160"
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 2 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 22 out of 160"
[1] "test complies with # rules: 22 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 1 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 17 out of 160"
[1] "test complies with # rules: 17 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 14 out of 160"
     len  freq   err
8602   6 0.128 0.071
                                                                                                                                                                                               condition
8602 SER_inside>0.838500730200154 & irregularity>0.940117684730477 & V5>3.49086708315126 & dce2SE7>0.625778063035859 & lateSE2>0.806322857503407 & T2texture_inversediffmoment_nondir<=0.201898180625362
     pred                                                             prob sumtempProb
8602    C 0.6746109, 0.4492241, 0.4756217, 0.4966826, 0.4517067, 0.6258111   0.5289428
[1] "test complies with # rules: 14 out of 160"
[1] "test complies with # rules: 1893 out of 160"
[1] "test complies with # rules: 3 out of 160"
[1] "test complies with # rules: 15 out of 160"
[1] "test complies with # rules: 15 out of 160"
[1] "test complies with # rules: 30 out of 160"
[1] "test complies with # rules: 30 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 4 out of 160"
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
4965   5 0.084 0.109
                                                                                                                                                                                                                                   condition
4965 UptakeRate_inside<=0.331212825563122 & texture_variance_nondir_post1<=188.263791397188 & texture_inversediffmoment_nondir_post1>0.208643298985366 & dce2SE17<=0.837805700721211 & T2texture_inversediffmoment_nondir<=0.221325574807973
     pred                                                  prob sumtempProb
4965   NC 0.3571228, 0.4603757, 0.3649369, 0.5094863, 0.4878860   0.4359615
[1] "test complies with # rules: 5 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 4 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 5 out of 160"
[1] "test complies with # rules: 5 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
3267   5 0.186 0.137
                                                                                                                                                                   condition
3267 SER_inside<=0.794461756516995 & texture_sumaverage_nondir_post2<=250.578758618558 & lateSE9>0.770890728378905 & lateSE9>1.24938940369975 & T2_lesionSI>54.1854160477181
     pred                                                  prob sumtempProb
3267   NC 0.3634066, 0.4839432, 0.4559588, 0.4990479, 0.5065191   0.4617751
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 9 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
      len  freq   err
11190   3 0.086 0.128
                                                                                                           condition
11190 texture_sumvariance_nondir_post1<=415.770720390637 & earlySE14>1.17045467703765 & T2_lesionSI>53.2834535991432
      pred                            prob sumtempProb
11190   NC 0.4248698, 0.3256649, 0.4487710   0.3997685
[1] "test complies with # rules: 4 out of 160"
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 11 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
2924   4 0.133 0.137
                                                                                                                                                                    condition
2924 texture_variance_nondir_post2<=180.876506362501 & texture_diffvariance_nondir_post3<=338.678697077251 & T2grad_margin_var>4935.73456149471 & T2RGH_mean>23.9432421868704
     pred                                       prob sumtempProb
2924   NC 0.5609790, 0.5372390, 0.3844625, 0.4686015   0.4878205
[1] "test complies with # rules: 4 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 2000 out of 160"
     len  freq   err
5857   3 0.035 0.105
                                                                                                                        condition
5857 skew_F_r_i>-0.0668608838335453 & texture_inversediffmoment_nondir_post4<=0.0913900931574419 & T2skew_F_r_i>0.134325109470664
     pred                            prob sumtempProb
5857   NC 0.4848748, 0.3959785, 0.4346656   0.4385063
[1] "test complies with # rules: 8 out of 160"
[1] "test complies with # rules: 13 out of 160"
[1] "test complies with # rules: 13 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 2 out of 160"
[1] "test complies with # rules: 7 out of 160"
[1] "test complies with # rules: 7 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 2000 out of 160"
[1] "test complies with # rules: 4 out of 160"
[1] "test complies with # rules: 18 out of 160"
[1] "test complies with # rules: 18 out of 160"
[1] "test complies with # rules: 20 out of 160"
     len  freq                err
3964   4 0.091 0.0600000000000001
                                                                                                                                        condition
3964 texture_energy_nondir_post3<=0.00710250947915509 & V12<=4.03021189922947 & earlySE10>0.412776483160783 & T2grad_margin_var<=13087.2498193737
     pred                                       prob sumtempProb
3964   NC 0.5148103, 0.4888280, 0.5169547, 0.4973049   0.5044745
[1] "test complies with # rules: 20 out of 160"
[1] "test complies with # rules: 23 out of 160"
     len  freq                err
6132   4 0.058 0.0620000000000001
                                                                                                                  condition
6132 mean_F_r_i>811.726530585355 & earlySE4>0.944321360871241 & earlySE8>0.962940140610659 & T2grad_margin>31.6517431192661
     pred                                       prob sumtempProb
6132    C 0.4419462, 0.6174297, 0.4525536, 0.5221466    0.508519
[1] "test complies with # rules: 23 out of 160"
[1] "test complies with # rules: 10 out of 160"
     len  freq   err
5658   3 0.177 0.093
                                                                                                                                 condition
5658 SER_countor<=0.675838603601509 & texture_diffvariance_nondir_post1<=52.6220679340441 & T2texture_correlation_nondir>0.107935297719362
     pred                            prob sumtempProb
5658   NC 0.5505793, 0.3818076, 0.5017071   0.4780313
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 12 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 10 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 6 out of 160"
[1] "test complies with # rules: 16 out of 160"
[1] "test complies with # rules: 16 out of 160"
[1] "Reading cv-fold: 8 left out cases"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 2000 out of 168"
    len  freq err
505   3 0.055 0.1
                                                                                                       condition
505 SER_countor<=0.707326058134249 & var_F_r_i<=19659.5152094893 & T2texture_energy_nondir<=0.000632018983719888
    pred                            prob sumtempProb
505   NC 0.4811075, 0.5600435, 0.3258341   0.4556617
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 8 out of 168"
      len  freq   err
13497   4 0.059 0.094
                                                                                                                                                         condition
13497 texture_sumaverage_nondir_post2<=226.942106218264 & lateSE8<=0.920961171375684 & T2RGH_var<=521.991169084647 & T2texture_energy_nondir<=0.000898111923071435
      pred                                       prob sumtempProb
13497   NC 0.4324431, 0.5009374, 0.2695222, 0.3900091   0.3982279
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 23 out of 168"
[1] "test complies with # rules: 23 out of 168"
[1] "test complies with # rules: 27 out of 168"
[1] "test complies with # rules: 27 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
4377   3 0.053 0.103
                                                                                                          condition
4377 texture_variance_nondir_post3<=223.487783416407 & dce2SE1>0.868363061127592 & T2kurt_F_r_i<=-0.520345577574295
     pred                            prob sumtempProb
4377   NC 0.5659013, 0.6675191, 0.4893800   0.5742668
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 2000 out of 168"
      len  freq   err
12729   5 0.142 0.104
                                                                                                                                                                  condition
12729 maxVr_inside>0.10451017124894 & maxVr_countor>0.0568363361355432 & texture_sumaverage_nondir_post3<=265.972128929092 & earlySE8<=0.967348276913911 & T2max_F_r_i<=406
      pred                                                  prob sumtempProb
12729   NC 0.5068127, 0.4770168, 0.4879522, 0.4684545, 0.3920593   0.4664591
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 16 out of 168"
     len  freq   err
7178   6 0.086 0.085
                                                                                                                                                                                         condition
7178 mean_F_r_i<=1117.30129009231 & irregularity<=0.974279294868243 & texture_sumaverage_nondir_post1>126.361051313361 & V0<=11.82706363796 & V12<=4.67723478611074 & T2_lesionSI>43.6236326923077
     pred                                                             prob sumtempProb
7178   NC 0.5334060, 0.5257874, 0.4944856, 0.4224559, 0.5227059, 0.5457399   0.5074301
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 4 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 16 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 22 out of 168"
[1] "test complies with # rules: 22 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq  err
6133   3 0.138 0.12
                                                                            condition pred
6133 A_inside<=226.058984902916 & V2<=5.23315348512354 & T2RGH_mean<=54.5119132719293   NC
                                prob sumtempProb
6133 0.5390625, 0.5116299, 0.5106702   0.5204542
[1] "test complies with # rules: 3 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 14 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
8099    C 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 7 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
8099    C 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq err
6891   3 0.055 0.1
                                                                                         condition
6891 ivVariance>0.0331643163273542 & circularity<=0.855734293227423 & T2RGH_mean<=54.5391564448446
     pred                            prob sumtempProb
6891   NC 0.4852294, 0.5920242, 0.5623309   0.5465282
[1] "test complies with # rules: 6 out of 168"
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 24 out of 168"
[1] "test complies with # rules: 24 out of 168"
[1] "test complies with # rules: 2000 out of 168"
      len  freq   err
12729   5 0.142 0.104
                                                                                                                                                                  condition
12729 maxVr_inside>0.10451017124894 & maxVr_countor>0.0568363361355432 & texture_sumaverage_nondir_post3<=265.972128929092 & earlySE8<=0.967348276913911 & T2max_F_r_i<=406
      pred                                                  prob sumtempProb
12729   NC 0.5068127, 0.4770168, 0.4879522, 0.4684545, 0.3920593   0.4664591
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 3 out of 168"
[1] "test complies with # rules: 15 out of 168"
[1] "test complies with # rules: 15 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
1833   3 0.088 0.104
                                                                                              condition
1833 irregularity<=0.984466492362561 & T2kurt_F_r_i>4.34329574043788 & LMSIR_predicted>1.98323389316333
     pred                            prob sumtempProb
1833   NC 0.4854764, 0.4002918, 0.4276336   0.4378006
[1] "test complies with # rules: 4 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 3 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 25 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 19 out of 168"
[1] "test complies with # rules: 8 out of 168"
     len  freq   err
1934   3 0.156 0.082
                                                                                   condition
1934 dce2SE12>0.417337428870872 & dce2SE15<=0.905346458673219 & T2RGH_mean<=27.5207599962328
     pred                            prob sumtempProb
1934   NC 0.3281155, 0.4799466, 0.2881014   0.3653878
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 2 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 17 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 4 out of 168"
[1] "test complies with # rules: 2000 out of 168"
     len  freq   err
1273   5 0.105 0.105
                                                                                                                                                                                 condition
1273 max_RGH_var>0.0733863679263791 & texture_sumvariance_nondir_post1>261.264914441484 & V16>7.87099690126964 & lateSE2>0.817789299079648 & T2texture_energy_nondir<=0.000543922788822772
     pred                                                  prob sumtempProb
1273    C 0.6443449, 0.4816719, 0.6141464, 0.6863575, 0.7402640   0.6333569
[1] "test complies with # rules: 8 out of 168"
[1] "test complies with # rules: 12 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
8099    C 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 12 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 18 out of 168"
[1] "test complies with # rules: 10 out of 168"
     len  freq   err
8099   3 0.094 0.078
                                                                                   condition
8099 mean_F_r_i>775.661788827359 & earlySE1>0.956285359163381 & T2_lesionSI>43.6236326923077
     pred                            prob sumtempProb
8099    C 0.4485668, 0.5767493, 0.4689847   0.4981003
[1] "test complies with # rules: 10 out of 168"
[1] "test complies with # rules: 11 out of 168"
     len  freq   err
5842   3 0.066 0.083
                                                                                                                       condition
5842 ivVariance<=0.0375174805963147 & T2texture_entropy_nondir<=3.48878658358655 & T2texture_diffentropy_nondir>1.70774926189513
     pred                            prob sumtempProb
5842   NC 0.4992823, 0.5380998, 0.4868618   0.5080813
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 7 out of 168"
     len  freq   err
8358   5 0.039 0.095
                                                                                                                                                                                 condition
8358 texture_sumaverage_nondir_post3>225.1975459895 & lateSE0>0.871147809637137 & lateSE0>0.987704612189311 & T2_lesionSI>98.472330941509 & T2texture_correlation_nondir>0.263195249129851
     pred                                                  prob sumtempProb
8358    C 0.4192649, 0.5621379, 0.5517737, 0.6701021, 0.4812520   0.5369061
[1] "test complies with # rules: 7 out of 168"
[1] "test complies with # rules: 13 out of 168"
[1] "test complies with # rules: 13 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 9 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 14 out of 168"
[1] "test complies with # rules: 26 out of 168"
[1] "test complies with # rules: 26 out of 168"
[1] "test complies with # rules: 11 out of 168"
     len  freq   err
5842   3 0.066 0.083
                                                                                                                       condition
5842 ivVariance<=0.0375174805963147 & T2texture_entropy_nondir<=3.48878658358655 & T2texture_diffentropy_nondir>1.70774926189513
     pred                            prob sumtempProb
5842   NC 0.4992823, 0.5380998, 0.4868618   0.5080813
[1] "test complies with # rules: 11 out of 168"
[1] "test complies with # rules: 2000 out of 168"
[1] "test complies with # rules: 3 out of 168"
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 5 out of 168"
[1] "test complies with # rules: 6 out of 168"
[1] "test complies with # rules: 6 out of 168"
[1] "Reading cv-fold: 9 left out cases"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 1 out of 47"
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
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4063   4 0.128 0.116
                                                                                                     condition
4063 min_F_r_i>389 & mean_F_r_i<=996.57033595173 & lateSE18<=2.01081577015833 & T2kurt_F_r_i>0.588717905867016
     pred                                       prob sumtempProb
4063   NC 0.5037287, 0.4757571, 0.4159564, 0.3816225   0.4442662
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9774   4 0.165 0.146
                                                                                                                                                               condition
9774 Tpeak_inside>2.36694862113409 & texture_diffvariance_nondir_post1<=120.496206781367 & texture_contrast_nondir_post2>194.701661779236 & T2_lesionSI>44.0377680747215
     pred                                       prob sumtempProb
9774   NC 0.5050244, 0.5244491, 0.4491952, 0.5204948   0.4997909
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9172   3 0.065 0.114
                                                                      condition pred
9172 dce2SE8>0.526488212923419 & lateSE13<=0.782897303276147 & T2min_F_r_i<=8.5   NC
                                prob sumtempProb
9172 0.5711305, 0.5520445, 0.4931413   0.5387721
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
9914   4 0.028 0.133
                                                                                                                         condition
9914 max_RGH_var>0.0719122415437442 & V13>12.797357398796 & V18<=6.75866036727473 & T2texture_diffvariance_nondir>28.0656417109182
     pred                                       prob sumtempProb
9914    C 0.4789980, 0.5192834, 0.6233938, 0.5384807    0.540039
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 2 out of 47"
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
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
7549   3 0.083 0.156
                                                                                                   condition
7549 iMax_Variance_uptake<=6.70552353229436 & lateSE6>0.586002782390577 & T2grad_margin_var>6372.15933273777
     pred                            prob sumtempProb
7549   NC 0.3006064, 0.4980637, 0.4829544   0.4272082
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
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
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 7 out of 47"
[1] "test complies with # rules: 7 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
2249   4 0.102 0.145
                                                                                                                                                                                         condition
2249 Vr_increasingRate_countor>0.0954987823275664 & Vr_decreasingRate_countor<=0.0636349491408924 & texture_sumvariance_nondir_post1>381.737085837675 & T2texture_contrast_nondir>338.784043932143
     pred                                       prob sumtempProb
2249    C 0.5529671, 0.5149765, 0.4829129, 0.6764732   0.5568325
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4635   3 0.052 0.107
                                                                                       condition
4635 edge_sharp_mean>0.671502660560278 & dce3SE1>1.22681788388568 & T2_lesionSI>50.6292782580445
     pred                            prob sumtempProb
4635   NC 0.4823546, 0.2729908, 0.4978068   0.4177174
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
5111   3 0.166 0.133
                                                                                   condition
5111 SER_inside<=0.682123129359999 & lateSE9<=1.19076743320431 & T2RGH_var<=575.314918315023
     pred                            prob sumtempProb
5111   NC 0.5941872, 0.5023692, 0.3768535   0.4911366
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 5 out of 47"
[1] "test complies with # rules: 5 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
5111   3 0.166 0.133
                                                                                   condition
5111 SER_inside<=0.682123129359999 & lateSE9<=1.19076743320431 & T2RGH_var<=575.314918315023
     pred                            prob sumtempProb
5111   NC 0.5941872, 0.5023692, 0.3768535   0.4911366
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4635   3 0.052 0.107
                                                                                       condition
4635 edge_sharp_mean>0.671502660560278 & dce3SE1>1.22681788388568 & T2_lesionSI>50.6292782580445
     pred                            prob sumtempProb
4635   NC 0.4823546, 0.2729908, 0.4978068   0.4177174
[1] "test complies with # rules: 3 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
1153   3 0.207 0.125
                                                                                                          condition
1153 Tpeak_inside>7.90722296680836 & texture_energy_nondir_post4<=0.00618966000712137 & T2RGH_var<=349.491222787327
     pred                            prob sumtempProb
1153   NC 0.4795554, 0.3563718, 0.3132517   0.3830596
[1] "test complies with # rules: 1 out of 47"
[1] "test complies with # rules: 2000 out of 47"
     len  freq   err
4294   4 0.109 0.119
                                                                                                                                 condition
4294 Tpeak_inside>2.54133465388862 & iiMin_change_Variance_uptake>0.395737622135059 & V17<=25.4632748734924 & T2_lesionSI>52.9125674273991
     pred                                       prob sumtempProb
4294   NC 0.4237933, 0.5097312, 0.4460827, 0.4451928      0.4562
[1] "test complies with # rules: 4 out of 47"
[1] "test complies with # rules: 2000 out of 47"
[1] "Reading cv-fold: 10 left out cases"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5012   3 0.176 0.167
                                                                                                                                  condition
5012 texture_diffvariance_nondir_post1<=60.9797692220328 & lateSE9>1.16541692526887 & T2texture_inversediffmoment_nondir<=0.246115364253849
     pred                            prob sumtempProb
5012   NC 0.4782870, 0.3576450, 0.4992456   0.4450592
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5637   4 0.104 0.123
                                                                                                                          condition
5637 A_inside<=1.5528673159683 & V9<=16.7536361430883 & earlySE6<=0.819922077850664 & T2texture_sumvariance_nondir>185.022460304943
     pred                                       prob sumtempProb
5637   NC 0.4832756, 0.5288236, 0.2834297, 0.5989886   0.4736294
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
2001   3 0.165 0.111
                                                                                                           condition
2001 SER_inside<=1.00080283237816 & texture_sumaverage_nondir_post2<=249.582981542509 & T2RGH_mean<=22.5838715481632
     pred                            prob sumtempProb
2001   NC 0.4448360, 0.3843497, 0.3452982   0.3914946
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 8 out of 51"
[1] "test complies with # rules: 8 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len freq   err
813   3 0.11 0.133
                                                                                                                                    condition
813 SER_inside<=0.806001166108312 & texture_inversediffmoment_nondir_post3>0.147244327606865 & T2texture_sumvariance_nondir<=381.387338187838
    pred                            prob sumtempProb
813   NC 0.6043602, 0.4880339, 0.4552047   0.5158663
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 7 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3657   3 0.209 0.123
                                                                                          condition
3657 SER_inside<=0.683454887275891 & skew_F_r_i<=0.724263044011434 & T2grad_margin>30.5173903603995
     pred                            prob sumtempProb
3657   NC 0.5138256, 0.2956922, 0.4628214   0.4241131
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq  err
5742   3 0.046 0.16
                                                                                                                      condition
5742 texture_diffvariance_nondir_post2<=24.0599300089026 & V15<=44.4612422149372 & T2texture_sumaverage_nondir>22.9708102937136
     pred                            prob sumtempProb
5742   NC 0.4570783, 0.5057757, 0.4847688   0.4825409
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5551   4 0.273 0.154
                                                                                                                                                condition
5551 A_inside<=312.396375500685 & edge_sharp_std<=0.991727213591119 & T2RGH_var<=348.526546275196 & T2texture_inversediffmoment_nondir<=0.254036393444923
     pred                                       prob sumtempProb
5551   NC 0.4330791, 0.2267610, 0.4129190, 0.4568792   0.3824096
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
5778   4 0.068 0.135
                                                                                                                                                                           condition
5778 Vr_post_1_inside>0.114510008249947 & texture_sumaverage_nondir_post1<=194.065626177225 & texture_inversediffmoment_nondir_post3<=0.104179970009203 & T2RGH_mean>21.889655177071
     pred                                       prob sumtempProb
5778   NC 0.4829825, 0.5148890, 0.4299427, 0.4467408   0.4686388
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6440   2 0.082 0.133
                                                                          condition pred
6440 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052    C
                     prob sumtempProb
6440 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6543   4 0.207 0.159
                                                                                                                             condition
6543 alpha_inside<=1.57751330223065 & circularity<=0.855732884884821 & V5>4.51415360177622 & T2texture_entropy_nondir>2.90485367714144
     pred                                       prob sumtempProb
6543   NC 0.4572355, 0.4567487, 0.4936095, 0.4366377   0.4610578
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6440   2 0.082 0.133
                                                                          condition pred
6440 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052    C
                     prob sumtempProb
6440 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6440   2 0.082 0.133
                                                                          condition pred
6440 texture_variance_nondir_post1>169.987422892515 & T2kurt_F_r_i>1.18418857120052    C
                     prob sumtempProb
6440 0.6384272, 0.5223976   0.5804124
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3408   3 0.051 0.143
                                                                                                      condition
3408 A_inside<=1.23329492688817 & T2RGH_mean<=52.7500988159803 & T2texture_correlation_nondir>0.255954782848448
     pred                            prob sumtempProb
3408   NC 0.7046357, 0.5281514, 0.5599449   0.5975773
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1024   3 0.313 0.164
                                                                                                                                         condition
1024 texture_variance_nondir_post2<=179.048343064544 & T2texture_energy_nondir>0.000478831110531278 & T2texture_sumaverage_nondir>40.2571158312255
     pred                            prob sumtempProb
1024   NC 0.4499608, 0.2332170, 0.3358754   0.3396844
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 9 out of 51"
[1] "test complies with # rules: 9 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len freq   err
813   3 0.11 0.133
                                                                                                                                    condition
813 SER_inside<=0.806001166108312 & texture_inversediffmoment_nondir_post3>0.147244327606865 & T2texture_sumvariance_nondir<=381.387338187838
    pred                            prob sumtempProb
813   NC 0.6043602, 0.4880339, 0.4552047   0.5158663
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
3675   3 0.059 0.125
                                                                                    condition
3675 max_RGH_mean>0.524313219791607 & lateSE4>0.953138290625636 & T2RGH_mean>47.5488747094981
     pred                            prob sumtempProb
3675   NC 0.5694536, 0.5043591, 0.5279330   0.5339152
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
    len  freq   err
550   3 0.148 0.198
                                                                                                                                     condition
550 texture_variance_nondir_post1>167.08364168864 & texture_inversediffmoment_nondir_post4>0.11508889696368 & LMSIR_predicted>1.76548498403474
    pred                            prob sumtempProb
550    C 0.4777159, 0.6997331, 0.5019843   0.5598111
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 7 out of 51"
[1] "test complies with # rules: 7 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
6277   4 0.134 0.137
                                                                                                                condition
6277 V17<=38.6058204388241 & dce2SE4>0.633227549675858 & dce3SE4<=0.950788958274215 & T2grad_margin_var<=6402.45036359788
     pred                                       prob sumtempProb
6277   NC 0.5104606, 0.5866031, 0.5544229, 0.5395111   0.5477494
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 10 out of 51"
     len  freq   err
1758   3 0.115 0.095
                                                                                                            condition
1758 circularity<=0.855726514911507 & texture_sumaverage_nondir_post2<=249.928310833517 & T2RGH_mean>26.4497935000351
     pred                            prob sumtempProb
1758   NC 0.5552296, 0.4969470, 0.4644038   0.5055268
[1] "test complies with # rules: 10 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 4 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1976   3 0.156 0.129
                                                                                  condition
1976 earlySE11<=0.882260181441451 & lateSE3<=1.10244673946561 & T2RGH_var<=348.899261925448
     pred                            prob sumtempProb
1976   NC 0.5758153, 0.4910063, 0.5130602   0.5266273
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
     len  freq   err
1462   3 0.048 0.115
                                                                                               condition
1462 peakVr_countor<=1.5 & mean_F_r_i>793.188105082232 & T2texture_correlation_nondir<=0.319617514478464
     pred                            prob sumtempProb
1462    C 0.4470454, 0.4989303, 0.6297145   0.5252301
[1] "test complies with # rules: 6 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 3 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 5 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 2 out of 51"
[1] "test complies with # rules: 1500 out of 51"
[1] "test complies with # rules: 1 out of 51"
```

## First Analysis of error and rules extracted:

```r
# 1) summary of all cases top rules
summary(avererrorTop5)
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00540 0.05520 0.07440 0.07578 0.09460 0.16320 
```

```r
summary(avererrorTop5 <= 0.1)
```

```
   Mode   FALSE    TRUE    NA's 
logical     129     498       0 
```

```r
as.numeric(summary(avererrorTop5 <= 0.1)[2:4])/length(avererrorTop5)
```

```
[1] 0.2057416 0.7942584 0.0000000
```

```r
## plot the error
cvstats = cbind(data.frame(avererrorTop5 = avererrorTop5), data.frame(T2wfeatureflag = T2wfeatureflag))
p <- ggplot(cvstats, aes(factor(T2wfeatureflag), avererrorTop5))
p + geom_boxplot(aes(fill = factor(T2wfeatureflag)))
```

![](resultsRulextraction_Section2_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
summary(subset(cvstats, T2wfeatureflag == FALSE))
```

```
 avererrorTop5     T2wfeatureflag 
 Min.   :0.02040   Mode :logical  
 1st Qu.:0.06340   FALSE:177      
 Median :0.07960   NA's :0        
 Mean   :0.08335                  
 3rd Qu.:0.10540                  
 Max.   :0.16320                  
```

```r
summary(subset(cvstats, T2wfeatureflag == TRUE))
```

```
 avererrorTop5     T2wfeatureflag
 Min.   :0.00540   Mode:logical  
 1st Qu.:0.05380   TRUE:450      
 Median :0.07170   NA's:0        
 Mean   :0.07281                 
 3rd Qu.:0.09120                 
 Max.   :0.15940                 
```

```r
# 2) summary of allcases use of T2w features
summary(T2wfeatureflag)
```

```
   Mode   FALSE    TRUE    NA's 
logical     177     450       0 
```

```r
as.numeric(summary(T2wfeatureflag)[2:4])/length(T2wfeatureflag)
```

```
[1] 0.2822967 0.7177033 0.0000000
```

```r
# 4) T2w features were among the top 5 scoring rules in XX lesions (\%: XX benign and XX
# malignant)
T2wamongTop5 = allTestinfo[T2wfeatureflag, -c(20)]
summary(as.factor(T2wamongTop5$lesion_label))
```

```
   massB    massM nonmassB nonmassM 
     188      106      104       52 
```

```r
# Top T1w only features were among the top 5 scoring rules in XX lesions (\%: XX benign and XX
# malignant)
T1wamongTop5 = allTestinfo[!T2wfeatureflag, -c(20)]
summary(as.factor(T1wamongTop5$lesion_label))
```

```
   massB    massM nonmassB nonmassM 
      52       61       38       26 
```

```r
# Top 5 rules accuracy
confusionMatrix(classTop5, allTraininfo$lesion_label)
```

```
Confusion Matrix and Statistics

          Reference
Prediction   C  NC
        C  176  93
        NC  69 289
                                          
               Accuracy : 0.7416          
                 95% CI : (0.7055, 0.7755)
    No Information Rate : 0.6093          
    P-Value [Acc > NIR] : 1.961e-12       
                                          
                  Kappa : 0.4667          
 Mcnemar's Test P-Value : 0.07075         
                                          
            Sensitivity : 0.7184          
            Specificity : 0.7565          
         Pos Pred Value : 0.6543          
         Neg Pred Value : 0.8073          
             Prevalence : 0.3907          
         Detection Rate : 0.2807          
   Detection Prevalence : 0.4290          
      Balanced Accuracy : 0.7375          
                                          
       'Positive' Class : C               
                                          
```

## Second Analysis: Rules extracted by Malignant/benign

```r
# per fold
topRulesIDC = c()
topRulesISDC = c()
topRulesFA = c()
topRulesFI = c()

for (k in 1:10) {
    cat(paste0("\n============== Reading cv-fold: ", k, " left out cases\n"))
    info = alltopRules[[k]]$info
    data = alltopRules[[k]]$data
    topRules = alltopRules[[k]]$topRules
    eachtestTopRules = alltopRules[[k]]$eachtestTopRules
    seltoprules = alltopRules[[k]]$allseltoprules
    # form feature dictionary
    fdict = feature_dictionary(imgT2pLMSIRtest)
    fnames = fdict$fnnames
    
    df = data.frame(cbind(seltoprules, info$lesion_id, info$lesion_label, info$lesion_diagnosis))
    colnames(df) <- c(as.character(1:nrow(topRules)), "lesion_id", "lesion_label", "lesion_diagnosis")
    cols = 1:nrow(topRules)
    df[, cols] = apply(df[, cols], 2, function(x) as.numeric(x))
    
    IDC = subset(df, lesion_diagnosis == "InvasiveDuctal")
    ISDC = subset(df, lesion_diagnosis == "InsituDuctal")
    FA = subset(df, lesion_diagnosis == "FIBROADENOMA")
    FI = subset(df, lesion_diagnosis == "FIBROCYSTIC")
    
    cat(paste0("\n============== Top five explaining rules for InvasiveDuctal (IDC): n = ", nrow(IDC)))
    topIDC = sort(colSums(IDC[, 1:nrow(topRules)]), decreasing = TRUE)
    print(topIDC[1:5])
    rulesTopIDC = as.numeric(names(topIDC[1:5]))
    print(mypresentRules(topRules[rulesTopIDC, ], colnames(data[, 2:ncol(data)]), fnames))
    # append
    topRulesIDC = rbind(topRulesIDC, topRules[rulesTopIDC, ])
    # display a case that meets them all
    casesIDC = IDC[, c(rulesTopIDC, nrow(topRules) + 1, nrow(topRules) + 2, nrow(topRules) + 3)]
    idxIDC = sort(rowSums(casesIDC[, 1:5]), decreasing = TRUE)[1]
    print("Top case meeting rules for InvasiveDuctal (IDC):")
    topcasesIDC = casesIDC[row.names(casesIDC) == names(idxIDC), ]
    print(topcasesIDC)
    
    cat(paste0("\n============== Top five explaining rules for InsituDuctal (ISDC): n = ", nrow(ISDC)))
    topISDC = sort(colSums(ISDC[, 1:nrow(topRules)]), decreasing = TRUE)
    print(topISDC[1:5])
    rulesTopISDC = as.numeric(names(topISDC[1:5]))
    print(mypresentRules(topRules[rulesTopISDC, ], colnames(data[, 2:ncol(data)]), fnames))
    # append
    topRulesISDC = rbind(topRulesISDC, topRules[rulesTopISDC, ])
    # display a case that meets them all
    casesISDC = ISDC[, c(rulesTopISDC, nrow(topRules) + 1, nrow(topRules) + 2, nrow(topRules) + 
        3)]
    idxISDC = sort(rowSums(casesISDC[, 1:5]), decreasing = TRUE)[1]
    print("Top case meeting rules for InsituDuctal (ISDC):")
    topcasesISDC = casesISDC[row.names(casesISDC) == names(idxISDC), ]
    print(topcasesISDC)
    
    cat(paste0("\n============== Top five explaining rules for FIBROADENOMA (FA): n = ", nrow(FA)))
    topFA = sort(colSums(FA[, 1:nrow(topRules)]), decreasing = TRUE)
    print(topFA[1:5])
    rulesTopFA = as.numeric(names(topFA[1:5]))
    print(mypresentRules(topRules[rulesTopFA, ], colnames(data[, 2:ncol(data)]), fnames))
    # append
    topRulesFA = rbind(topRulesFA, topRules[rulesTopFA, ])
    # display a case that meets them all
    casesFA = FA[, c(rulesTopFA, nrow(topRules) + 1, nrow(topRules) + 2, nrow(topRules) + 3)]
    idxFA = sort(rowSums(casesFA[, 1:5]), decreasing = TRUE)[1]
    print("Top case meeting rules for FIBROADENOMA (FA):")
    topcasesFA = casesFA[row.names(casesFA) == names(idxFA), ]
    print(topcasesFA)
    
    cat(paste0("\n============== Top five explaining rules for FIBROCYSTIC (FI): n = ", nrow(FI)))
    topFI = sort(colSums(FI[, 1:nrow(topRules)]), decreasing = TRUE)
    print(topFI[1:5])
    rulesTopFI = as.numeric(names(topFI[1:5]))
    print(mypresentRules(topRules[rulesTopFI, ], colnames(data[, 2:ncol(data)]), fnames))
    # append
    topRulesFI = rbind(topRulesFI, topRules[rulesTopFI, ])
    # display a case that meets them all
    casesFI = FI[, c(rulesTopFI, nrow(topRules) + 1, nrow(topRules) + 2, nrow(topRules) + 3)]
    print(casesFI)
    idxFI = sort(rowSums(casesFI[, 1:5]), decreasing = TRUE)[1]
    print("Top case meeting rules for FIBROCYSTIC (FI):")
    topcasesFI = casesFI[row.names(casesFI) == names(idxFI), ]
    print(topcasesFI)
    
    cat(paste0("\n============== Done cv-fold: ", k, " ==============\n"))
    
}
```

```

============== Reading cv-fold: 1 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 1363 88 93 96 53 
 5  5  5  5  4 
     len  freq   err
2963   3 0.091  0.08
615    2  0.14 0.091
85     2 0.138 0.092
155    2 0.155 0.094
310    2 0.144 0.076
                                                                          condition pred
2963 uptake average high  & Energy ptime4 very low to high  & last post-SE s2 high     C
615                            irregularity high  & 1st post-SE s15 medium to high     C
85                      Variance ptime1 very high  & dispersion s16 medium to high     C
155                     Variance ptime1 very high  & dispersion s16 medium to high     C
310                                    SER(in) medium to high  & irregularity high     C
                                prob sumtempProb                                  Ruletype
2963 0.4627302, 0.4869596, 0.5895322    0.513074 T1wdynamic & T1wtexture & single-time-Enh
615             0.4819968, 0.7471799   0.6145883           T1wmorphology & single-time-Enh
85              0.4344588, 0.7756773   0.6050681                   T1wtexture & dispersion
155             0.4291609, 0.7692691    0.599215                   T1wtexture & dispersion
310             0.5006988, 0.6626664   0.5816826                T1wdynamic & T1wmorphology
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   63 88 93 96 53 lesion_id lesion_label lesion_diagnosis
49  1  1  1  1  1       515     nonmassM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 7105  16  33  38  41 
  3   2   2   2   2 
     len  freq                err
593    2 0.248              0.096
6488   3 0.102              0.036
2897   3 0.031 0.0590000000000001
1666   2 0.055 0.0669999999999999
90     2 0.135 0.0679999999999999
                                                                                                                              condition
593                                                                                  SER(in) low  & max Radial gradient medium to high 
6488                                                   Time-to-peak(in) medium to high  & uptake skewness low  & irregularity very low 
2897 enhancement-variance decreasing rate(in) low to medium  & Rate of signal increase(rim) very low  & dispersion s9 low to very high 
1666                                                                    Curvature at Time-to-peak(rim) low  & 1st post-SE s16 very low 
90                                                                                       irregularity high  & Sum variance ptime1 high 
     pred                            prob sumtempProb                     Ruletype
593    NC            0.5088292, 0.2628051   0.3858172   T1wdynamic & T1wmorphology
6488   NC 0.4742468, 0.3732442, 0.4938563   0.4471158   T1wdynamic & T1wmorphology
2897   NC 0.4995175, 0.6665500, 0.6105294    0.592199      T1wdynamic & dispersion
1666   NC            0.4607474, 0.5917962   0.5262718 T1wdynamic & single-time-Enh
90      C            0.5087730, 0.8044078   0.6565904   T1wmorphology & T1wtexture
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   105 16 33 38 41 lesion_id lesion_label lesion_diagnosis
35   1  0  1  1  1       416     nonmassM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 616 24 44 50 65 
 2  2  2  2  2 
     len  freq   err
6488   3 0.102 0.036
9023   4 0.038 0.048
1996   2 0.184 0.069
411    2 0.228 0.072
336    2 0.246 0.081
                                                                                                                                               condition
6488                                                                    Time-to-peak(in) medium to high  & uptake skewness low  & irregularity very low 
9023 Inverse difference moment ptime1 low to medium  & dispersion s12 very low  & 3rd post-SE s17 low to very high  & T2w Correlation  low to very high 
1996                                                                                                    irregularity low  & Energy ptime4 low to medium 
411                                                                                                           SER(in) low to medium  & irregularity low 
336                                                                                                           SER(in) low to medium  & irregularity low 
     pred                                       prob sumtempProb
6488   NC            0.4742468, 0.3732442, 0.4938563   0.4471158
9023   NC 0.5515969, 0.6704627, 0.4683166, 0.4878620   0.5445595
1996   NC                       0.5007048, 0.3009807   0.4008427
411    NC                       0.4405761, 0.3043518    0.372464
336    NC                       0.5053493, 0.3706932   0.4380213
                                            Ruletype
6488                      T1wdynamic & T1wmorphology
9023 T1wtexture & dispersion & single-time-Enh & T2w
1996                      T1wmorphology & T1wtexture
411                       T1wdynamic & T1wmorphology
336                       T1wdynamic & T1wmorphology
[1] "Top case meeting rules for FIBROADENOMA (FA):"
  16 24 44 50 65 lesion_id lesion_label lesion_diagnosis
7  1  1  1  1  1       104        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 534 41 42 44 62 
 3  2  2  2  2 
     len  freq                err
1073   3  0.06 0.0610000000000001
90     2 0.135 0.0679999999999999
6000   3 0.107 0.0679999999999999
1996   2 0.184              0.069
315    2  0.14              0.078
                                                                                                            condition
1073 circularity high & T2w radial gradient variance medium to high  & T2w radial gradient variance very low to high 
90                                                                     irregularity high  & Sum variance ptime1 high 
6000         irregularity medium to high  & Max Radial gradient variance medium to high  & Variance ptime1 very high 
1996                                                                 irregularity low  & Energy ptime4 low to medium 
315                                                                    irregularity high  & Sum variance ptime1 high 
     pred                            prob sumtempProb                   Ruletype
1073   NC 0.6781025, 0.4757408, 0.5921319   0.5819917        T1wmorphology & T2w
90      C            0.5087730, 0.8044078   0.6565904 T1wmorphology & T1wtexture
6000    C 0.4969203, 0.6008289, 0.5469925   0.5482472 T1wmorphology & T1wtexture
1996   NC            0.5007048, 0.3009807   0.4008427 T1wmorphology & T1wtexture
315     C            0.4957177, 0.7767643    0.636241 T1wmorphology & T1wtexture
   34 41 42 44 62 lesion_id lesion_label lesion_diagnosis
13  0  1  1  0  1       165        massB      FIBROCYSTIC
19  1  0  0  1  0       228        massB      FIBROCYSTIC
25  1  0  0  1  0       301     nonmassB      FIBROCYSTIC
34  0  1  1  0  1       415     nonmassB      FIBROCYSTIC
47  1  0  0  0  0       503        massB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   34 41 42 44 62 lesion_id lesion_label lesion_diagnosis
13  0  1  1  0  1       165        massB      FIBROCYSTIC

============== Done cv-fold: 1 ==============

============== Reading cv-fold: 2 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 17 24 104 137 162  49 
  9   5   5   5   4 
      len  freq                err
322     2 0.113              0.049
300     3 0.115              0.081
13991   4 0.106              0.088
7286    4 0.159              0.093
11396   3 0.063 0.0590000000000001
                                                                                                                                                 condition
322                                                                                                             irregularity high  & Variance ptime1 high 
300                                                                          SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
13991                                  uptake average high  & dispersion s19 low to very high  & 3rd post-SE s5 very low  & skewness T2w SI low to medium 
7286  irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
11396                                                              Washout rate(in) very high  & Sum Entropy ptime1 very high  & Entropy ptime1 very high 
      pred                                       prob sumtempProb
322      C                       0.4675507, 0.7419192    0.604735
300      C            0.5949279, 0.7518615, 0.4569649   0.6012514
13991   NC 0.4503873, 0.5643563, 0.5306730, 0.2458475   0.4478161
7286     C 0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
11396    C            0.5191539, 0.4959777, 0.3986381   0.4712566
                                             Ruletype
322                        T1wmorphology & T1wtexture
300                        T1wdynamic & T1wmorphology
13991 T1wdynamic & dispersion & single-time-Enh & T2w
7286                 T1wmorphology & T1wtexture & T2w
11396                         T1wdynamic & T1wtexture
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   24 104 137 162 49 lesion_id lesion_label lesion_diagnosis
19  1   1   0   1  1       235        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 13162 139  30 104 106 
  8   6   5   5   5 
     len  freq   err
7286   4 0.159 0.093
4150   5 0.104 0.089
8280   4 0.035 0.053
300    3 0.115 0.081
4592   2 0.091 0.082
                                                                                                                                                                      condition
7286                       irregularity medium to high  & Variance ptime1 very high  & Difference variance ptime4 very low to high  & T2w Difference entropy  low to very high 
4150 uptake average high  & Max Radial gradient variance low to very high  & Energy ptime4 very low to high  & dispersion s19 medium to high  & 1st post-SE s12 medium to high 
8280                                                          Rate of signal decrease(in) high & uptake variance very high  & circularity low & last post-SE s8 medium to high 
300                                                                                               SER(in) medium to high  & uptake skewness low  & irregularity medium to high 
4592                                                                                                                                uptake average high  & 2nd post-SE s2 high 
     pred                                                  prob sumtempProb
7286    C            0.7622081, 0.4384610, 0.4232979, 0.4466687   0.5176589
4150    C 0.5303059, 0.5125338, 0.5635975, 0.6702745, 0.6316332    0.581669
8280    C            0.4768186, 0.4505621, 0.7410536, 0.6456018    0.578509
300     C                       0.5949279, 0.7518615, 0.4569649   0.6012514
4592    C                                  0.4996731, 0.5976669     0.54867
                                                                   Ruletype
7286                                       T1wmorphology & T1wtexture & T2w
4150 T1wdynamic & T1wmorphology & T1wtexture & dispersion & single-time-Enh
8280                           T1wdynamic & T1wmorphology & single-time-Enh
300                                              T1wdynamic & T1wmorphology
4592                                           T1wdynamic & single-time-Enh
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   162 139 30 104 106 lesion_id lesion_label lesion_diagnosis
64   1   1  1   1   1       574        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 1390 29 80  7 77 
 6  5  5  4  4 
     len  freq   err
1338   2 0.193 0.077
519    2 0.176 0.053
3987   5 0.078 0.071
6093   6 0.031     0
331    3 0.106  0.07
                                                                                                                                                                                                                                              condition
1338                                                                                                                                                                                                   irregularity low  & Energy ptime4 low to medium 
519                                                                                                                                                                                                    irregularity low  & Energy ptime4 low to medium 
3987                                                             Inverse difference moment ptime1 very low to high  & dispersion s5 low to very high  & dispersion s7 low to very high  & 1st post-SE s19 low to very high  & 3rd post-SE s12 very low 
6093 SER(in) low to medium  & Variance of spatial Margin Gradient low to very high  & std 3D Sharpness of lesion margin  low to medium  & Inverse difference moment ptime4 medium to high  & 2nd post-SE s12 very low  & 2nd post-SE s17 low to medium 
331                                                                                                                                                                          SER(in) medium to high  & irregularity low  & Energy ptime4 low to medium 
     pred                                                             prob sumtempProb
1338   NC                                             0.4991863, 0.2996197    0.399403
519    NC                                             0.4703946, 0.2461935    0.358294
3987   NC            0.5232709, 0.5661026, 0.5040197, 0.5410826, 0.4760760   0.5221104
6093   NC 0.6855391, 0.5934397, 0.5476637, 0.4856720, 0.6810844, 0.5714698   0.5941448
331    NC                                  0.4648444, 0.5897234, 0.3071130   0.4538936
                                                      Ruletype
1338                                T1wmorphology & T1wtexture
519                                 T1wmorphology & T1wtexture
3987                 T1wtexture & dispersion & single-time-Enh
6093 T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
331                    T1wdynamic & T1wmorphology & T1wtexture
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   90 29 80 7 77 lesion_id lesion_label lesion_diagnosis
48  1  1  1 1  1       436        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 10130 160  90 107 142 
  6   6   5   5   5 
     len  freq   err
5547   5  0.17 0.087
6375   4 0.202 0.092
1338   2 0.193 0.077
7750   2 0.091 0.082
8424   3 0.185  0.09
                                                                                                                                                      condition
5547 circularity high & Sum average ptime2 very low to high  & dispersion s0 very low to high  & dispersion s5 low to very high  & T2w std SI low to very high 
6375                                 irregularity low  & dispersion s5 low to very high  & dispersion s10 low to very high  & predicted LMSIR low to very high 
1338                                                                                                           irregularity low  & Energy ptime4 low to medium 
7750                                                                                                                 circularity high & Sum average ptime1 low 
8424                                                       Curvature at Time-to-peak(in) low to very high  & uptake skewness low to medium  & irregularity low 
     pred                                                  prob sumtempProb
5547   NC 0.3789780, 0.4486299, 0.4679563, 0.4074943, 0.4236785   0.4253474
6375   NC            0.4029223, 0.4694528, 0.4305508, 0.4448463   0.4369431
1338   NC                                  0.4991863, 0.2996197    0.399403
7750   NC                                  0.5779883, 0.4539037    0.515946
8424   NC                       0.2901904, 0.4699588, 0.5212986   0.4271493
                                          Ruletype
5547 T1wmorphology & T1wtexture & dispersion & T2w
6375              T1wmorphology & dispersion & T2w
1338                    T1wmorphology & T1wtexture
7750                    T1wmorphology & T1wtexture
8424                    T1wdynamic & T1wmorphology
   130 160 90 107 142 lesion_id lesion_label lesion_diagnosis
5    0   0  0   0   0        57     nonmassB      FIBROCYSTIC
6    0   0  0   0   0        58     nonmassB      FIBROCYSTIC
7    1   1  0   1   1        66     nonmassB      FIBROCYSTIC
8    1   1  1   0   1        74        massB      FIBROCYSTIC
14   0   0  1   0   0       123        massB      FIBROCYSTIC
23   0   0  0   0   0       287        massB      FIBROCYSTIC
24   1   1  0   1   0       288        massB      FIBROCYSTIC
41   1   1  1   1   1       383     nonmassB      FIBROCYSTIC
45   1   1  1   1   1       411     nonmassB      FIBROCYSTIC
46   1   1  1   1   1       412     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   130 160 90 107 142 lesion_id lesion_label lesion_diagnosis
41   1   1  1   1   1       383     nonmassB      FIBROCYSTIC

============== Done cv-fold: 2 ==============

============== Reading cv-fold: 3 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 10 40 127 139   6  31 
  3   3   3   2   2 
      len  freq   err
420     3 0.148 0.049
182     2 0.123 0.088
10681   2  0.08 0.091
3499    4 0.031     0
4984    2 0.042 0.043
                                                                                                                                  condition
420                                          Variance ptime1 very high  & dispersion s16 low to very high  & 3rd post-SE s6 medium to high 
182                                                                              Variance ptime1 very high  & dispersion s9 medium to high 
10681                                                  Difference variance ptime1 very high  & Difference variance ptime2 very low to high 
3499  max uptake very high  & Sum average ptime2 very low to high  & last post-SE s8 low to very high  & T2w Sum average  low to very high 
4984                                                                              max uptake very high  & last post-SE s8 low to very high 
      pred                                       prob sumtempProb
420      C            0.4705374, 0.8140437, 0.7544592   0.6796801
182      C                       0.4279656, 0.7763369   0.6021512
10681    C                       0.4832205, 0.4638078   0.4735142
3499     C 0.4648681, 0.4791564, 0.4145251, 0.4431299   0.4504199
4984     C                       0.5415453, 0.4910949   0.5163201
                                             Ruletype
420         T1wtexture & dispersion & single-time-Enh
182                           T1wtexture & dispersion
10681                                      T1wtexture
3499  T1wdynamic & T1wtexture & single-time-Enh & T2w
4984                     T1wdynamic & single-time-Enh
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   40 127 139 6 31 lesion_id lesion_label lesion_diagnosis
37  0   1   1 1  1       404        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 6 29  30  76  81 103 
  2   2   2   2   2 
      len  freq                err
3314    5 0.042              0.043
4258    2 0.085              0.043
12778   3 0.107 0.0679999999999999
11518   4 0.051              0.071
9510    5  0.09               0.08
                                                                                                                                                                             condition
3314  Area Under Uptake curve(rim) high  & Difference variance ptime2 very low to high  & dispersion s17 low to medium  & 1st post-SE s13 low to very high  & 3rd post-SE s5 very low 
4258                                                                                                                      2nd post-SE s14 low to very high  & 3rd post-SE s9 very low 
12778                                                                                   irregularity very low  & 1st post-SE s14 very low to high  & kurtosis T2w SI low to very high 
11518                            Inverse difference moment ptime4 low to very high  & 2nd post-SE s7 low to very high  & 3rd post-SE s9 very low  & T2w Correlation  very low to high 
9510               Energy ptime4 very low to high  & Contrast ptime4 medium to high  & dispersion s0 very low to high  & 3rd post-SE s19 low to very high  & last post-SE s5 very low 
      pred                                                  prob sumtempProb
3314    NC 0.5324878, 0.5164865, 0.6178433, 0.6869821, 0.6636837   0.6034967
4258    NC                                  0.4849738, 0.4662653   0.4756196
12778   NC                       0.4735779, 0.4745754, 0.5102290   0.4861274
11518   NC            0.4891468, 0.5655790, 0.5362972, 0.5393142   0.5325843
9510    NC 0.5528926, 0.4808834, 0.5060270, 0.5276022, 0.5374895   0.5209789
                                                    Ruletype
3314  T1wdynamic & T1wtexture & dispersion & single-time-Enh
4258                                         single-time-Enh
12778                  T1wmorphology & single-time-Enh & T2w
11518                     T1wtexture & single-time-Enh & T2w
9510               T1wtexture & dispersion & single-time-Enh
[1] "Top case meeting rules for InsituDuctal (ISDC):"
  29 30 76 81 103 lesion_id lesion_label lesion_diagnosis
1  1  1  0  1   1        33        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 8103  86  88 142  21 
  5   3   3   3   2 
      len  freq   err
9510    5  0.09  0.08
9538    3 0.074 0.073
515     4 0.195 0.074
10600   2 0.098 0.093
11672   4 0.067 0.027
                                                                                                                                                                condition
9510  Energy ptime4 very low to high  & Contrast ptime4 medium to high  & dispersion s0 very low to high  & 3rd post-SE s19 low to very high  & last post-SE s5 very low 
9538                                                          Max Radial gradient variance very low  & Energy ptime3 medium to high  & Sum variance ptime3 low to medium 
515                                        SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
10600                                                                                                        2nd post-SE s14 low to very high  & 3rd post-SE s9 very low 
11672                       Initial Uptake slope(in) very low to high  & 1st post-SE s10 low to very high  & 3rd post-SE s9 very low  & last post-SE s4 low to very high 
      pred                                                  prob sumtempProb
9510    NC 0.5528926, 0.4808834, 0.5060270, 0.5276022, 0.5374895   0.5209789
9538    NC                       0.5114794, 0.4944816, 0.5861188   0.5306933
515     NC            0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
10600   NC                                  0.4828865, 0.4670870   0.4749867
11672   NC            0.5048007, 0.5222387, 0.5252963, 0.4859425   0.5095695
                                                       Ruletype
9510                  T1wtexture & dispersion & single-time-Enh
9538                                 T1wmorphology & T1wtexture
515   T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
10600                                           single-time-Enh
11672                              T1wdynamic & single-time-Enh
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   103 86 88 142 21 lesion_id lesion_label lesion_diagnosis
12   1  1  1   0  1       127        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 10 85  88 146 107 131 
  4   4   4   3   3 
     len  freq   err
883    2 0.197 0.073
515    4 0.195 0.074
750    2  0.25 0.094
723    3 0.108 0.083
8534   4 0.103 0.088
                                                                                                                          condition
883                                                                                irregularity low  & Energy ptime4 low to medium 
515  SER(in) low  & max Radial gradient medium to high  & Sum variance ptime2 very low to high  & 1st post-SE s10 low to very high 
750                                                             irregularity low  & Inverse difference moment ptime4 low to medium 
723                                          irregularity low to medium  & 3rd post-SE s10 low to very high  & 3rd post-SE s16 low 
8534                              SER(rim) very low to high  & dispersion s1 very low  & dispersion s14 low  & 2nd post-SE s16 high
     pred                                       prob sumtempProb
883    NC                       0.4694149, 0.2405692   0.3549921
515    NC 0.2902278, 0.4509106, 0.3336921, 0.3700841   0.3612286
750    NC                       0.4236079, 0.2408969   0.3322524
723    NC            0.5052503, 0.4047316, 0.3668232   0.4256017
8534   NC 0.4581848, 0.5217262, 0.6139028, 0.5050821    0.524724
                                                      Ruletype
883                                 T1wmorphology & T1wtexture
515  T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
750                                 T1wmorphology & T1wtexture
723                            T1wmorphology & single-time-Enh
8534                 T1wdynamic & dispersion & single-time-Enh
   85 88 146 107 131 lesion_id lesion_label lesion_diagnosis
5   1  0   1   1   1        93        massB      FIBROCYSTIC
14  0  0   0   0   0       146     nonmassB      FIBROCYSTIC
19  1  0   1   0   1       212        massB      FIBROCYSTIC
20  1  1   1   1   0       213        massB      FIBROCYSTIC
23  0  1   0   0   0       233        massB      FIBROCYSTIC
27  1  1   1   1   0       311        massB      FIBROCYSTIC
34  0  0   0   0   0       339     nonmassB      FIBROCYSTIC
40  0  0   0   0   0       419     nonmassB      FIBROCYSTIC
41  0  0   0   0   0       461     nonmassB      FIBROCYSTIC
42  0  1   0   0   1       462        massB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
  85 88 146 107 131 lesion_id lesion_label lesion_diagnosis
5  1  0   1   1   1        93        massB      FIBROCYSTIC

============== Done cv-fold: 3 ==============

============== Reading cv-fold: 4 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 14147 121 148 183  92 
  8   7   7   7   6 
      len  freq   err
62      3  0.16  0.09
455     2 0.128 0.085
588     2 0.121  0.09
511     2  0.13 0.097
12068   3 0.094 0.077
                                                                                                                             condition
62    max relative signal enhancement(rim) medium to high  & Variance ptime1 very high  & Difference variance ptime2 very low to high 
455                                                                                       SER(in) medium to high  & irregularity high 
588                                                                                                 SER(in) high  & irregularity high 
511                                                                                       SER(in) medium to high  & irregularity high 
12068                                 irregularity low to very high  & Sum variance ptime1 very high  & dispersion s18 medium to high 
      pred                            prob sumtempProb                                Ruletype
62       C 0.8237662, 0.4610768, 0.7663881   0.6837437                 T1wdynamic & T1wtexture
455      C            0.7662063, 0.4502787   0.6082425              T1wdynamic & T1wmorphology
588      C            0.4565346, 0.6397499   0.5481423              T1wdynamic & T1wmorphology
511      C            0.7969179, 0.5002271   0.6485725              T1wdynamic & T1wmorphology
12068    C 0.4815758, 0.4739428, 0.5346694   0.4967293 T1wmorphology & T1wtexture & dispersion
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   147 121 148 183 92 lesion_id lesion_label lesion_diagnosis
14   1   1   1   1  1       168        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 6159  16  54  57  71 
  3   2   2   2   2 
      len  freq                err
12921   3 0.079              0.091
11786   5 0.049                  0
4983    3 0.059 0.0610000000000001
2013    4 0.029 0.0620000000000001
7718    6 0.135 0.0669999999999999
                                                                                                                                                                                               condition
12921                                                                                                           uptake kurtosis high  & 1st post-SE s16 low to medium  & 3rd post-SE s18 medium to high 
11786                            Curvature at Time-to-peak(in) high  & Time-to-peak(rim) medium to high  & Variance ptime1 very low  & dispersion s5 low to very high  & dispersion s9 low to very high 
4983                                                                                                                     irregularity high  & last post-SE s12 very high  & max T2w SI very low to high 
2013                                                                                                     max uptake low  & irregularity high  & Variance ptime2 high  & dispersion s15 low to very high 
7718  Sum variance ptime1 low  & Energy ptime4 very low to high  & dispersion s0 medium to high  & 2nd post-SE s3 low to medium  & 2nd post-SE s12 low to very high  & 2nd post-SE s12 low to very high 
      pred                                                             prob sumtempProb
12921   NC                                  0.3424897, 0.4924829, 0.4497974   0.4282567
11786   NC            0.5062743, 0.5911322, 0.5236123, 0.5251700, 0.4948555   0.5282089
4983     C                                  0.3542178, 0.5605006, 0.5435639   0.4860941
2013     C                       0.5774260, 0.3148476, 0.4821741, 0.6093345   0.4959456
7718    NC 0.5216339, 0.5543512, 0.4942104, 0.2251983, 0.4688597, 0.4517625   0.4526693
                                                  Ruletype
12921                         T1wdynamic & single-time-Enh
11786                 T1wdynamic & T1wtexture & dispersion
4983                 T1wmorphology & single-time-Enh & T2w
2013  T1wdynamic & T1wmorphology & T1wtexture & dispersion
7718             T1wtexture & dispersion & single-time-Enh
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   159 16 54 57 71 lesion_id lesion_label lesion_diagnosis
21   1  1  0  0  1       195     nonmassM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 6119 174  18  38  55 
  4   3   2   2   2 
      len  freq                err
13980   3 0.173              0.083
2209    3 0.076              0.095
12940   4 0.049                  0
10957   4 0.038              0.048
7246    2 0.059 0.0610000000000001
                                                                                                                                                      condition
13980                                                                          irregularity low  & Energy ptime4 low to medium  & dispersion s17 low to medium 
2209                                                enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
12940                           Uptake(in) curve amplitude medium to high  & SER(in) low to medium  & Variance ptime3 low to medium  & dispersion s11 very low 
10957 Rate of signal increase(rim) low to medium  & Difference variance ptime2 very low to high  & dispersion s7 very low  & last post-SE s15 low to very high 
7246                                                                                                         SER(in) low to medium  & last post-SE s2 very low 
      pred                                       prob sumtempProb
13980   NC            0.5332865, 0.3120332, 0.4861676   0.4438291
2209    NC            0.3132262, 0.5008952, 0.5531483   0.4557566
12940   NC 0.4874018, 0.5424656, 0.3437655, 0.4102728   0.4459764
10957   NC 0.4644762, 0.4463054, 0.3475103, 0.3686232   0.4067288
7246    NC                       0.5333206, 0.4317962   0.4825584
                                                    Ruletype
13980                T1wmorphology & T1wtexture & dispersion
2209               T1wdynamic & dispersion & single-time-Enh
12940                   T1wdynamic & T1wtexture & dispersion
10957 T1wdynamic & T1wtexture & dispersion & single-time-Enh
7246                            T1wdynamic & single-time-Enh
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   119 174 18 38 55 lesion_id lesion_label lesion_diagnosis
34   1   0  1  0  1       377     nonmassB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 7 82 174  36  85  89 
  4   3   2   2   2 
     len  freq   err
5316   5 0.168 0.075
2209   3 0.076 0.095
4833   5 0.038 0.048
3314   5  0.07 0.077
7111   4 0.094 0.077
                                                                                                                                                                             condition
5316 Initial Uptake slope(in) low  & uptake average low to very high  & 1st post-SE s10 low to very high  & T2w std SI low to medium  & T2w margin gradient variance very low to high 
2209                                                                       enhancement-variance increasing rate(rim) low to medium  & dispersion s16 very low  & last post-SE s8 high 
4833          Max Radial gradient variance very low  & Sum average ptime3 medium to high  & Entropy ptime4 very low to high  & dispersion s18 low  & 2nd post-SE s15 very low to high 
3314                                   SER(rim) very low to high  & dispersion s14 very low  & dispersion s15 very low to high  & 1st post-SE s4 low to very high  & T2w Entropy  low 
7111                              Initial Uptake slope(in) low to medium  & Curvature at Time-to-peak(in) low to very high  & dispersion s10 low to very high  & 2nd post-SE s11 high 
     pred                                                  prob sumtempProb
5316   NC 0.4508790, 0.4713496, 0.5509649, 0.5310208, 0.5073532   0.5023135
2209   NC                       0.3132262, 0.5008952, 0.5531483   0.4557566
4833   NC 0.4289456, 0.6226296, 0.5110246, 0.5385531, 0.4767274   0.5155761
3314   NC 0.4499065, 0.4352475, 0.4759258, 0.4943733, 0.5792785   0.4869463
7111   NC            0.4452511, 0.5362599, 0.3412410, 0.4844197   0.4517929
                                                      Ruletype
5316                        T1wdynamic & single-time-Enh & T2w
2209                 T1wdynamic & dispersion & single-time-Enh
4833 T1wmorphology & T1wtexture & dispersion & single-time-Enh
3314           T1wdynamic & dispersion & single-time-Enh & T2w
7111                 T1wdynamic & dispersion & single-time-Enh
   82 174 36 85 89 lesion_id lesion_label lesion_diagnosis
9   1   1  1  0  1       135        massB      FIBROCYSTIC
26  1   1  1  1  0       302        massB      FIBROCYSTIC
28  0   1  0  1  0       317     nonmassB      FIBROCYSTIC
31  0   0  0  0  0       338     nonmassB      FIBROCYSTIC
36  0   0  0  0  0       396     nonmassB      FIBROCYSTIC
47  1   0  0  0  0       558        massB      FIBROCYSTIC
53  1   0  0  0  1       613     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
  82 174 36 85 89 lesion_id lesion_label lesion_diagnosis
9  1   1  1  0  1       135        massB      FIBROCYSTIC

============== Done cv-fold: 4 ==============

============== Reading cv-fold: 5 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 13 2  7 10 15 16 
 3  2  2  2  2 
     len  freq   err
1342   3 0.053 0.034
3290   2 0.024 0.077
1063   3 0.067 0.081
858    2 0.062 0.088
3902   2 0.124 0.088
                                                                                                          condition
1342 enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
3290                                                Washout rate(in) very low to high  & 3rd post-SE s11 very high 
1063                    irregularity high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
858                                                SER(in) medium to high  & Max Radial gradient variance very low 
3902                                                circularity high & T2w radial gradient variance medium to high 
     pred                            prob sumtempProb
1342   NC 0.4848342, 0.5546800, 0.6219773   0.5538305
3290   NC            0.5027937, 0.4296802    0.466237
1063    C 0.5128461, 0.6452170, 0.4944342   0.5508324
858    NC            0.4812519, 0.5830327   0.5321423
3902   NC            0.5514954, 0.4847561   0.5181258
                                         Ruletype
1342             T1wdynamic & T1wmorphology & T2w
3290                 T1wdynamic & single-time-Enh
1063 T1wmorphology & T1wtexture & single-time-Enh
858                    T1wdynamic & T1wmorphology
3902                          T1wmorphology & T2w
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   2 7 10 15 16 lesion_id lesion_label lesion_diagnosis
15 0 1  1  0  0       200        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 9 5  6  7 21 25 
 2  2  2  2  2 
     len  freq                err
715    3  0.08 0.0679999999999999
2980   3 0.076              0.071
3290   2 0.024              0.077
478    3 0.156              0.093
2203   3 0.093              0.098
                                                                                                                    condition
715                             uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
2980                                                     uptake average high  & 1st post-SE s8 high  & last post-SE s11 high 
3290                                                          Washout rate(in) very low to high  & 3rd post-SE s11 very high 
478  irregularity low  & std 3D Sharpness of lesion margin  medium to high  & Inverse difference moment ptime3 low to medium 
2203                            uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
     pred                            prob sumtempProb
715     C 0.5436611, 0.6397459, 0.5003391   0.5612487
2980    C 0.5757935, 0.5294741, 0.6741425   0.5931367
3290   NC            0.5027937, 0.4296802    0.466237
478    NC 0.5702187, 0.4669343, 0.3339196   0.4570242
2203    C 0.5260448, 0.6336959, 0.5132996   0.5576801
                                      Ruletype
715  T1wdynamic & T1wtexture & single-time-Enh
2980              T1wdynamic & single-time-Enh
3290              T1wdynamic & single-time-Enh
478                 T1wmorphology & T1wtexture
2203 T1wdynamic & T1wtexture & single-time-Enh
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   5 6 7 21 25 lesion_id lesion_label lesion_diagnosis
51 1 1 1  0  1       549        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 222  3 18  1  2 
 2  1  1  0  0 
     len  freq                err
726    2 0.215              0.093
1581   3 0.096 0.0570000000000001
2249   2  0.04              0.091
3843   2 0.047                  0
1342   3 0.053              0.034
                                                                                                          condition
726                                                                      SER(in) low  & T2w Energy  medium to high 
1581                     Difference variance ptime1 low  & dispersion s16 low to medium  & T2w margin gradient low 
2249                                                  1st post-SE s19 very low  & T2w radial gradient variance low 
3843                                                Time-to-peak(in) very high  & T2w radial gradient variance low 
1342 enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
     pred                            prob sumtempProb                         Ruletype
726    NC            0.4881946, 0.2745578   0.3813762                 T1wdynamic & T2w
1581   NC 0.5724939, 0.3448565, 0.4659352   0.4610952    T1wtexture & dispersion & T2w
2249   NC            0.4717060, 0.6599352   0.5658206            single-time-Enh & T2w
3843   NC            0.4480160, 0.5904218   0.5192189                 T1wdynamic & T2w
1342   NC 0.4848342, 0.5546800, 0.6219773   0.5538305 T1wdynamic & T1wmorphology & T2w
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   22 3 18 1 2 lesion_id lesion_label lesion_diagnosis
31  1 0  1 0 0       345        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 5 2 16 22  5  8 
 2  2  2  1  1 
     len  freq                err
1342   3 0.053              0.034
3902   2 0.124              0.088
726    2 0.215              0.093
715    3  0.08 0.0679999999999999
1977   2 0.045               0.08
                                                                                                          condition
1342 enhancement-variance at first time-point(in) very low  & irregularity very low  & max T2w SI low to very high 
3902                                                circularity high & T2w radial gradient variance medium to high 
726                                                                      SER(in) low  & T2w Energy  medium to high 
715                   uptake average high  & Inverse difference moment ptime2 low to medium  & 1st post-SE s8 high 
1977                                                                 dispersion s5 very low  & 1st post-SE s11 low 
     pred                            prob sumtempProb
1342   NC 0.4848342, 0.5546800, 0.6219773   0.5538305
3902   NC            0.5514954, 0.4847561   0.5181258
726    NC            0.4881946, 0.2745578   0.3813762
715     C 0.5436611, 0.6397459, 0.5003391   0.5612487
1977   NC            0.4735013, 0.6948949   0.5841981
                                      Ruletype
1342          T1wdynamic & T1wmorphology & T2w
3902                       T1wmorphology & T2w
726                           T1wdynamic & T2w
715  T1wdynamic & T1wtexture & single-time-Enh
1977              dispersion & single-time-Enh
   2 16 22 5 8 lesion_id lesion_label lesion_diagnosis
1  0  0  0 1 0        10        massB      FIBROCYSTIC
25 0  0  0 0 0       255        massB      FIBROCYSTIC
26 1  1  0 0 0       256        massB      FIBROCYSTIC
27 1  1  1 0 1       274        massB      FIBROCYSTIC
59 0  0  1 0 0       636     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   2 16 22 5 8 lesion_id lesion_label lesion_diagnosis
27 1  1  1 0 1       274        massB      FIBROCYSTIC

============== Done cv-fold: 5 ==============

============== Reading cv-fold: 6 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 1134 15 17 18 21 
 5  2  2  2  2 
     len  freq                err
184    2 0.156              0.094
2275   2 0.059 0.0620000000000001
5282   2 0.027 0.0669999999999999
40     2 0.132              0.069
6838   2 0.046               0.08
                                                               condition pred
184                         SER(in) medium to high  & irregularity high     C
2275                    dispersion s11 high  & 1st post-SE s18 very low    NC
5282 Sum average ptime1 very low  & Sum average ptime2 low to very high    NC
40                            irregularity high  & Variance ptime1 high     C
6838         Sum average ptime1 very low  & T2w Entropy  medium to high    NC
                     prob sumtempProb                     Ruletype
184  0.7579719, 0.4714506   0.6147112   T1wdynamic & T1wmorphology
2275 0.5553819, 0.4554229   0.5054024 dispersion & single-time-Enh
5282 0.5302071, 0.5473805   0.5387938                   T1wtexture
40   0.5017755, 0.7615645     0.63167   T1wmorphology & T1wtexture
6838 0.5821745, 0.4920792   0.5371268             T1wtexture & T2w
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   34 15 17 18 21 lesion_id lesion_label lesion_diagnosis
40  1  1  1  1  1       341        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 727 29 35  7 10 
 2  2  2  1  1 
     len  freq   err
7002   2 0.084 0.087
972    3 0.101 0.091
3648   2 0.059 0.094
7014   3 0.092  0.04
2382   3 0.035 0.053
                                                                                                   condition
7002                                            2nd post-SE s4 low to very high  & last post-SE s4 very low 
972                           SER(in) medium to high  & irregularity high  & Difference entropy ptime3 high 
3648                                             dispersion s1 low to very high  & 1st post-SE s8 very high 
7014           SER(in) low  & Energy ptime1 medium to high  & T2w radial gradient variance very low to high 
2382 uptake variance very low to high  & 2nd post-SE s10 very low to high  & T2w average radial gradient low
     pred                            prob sumtempProb                                Ruletype
7002   NC            0.4880460, 0.4469655   0.4675058                         single-time-Enh
972     C 0.4542278, 0.3886270, 0.5466488   0.4631679 T1wdynamic & T1wmorphology & T1wtexture
3648    C            0.4857548, 0.5186936   0.5022242            dispersion & single-time-Enh
7014   NC 0.6349908, 0.5247829, 0.4750272   0.5449336           T1wdynamic & T1wtexture & T2w
2382   NC 0.5287212, 0.4843964, 0.4342890   0.4824689      T1wdynamic & single-time-Enh & T2w
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   27 29 35 7 10 lesion_id lesion_label lesion_diagnosis
16  1  0  0 0  1       131        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 1016 33  2  8 12 
 3  3  2  2  2 
     len  freq                err
1578   3 0.057 0.0649999999999999
2330   3 0.216              0.093
3243   3 0.018                  0
7103   2 0.044              0.042
459    2 0.064 0.0570000000000001
                                                                                                                condition
1578 irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
2330                                   irregularity low  & Energy ptime3 low to medium  & dispersion s5 low to very high 
3243               Contrast ptime3 low  & Difference variance ptime3 medium to high  & T2w Sum average  low to very high 
7103                                                irregularity very low  & T2w radial gradient variance medium to high 
459                                                                irregularity very low  & dispersion s9 medium to high 
     pred                            prob sumtempProb                                Ruletype
1578   NC 0.5914519, 0.4717202, 0.4894149    0.517529        T1wmorphology & T1wtexture & T2w
2330   NC 0.4678680, 0.2644552, 0.4902931   0.4075388 T1wmorphology & T1wtexture & dispersion
3243    C 0.5483842, 0.4931628, 0.5762954   0.5392808                        T1wtexture & T2w
7103   NC            0.5748740, 0.4661645   0.5205192                     T1wmorphology & T2w
459    NC            0.6755135, 0.4776196   0.5765665              T1wmorphology & dispersion
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   16 33 2 8 12 lesion_id lesion_label lesion_diagnosis
44  1  1 0 1  1       366        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 4 4  8 16 36 41 
 1  1  1  1  1 
     len  freq                err
178    3 0.048              0.038
7103   2 0.044              0.042
1578   3 0.057 0.0649999999999999
2931   3 0.038              0.095
190    3 0.057              0.097
                                                                                                                condition
178                                               SER(in) low  & irregularity medium to high  & last post-SE s2 very low 
7103                                                irregularity very low  & T2w radial gradient variance medium to high 
1578 irregularity very low  & Difference variance ptime3 very low to high  & T2w radial gradient variance medium to high 
2931            Rate of signal increase(rim) very low  & dispersion s17 low to very high  & T2w mean SI low to very high 
190                                    SER(rim) low to medium  & irregularity medium to high  & last post-SE s2 very low 
     pred                            prob sumtempProb
178    NC 0.6131108, 0.4701726, 0.3805070   0.4879301
7103   NC            0.5748740, 0.4661645   0.5205192
1578   NC 0.5914519, 0.4717202, 0.4894149    0.517529
2931   NC 0.5582364, 0.5165599, 0.5742535   0.5496832
190    NC 0.5743398, 0.4464800, 0.3780448   0.4662882
                                         Ruletype
178  T1wdynamic & T1wmorphology & single-time-Enh
7103                          T1wmorphology & T2w
1578             T1wmorphology & T1wtexture & T2w
2931                T1wdynamic & dispersion & T2w
190  T1wdynamic & T1wmorphology & single-time-Enh
   4 8 16 36 41 lesion_id lesion_label lesion_diagnosis
9  1 0  0  0  1        75        massB      FIBROCYSTIC
10 0 1  1  1  0        76        massB      FIBROCYSTIC
31 0 0  0  0  0       236     nonmassB      FIBROCYSTIC
58 0 0  0  0  0       524        massB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   4 8 16 36 41 lesion_id lesion_label lesion_diagnosis
10 0 1  1  1  0        76        massB      FIBROCYSTIC

============== Done cv-fold: 6 ==============

============== Reading cv-fold: 7 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 1831 10 38 58 95 
 7  6  6  6  6 
      len  freq                err
68      3 0.111              0.049
10957   5 0.051                  0
9295    3 0.073               0.05
112     2 0.089 0.0610000000000001
2791    3 0.117              0.078
                                                                                                                                                               condition
68                                                                                  Initial Uptake slope(in) high  & uptake skewness low  & irregularity medium to high 
10957 Time-to-peak(in) low  & Variance of spatial Margin Gradient very low to high  & circularity low & dispersion s12 very low to high  & dispersion s17 low to medium 
9295                                                                                                  Time-to-peak(in) low  & irregularity high  & 1st post-SE s10 high 
112                                                                                                                  irregularity very high  & Sum variance ptime1 high 
2791                                                                            Time-to-peak(in) low  & circularity low & Max Radial gradient variance low to very high 
      pred                                                  prob sumtempProb
68       C                       0.4479334, 0.7561858, 0.6379134   0.6140108
10957    C 0.5682319, 0.6838746, 0.4879725, 0.4530955, 0.4300886   0.5246526
9295     C                       0.7741849, 0.4788727, 0.7091933   0.6540836
112      C                                  0.7075908, 0.4421720   0.5748814
2791     C                       0.6633012, 0.5110174, 0.6863620   0.6202269
                                          Ruletype
68                      T1wdynamic & T1wmorphology
10957      T1wdynamic & T1wmorphology & dispersion
9295  T1wdynamic & T1wmorphology & single-time-Enh
112                     T1wmorphology & T1wtexture
2791                    T1wdynamic & T1wmorphology
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   31 10 38 58 95 lesion_id lesion_label lesion_diagnosis
15  1  1  1  1  1       154        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 4 96 111   5   6  29 
  2   2   1   1   1 
      len  freq   err
841     2 0.139 0.079
12938   4 0.151 0.084
4715    4 0.027     0
5887    4  0.04     0
9188    5 0.038 0.048
                                                                                                                                                                    condition
841                                                                                                                         SER(rim) very low  & dispersion s9 low to medium 
12938                                                  Curvature at Time-to-peak(rim) high & SER(rim) very low  & max uptake very low to high  & max T2w SI very low to high 
4715                                                 2nd post-SE s18 low to very high  & 3rd post-SE s4 high  & last post-SE s8 low to very high  & predicted LMSIR very low 
5887          enhancement-variance increasing rate(in) low to very high  & max uptake very high  & max Radial gradient low to very high  & last post-SE s11 low to very high 
9188  max uptake medium to high  & Entropy ptime2 low to very high  & 1st post-SE s16 low to medium  & 3rd post-SE s10 low to very high  & T2w Difference variance  very low 
      pred                                                  prob sumtempProb
841     NC                                  0.4450006, 0.5071265   0.4760636
12938   NC            0.4801385, 0.4899903, 0.5160871, 0.4976432   0.4959648
4715     C            0.5430140, 0.7717451, 0.5008960, 0.5124503   0.5820264
5887     C            0.4887420, 0.4518504, 0.5100942, 0.4596818   0.4775921
9188    NC 0.4555243, 0.4955115, 0.5378557, 0.4726280, 0.4208891   0.4764817
                                             Ruletype
841                           T1wdynamic & dispersion
12938                                T1wdynamic & T2w
4715                            single-time-Enh & T2w
5887     T1wdynamic & T1wmorphology & single-time-Enh
9188  T1wdynamic & T1wtexture & single-time-Enh & T2w
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   96 111 5 6 29 lesion_id lesion_label lesion_diagnosis
13  1   1 0 0  0       113     nonmassM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 3 32  67 112 103 107 
  3   3   3   2   2 
      len  freq                err
3013    3 0.075              0.049
4569    4 0.109 0.0669999999999999
1128    2 0.215              0.085
11043   3 0.111              0.082
4925    2 0.066              0.083
                                                                                                          condition
3013        irregularity very low  & dispersion s3 low to very high  & T2w radial gradient variance medium to high 
4569  uptake skewness low  & irregularity low  & Entropy ptime3 low to very high  & dispersion s5 low to very high 
1128                                                               irregularity low  & Energy ptime4 low to medium 
11043                    Sum Entropy ptime1 low to medium  & 1st post-SE s8 low to medium  & last post-SE s17 high 
4925                                          irregularity very low  & T2w radial gradient variance medium to high 
      pred                                       prob sumtempProb
3013    NC            0.5124585, 0.5480302, 0.4529533   0.5044807
4569    NC 0.4028134, 0.5669529, 0.4266298, 0.4695194   0.4664789
1128    NC                       0.4553747, 0.2413562   0.3483655
11043   NC            0.4335317, 0.4984163, 0.3944782   0.4421421
4925    NC                       0.5857864, 0.4576238   0.5217051
                                                  Ruletype
3013                      T1wmorphology & dispersion & T2w
4569  T1wdynamic & T1wmorphology & T1wtexture & dispersion
1128                            T1wmorphology & T1wtexture
11043                         T1wtexture & single-time-Enh
4925                                   T1wmorphology & T2w
[1] "Top case meeting rules for FIBROADENOMA (FA):"
  32 67 112 103 107 lesion_id lesion_label lesion_diagnosis
9  1  1   1   1   1        78        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 7133  13  19 155  15 
  4   3   3   3   2 
      len  freq   err
5658    3 0.177 0.093
1233    4 0.066 0.028
1593    6 0.047 0.038
9527    2 0.057 0.097
10203   4 0.057 0.032
                                                                                                                                                                               condition
5658                                                                                                SER(rim) low  & Difference variance ptime1 low  & T2w Correlation  low to very high 
1233                                    enhancement-variance increasing rate(rim) low to medium  & Max Margin Gradient low  & Entropy ptime1 low  & Sum average ptime2 very low to high 
1593  circularity high & dispersion s5 low to very high  & dispersion s10 low to very high  & 3rd post-SE s7 very low to high  & 3rd post-SE s12 medium to high  & last post-SE s17 low 
9527                                                                     max relative signal enhancement(rim) medium to high  & change in Variance of spatial Margin Gradient very high 
10203                 Max Radial gradient variance very low  & Sum average ptime1 medium to high  & 1st post-SE s16 low to very high  & T2w Inverse difference moment  low to very high 
      pred                                                             prob sumtempProb
5658    NC                                  0.5505793, 0.3818076, 0.5017071   0.4780313
1233    NC                       0.4810853, 0.2872408, 0.4052440, 0.3282516   0.3754554
1593    NC 0.6372325, 0.7199635, 0.4684946, 0.4917864, 0.5568750, 0.4652017   0.5565923
9527    NC                                             0.6236805, 0.4747261   0.5492033
10203   NC                       0.6231637, 0.2928371, 0.5364463, 0.5693335   0.5054452
                                                Ruletype
5658                       T1wdynamic & T1wtexture & T2w
1233             T1wdynamic & T1wmorphology & T1wtexture
1593        T1wmorphology & dispersion & single-time-Enh
9527                          T1wdynamic & T1wmorphology
10203 T1wmorphology & T1wtexture & single-time-Enh & T2w
   133 13 19 155 15 lesion_id lesion_label lesion_diagnosis
10   0  0  1   0  0        79     nonmassB      FIBROCYSTIC
14   0  0  0   1  1       137        massB      FIBROCYSTIC
19   1  1  0   1  0       166        massB      FIBROCYSTIC
24   1  1  0   0  1       193     nonmassB      FIBROCYSTIC
36   1  0  0   1  0       348     nonmassB      FIBROCYSTIC
39   0  0  1   0  0       356     nonmassB      FIBROCYSTIC
58   1  1  1   0  0       579     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   133 13 19 155 15 lesion_id lesion_label lesion_diagnosis
19   1  1  0   1  0       166        massB      FIBROCYSTIC

============== Done cv-fold: 7 ==============

============== Reading cv-fold: 8 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 16 59 162  87 117  66 
  7   6   5   5   4 
      len  freq                err
11971   3 0.136 0.0679999999999999
454     3 0.132              0.097
6241    3 0.072              0.077
3487    3 0.153              0.084
442     3 0.237               0.07
                                                                                       condition
11971           Correlation ptime2 low to medium  & dispersion s2 low  & dispersion s2 very low 
454       Time-to-peak(in) low  & irregularity medium to high  & 1st post-SE s12 medium to high 
6241  uptake skewness low  & max Radial gradient very high  & T2w Sum Entropy  low to very high 
3487           circularity high & Variance ptime1 low to medium  & Entropy ptime2 low to medium 
442         irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
      pred                            prob sumtempProb
11971   NC 0.3216921, 0.4685461, 0.5324010   0.4408798
454      C 0.7569990, 0.5009352, 0.6581838    0.638706
6241    NC 0.4645884, 0.6111328, 0.4870587   0.5209266
3487    NC 0.4171425, 0.4806995, 0.3749385   0.4242602
442     NC 0.4378788, 0.2352466, 0.2721108   0.3150787
                                          Ruletype
11971                      T1wtexture & dispersion
454   T1wdynamic & T1wmorphology & single-time-Enh
6241              T1wdynamic & T1wmorphology & T2w
3487                    T1wmorphology & T1wtexture
442                     T1wmorphology & T1wtexture
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   59 162 87 117 66 lesion_id lesion_label lesion_diagnosis
45  1   0  1   1  1       444        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 8156  41  59  66  70 
  4   3   3   3   3 
      len  freq                err
7118    2 0.039              0.095
13591   4  0.09 0.0610000000000001
11971   3 0.136 0.0679999999999999
442     3 0.237               0.07
2220    2 0.051              0.071
                                                                                                                                           condition
7118                                                                                             Sum average ptime2 very high  & dispersion s17 low 
13591 enhancement-variance decreasing rate(rim) very low to high  & irregularity low  & Sum average ptime1 high  & 2nd post-SE s14 very low to high 
11971                                                               Correlation ptime2 low to medium  & dispersion s2 low  & dispersion s2 very low 
442                                                             irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
2220                                                                                              dispersion s5 very low  & dispersion s19 very low 
      pred                                       prob sumtempProb
7118    NC                       0.4785108, 0.7332461   0.6058785
13591   NC 0.4442097, 0.5323714, 0.4614969, 0.5905932   0.5071678
11971   NC            0.3216921, 0.4685461, 0.5324010   0.4408798
442     NC            0.4378788, 0.2352466, 0.2721108   0.3150787
2220    NC                       0.4927936, 0.6679542   0.5803739
                                                       Ruletype
7118                                    T1wtexture & dispersion
13591 T1wdynamic & T1wmorphology & T1wtexture & single-time-Enh
11971                                   T1wtexture & dispersion
442                                  T1wmorphology & T1wtexture
2220                                                 dispersion
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   156 41 59 66 70 lesion_id lesion_label lesion_diagnosis
31   1  1  1  1  0       277        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 5 59  66 134  37  96 
  3   3   3   2   2 
      len  freq                err
11971   3 0.136 0.0679999999999999
442     3 0.237               0.07
5343    4 0.142              0.091
1541    2 0.162 0.0570000000000001
183     2 0.232              0.079
                                                                                                                     condition
11971                                         Correlation ptime2 low to medium  & dispersion s2 low  & dispersion s2 very low 
442                                       irregularity low  & Variance ptime2 low to medium  & Energy ptime4 very low to high 
5343  last post-SE s1 low  & last post-SE s9 low to medium  & T2w mean SI low to very high  & T2w average radial gradient low 
1541                                                                          irregularity low  & Energy ptime4 low to medium 
183                                                        irregularity low  & Inverse difference moment ptime4 low to medium 
      pred                                       prob sumtempProb                   Ruletype
11971   NC            0.3216921, 0.4685461, 0.5324010   0.4408798    T1wtexture & dispersion
442     NC            0.4378788, 0.2352466, 0.2721108   0.3150787 T1wmorphology & T1wtexture
5343    NC 0.5256552, 0.4322540, 0.4772128, 0.4044395   0.4598904      single-time-Enh & T2w
1541    NC                       0.4792871, 0.2735542   0.3764207 T1wmorphology & T1wtexture
183     NC                       0.4863329, 0.2652150   0.3757739 T1wmorphology & T1wtexture
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   59 66 134 37 96 lesion_id lesion_label lesion_diagnosis
29  0  1   1  1  1       272        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 220 23 37 38 40 
 1  1  1  1  1 
      len  freq                err
10158   2 0.048              0.038
12721   3 0.127              0.043
1541    2 0.162 0.0570000000000001
12632   6 0.096 0.0580000000000001
6591    4 0.031 0.0590000000000001
                                                                                                                                                                                                                  condition
10158                                                                                                                                                            Energy ptime1 very low to high  & Entropy ptime2 very low 
12721                                                                                                                                       uptake skewness low  & irregularity very low  & dispersion s1 low to very high 
1541                                                                                                                                                                       irregularity low  & Energy ptime4 low to medium 
12632 Area Under Uptake curve(in) very low to high  & Variance ptime1 very low  & Energy ptime4 very low to high  & dispersion s16 low to very high  & 2nd post-SE s17 low to very high  & 3rd post-SE s6 low to very high 
6591                                                                 Rate of signal increase(in) low  & Rate of signal increase(in) medium to high  & 1st post-SE s15 low to very high  & last post-SE s8 low to very high 
      pred                                                             prob sumtempProb
10158   NC                                             0.5127610, 0.4835754   0.4981682
12721   NC                                  0.4677355, 0.5728223, 0.2321752   0.4242443
1541    NC                                             0.4792871, 0.2735542   0.3764207
12632   NC 0.4755872, 0.5433402, 0.4788114, 0.4381729, 0.4557238, 0.5067061   0.4830569
6591    NC                       0.5271721, 0.3903370, 0.5439983, 0.4908403    0.488087
                                                    Ruletype
10158                                             T1wtexture
12721                T1wdynamic & T1wmorphology & dispersion
1541                              T1wmorphology & T1wtexture
12632 T1wdynamic & T1wtexture & dispersion & single-time-Enh
6591                            T1wdynamic & single-time-Enh
   20 23 37 38 40 lesion_id lesion_label lesion_diagnosis
19  1  0  0  1  0       211        massB      FIBROCYSTIC
67  0  1  1  0  1       627        massB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   20 23 37 38 40 lesion_id lesion_label lesion_diagnosis
67  0  1  1  0  1       627        massB      FIBROCYSTIC

============== Done cv-fold: 8 ==============

============== Reading cv-fold: 9 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 618 32 45  1  3 
 2  2  2  1  1 
     len  freq                err
3868   3 0.031 0.0590000000000001
3524   2 0.043              0.087
124    3 0.133              0.097
1193   3  0.03                  0
7946   2  0.03                  0
                                                                                                     condition
3868                       dispersion s15 low to very high  & 3rd post-SE s4 high  & predicted LMSIR very low 
3524                            change in Variance of spatial Margin Gradient low  & dispersion s17 very high 
124        Time-to-peak(in) low  & Rate of signal increase(rim) medium to high  & irregularity medium to high 
1193          Time-to-peak(rim) medium to high  & max Radial gradient low  & T2w radial gradient variance low 
7946 change in Variance of spatial Margin Gradient low & Variance of spatial Margin Gradient very low to high 
     pred                            prob sumtempProb                           Ruletype
3868    C 0.4853810, 0.6962119, 0.4631950   0.5482627 dispersion & single-time-Enh & T2w
3524    C            0.3885357, 0.5447962    0.466666         T1wmorphology & dispersion
124     C 0.7280304, 0.4252003, 0.5877582   0.5803296         T1wdynamic & T1wmorphology
1193   NC 0.4662758, 0.3631702, 0.5977509   0.4757323   T1wdynamic & T1wmorphology & T2w
7946   NC            0.5439387, 0.5092900   0.5266144                      T1wmorphology
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
  18 32 45 1 3 lesion_id lesion_label lesion_diagnosis
4  1  1  0 1 0         9     nonmassM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 12 2  5 30 47  6 
 2  2  2  2  1 
     len  freq   err
7536   3 0.031     0
4971   4 0.055 0.033
6675   3 0.087 0.085
315    2 0.152 0.098
5728   3  0.05 0.037
                                                                                                                             condition
7536                                                         SER(rim) very low  & 1st post-SE s6 very low  & last post-SE s7 very low 
4971 min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
6675                      Sum variance ptime1 high  & Inverse difference moment ptime4 medium to high  & dispersion s2 medium to high 
315                                                                                         irregularity high  & Variance ptime1 high 
5728                                                 uptake variance very high  & last post-SE s6 low to very high  & T2w Energy  low 
     pred                                       prob sumtempProb
7536   NC            0.6420296, 0.5124139, 0.4394809   0.5313081
4971   NC 0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
6675    C            0.3978654, 0.7406423, 0.5150453   0.5511844
315     C                       0.4852856, 0.7650575   0.6251715
5728    C            0.7332721, 0.4739669, 0.5252209   0.5774866
                                                 Ruletype
7536                         T1wdynamic & single-time-Enh
4971 T1wdynamic & T1wmorphology & T1wtexture & dispersion
6675                              T1wtexture & dispersion
315                            T1wmorphology & T1wtexture
5728                   T1wdynamic & single-time-Enh & T2w
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   2 5 30 47 6 lesion_id lesion_label lesion_diagnosis
37 1 0  0  0 1       298     nonmassM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 8 7 13 20 31 35 
 2  2  2  2  2 
     len  freq                err
5832   3 0.048              0.038
3144   3 0.105              0.053
9672   3 0.061 0.0610000000000001
1538   4 0.085              0.087
843    2 0.063              0.088
                                                                                                                                   condition
5832                                  Initial Uptake slope(in) very low to high  & dispersion s6 very low  & last post-SE s8 medium to high 
3144                                               irregularity very low  & max T2w SI low to very high  & kurtosis T2w SI low to very high 
9672                                         dispersion s0 very low to high  & 1st post-SE s16 very low  & last post-SE s4 low to very high 
1538 Uptake(in) curve amplitude medium to high  & Sum Entropy ptime4 low to medium  & dispersion s17 low to very high  & 3rd post-SE s4 low 
843                                                                    irregularity very low  & T2w radial gradient variance medium to high 
     pred                                       prob sumtempProb
5832   NC            0.3299171, 0.5255908, 0.5713447   0.4756175
3144   NC            0.4195477, 0.4654747, 0.4501588   0.4450604
9672   NC            0.5325642, 0.5042174, 0.4817089   0.5061635
1538   NC 0.6057217, 0.5218727, 0.4647596, 0.4578509   0.5125512
843    NC                       0.6537632, 0.5159319   0.5848476
                                                   Ruletype
5832              T1wdynamic & dispersion & single-time-Enh
3144                                    T1wmorphology & T2w
9672                           dispersion & single-time-Enh
1538 T1wdynamic & T1wtexture & dispersion & single-time-Enh
843                                     T1wmorphology & T2w
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   7 13 20 31 35 lesion_id lesion_label lesion_diagnosis
10 1  1  0  0  1        49        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 1244  4 27  2  5 
 3  2  2  1  1 
     len  freq   err
5386   3 0.116 0.095
9109   4 0.041     0
2076   2 0.052 0.071
7536   3 0.031     0
4971   4 0.055 0.033
                                                                                                                                        condition
5386 SER(in) very low  & enhancement-variance increasing rate(in) very low to high  & enhancement-variance decreasing rate(rim) very low to high 
9109                 Initial Uptake slope(in) low  & uptake average high  & 1st post-SE s15 low to very high  & 1st post-SE s16 low to very high 
2076                                                                            SER(in) low  & enhancement-variance increasing rate(in) very low 
7536                                                                    SER(rim) very low  & 1st post-SE s6 very low  & last post-SE s7 very low 
4971            min uptake very low to high  & mean 3D Sharpness of lesion margin  very low  & Energy ptime4 low to medium  & dispersion s16 low 
     pred                                       prob sumtempProb
5386   NC            0.5192031, 0.5443866, 0.4975629   0.5203842
9109   NC 0.6920559, 0.4833696, 0.7240647, 0.5912338    0.622681
2076   NC                       0.6114913, 0.4662233   0.5388573
7536   NC            0.6420296, 0.5124139, 0.4394809   0.5313081
4971   NC 0.4686661, 0.3111970, 0.4932763, 0.4507027   0.4309605
                                                 Ruletype
5386                                           T1wdynamic
9109                         T1wdynamic & single-time-Enh
2076                                           T1wdynamic
7536                         T1wdynamic & single-time-Enh
4971 T1wdynamic & T1wmorphology & T1wtexture & dispersion
   44 4 27 2 5 lesion_id lesion_label lesion_diagnosis
5   0 0  0 0 0        26        massB      FIBROCYSTIC
6   0 0  0 0 0        27     nonmassB      FIBROCYSTIC
8   0 1  1 0 0        31        massB      FIBROCYSTIC
9   0 0  0 0 1        48        massB      FIBROCYSTIC
29  0 0  0 0 0       263        massB      FIBROCYSTIC
30  0 0  0 0 0       266     nonmassB      FIBROCYSTIC
31  0 0  0 0 0       267     nonmassB      FIBROCYSTIC
39  0 0  0 0 0       307        massB      FIBROCYSTIC
48  0 0  0 0 0       426        massB      FIBROCYSTIC
50  1 0  0 0 0       432        massB      FIBROCYSTIC
59  1 0  0 1 0       537     nonmassB      FIBROCYSTIC
65  1 1  1 0 0       590     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   44 4 27 2 5 lesion_id lesion_label lesion_diagnosis
65  1 1  1 0 0       590     nonmassB      FIBROCYSTIC

============== Done cv-fold: 9 ==============

============== Reading cv-fold: 10 left out cases

============== Top five explaining rules for InvasiveDuctal (IDC): n = 1642 20 39 48 17 
 9  6  5  5  4 
     len  freq                err
135    2 0.115              0.095
670    3 0.079               0.07
7460   3 0.062              0.088
2589   2 0.095              0.096
4999   3 0.059 0.0620000000000001
                                                                                              condition
135                                                        SER(in) medium to high  & irregularity high 
670                          irregularity high  & Variance ptime2 high  & Sum Entropy ptime2 very high 
7460 uptake average high  & std 3D Sharpness of lesion margin  high  & last post-SE s18 medium to high 
2589                          Variance ptime2 very high  & Difference variance ptime2 very low to high 
4999                                uptake average high  & 1st post-SE s8 high  & T2w Correlation  low 
     pred                            prob sumtempProb
135     C            0.7467186, 0.4575242   0.6021214
670     C 0.5627317, 0.5135349, 0.6720011   0.5827559
7460    C 0.2956963, 0.4250158, 0.5423104   0.4210075
2589    C            0.4973411, 0.4768212   0.4870812
4999    C 0.3727659, 0.5435277, 0.4635614   0.4599517
                                         Ruletype
135                    T1wdynamic & T1wmorphology
670                    T1wmorphology & T1wtexture
7460 T1wdynamic & T1wmorphology & single-time-Enh
2589                                   T1wtexture
4999           T1wdynamic & single-time-Enh & T2w
[1] "Top case meeting rules for InvasiveDuctal (IDC):"
   42 20 39 48 17 lesion_id lesion_label lesion_diagnosis
12  1  1  1  1  1       160        massM   InvasiveDuctal

============== Top five explaining rules for InsituDuctal (ISDC): n = 842 21 13 20 48 
 4  3  2  2  2 
     len  freq                err
135    2 0.115              0.095
6220   3 0.104               0.07
138    3  0.09 0.0610000000000001
670    3 0.079               0.07
2589   2 0.095              0.096
                                                                                condition pred
135                                          SER(in) medium to high  & irregularity high     C
6220 Variance ptime1 very high  & dispersion s16 medium to high  & last post-SE s11 high     C
138                      SER(in) medium to high  & uptake skewness low  & circularity low    C
670            irregularity high  & Variance ptime2 high  & Sum Entropy ptime2 very high     C
2589            Variance ptime2 very high  & Difference variance ptime2 very low to high     C
                                prob sumtempProb                                  Ruletype
135             0.7467186, 0.4575242   0.6021214                T1wdynamic & T1wmorphology
6220 0.3621254, 0.7501070, 0.4513578   0.5211968 T1wtexture & dispersion & single-time-Enh
138  0.4310849, 0.6829253, 0.5595104   0.5578402                T1wdynamic & T1wmorphology
670  0.5627317, 0.5135349, 0.6720011   0.5827559                T1wmorphology & T1wtexture
2589            0.4973411, 0.4768212   0.4870812                                T1wtexture
[1] "Top case meeting rules for InsituDuctal (ISDC):"
   42 21 13 20 48 lesion_id lesion_label lesion_diagnosis
59  1  1  1  1  1       608        massM     InsituDuctal

============== Top five explaining rules for FIBROADENOMA (FA): n = 838 43 45  1 25 
 3  3  3  2  2 
     len  freq   err
7336   4 0.104 0.088
722    4 0.136 0.095
1758   3 0.115 0.095
2598   3 0.066 0.028
766    2   0.2 0.073
                                                                                                                                                              condition
7336 max relative signal enhancement(rim) very low to high  & Difference variance ptime1 very low to high  & Inverse difference moment ptime2 low  & dispersion s6 low 
722                                              uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
1758                                                              circularity high & Sum average ptime2 very low to high  & T2w average radial gradient medium to high 
2598                                                            irregularity very low  & dispersion s19 low to very high  & T2w average radial gradient medium to high 
766                                                                                                                    irregularity low  & Energy ptime4 low to medium 
     pred                                       prob sumtempProb
7336   NC 0.4074622, 0.5043908, 0.6186548, 0.4835740   0.5035205
722    NC 0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
1758   NC            0.5552296, 0.4969470, 0.4644038   0.5055268
2598   NC            0.4959181, 0.4310173, 0.4138270   0.4469208
766    NC                       0.4867494, 0.2974405    0.392095
                                               Ruletype
7336               T1wdynamic & T1wtexture & dispersion
722  T1wdynamic & T1wmorphology & single-time-Enh & T2w
1758                   T1wmorphology & T1wtexture & T2w
2598                   T1wmorphology & dispersion & T2w
766                          T1wmorphology & T1wtexture
[1] "Top case meeting rules for FIBROADENOMA (FA):"
   38 43 45 1 25 lesion_id lesion_label lesion_diagnosis
24  1  1  0 0  1       310        massB     FIBROADENOMA

============== Top five explaining rules for FIBROCYSTIC (FI): n = 843 32 50 10 49 
 4  3  3  2  2 
     len  freq   err
722    4 0.136 0.095
4749   4 0.066 0.083
1308   4 0.057 0.097
2822   2 0.071 0.051
4809   3 0.095 0.096
                                                                                                                                                     condition
722                                     uptake average high  & irregularity low to medium  & 1st post-SE s8 low to medium  & T2w radial gradient variance low 
4749                uptake skewness medium to high  & irregularity medium to high  & kurtosis T2w SI low to very high  & T2w average radial gradient very low 
1308 enhancement-variance at first time-point(rim) low to medium  & uptake variance medium to high  & max Radial gradient low  & dispersion s2 medium to high 
2822                                                                                            3rd post-SE s10 low to very high  & last post-SE s11 very low 
4809                                   enhancement-variance decreasing rate(rim) high  & Sum variance ptime2 low to medium  & last post-SE s10 medium to high 
     pred                                       prob sumtempProb
722    NC 0.4364043, 0.5007630, 0.5571853, 0.4727589   0.4917779
4749   NC 0.6209921, 0.4870514, 0.5075204, 0.4635613   0.5197813
1308   NC 0.5904771, 0.6088937, 0.5063315, 0.6877946   0.5983742
2822   NC                       0.5144900, 0.4830707   0.4987804
4809   NC            0.4806604, 0.4285987, 0.5949174   0.5013922
                                               Ruletype
722  T1wdynamic & T1wmorphology & single-time-Enh & T2w
4749                   T1wdynamic & T1wmorphology & T2w
1308            T1wdynamic & T1wmorphology & dispersion
2822                                    single-time-Enh
4809          T1wdynamic & T1wtexture & single-time-Enh
   43 32 50 10 49 lesion_id lesion_label lesion_diagnosis
3   0  0  0  0  1        25        massB      FIBROCYSTIC
4   0  0  0  0  0        29        massB      FIBROCYSTIC
7   0  0  0  0  1        56     nonmassB      FIBROCYSTIC
10  1  1  1  0  0       141        massB      FIBROCYSTIC
18  1  1  1  0  0       231     nonmassB      FIBROCYSTIC
28  0  0  0  1  0       330     nonmassB      FIBROCYSTIC
42  1  0  1  0  0       468        massB      FIBROCYSTIC
43  1  1  0  1  0       469     nonmassB      FIBROCYSTIC
[1] "Top case meeting rules for FIBROCYSTIC (FI):"
   43 32 50 10 49 lesion_id lesion_label lesion_diagnosis
10  1  1  1  0  0       141        massB      FIBROCYSTIC

============== Done cv-fold: 10 ==============
```

  
  

```r
# save current state k patient out
save.image("Outputs/allcvRules_noT2SIpredLMSIR_boost_addeddiagvalue.RData")
```
