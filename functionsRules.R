library(inTrees)
library(stringr)

identifyT2wrule<- function(RulesIx, colN, fnames) 
{
  # colN = colnames(X)
  RulesIxcond = presentRules(RulesIx, colnames(X))$condition
  
  rules = c()
  nameRules = c()
  for(ixr in RulesIxcond){
    # print(ixr)
    C = unlist(strsplit(ixr, split = " & "))
    lnameRule = strsplit(C[grepl("<", C)], split = "<")
    ltnameRule = strsplit(C[grepl("<=", C)], split = "<=")
    gnameRule = strsplit(C[grepl(">", C)], split = ">")
    gtnameRule = strsplit(C[grepl(">=", C)], split = ">=")
    
    if(!is.null(ltnameRule)){
      nameRules = c(nameRules, unlist( lapply(ltnameRule, function(x) strsplit(paste(x, collapse = " "), split = " ")[[1]][1]) ))
    }
    if(!is.null(lnameRule)){
      nameRules = c(nameRules, unlist( lapply(lnameRule, function(x) strsplit(paste(x, collapse = " "), split = " ")[[1]][1]) ))
    }
    if(!is.null(gnameRule)){
      nameRules = c(nameRules, unlist( lapply(gnameRule, function(x) strsplit(paste(x, collapse = " "), split = " ")[[1]][1]) ))
    }
    if(!is.null(gtnameRule)){
      nameRules = c(nameRules, unlist( lapply(gtnameRule, function(x) strsplit(paste(x, collapse = " "), split = " ")[[1]][1]) ))
    }
  }
  
  ## check for T2w in rules
  identifyT2wruleNames = "T2w" %in% unlist(lapply(nameRules, function(f) fnames$type[f == fnames$f]))
  return(identifyT2wruleNames)
}


# getAnywhere(applyLearner)
# A single object matching applyLearner was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value
mynewapplyLearnerxrules <- function(alearner, X, y, minerr, minfrq, classes, gbmModel) 
{
  out = list()
  # alearner=tempruleMetric
  Rules <- alearner[order(alearner[, "err"], decreasing = FALSE),]
  topRules = Rules[Rules[,"err"]<minerr,]
  topRules = topRules[topRules[,"freq"]>=minfrq,]
  out[[1]] = topRules
  
  ## find the rules that the test case comply
  predY <- c()
  RulesIx <- NULL
  seltoprules <- c()
  for (i in 1:nrow(topRules)) {
    ixMatch <- eval(parse(text = paste("which(", topRules[i,"condition"], ")")))
    if(length(ixMatch) > 0) {
      predY <- c(predY, topRules[i, "pred"])
      RulesIx <- rbind(RulesIx, topRules[i,])
      seltoprules <- c(seltoprules,1)
    }else{
      seltoprules <- c(seltoprules,0)
    }
  }
  
  ## identify if T2w rule used
  if(!is.null(RulesIx)){
    ## print top rules
    print(paste0("test complies with # rules: ",nrow(RulesIx)," out of ",nrow(topRules)))
    nr = nrow(RulesIx)
  
    ### organize and present result s
    ## top 5 rule
    resPredY = which.max(table(predY[1:nr]))
    resPred = names(resPredY)  
    #final
    out[[2]] = RulesIx
    out[[3]] = seltoprules
  }
  return(out)
}



myapplyLearner <- function(alearner, X, y, minerr, minfrq, classes, gbmModel) 
{
  # alearner=tempruleMetric
  Rules <- alearner[order(alearner[, "err"], decreasing = FALSE),]
  topRules = Rules[Rules[,"err"]<minerr,]
  topRules = topRules[topRules[,"freq"]>=minfrq,]
  
  ## find the rules that the test case comply
  predY <- c()
  RulesIx <- NULL
  for (i in 1:nrow(topRules)) {
    ixMatch <- eval(parse(text = paste("which(", topRules[i,"condition"], ")")))
    if (length(ixMatch) > 0) {
      predY <- c(predY, topRules[i, "pred"])
      RulesIx <- rbind(RulesIx, topRules[i,])
    }
  }
  
  ## identify if T2w rule used
  if(is.null(RulesIx) || nrow(RulesIx)<5 || !identifyT2wrule(RulesIx, colnames(X), fnames) ){
    predY <- c()
    RulesIx <- NULL
    for (i in 1:nrow(Rules)) {
      ixMatch <- eval(parse(text = paste("which(", Rules[i,"condition"], ")")))
      if (length(ixMatch) > 0) {
        predY <- c(predY, Rules[i, "pred"])
        RulesIx <- rbind(RulesIx, Rules[i,])
      }
    }
  }
    
  ## print top rules
  print(paste0("test complies with # rules: ",nrow(RulesIx)," out of ",nrow(topRules)))
  nr = c(5,nrow(RulesIx))[ which.min(c(5,nrow(RulesIx))) ]
  selRuleswT2Ix=NULL
  if(nr > 1){
    # select top 5
    selRulesIx = presentRules(RulesIx[1:nr,],colnames(X))
    # check if any T2w in top 5
    if(identifyT2wrule(selRulesIx, colnames(X), fnames)){
      #print( selRulesIx[,c(1,2,3,5,4)] )
    }else{
      aflag = FALSE
      aT2wrule = NULL
      while(!aflag && nr<=nrow(RulesIx)-1){
        new_nr = nr+1
        aT2wrule = presentRules(RulesIx[new_nr,], colnames(X))
        aflag = identifyT2wrule(aT2wrule, colnames(X), fnames)
        nr = new_nr
      }
      if(!is.null(aT2wrule)){
        selRuleswT2Ix = rbind(selRulesIx, aT2wrule)
        print( aT2wrule )
      }
      
    }
  }
  
  probRules = sum(unlist(selRulesIx$sumtempProb))/nrow(selRulesIx)
  #print(probRules)
  probModel = predict(object=gbmModel$finalModel, X, type='response', n.trees = gbmModel$finalModel$n.trees)
  #print(probModel)
  
  ### find extracted rules that test case comply that contain "LMSIR_predicted" 
  selrules = c() #"LMSIR_predicted" == 'X[,241]' == ',241' (for grlp)
  for(k in 1:length(RulesIx[,"freq"])){
    if(grepl(',241]',RulesIx[k, "condition"])){
      selrules = rbind(selrules, RulesIx[k,])
    }
  }
  if(!is.null(selrules)){
    wLMSIR = presentRules(selrules,colnames(X))[,"condition"]
    #print(wLMSIR) 
  }
  if(is.null(selrules)){
    wLMSIR = "0"
  }
    
  ### organize and present result s
  ## top 1 rule
  # if(!is.na(as.numeric(predY[1]))){
  #   resPredY = ifelse(as.numeric(predY[1])==1,"C","NC")
  #   resPred = classes[classes==resPredY] #which.max(table(resPredY[1:2]))
  # }else{
  #   resPred = classes[classes==predY[1]]   # which.max(table(predY[1:2]))]
  # }
  ## top 5 rule
  if(!is.na(as.numeric(predY[1]))){
    resPredY = which.max(table(predY[1:nr]))
    resPred = ifelse(names(resPredY)=="1","C","NC")
  }else{
    resPredY = which.max(table(predY[1:nr]))
    resPred = names(resPredY)   # which.max(table(predY[1:2]))]
  }
  
  
  thistestcases <- c()
  thistestcases[c("obsR", "predR", "probRules", "probModel",
                  "id","wLMSIR")] <- c(as.character(y), resPred, probRules, probModel,
                                                   rownames(X), paste(wLMSIR, collapse = ' // '))
  #print(thistestcases)

  out = list()
  out[[1]] = selRulesIx
  out[[2]] = thistestcases
  out[[3]] = selRuleswT2Ix
  
  return(out)
}


# getAnywhere(extractRules)
# A single object matching ???getRuleMetric??? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value
myextractRules <- function(treeList, X, ntree = 1000, maxdepth = 6, random = FALSE) 
{
  levelX = list()
  for (iX in 1:ncol(X)) levelX <- c(levelX, list(levels(X[,iX])))
  ntree = min(treeList$ntree, ntree)
  allRulesList = list()
  for (iTree in 1:ntree) {
    if (random == TRUE) {
      max_length = sample(1:maxdepth, 1, replace = FALSE)
    }
    else {
      max_length = maxdepth
    }
    rule = list()
    count = 0
    rowIx = 1
    tree <- treeList$list[[iTree]]
    ruleSet = vector("list", length(which(tree[, "status"] == -1)))
    myres = mytreeVisit(tree, rowIx = rowIx, count, ruleSet, 
                    rule, levelX, length = 0, max_length = max_length)
    allRulesList = c(allRulesList, myres$ruleSet)
  }
  
  # finish all trees
  allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
  cat(paste(length(allRulesList), " rules (length<=", max_length, 
            ") were extracted from the first ", ntree, " trees.", 
            "\n", sep = ""))
  
  # parse probabilities
  rulesExecnProb <- myruleList2Exec(X, allRulesList)
  
  return(rulesExecnProb)
}

mytreeVisit <- function (tree, rowIx, count, ruleSet, rule, levelX, length, max_length) 
  {
    if (tree[rowIx, "status"] == -1 | length == max_length) {
      count = count + 1
      ruleSet[[count]] = rule
      return(list(ruleSet = ruleSet, count = count))
    }
    ####
    xIx <- tree[rowIx, "split var"]
    xValue <- tree[rowIx, "split point"]
    if (is.null(levelX[[xIx]])) {
      lValue <- paste(paste("X[,", xIx, "]<=", xValue, sep = ""), tree[rowIx, "Prediction"], sep = ";")
      rValue <- paste(paste("X[,", xIx, "]>", xValue, sep = ""), tree[rowIx, "Prediction"], sep = ";")
    }
    else {
      xValue <- which(as.integer(intToBits(as.integer(xValue))) > 0)
      lValue <- levelX[[xIx]][xValue]
      rValue <- setdiff(levelX[[xIx]], lValue)
    }
    
    ####
    xValue <- NULL
    ruleleft <- rule
    if (length(ruleleft) == 0) {
      ruleleft[[as.character(xIx)]] <- lValue
    }
    else {
      if (as.character(xIx) %in% ls(ruleleft)) {
        if (!is.null(levelX[[xIx]])) {
          lValue <- intersect(ruleleft[[as.character(xIx)]], lValue)
          ruleleft[[as.character(xIx)]] <- lValue
        }
        else {
          ruleleft[[as.character(xIx)]] <- paste(ruleleft[[as.character(xIx)]], 
                                                 "&", lValue)
        }
      }
      else {
        ruleleft[[as.character(xIx)]] <- lValue
      }
    }
    
    ####
    ruleright <- rule
    if (length(ruleright) == 0) {
      ruleright[[as.character(xIx)]] <- rValue
    }
    else {
      if (as.character(xIx) %in% ls(ruleright)) {
        if (!is.null(levelX[[xIx]])) {
          rValue <- intersect(ruleright[[as.character(xIx)]], rValue)
          ruleright[[as.character(xIx)]] <- rValue
        }
        else {
          ruleright[[as.character(xIx)]] <- paste(ruleright[[as.character(xIx)]], 
                                                  "&", rValue)
        }
      }
      else {
        ruleright[[as.character(xIx)]] <- rValue
      }
    }
    
    ###### 
    thisList = mytreeVisit(tree, tree[rowIx, "left daughter"], 
                         count, ruleSet, ruleleft, levelX, length + 1, max_length)
    ruleSet = thisList$ruleSet
    count = thisList$count
    
    ###### 
    thisList = mytreeVisit(tree, tree[rowIx, "right daughter"], 
                         count, ruleSet, ruleright, levelX, length + 1, max_length)
    ruleSet = thisList$ruleSet
    count = thisList$count
    
    ###### 
    return(list(ruleSet = ruleSet, count = count))
}

myruleList2Exec <- function(X, allRulesList) 
{
  typeX = getTypeX(X)
  ruleExec <- unique(t(sapply(allRulesList, singleRuleList2Exec, typeX = typeX)))
  ruleExec <- t(ruleExec)
  colnames(ruleExec) <- "condition"
 
  # format to remove probabilies
  rules <- unique(t(sapply(ruleExec, myexecRule)))
  rules <- t(rules)
  colnames(rules) <- "condition"
  rownames(rules) <- 1:nrow(rules)
  
  # format to add probabilies
  rulesProb <- t(sapply(ruleExec[, "condition", drop = FALSE], myprobabilityRule))
  
  return(list(ruleExec = ruleExec, rules=rules, rulesProb = rulesProb))
}

myprobabilityRule<- function(ruleExec)
{
  rulesProb = c()
  splitRules <- unlist(strsplit(ruleExec, split = " & "))
  for(ixr in seq_along(splitRules)){
    rulesProb <- c(rulesProb, unlist(strsplit(splitRules[ixr], split = ";"))[2] )
  }
  return(rulesProb)
}

myexecRule<- function(aruleExec)
{
  rules = c()
  splitRules <- unlist(strsplit(aruleExec, split = " & "))
  for(ixr in seq_along(splitRules)){
    rules = c(rules, unlist(strsplit(splitRules[ixr], split = ";"))[1])
  }
  rules = paste(rules, collapse=' & ')
  return(rules)
}


# getAnywhere(getRuleMetric)
# A single object matching ???getRuleMetric??? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value
mygetRuleMetric <- function(rules, rulesProb, X, target, gbmModel) 
{
  ruleMetric <- t(sapply(rules[, "condition", drop = FALSE], mymeasureRule, X, target))
  rownames(ruleMetric) = NULL
  colnames(ruleMetric) <- c("len", "freq", "err", "condition", "pred")
  
  tempProb <- lapply(sapply(rulesProb, unlist), function(x) 
              plogis(2*as.numeric(x), scale=gbmModel$finalModel$shrinkage))
  sumtempProb <- lapply(tempProb, function(x) sum(x)/length(x))
  tempruleMetric <- as.data.frame(ruleMetric, stringsAsFactors = FALSE)
  tempruleMetric$prob <- tempProb
  tempruleMetric$sumtempProb <- sumtempProb
  
  dIx <- which(tempruleMetric[, "len"] == "-1")
  if (length(dIx) > 0) {
    tempruleMetric <- tempruleMetric[-dIx, ]
    print(paste(length(dIx), " paths are ignored.", sep = ""))
  }
  return(tempruleMetric)
}

# getAnywhere(measureRule)
# A single object matching ??????buildLearner???? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value
mymeasureRule <- function (ruleExec, X, target, pred = NULL, regMethod = "mean") 
{
  len <- length(unlist(strsplit(ruleExec, split = " & ")))
  origRule <- ruleExec
  # find the rulesExec that comply with cases in X
  ruleExec <- paste("which(", ruleExec, ")")
  ixMatch <- eval(parse(text = ruleExec))
  if (length(ixMatch) == 0) {
    v <- c("-1", "-1", "-1", "", "")
    names(v) <- c("len", "freq", "err", "condition", "pred")
    return(v)
  }
  # find the labels of ixMatch cases that comply with ruleExec
  ys <- target[ixMatch]
  # find the freq of ruleExec based on the fraction of (ixMatch cases/ total X cases)
  freq <- round(length(ys)/nrow(X), digits = 3)
  # find the pred of ruleExec based on which.max(table(ys))
  if (is.numeric(target)) {
    if (regMethod == "median") {
      ysMost = median(ys)
    }
    else {
      ysMost <- mean(ys)
    }
    err <- sum((ysMost - ys)^2)/length(ys)
  }
  else {
    if (length(pred) > 0) {
      ysMost = as.character(pred)
    }
    else {
      ysMost <- names(which.max(table(ys)))
    }
    # find the error of ruleExec based on the fraction of (I(y'==y) cases/ total ixMatch cases)
    ly <- sum(as.character(ys) == ysMost)
    conf <- round(ly/length(ys), digits = 3)
    err <- 1 - conf
  }
  rule <- origRule
  v <- c(len, freq, err, rule, ysMost)
  names(v) <- c("len", "freq", "err", "condition", "pred")
  return(v)
}


# getAnywhere(buildLearner)
# A single object matching ??????buildLearner???? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value

mybuildLearner <- function(ruleMetric, X, target, minFreq = 0.01) 
{
  ruleMetric <- ruleMetric[, c("len", "freq", "err", "condition", 
                               "pred"), drop = FALSE]
  learner <- NULL
  listIxInst <- vector("list", nrow(ruleMetric))
  for (i in 1:nrow(ruleMetric)) {
    ixMatch <- eval(parse(text = paste("which(", ruleMetric[i,"condition"], ")")))
    if (length(ixMatch) == 0) 
      next
    listIxInst[[i]] = ixMatch
  }
  ixInstLeft <- 1:length(target)
  while (TRUE) {
    infor = NULL
    restErr <- 1 - max(table(target[ixInstLeft]))/length(target[ixInstLeft])
    for (i in 1:length(listIxInst)) {
      thisInfor <- computeRuleInfor(listIxInst[[i]], ruleMetric[i,"pred"], target)
      infor <- rbind(infor, c(thisInfor, len = as.numeric(ruleMetric[i,"len"])))
    }
    topIx <- order(infor[, "err"], -infor[, "freq"], infor[,"len"], decreasing = FALSE)
    minSupIx <- which(infor[, "freq"] < minFreq)
    if (length(minSupIx) > 0) 
      topIx <- setdiff(topIx, minSupIx)
    if (length(topIx) > 0) 
      topIx <- topIx[1]
    if (length(topIx) == 0) {
      restCondition <- paste("X[,1]==X[,1]")
      restPred <- names(table(target[ixInstLeft]))[which.max(table(target[ixInstLeft]))]
      restSup <- length(ixInstLeft)/length(target)
      thisRuleMetric <- c(len = 1, freq = restSup, err = restErr, 
                          condition = restCondition, pred = restPred)
      learner <- rbind(learner, thisRuleMetric)
      break
    }
    else if (infor[topIx, "err"] >= restErr) {
      restCondition <- paste("X[,1]==X[,1]")
      restPred <- names(table(target[ixInstLeft]))[which.max(table(target[ixInstLeft]))]
      restSup <- length(ixInstLeft)/length(target)
      thisRuleMetric <- c(len = 1, freq = restSup, err = restErr, 
                          condition = restCondition, pred = restPred)
      learner <- rbind(learner, thisRuleMetric)
      break
    }
    thisRuleMetric <- ruleMetric[topIx, , drop = FALSE]
    thisRuleMetric[, c("freq", "err", "len")] <- infor[topIx, 
                                                       c("freq", "err", "len")]
    learner <- rbind(learner, thisRuleMetric)
    ixInstLeft <- setdiff(ixInstLeft, listIxInst[[topIx]])
    listIxInst <- sapply(listIxInst, setdiff, listIxInst[[topIx]])
    if (length(ixInstLeft) == 0) {
      restCondition <- paste("X[,1]==X[,1]")
      restPred <- names(table(target))[which.max(table(target))]
      restSup <- 0
      restErr <- 0
      thisRuleMetric <- c(len = 1, freq = restSup, err = restErr, 
                          condition = restCondition, pred = restPred)
      learner <- rbind(learner, thisRuleMetric)
      break
    }
  }
  rownames(learner) <- NULL
  return(learner)
}

# getAnywhere(computeRuleInfor)
# A single object matching ???computeRuleInfor??? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value

mycomputeRuleInfor <- function(instIx, pred, target) 
{
  trueCls <- as.character(target[instIx])
  err <- 1 - length(which(trueCls == pred))/length(trueCls)
  return(c(err = err, freq = length(instIx)/length(target)))
}

###########################
## GBM rules STEL with probabilities
############################
# getAnywhere(myGBM2List)
# A single object matching ???RF2List??? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value
myGBM2List <- function (object, X) 
{
  treeList <- NULL
  treeList$ntree <- length(object$trees)
  treeList$list <- vector("list", length(object$trees))
  for (i in 1:treeList$ntree) {
    treeList$list[[i]] <-  pretty.gbm.tree(object, i.tree = i)
  }
  v2int <- function(v) {
    sum((-v + 1)/2 * 2^seq(0, (length(v) - 1), 1))
  }
  splitBin = sapply(object$c.splits, v2int)
  
  ## final format
  mytreeList = myformatGBM(treeList, splitBin, X)
  
  return(mytreeList)
}

myformatGBM <- function(gbmList, splitBin, X) 
{
  for (j in 1:length(gbmList$list)) {
    a <- gbmList$list[[j]]
    rownames(a) <- 1:nrow(a)
    a$status <- a$SplitVar
    a <- a[, c("LeftNode", "RightNode", "MissingNode", "SplitVar","SplitCodePred","Prediction","status")]
    a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode", "RightNode", "MissingNode")] <- 
      a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode", "RightNode", "MissingNode")] + 1
    ix <- a$MissingNode[which(a$MissingNode > 0)]
    if (length(ix) > 0) 
      a$status[ix] <- 10
    a <- a[, c("LeftNode","RightNode","SplitVar","SplitCodePred","Prediction","status")]
    
    # deal with categorical or ordered features
    cat <- which(sapply(X, is.factor) & !sapply(X, is.ordered))
    ix <- which(a[, "SplitVar"] %in% cat)
    for (i in ix) a[i, "SplitCodePred"] <- splitBin[a[i, "SplitCodePred"] + 1]
    
    ## finalize results of gbm tree with logit one.sided lower tail probabilities are P[X ??? x] for prediction of rule
    #a$Prediction = plogis(2*a$Prediction, log=FALSE)
    colnames(a) <- c("left daughter", "right daughter", "split var","split point", "Prediction","status")
    gbmList$list[[j]] <- a
  }
  return(gbmList)
}



###########################
## adaboost rules STEL with probabilities
############################
adaboost2List <- function (object, X) 
{
  # object = rpart_adaboost
  treeList <- NULL
  treeList$ntree <- length(object$trees)
  treeList$list <- vector("list", length(object$trees))
  for (i in 1:treeList$ntree) {
    treeList$list[[i]] <- pretty.adaboost.tree(object, i.tree = i)
  }
  mytreeList = myformatadaboost(treeList, X)
  return(mytreeList)
}


myformatadaboost <- function (gbmList, X) 
{
  for (j in 1:length(gbmList$list)) {
    a <- as.data.frame(gbmList$list[[j]])
    rownames(a) <- 1:nrow(a)
    a$status <- a$SplitVar
    a <- a[, c("LeftNode", "RightNode", "MissingNode", "SplitVar", 
               "SplitCodePred", "status")]
    a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode", 
                                     "RightNode", "MissingNode")] <- 
      a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode", "RightNode", "MissingNode")] + 1
    ix <- a$MissingNode[which(a$MissingNode > 0)]
    if (length(ix) > 0) 
      a$status[ix] <- 10
    a <- a[, c("LeftNode", "RightNode", "SplitVar", "SplitCodePred", 
               "status")]
    cat <- which(sapply(X, is.factor) & !sapply(X, is.ordered))
    ix <- which(a[, "SplitVar"] %in% cat)
    for (i in ix) a[i, "SplitCodePred"] <- splitBin[a[i,"SplitCodePred"] + 1]
    colnames(a) <- c("left daughter", "right daughter", "split var", 
                     "split point", "status")
    gbmList$list[[j]] <- a
  }
  return(gbmList)
}


pretty.adaboost.tree <- function (object, i.tree = 1) 
{
  # pretty.gbm.tree returns a data frame. Each row corresponds to a node in the tree. Columns indicate
  # 
  # SplitVar:   index of which variable is used to split. -1 indicates a terminal node.
  # SplitCodePred: if the split variable is continuous then this component is the split point. If the split variable is categorical then this component contains the index of object$c.split that describes the categorical split. If the node is a terminal node then this is the prediction.
  # LeftNode:   the index of the row corresponding to the left node.
  # RightNode:   the index of the row corresponding to the right node.
  # ErrorReduction:  the reduction in the loss function as a result of splitting this node.
  # Weight:   the total weight of observations in the node. If weights are all equal to 1 then this is the number of observations in the node.
  # 
  # object = rpart_adaboost
  library(tree)
  t1<-object$trees[[i.tree]]
  # plot(t1)
  # text(t1,pretty=0)
  var.names = attr(object$terms, "term.labels")[2:length(attr(object$terms, "term.labels"))]
  vars = 1:length( var.names)-1
  #### ff
  ff <- t1$frame
  ylevels <- attr(t1, "ylevels")
  is.leaf <- ff$var == "<leaf>"
  ## 
  splitname <- rownames(t1$splits)
  # struct
  prettyboostT <- matrix("", nrow = nrow(ff), ncol = 8)
  colnames(prettyboostT) <- c("SplitVar", "SplitCodePred", "LeftNode", 
                              "RightNode", "MissingNode", "ErrorReduction", "Weight", 
                              "Prediction")
  row.names(prettyboostT) <- 0:(nrow(prettyboostT)-1)
  ln=1; rn=2;
  for (i in 1:nrow(ff)) {
    if(!is.leaf[i]){
      prettyboostT[i,"SplitVar"] <- vars[var.names==splitname[i]]+1
      prettyboostT[i,"SplitCodePred"] <- format(signif(t1$splits[i, 4L], 3))
      prettyboostT[i,"LeftNode"] <- ln
      prettyboostT[i,"RightNode"] <- rn
      prettyboostT[i,"MissingNode"] <- -1
      ln=rn+1; rn=ln+1;
    } else {
      prettyboostT[i,"SplitVar"] <- -1
      prettyboostT[i,"LeftNode"] <- -1
      prettyboostT[i,"RightNode"] <- -1
      prettyboostT[i,"MissingNode"] <- -1
    }
  }
  prettyboostT[,"ErrorReduction"] = ff$dev
  prettyboostT[,"Weight"] = ff$wt
  prettyboostT[,"Prediction"] = paste0(ff$yval,',',ff$yval2[,4],',',ylevels[ff$yval])
  
  return(prettyboostT)
}




###########################
## Rule extraction and presentation
############################
feature_dictionary<- function(features){
  # create dictionary (34 dynamic features (inside=in) and from (contour=rim))
  dict_dyn = data.frame(f="A_inside", nname="Uptake(in) curve amplitude", stringsAsFactors = FALSE)
  dict_dyn = rbind(dict_dyn, c("alpha_inside","Rate of signal increase(in)" ))   
  dict_dyn = rbind(dict_dyn, c("beta_inside","Rate of signal decrease(in)" ))  
  dict_dyn = rbind(dict_dyn, c("iAUC1_inside","Area Under Uptake curve(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Slope_ini_inside","Initial Uptake slope(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Tpeak_inside","Time-to-peak(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Kpeak_inside","Curvature at Time-to-peak(in)" ))  
  dict_dyn = rbind(dict_dyn, c("SER_inside","SER(in)" ))  
  dict_dyn = rbind(dict_dyn, c("maxCr_inside","max relative signal enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("peakCr_inside","time to Max relative signal enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("UptakeRate_inside","Uptake rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("washoutRate_inside","Washout rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("maxVr_inside","Max spatial variance of enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("peakVr_inside","time to Max spatial variance of enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_increasingRate_inside","enhancement-variance increasing rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_decreasingRate_inside","enhancement-variance decreasing rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_post_1_inside","enhancement-variance at first time-point(in)" ))  
  
  # for contour or rim
  dict_dyn = rbind(dict_dyn, c("A_countor", "Uptake(rim) curve amplitude") )
  dict_dyn = rbind(dict_dyn, c("alpha_countor","Rate of signal increase(rim)" ))   
  dict_dyn = rbind(dict_dyn, c("beta_countor","Rate of signal decrease(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("iAUC1_countor","Area Under Uptake curve(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Slope_ini_countor","Initial Uptake slope(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Tpeak_countor","Time-to-peak(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Kpeak_countor","Curvature at Time-to-peak(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("SER_countor","SER(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("maxCr_countor","max relative signal enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("peakCr_countor","time to Max relative signal enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("UptakeRate_countor","Uptake rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("washoutRate_countor","Washout rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("maxVr_countor","Max spatial variance of enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("peakVr_countor","time to Max spatial variance of enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_increasingRate_countor","enhancement-variance increasing rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_decreasingRate_countor","enhancement-variance decreasing rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_post_1_countor","enhancement-variance at first time-point(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("min_F_r_i","min uptake" ))  
  dict_dyn = rbind(dict_dyn, c("max_F_r_i","max uptake" ))  
  dict_dyn = rbind(dict_dyn, c("mean_F_r_i","uptake average" ))  
  dict_dyn = rbind(dict_dyn, c("var_F_r_i","uptake variance" ))  
  dict_dyn = rbind(dict_dyn, c("skew_F_r_i","uptake skewness" ))  
  dict_dyn = rbind(dict_dyn, c("kurt_F_r_i","uptake kurtosis" )) 
  dict_dyn$type = "T1wdynamic"
  
  # for morphology
  dict_mor = data.frame(f="ivVariance", nname="Variance of spatial Margin Gradient",  stringsAsFactors = FALSE)
  dict_mor = rbind(dict_mor, c("iMax_Variance_uptake", "max Variance of spatial Margin Gradient"))
  dict_mor = rbind(dict_mor, c("iiMin_change_Variance_uptake","change in Variance of spatial Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("iiiMax_Margin_Gradient","Max Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("k_Max_Margin_Grad","time to Max Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("circularity","circularity" )) 
  dict_mor = rbind(dict_mor, c("irregularity","irregularity" )) 
  dict_mor = rbind(dict_mor, c("irregularity","irregularity" )) 
  dict_mor = rbind(dict_mor, c("edge_sharp_mean","mean 3D Sharpness of lesion margin " )) 
  dict_mor = rbind(dict_mor, c("edge_sharp_std","std 3D Sharpness of lesion margin " )) 
  dict_mor = rbind(dict_mor, c("max_RGH_mean","max Radial gradient" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_mean_k","Time to max Radial gradient" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_var","Max Radial gradient variance" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_var_k","Time to max variance Radial gradient" )) 
  dict_mor$type = "T1wmorphology"
  
  # for T1w texture
  dict_tex = data.frame(f="texture_energy_nondir_post1", nname="Energy ptime1",  stringsAsFactors = FALSE)
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post1", "Contrast ptime1"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post1","Correlation ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post1","Variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post1","Inverse difference moment ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post1","Sum average ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post1","Sum variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post1","Sum Entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post1","Entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post1","Difference variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post1","Difference entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post2","Energy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post2", "Contrast ptime2"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post2","Correlation ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post2","Variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post2","Inverse difference moment ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post2","Sum average ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post2","Sum variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post2","Sum Entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post2","Entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post2","Difference variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post2","Difference entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post3","Energy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post3", "Contrast ptime3"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post3","Correlation ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post3","Variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post3","Inverse difference moment ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post3","Sum average ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post3","Sum variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post3","Sum Entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post3","Entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post3","Difference variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post3","Difference entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post4","Energy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post4", "Contrast ptime4"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post4","Correlation ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post4","Variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post4","Inverse difference moment ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post4","Sum average ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post4","Sum variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post4","Sum Entropy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post4","Entropy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post4","Difference variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post4","Difference entropy ptime4" )) 
  dict_tex$type = "T1wtexture"
  
  # for dispersion
  dict_dis = data.frame(f="V0", nname="dispersion s0",  stringsAsFactors = FALSE)
  dict_dis = rbind(dict_dis, c("V1","dispersion s1" )) 
  dict_dis = rbind(dict_dis, c("V2","dispersion s2" )) 
  dict_dis = rbind(dict_dis, c("V3","dispersion s3" )) 
  dict_dis = rbind(dict_dis, c("V4","dispersion s4" )) 
  dict_dis = rbind(dict_dis, c("V5","dispersion s5" )) 
  dict_dis = rbind(dict_dis, c("V6","dispersion s6" )) 
  dict_dis = rbind(dict_dis, c("V7","dispersion s7" )) 
  dict_dis = rbind(dict_dis, c("V8","dispersion s8" )) 
  dict_dis = rbind(dict_dis, c("V9","dispersion s9" )) 
  dict_dis = rbind(dict_dis, c("V10","dispersion s10" )) 
  dict_dis = rbind(dict_dis, c("V11","dispersion s11" )) 
  dict_dis = rbind(dict_dis, c("V12","dispersion s12" )) 
  dict_dis = rbind(dict_dis, c("V13","dispersion s13" )) 
  dict_dis = rbind(dict_dis, c("V14","dispersion s14" )) 
  dict_dis = rbind(dict_dis, c("V15","dispersion s15" )) 
  dict_dis = rbind(dict_dis, c("V16","dispersion s16" )) 
  dict_dis = rbind(dict_dis, c("V17","dispersion s17" )) 
  dict_dis = rbind(dict_dis, c("V18","dispersion s18" )) 
  dict_dis = rbind(dict_dis, c("V19","dispersion s19" )) 
  dict_dis$type = "dispersion"
  
  # for single enhancement time point samples
  dict_timenh = data.frame(f="earlySE0", nname="1st post-SE s0",  stringsAsFactors = FALSE)
  dict_timenh = rbind(dict_timenh, c("earlySE1","1st post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE2","1st post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE3","1st post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE4","1st post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE5","1st post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE6","1st post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE7","1st post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE8","1st post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE9","1st post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE10","1st post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE11","1st post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE12","1st post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE13","1st post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE14","1st post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE15","1st post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE16","1st post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE17","1st post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE18","1st post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE19","1st post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE0","2nd post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE1","2nd post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE2","2nd post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE3","2nd post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE4","2nd post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE5","2nd post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE6","2nd post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE7","2nd post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE8","2nd post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE9","2nd post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE10","2nd post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE11","2nd post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE12","2nd post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE13","2nd post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE14","2nd post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE15","2nd post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE16","2nd post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE17","2nd post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE18","2nd post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE19","2nd post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE0","3rd post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE1","3rd post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE2","3rd post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE3","3rd post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE4","3rd post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE5","3rd post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE6","3rd post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE7","3rd post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE8","3rd post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE9","3rd post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE10","3rd post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE11","3rd post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE12","3rd post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE13","3rd post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE14","3rd post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE15","3rd post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE16","3rd post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE17","3rd post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE18","3rd post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE19","3rd post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE0","last post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE1","last post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE2","last post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE3","last post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE4","last post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE5","last post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE6","last post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE7","last post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE8","last post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE9","last post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE10","last post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE11","last post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE12","last post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE13","last post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE14","last post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE15","last post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE16","last post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE17","last post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE18","last post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE19","last post-SE s19" )) 
  dict_timenh$type = "single-time-Enh"
  
  # create dictionaries ( T2 )
  dict_T2w = data.frame(f="find_t2_signal_int", nname="BIRADS T2w SI cat",  stringsAsFactors = FALSE)
  dict_T2w = rbind(dict_T2w, c("T2_lesionSI","T2w mean SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2_lesionSIstd","T2w std SI" )) 
  dict_T2w = rbind(dict_T2w, c("LMSIR","measured LMSIR" )) 
  dict_T2w = rbind(dict_T2w, c("T2min_F_r_i","min T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2max_F_r_i","max T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2mean_F_r_i","mean T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2var_F_r_i","variance T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2skew_F_r_i","skewness T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2kurt_F_r_i","kurtosis T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2grad_margin","T2w margin gradient" ))
  dict_T2w = rbind(dict_T2w, c("T2grad_margin_var","T2w margin gradient variance" )) 
  dict_T2w = rbind(dict_T2w, c("T2RGH_mean","T2w average radial gradient" )) 
  dict_T2w = rbind(dict_T2w, c("T2RGH_var","T2w radial gradient variance" )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_energy_nondir","T2w Energy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_contrast_nondir", "T2w Contrast "))
  dict_T2w = rbind(dict_T2w, c("T2texture_correlation_nondir","T2w Correlation " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_variance_nondir","T2w Variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_inversediffmoment_nondir","T2w Inverse difference moment " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumaverage_nondir","T2w Sum average " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumvariance_nondir","T2w Sum variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumentropy_nondir","T2w Sum Entropy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_entropy_nondir","T2w Entropy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_diffvariance_nondir","T2w Difference variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_diffentropy_nondir","T2w Difference entropy " )) 
  ## add predictive features
  dict_T2w = rbind(dict_T2w, c("LMSIR_predicted","predicted LMSIR" )) 
  dict_T2w = rbind(dict_T2w, c("T2wSI_predicted","predicted BIRADS T2w SI cat" )) 
  ## add predictive features
  dict_T2w = rbind(dict_T2w, c("ave_T20","average T2w SI")) 
  dict_T2w = rbind(dict_T2w, c("ave_T21","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T22","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T23","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T24","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T25","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T26","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T27","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T28","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T29","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T210","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T211","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T212","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T213","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T214","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T215","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T216","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T217","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T218","average T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("ave_T219","average T2w SI" )) 
  dict_T2w$type = "T2w"
  
  alldicts = list(dict_dyn, dict_mor, dict_tex, dict_dis, dict_timenh, dict_T2w)
  
  library(RColorBrewer) # hist(discoveries, col = colors)
  n=length(alldicts) 
  colors <- brewer.pal(n, "Dark2")
  
  # based on number of dictionaries asign colors
  for(k in 1:n){
    alldicts[[k]]$color = colors[k]
  }
  
  # translated names of features
  namesf = colnames(features)
  fnnames = data.frame()
  for(k in 2:ncol(features)){
    for(d in 1:length(alldicts)){
      if(namesf[k] %in% alldicts[[d]]$f){
        namentype = alldicts[[d]][namesf[k] == alldicts[[d]]$f,]
        fnnames = rbind(fnnames, namentype )
        #print(namentype$nname)
        break;
      }
    }
  }
  output = list(alldicts=alldicts, fnnames=fnnames)
  
  return(output)
}

# getAnywhere(presentRules)
# A single object matching ???computeRuleInfor??? was found
# It was found in the following places
# package:inTrees
# namespace:inTrees
# with value

mypresentRules <- function (rules, colN, fnames) 
{
  # colN = colnames(X)
  # rules = selRulesIx[,c(1:5,7)]
  namedrules = presentRules(rules, colN)
  
  # per rule
  for (i in 1:nrow(namedrules[, "condition", drop = FALSE])) {
    origRule = unlist(strsplit(namedrules[i, "condition"], split = " & "))
    
    ## extract rule values
    newRules = c()
    types = c()
    for(C in origRule){
      if(!is.na(str_match(C,">"))){
        if(!is.na(str_match(C,">="))){
          resC = str_split(C, ">=")
          dirC = ">="
        }
        else{
          resC = str_split(C, ">")
          dirC = ">"
        }
      }
      
      if(!is.na(str_match(C,"<"))){
        if(!is.na(str_match(C,"<="))) {
          resC = str_split(C, "<=")
          dirC = "<="
        }
        else{
          resC = str_split(C, "<")
          dirC = "<"
        }
      }
      
      # get rule value
      #print(C)
      Cval = as.numeric(resC[[1]][2])
      Cname = resC[[1]][1]
      Cnewname = fnames$nname[Cname == fnames$f][1]
      Ctype = fnames$type[Cname == fnames$f][1]
      
      # find dist
      fqrts = summary(na.omit(imgT2pLMSIRtest[,Cname]))[-4] # exclude mean from quartiles
      #cat(Cname,":==========\n")
      #print(fqrts)
    
      ####### flag quartile
      if (is.infinite(fqrts[1])){  # is lower than
        fqrts[1] = Cval
      }
      if (is.infinite(fqrts[5])){ 
        fqrts[5] = Cval
      }
      
      if(!is.infinite(fqrts[1]) && !is.infinite(fqrts[5])){
        if(dirC=='>='){
          flag = Cval >= fqrts # find ranges for threshold
        }
        if(dirC=='>'){
          flag = Cval > fqrts # find ranges for threshold
        }
        if(dirC=='<='){
          flag = Cval <= fqrts # find ranges for threshold
        }
        if(dirC=='<'){
          flag = Cval < fqrts
        }
        flagqrt = !flag
      }
      
      ####### flag quartile
      #print(fqrts[flagqrt])
      nqrts = length(names(fqrts[flagqrt]))
      if(nqrts==0){
        textvar = "low"
      }
      else if(nqrts==1){
        if( all(names(fqrts[flagqrt]) == c("Min."))){textvar = "very low "}
        if( all(names(fqrts[flagqrt]) == c("1st Qu."))){textvar = "low"}
        if( all(names(fqrts[flagqrt]) == c("Median"))){textvar = "Median"}
        if( all(names(fqrts[flagqrt]) == c("3rd Qu."))){textvar = "high"}
        if( all(names(fqrts[flagqrt]) == c("Max."))){textvar = "very high "}
      }
      else if(nqrts==2){
        if( all(names(fqrts[flagqrt]) == c("Min.", "1st Qu."))){textvar = "low "}
        if( all(names(fqrts[flagqrt]) == c("1st Qu.", "Median"))){textvar = "slightly low to medium"}
        if( all(names(fqrts[flagqrt]) == c("Median","3rd Qu."))){textvar = "medium to slightly high"}
        if( all(names(fqrts[flagqrt]) == c("3rd Qu.", "Max."))){textvar = "high "}
      }
      else if(nqrts==3){
        if( all(names(fqrts[flagqrt]) == c("Min.", "1st Qu.", "Median"))){textvar = "low to medium "}
        if( all(names(fqrts[flagqrt]) == c("1st Qu.", "Median", "3rd Qu."))){textvar = "medium"}
        if( all(names(fqrts[flagqrt]) == c("Median", "3rd Qu.","Max."))){textvar = "medium to high "}
      }
      else if(nqrts==4){
        if( all(names(fqrts[flagqrt]) == c("Min.", "1st Qu.", "Median", "3rd Qu."))){textvar = "very low to high "}
        if( all(names(fqrts[flagqrt]) == c("1st Qu.", "Median", "3rd Qu.","Max."))){textvar = "low to very high "}
      }
      else if(nqrts==5){
        textvar = "high"
      }
      else{
        textvar = ""
      }
      newRules <- c(newRules, paste(Cnewname, textvar, sep = " "))
      types <- c(types, Ctype)
    }
    
    # append per conndiction
    newconditions <- paste(newRules, sep="", collapse=" & ")
    newtypes <- paste(unique(types), sep="", collapse=" & ")
    #print(newconditions)
    
    # replace new formated condition in rule
    namedrules[i, "condition"] <- newconditions
    namedrules[i, "Ruletype"]  <- newtypes
  }
  
  return(namedrules)
}

  
#################
# library(gbm)
# library(randomForest)
# library(inTrees)
# library(reprtree)
# data(iris)
# iris = subset(iris, Species!="setosa") # for binary
# X <- iris[,1:(ncol(iris)-1)]
# target <- iris[,"Species"]
# 
# ## RF
# # rf <- randomForest(Species ~ ., data=iris, importance=TRUE,  proximity=TRUE)
# # tree <- getTree(rf, k=1, labelVar=TRUE)
# # realtree <- reprtree:::as.tree(a, rf)
# # gTree = tree
# # rforest = rf
# # reprtree:::plot.getTree(rforest, k=3, depth=4)
# 
# ## GBM
# iris$Species = ifelse(iris$Species=="versicolor",1,0)
# gbmT <- gbm(Species ~ ., data=iris, dist="adaboost", n.trees = 5, interaction.depth = 2,
#                     shrinkage = 1)
# treeList <- GBM2List(gbmT,X)
# ruleExec = extractRules(treeList,X)
# 
# # ## The adaboost method gives the predictions on logit scale. You can convert it to the 0-1 output:
# # gbm_predicted<-plogis(2*gbm_predicted)
# # ## You can also directly obtain the probabilities from the predict.gbm function;
# # predict(gbm_algorithm, test_dataset, n.trees = 5000, type = 'response')
# library(adabag)
# library(rpart)
# iris$Species = as.factor(ifelse(iris$Species==1,"versicolor","virginica"))
# rpart_adaboost <- boosting(Species ~ ., data=iris, mfinal=5, coeflearn="Freund",
#                            boos = FALSE,
#                            control = rpart.control(maxdepth = 2))
# 
# mytreeList  <- adaboost2List(rpart_adaboost,X) 
# myruleExec = extractRules(mytreeList,X)



#################
# k=1
# load(paste0("Z:/Cristina/Section2/finalTreebased-boosting/Outputs/T2SIvspredLMSIR_boost_addeddiagvalue_cv",k,".RData"))
# 
# ##################
# # Define datasets
# ##################
# imgT1train = T1train[,-c(ncol(T1train))]
# imgT2mLMSIR = T1T2train[,-c(ncol(T1T2train))] #"T1w+T2wSI_measureLMSIR"
# imgT2pLMSIR = cbind(imgT2mLMSIR[,-c(201)], trainpredT2LMSIR[1]) #"T1w+T2wSI_predictedLMSIR"
# 
# # test set
# imgT1test = T1test[,-c(ncol(T1test))]
# T1T2test = allfT1T2[[2]][-c(201,ncol(allfT1T2[[2]]))]
# imgT2mLMSIRtest  = cbind(T1T2test,  testpredT2LMSIR[2])
# imgT2pLMSIRtest  = cbind(T1T2test,  testpredT2LMSIR[1])
# classes = levels(imgT2mLMSIRtest[,"lesion_label"])
# 
# ############
# # Build a train trees
# imgT2pLMSIR$lesion_label = as.factor(ifelse(imgT2pLMSIR$lesion_label=="C",1,0))
# gbmT2 <- gbm(lesion_label ~ .,
#              data=imgT2pLMSIR,
#              dist="adaboost", n.tree = 1000, interaction.depth = 6, shrinkage = 0.01, cv.folds = 10)
# 
# # check performance using 5-fold cross-validation
# best.iter <- gbm.perf(gbmT2,method="cv")
# print(best.iter)
# 
# 
# objControl <- trainControl(method='cv', number=10, 
#                            returnResamp='final', 
#                            savePredictions='final', 
#                            summaryFunction = twoClassSummary, 
#                            classProbs = TRUE)
# 
# gbmGrid = expand.grid(.n.trees = c(100,250,500,1000,1500,2000),
#                       .interaction.depth = c(1,3,6),
#                       .shrinkage = c(0.1,0.01,0.001),
#                       .n.minobsinnode = c(5))
# 
# imgT2pLMSIR = na.omit(imgT2pLMSIR)
# imgT2pLMSIR$lesion_label = as.factor(ifelse(imgT2pLMSIR$lesion_label=="C",1,0))
# 
# gbmModel <- train(lesion_label ~ .,
#                   data=imgT2pLMSIR, 
#                   method='gbm', 
#                   distribution='adaboost',
#                   trControl=objControl,
#                   tuneGrid=gbmGrid,
#                   verbose=TRUE,
#                   metric = "ROC")
# 
# #predictions <- predict(object=gbmModel, X, type='prob')
# # summarize the model  # plot the effect of parameters on accuracy
# print(gbmModel)
# plot(gbmModel)
# print(gbmModel$results)
# 
# # transform rf object to an inTrees' format
# X <- imgT2pLMSIR[,2:(ncol(imgT2pLMSIR))]
# target <- imgT2pLMSIR[,"lesion_label"]
# # treeList <- adaboost2List(treedata_T2wpLMSIR$forest,X)
# # treeList <- GBM2List(gbmT2,X) #gbmT2,X)
# # treeList <- myGBM2List(gbmT2,X) #gbmT2,X)
# treeList <- GBM2List(gbmModel$finalModel,X)
# 
# # extract rules
# ruleExec = extractRules(treeList, X, gbmModel$finalModel$n.trees)
# myruleExec = myextractRules(treeList, X, ntree = 1000) # treedata_T2wpLMSIR$mfinal)
# ruleExec = myruleExec$ruleExec
# rules = myruleExec$rules
# rulesProb = myruleExec$rulesProb
# 
# # how good are these rules or what classification we obtain by applying them. 
# # This can be done with the handy function getRuleMetric. For the first 10 extracted rules:
# ruleMetric <- getRuleMetric(ruleExec,X,target)
# ruleMetricc2 = pruneRule(ruleMetric, X, target)
# ruleMetric <- mygetRuleMetric(rules,rulesProb,X,target)  # get rule metrics
# ruleMetric[1:10,]
# 
# # a matrix including the conditions, prediction, and metrics, ordered by priority.
# alearner <- ruleMetricc2
# restest = c()
# for(i in 1:nrow(imgT2pLMSIRtest))
# {
#   X = imgT2pLMSIRtest[i,2:ncol(imgT2pLMSIRtest)]
#   y = imgT2pLMSIRtest[i,"lesion_label"]
#   rulesoutput = myapplyLearner(alearner, X, y, minerr=0.10, classes)
#   restest = rbind(restest, rulesoutput[[2]])
# }
# #perf_all = cbind(perf_T2wpLMSIR, restest)
# confusionMatrix(perf_T2wpLMSIR$pred, perf_T2wpLMSIR$obs)
# confusionMatrix(restest[,"predR"], restest[,"obsR"])
# 
# # v,s perf_T2wpLMSIR[1:10,]
# # v.s perf_imgT1[1:10,]
# # library(xtable)
# # print(xtable(RulesIx), include.rownames=FALSE)
# allrulesres = rbind(allrulesres, restest) 
# }

