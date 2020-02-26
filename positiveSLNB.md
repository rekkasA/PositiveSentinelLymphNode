---
title: "Positive SLNB"
output: 
  bookdown::html_document2:
    keep_md: true
    code_folding: hide
---





# Model development

The primary goal is to create a model that predicts recurrence using the available variables except the variable with additional positive lymph nodes.
The set of variables presented below was considered of importance for the following analyses: 


- `age`
- `sex`
- `simpleloc` (simple tumor location)
- `clark` level
- `histology_simple`
- `breslow` thickness
- `AJCC_sub_SLNB_8th`
- `ulceration`
- `SN_tumor_burden`
- `no_removed_SNs`
- `no_removed_nonSNs` (not used for the first set of analyses)
- `no_pos_SNs`
- `no_pos_SNs_cat`



```r
dd <- datadist(dat)
options(datadist = 'dd')
options(digits = 8)
miVariables <- c("center", "age", "sex", "simpleloc", "SN_tumor_burden_extended",
                 "histology_simple","breslow", "AJCC_sub_SLNB_8th", 
                 "no_pos_SNs", "no_removed_SNs", "EJC_groups",
                 "ulceration", "SLNBdate",
                 "SN_tumor_burden", "no_removed_nonSNs",
                 "status_RFS_FDA", "timeToRecurrence", 
                 "status_DMFS_FDA", "timeToDistant",
                 "status_OS_FDA", "timeToFollowup")
miData <- dat[miVariables] 
miData$timeToRecurrence <- log(miData$timeToRecurrence)
miData$timeToDistant <- log(miData$timeToDistant)
miData$timeToFollowup <- log(miData$timeToFollowup)
# multiple imputations
# miDataResult <- mice(data = miData, m = 20)
# saveRDS(miDataResult, file = "miDataResult.rds")
miDataResult <- readRDS("miDataResult.rds")
```

## Models for recurrence


```r
S <- Surv(dat$timeToRecurrence,
          dat$status_RFS_FDA)
form <- S ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
# tableCount <- tableCount + 1
kable(round(summary(f.mi)[, c(4, 6, 7)], 2), caption = "Multivariable Cox analysis of recurrence") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-3)Multivariable Cox analysis of recurrence</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Effect </th>
   <th style="text-align:right;"> Lower 0.95 </th>
   <th style="text-align:right;"> Upper 0.95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> no_removed_SNs </td>
   <td style="text-align:right;"> -0.13 </td>
   <td style="text-align:right;"> -0.30 </td>
   <td style="text-align:right;"> 0.03 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 0.87 </td>
   <td style="text-align:right;"> 0.74 </td>
   <td style="text-align:right;"> 1.03 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> no_pos_SNs </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> -0.37 </td>
   <td style="text-align:right;"> 0.99 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.37 </td>
   <td style="text-align:right;"> 0.69 </td>
   <td style="text-align:right;"> 2.69 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> breslow </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.40 </td>
   <td style="text-align:right;"> 1.16 </td>
   <td style="text-align:right;"> 1.68 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.69 </td>
   <td style="text-align:right;"> 1.39 </td>
   <td style="text-align:right;"> 2.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> age </td>
   <td style="text-align:right;"> 0.22 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 0.35 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.25 </td>
   <td style="text-align:right;"> 1.10 </td>
   <td style="text-align:right;"> 1.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sex - female:male </td>
   <td style="text-align:right;"> -0.22 </td>
   <td style="text-align:right;"> -0.39 </td>
   <td style="text-align:right;"> -0.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;"> 0.67 </td>
   <td style="text-align:right;"> 0.95 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ulceration - present:absent </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.19 </td>
   <td style="text-align:right;"> 0.57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.46 </td>
   <td style="text-align:right;"> 1.20 </td>
   <td style="text-align:right;"> 1.76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simpleloc - trunk:extremity </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:right;"> -0.14 </td>
   <td style="text-align:right;"> 0.23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.05 </td>
   <td style="text-align:right;"> 0.87 </td>
   <td style="text-align:right;"> 1.26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simpleloc - head&amp;neck:extremity </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> -0.35 </td>
   <td style="text-align:right;"> 0.52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.09 </td>
   <td style="text-align:right;"> 0.71 </td>
   <td style="text-align:right;"> 1.68 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> histology_simple - NM:SSM </td>
   <td style="text-align:right;"> -0.06 </td>
   <td style="text-align:right;"> -0.26 </td>
   <td style="text-align:right;"> 0.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;"> 0.77 </td>
   <td style="text-align:right;"> 1.15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> histology_simple - other:SSM </td>
   <td style="text-align:right;"> 0.28 </td>
   <td style="text-align:right;"> -0.09 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hazard Ratio </td>
   <td style="text-align:right;"> 1.32 </td>
   <td style="text-align:right;"> 0.92 </td>
   <td style="text-align:right;"> 1.90 </td>
  </tr>
</tbody>
</table>

We performed a backwards selection of to come up with the final model that included `ulceration`, `age`, `breslow` and `SN_tumor_burden` (Table \@ref(tab:finalRec)). Logarithmic transformations of the continuous covariates -i.e. `age`, `breslow`  and `SN_tumor_burden`- adequately represented their effects.


```r
form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
summaryFit <- summary(f.mi)
namesSummaryFit <- rownames(summaryFit)[seq(1, 7, 2)]
summaryFit <- summaryFit[seq(2, 8, 2), c(4, 6, 7)]
rownames(summaryFit) <- namesSummaryFit
colnames(summaryFit)[1] <- "Hazard ratio"
kable(
  round(summaryFit, 2), 
  # "latex",
  # booktabs = TRUE,
  caption =  "Final model for recurrence") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:finalRec)Final model for recurrence</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Hazard ratio </th>
   <th style="text-align:right;"> Lower 0.95 </th>
   <th style="text-align:right;"> Upper 0.95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> age </td>
   <td style="text-align:right;"> 1.27 </td>
   <td style="text-align:right;"> 1.13 </td>
   <td style="text-align:right;"> 1.43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden </td>
   <td style="text-align:right;"> 1.50 </td>
   <td style="text-align:right;"> 1.33 </td>
   <td style="text-align:right;"> 1.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> breslow </td>
   <td style="text-align:right;"> 1.36 </td>
   <td style="text-align:right;"> 1.20 </td>
   <td style="text-align:right;"> 1.55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ulceration - present:absent </td>
   <td style="text-align:right;"> 1.40 </td>
   <td style="text-align:right;"> 1.16 </td>
   <td style="text-align:right;"> 1.69 </td>
  </tr>
</tbody>
</table>

We used the variables in the final model to predict 5-year recurrence (Table \@ref(tab:finalRec5)).


```r
S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)
form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
summaryFit <- summary(f.mi)
namesSummaryFit <- rownames(summaryFit)[seq(1, 7, 2)]
summaryFit <- summaryFit[seq(2, 8, 2), c(4, 6, 7)]
rownames(summaryFit) <- namesSummaryFit
colnames(summaryFit)[1] <- "Hazard ratio"
kable(round(summaryFit, 2),
      # "latex",
      # booktabs = TRUE,
      caption = "Final model for 5-year recurrence") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:finalRec5)Final model for 5-year recurrence</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Hazard ratio </th>
   <th style="text-align:right;"> Lower 0.95 </th>
   <th style="text-align:right;"> Upper 0.95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> age </td>
   <td style="text-align:right;"> 1.28 </td>
   <td style="text-align:right;"> 1.12 </td>
   <td style="text-align:right;"> 1.45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden </td>
   <td style="text-align:right;"> 1.59 </td>
   <td style="text-align:right;"> 1.39 </td>
   <td style="text-align:right;"> 1.81 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> breslow </td>
   <td style="text-align:right;"> 1.40 </td>
   <td style="text-align:right;"> 1.23 </td>
   <td style="text-align:right;"> 1.61 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ulceration - present:absent </td>
   <td style="text-align:right;"> 1.43 </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> 1.74 </td>
  </tr>
</tbody>
</table>

The final 5-year recurrence model had c-index of 0.68 (95 percent c.i. 0.65 to 0.7). 

We assess calibration of the final model using a leave-one-center-out cross validation approach. The prediction model is built on 8 centers and calibration is evaluated on the 9th. That is performed recursively, each time leaving a different center out for validation. We separate multiple imputation in the training set and the test set to avoid using information of missingness in the training centers to the test center.

In general, we see quite adequate performance across centers, where confidence intervals include the diagonal. However, in smaller centers such as the one in Groningen there is substantial underestimation of risk.


```r
centers <- unique(dat$center)
discrimination <- rep(0, length(centers))
par(mfrow=c(3, 3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence,
                      datTestCenter$status_RFS_FDA)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  # miDataResultTrain <- mice(data = miDataTrain, m = 20)
  # saveRDS(miDataResultTrain, file.path(".", "mi",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataResultTrain <- readRDS(file.path(".", "mi",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  if(any(table(datTestCenter$AJCC_sub_SLNB_8th) == 0))
    miDataTest <- subset(miDataTest, select = -AJCC_sub_SLNB_8th)
  
  # miDataResultTest <- mice(data = miDataTest, m = 20)
  # saveRDS(miDataResultTest, file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  miDataResultTest <- readRDS(file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 20,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(f.miTrainCenters, TRUE)
  timePoint <- 60
  cIndex <- rep(0, 20)
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex[j] <- rc[1]
    
    linearPredictor <- cbind(linearPredictor, lpTestCenter)
    
  }
  
  discrimination[i] <- mean(cIndex)
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = baselineHazard,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToRecurrence60,
                                status = datTestCenter$status_RFS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalRec-1.png" alt="\label{fig:leave1outCalRec}Leave one center out cross validation for the prediction of 5-year recurrence"  />
<p class="caption">(\#fig:leave1outCalRec)\label{fig:leave1outCalRec}Leave one center out cross validation for the prediction of 5-year recurrence</p>
</div>

## Distant metastasis


```r
SDistant <- Surv(dat$timeToDistant60,
                 dat$status_DMFS_FDA60)
form <- SDistant ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)

f.miDistantFull <- fit.mult.impute(form,
                                   cph,
                                   xtrans = miDataResult,
                                   data = miData,
                                   n.impute = 20,
                                   pr = FALSE,
                                   fit.reps = TRUE,
                                   y = TRUE,
                                   x = TRUE,
                                   se.fit = TRUE)
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResult,
                           data = miData,
                           n.impute = 20,
                           pr = FALSE,
                           fit.reps = TRUE,
                           y = TRUE,
                           x = TRUE,
                           se.fit = TRUE)
linearPredictor <- f.miRec$linear
dd <- datadist(linearPredictor)
options(datadist = 'dd')
options(digits = 8)
fDistantCalibrated <- cph(SDistant ~ linearPredictor, 
                          y = TRUE,
                          x = TRUE)
cIndexCalibrated <- rcorr.cens(-linearPredictor*coef(fDistantCalibrated), SDistant)
```

For the assessment of distant metastasis we considered a calibrated version of the 5-year recurrence model. The association between distant metastasis and was of the same size (calibration slope of 1.01, 95 percent c.i. 0.87 to 1.16). We compare the performance of the considered model to that of multivariable Cox regression model including all 9 covariates of interest (`sex`, `ulceration`, `no_removed_SNs`, `no_pos_SNs`, `simpleloc`, `histology_simple`, `breslow`, `SN_tumor_burden` and `age`). The full model had a c-index of 0.7 (95 percent c.i. 0.68 to 0.73) while the calibrated model had a c-index of 0.7 (95 percent c.i. 0.67 to 0.72).

In terms of calibation, we first assess a leave-one-center-out cross validation, where the baseline hazard is estimated from the training set of 8 centers based on the cox regression model for 5-year recurrence and the calibration slope for the linear predictor is derived from a cox model predicting risk of distant metastasis form the linear predictor of the previous model (Figure \@ref(fig:leave1outCalDistant)). 


```r
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
centers <- unique(dat$center)
discrimination <- rep(0, length(centers))
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence60,
                      datTestCenter$status_RFS_FDA60)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  miDataResultTrain <- readRDS(file.path(".", "mi",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  miDataResultTest <- readRDS(file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 20,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  SDistantTrainCenter <- Surv(datTrainCenters$timeToDistant60,
                              datTrainCenters$status_DMFS_FDA60)
  lp <- f.miTrainCenters$linear
  fDistant <- cph(SDistantTrainCenter ~ lp, 
                  y = TRUE,
                  x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(fDistant, TRUE) 
  timePoint <- 60
  cIndex <- rep(0, 20)
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    lpTestCenter <- lpTestCenter*coef(fDistant)
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex[j] <- rc[1]
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  discrimination[i] <- mean(cIndex)
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  # baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = f.basehazTrainCenters,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToDistant60,
                                status = datTestCenter$status_DMFS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalDistant-1.png" alt="\label{fig:leave1outCal}Leave one center out cross validation of the calibrated model for 5-year distant metastasis"  />
<p class="caption">(\#fig:leave1outCalDistant)\label{fig:leave1outCal}Leave one center out cross validation of the calibrated model for 5-year distant metastasis</p>
</div>

## Overall mortality



```r
SDeath <- Surv(dat$timeToFollowup60,
               dat$status_OS_FDA60)
form <- SDeath ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)

f.miDeathFull <- fit.mult.impute(form,
                                 cph,
                                 xtrans = miDataResult,
                                 data = miData,
                                 n.impute = 20,
                                 pr = FALSE,
                                 fit.reps = TRUE,
                                 y = TRUE,
                                 x = TRUE,
                                 se.fit = TRUE)
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResult,
                           data = miData,
                           n.impute = 20,
                           pr = FALSE,
                           fit.reps = TRUE,
                           y = TRUE,
                           x = TRUE,
                           se.fit = TRUE)
linearPredictor <- f.miRec$linear
dd <- datadist(linearPredictor)
options(datadist = 'dd')
options(digits = 8)
fDeathCalibrated <- cph(SDeath ~ linearPredictor, 
                        y = TRUE,
                        x = TRUE)
cIndexCalibrated <- rcorr.cens(-linearPredictor*coef(fDeathCalibrated), SDeath)
```

For the assessment 5-year overall mortality we considered a calibrated version of the 5-year recurrence model. We compare the performance of the considered model to that of multivariable Cox regression model including all 9 covariates of interest. The association between recurrence and overall mortality was not different (calibration slope 1.04, 95 percent c.i. 0.88 to 1.20). The full model had a c-index of 0.71 (95 percent c.i. 0.68 to 0.73) while the calibrated model had a c-index of 0.7 (95 percent c.i. 0.68 to 0.73). We assess calibation as previously using leave-one-center-out cross validation (Figure \@ref(fig:leave1outCalDeath)).



```r
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
centers <- unique(dat$center)
discrimination <- rep(0, length(centers))
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence60,
                      datTestCenter$status_RFS_FDA60)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  miDataResultTrain <- readRDS(file.path(".", "mi",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  miDataResultTest <- readRDS(file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 20,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  SDeathTrainCenter <- Surv(datTrainCenters$timeToFollowup60,
                            datTrainCenters$status_OS_FDA60)
  lp <- f.miTrainCenters$linear
  fDeath <- cph(SDeathTrainCenter ~ lp, 
                y = TRUE,
                x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(fDeath, TRUE)
  timePoint <- 60
  cIndex <- rep(0, 20)
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    lpTestCenter <- lpTestCenter*coef(fDeath)
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex[j] <- rc[1]
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  discrimination[i] <- mean(cIndex)
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  # baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = f.basehazTrainCenters,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToFollowup60,
                                status = datTestCenter$status_OS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalDeath-1.png" alt="\label{fig:leave1outCalDeath}Leave one center out cross validation of the calibrated model for 5-year overall mortality"  />
<p class="caption">(\#fig:leave1outCalDeath)\label{fig:leave1outCalDeath}Leave one center out cross validation of the calibrated model for 5-year overall mortality</p>
</div>

## Nomogram

We developed a 4-item score assigning points to each prognostic factor based on the magnitude of their association with recurrence. The nomogram to calculate the score is given in \@ref(fig:nomogram). The calibrated results for 5-year distant recurrence and 5-year mortality can be found in Figure \@ref(fig:combinedPlot).


```r
breslow.class <- c(0.1, .15, .25, .4, .6, 1, 1.5, 2.5, 4, 6, 8, 10)
t.class <- c(.05,.1, .2, .4, .8, 1.5, 2.5, 4, 6, 8)
age.class <- c(25, 30, 40, 50, 60, 75)
nom <- nomogram(f.mi,
                lp = FALSE,
                maxscale = 10,
                age = age.class,
                SN_tumor_burden = t.class,
                breslow = breslow.class)
g <- Newlabels(f.mi, c(ulceration = "Ulceration",
                       age = "Age (years)",
                       SN_tumor_burden = "Tumor burden (mm)",
                       breslow = "Breslow (mm)"))
nom1 <- nomogram(g, lp = FALSE,
                 maxscale = 10,
                 age = age.class,
                 SN_tumor_burden = t.class,
                 breslow = breslow.class)
# png("nomogram.png", width = 900, height = 700)
plot(nom1,
     total.sep.page = FALSE,
     col.grid = gray(c(.7, .9)),
     vnames = "names")
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/nomogram-1.png" alt="\label{fig:nomogram}Nomogram for 5-year recurrence"  />
<p class="caption">(\#fig:nomogram)\label{fig:nomogram}Nomogram for 5-year recurrence</p>
</div>

```r
# dev.off()
```



```r
dat.complete <- complete(miDataResult, 1)
breslowTrunc <- dat.complete$breslow
breslowTrunc[dat.complete$breslow < .2] <- .2
breslowTrunc[dat.complete$breslow > 10] <- 10
SN_tumor_burdenTrunc <- dat.complete$SN_tumor_burden
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden < .2] <- .2
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden > 7] <- 7
ageTrunc <- dat.complete$age
ageTrunc[dat.complete$age < 25] <- 25
ageTrunc[dat.complete$age > 75] <- 75
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResult,
                           data = miData,
                           n.impute = 20,
                           pr = FALSE,
                           fit.reps = TRUE,
                           y = TRUE,
                           x = TRUE,
                           se.fit = TRUE)
SDistant <- Surv(dat$timeToDistant60,
                 dat$status_DMFS_FDA60)
lp <- f.miRec$linear
fDistant <- cph(SDistant ~ lp, 
                y = TRUE,
                x = TRUE)
SDeath <- Surv(dat$timeToFollowup60, 
               dat$status_OS_FDA60)
fDeath <- cph(SDeath ~ lp,
              y = TRUE, 
              x = TRUE)
baseHazRec <- basehaz(f.miRec)
baseHazDistant <- basehaz(fDistant)
baseHazDeath <- basehaz(fDeath)
baseSurv <- estimtateBaselineSurvival(baseHazRec, 60)
baseSurvDistant <- estimtateBaselineSurvival(baseHazDistant, 60)
baseSurvDeath <- estimtateBaselineSurvival(baseHazDeath, 60)
interceptNom <- attr(nom, "info")$Intercept
pointwiseIncrease <- 1/attr(nom, "info")$sc
points <- 0:23
lp.score<-predict(f.miRec,
                  newdata=data.frame(ulceration=dat.complete$ulceration,
                                     age=ageTrunc,
                                     SN_tumor_burden=SN_tumor_burdenTrunc,
                                     breslow = breslowTrunc))
score <- (lp.score - interceptNom)/pointwiseIncrease
score <- ifelse(score > 22, 22, score)
score <- round((lp.score - interceptNom)/pointwiseIncrease, 0)
score[score>22] = 22
h <- hist(score,
          plot = FALSE,
          breaks = 0:24,
          right = FALSE)
recPoints <- interceptNom + pointwiseIncrease*points
distantPoints <- coef(fDistant)*(interceptNom + pointwiseIncrease*points)
deathPoints <- coef(fDeath)*(interceptNom + pointwiseIncrease*points)
recProbs <- 1 - baseSurv^exp(recPoints)
distantProbs <- 1 - baseSurvDistant^exp(distantPoints)
deathProbs <- 1 - baseSurvDeath^exp(deathPoints)
pointsToProbs <- data.frame(points = points,
                            recProbs = recProbs*100,
                            distantProbs = distantProbs*100,
                            deathProbs = deathProbs*100,
                            hist = h$density*100)
pointsToProbs$riskLevel <- c(rep(1, 10),
                             rep(2, 6),
                             rep(3, 4),
                             rep(4, 4))
pointsToProbs$riskLevel <- factor(pointsToProbs$riskLevel, labels = c(" < 25%", " 25-50%", " 50-75%", " > 75%"))
# write.csv(pointsToProbs, "newPointConversion.csv")
```


```r
ggplot(data = pointsToProbs, aes(x = points)) +
  geom_bar(mapping = aes(x = points-.5, y = hist*5, fill = factor(riskLevel)),
           stat = "identity", 
           width = 1, 
           alpha = .7,
           position = position_nudge(x = 0.5),
           color = "black") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(sec.axis = sec_axis(~./5, name = "% of patients")) +
  scale_x_continuous(breaks = seq(0, 25, 5)) +
  geom_line(aes(x = points, y = recProbs, col = "Recurrence"), 
            size = 1.2) +
  geom_line(aes(x = points, y = distantProbs, col = "Distant metastasis"), 
            size = 1.2) +
  geom_line(aes(x = points, y = deathProbs, col = "Death"),
            size = 1.2) +
  ylab("5-year probability (%)") + 
  xlab("Risk score") +
  theme_classic() +
  theme(legend.position = c(0.15, 0.75),
        legend.title = element_blank())
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/combinedPlot-1.png" alt="\label{fig:combinedPlot}Absolute risks along with risk distribution"  />
<p class="caption">(\#fig:combinedPlot)\label{fig:combinedPlot}Absolute risks along with risk distribution</p>
</div>

```r
ggsave("newCombinedPlot.png")
```

```
## Saving 8.5 x 7 in image
```


# External validation{#externalVal}


```r
dd <- datadist(dat)
options(datadist = 'dd')
options(digits = 8)
S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)
form <- S ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
summaryFit <- summary(f.mi)
namesSummaryFit <- rownames(summaryFit)
namesSummaryFit <- namesSummaryFit[seq(1, 15, 2)]
summaryFit <- summaryFit[seq(2, 16, 2), c(4, 6, 7)]
rownames(summaryFit) <- namesSummaryFit
colnames(summaryFit)[1] <- "Hazard ratio"
```

The variable `SN_tumor_burden` is not available in the valaidation set, but `SN_tumor_burden_extended` is. For that reason we first need to assess the performance of a model derived in the training set, where we substitute `SN_tumor_burden`. The new model has an apparent c-index of 0.68 (95 percent c.i. 0.65 to 0.7). The model to be validated in the validation set is given in Table \@ref(tab:validRec)


```r
kable(
  round(summaryFit, 2), 
  # "latex",
  # booktabs = TRUE,
  caption =  "Model for recurrence to be validated") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:validRec)Model for recurrence to be validated</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Hazard ratio </th>
   <th style="text-align:right;"> Lower 0.95 </th>
   <th style="text-align:right;"> Upper 0.95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> age </td>
   <td style="text-align:right;"> 1.28 </td>
   <td style="text-align:right;"> 1.13 </td>
   <td style="text-align:right;"> 1.46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> breslow </td>
   <td style="text-align:right;"> 1.42 </td>
   <td style="text-align:right;"> 1.23 </td>
   <td style="text-align:right;"> 1.62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ulceration - present:absent </td>
   <td style="text-align:right;"> 1.43 </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> 1.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden_extended - single cells:0.5 - 1.0 mm </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.30 </td>
   <td style="text-align:right;"> 0.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden_extended - &lt;0.5 mm:0.5 - 1.0 mm </td>
   <td style="text-align:right;"> 0.84 </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 1.13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden_extended - &gt;1.0 - 2.0 mm:0.5 - 1.0 mm </td>
   <td style="text-align:right;"> 1.19 </td>
   <td style="text-align:right;"> 0.90 </td>
   <td style="text-align:right;"> 1.57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden_extended - &gt;2.0 - 5.0 mm:0.5 - 1.0 mm </td>
   <td style="text-align:right;"> 1.57 </td>
   <td style="text-align:right;"> 1.20 </td>
   <td style="text-align:right;"> 2.06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden_extended - &gt;5.0 mm:0.5 - 1.0 mm </td>
   <td style="text-align:right;"> 1.65 </td>
   <td style="text-align:right;"> 1.21 </td>
   <td style="text-align:right;"> 2.24 </td>
  </tr>
</tbody>
</table>







The altered model for recurrence gave very similar performance in the validation set, with a c-index of 0.7 (95 percent c.i. 0.66 to 0.73). From the calibration plot (Figure \@ref(fig:calibrationPlotRec)) we see that there may be slight under-estimation for higher risk patients.



```r
ggplot(data = quintileResults, aes(x = mean, y = surv)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2, size = 1, alpha = .4) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Predicted probability") +
  ylab("Observed frequency") +
  ggtitle("External validation") +
  theme_bw()
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/calibrationPlotRec-1.png" alt="\label{fig:calibrationPlotRec}Calibration plot of 5-year recurrence model"  />
<p class="caption">(\#fig:calibrationPlotRec)\label{fig:calibrationPlotRec}Calibration plot of 5-year recurrence model</p>
</div>




For the case of distant metastasis, the c-index of the calibrated model in external validation was 0.71 (95 percent c.i. 0.68 to 0.75). The model performed very well in terms of calibration as well (Figure \@ref(fig:calibrationPlotDistant)). 


```r
pp$calibration + ggtitle("Calibration plot for distant metastasis")
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/calibrationPlotDistant-1.png" alt="\label{fig:calibrationPlotDistant}Calibration plot of 5-year distant metastasis model"  />
<p class="caption">(\#fig:calibrationPlotDistant)\label{fig:calibrationPlotDistant}Calibration plot of 5-year distant metastasis model</p>
</div>



Similar conclusions can be drawn for the case of overall mortality. The c-index of the calibrated model in external validation was 0.74 (95 percent c.i. 0.7 to 0.78) while calibration to the test set was again very good plot except for the case of high risk patients where risks were slightly underestimated (Figure \@ref(fig:calibrationPlotDeath))


```r
pp$calibration + ggtitle("Calibration plot for overall mortality")
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/calibrationPlotDeath-1.png" alt="\label{fig:calibrationPlotDeath}Calibration plot of 5-year overall mortality model"  />
<p class="caption">(\#fig:calibrationPlotDeath)\label{fig:calibrationPlotDeath}Calibration plot of 5-year overall mortality model</p>
</div>

# Tertiary objective

The tertiary goal is to create a prediction model for recurrence (and DMFS and OS) using the extra variable positive additional positive nodes after CLND. This model can only be based on those patients who actually underwent a CLND of course. This would tell us something about how much this variable would add to the discrimination of the prediction model. 



More specifically, in the case of 5-year recurrence a model including `ulceration`, `age`, `SN_tumor_burden`, `breslow` AND `no_pos_nonSNs_cat` has a c-index of 0.69 with a 95 percent c.i. 0.67 to 0.72, which a slight increase from 0.68 (95 percent c.i. 0.65 to 0.7).

For 5-year distant metastasis the c-index is 0.72 with a 95 percent c.i. 0.69 to 0.74 compared to 0.7 (95 percent c.i. 0.67 to 0.72) of the origial model.



Finally, for 5-year overall mortality the c-index is 0.72 with a 95 percent c.i. 0.69 to 0.75 comapared to 0.70 (95 percent c.i. 0.67 to 0.72) of the origial model.

When we use the simpler EJC risk groups to assign patients to risk categories the performance drops significantly. In Table \@ref(tab:cIndexOverallModels) we give an overview of the perfomrmance regarding discrimination for the different modeling approaches considered.


```r
tableData <- data.frame(Models = c("EORTC prediction model", "EORTC prediction model with no nonSN's", "EORTC - simple EJC groups"),
                        Recurrence = c("0.68 (0.65 - 0.70)", "0.69 (0.67 - 0.72)", "0.61 (0.59 - 0.63)"),
                        Distant_Metastasis = c("0.70 (0.67 - 0.72)", "0.72 (0.69 - 0.74)", "0.61 (0.59 - 0.63)"),
                        Overall_Mortality = c("0.70 (0.67 - 0.73)", "0.72 (0.69 - 0.75)", "0.60 (0.58 - 0.63)"))
kable(
  tableData, 
  # "latex",
  caption =  "C-indices for the different prediction models considered internal discrimination",
  booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:cIndexOverallModels)C-indices for the different prediction models considered internal discrimination</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Models </th>
   <th style="text-align:left;"> Recurrence </th>
   <th style="text-align:left;"> Distant_Metastasis </th>
   <th style="text-align:left;"> Overall_Mortality </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> EORTC prediction model </td>
   <td style="text-align:left;"> 0.68 (0.65 - 0.70) </td>
   <td style="text-align:left;"> 0.70 (0.67 - 0.72) </td>
   <td style="text-align:left;"> 0.70 (0.67 - 0.73) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EORTC prediction model with no nonSN's </td>
   <td style="text-align:left;"> 0.69 (0.67 - 0.72) </td>
   <td style="text-align:left;"> 0.72 (0.69 - 0.74) </td>
   <td style="text-align:left;"> 0.72 (0.69 - 0.75) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EORTC - simple EJC groups </td>
   <td style="text-align:left;"> 0.61 (0.59 - 0.63) </td>
   <td style="text-align:left;"> 0.61 (0.59 - 0.63) </td>
   <td style="text-align:left;"> 0.60 (0.58 - 0.63) </td>
  </tr>
</tbody>
</table>

# Merging data sets

Due to the very good performance of the developed model in the test dataset, we combine the training with the test datasets to for the development of the final model.

## Recurrence




```r
S <- Surv(datMerged$timeToRecurrence,
          datMerged$status_RFS_FDA)
form <- S ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResultMerged,
                        data = miDataMerged,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
# tableCount <- tableCount + 1
# kable(round(summary(f.mi)[, c(4, 6, 7)], 2), caption = paste("Table", paste0(tableCount, ":"), "Multivariable Cox analysis of recurrence"))
```


```r
S <- Surv(datMerged$timeToRecurrence60,
          datMerged$status_RFS_FDA60)
form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResultMerged,
                        data = miDataMerged,
                        n.impute = 20,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
summaryFit <- summary(f.mi)
namesSummaryFit <- rownames(summaryFit)[seq(1, 7, 2)]
summaryFit <- summaryFit[seq(2, 8, 2), c(4, 6, 7)]
rownames(summaryFit) <- namesSummaryFit
colnames(summaryFit)[1] <- "Hazard ratio"
kable(round(summaryFit, 2),
      # "latex",
      # booktabs = TRUE,
      caption = "Final model for 5-year recurrence using the merged dataset") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:finalRec5Merged)Final model for 5-year recurrence using the merged dataset</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Hazard ratio </th>
   <th style="text-align:right;"> Lower 0.95 </th>
   <th style="text-align:right;"> Upper 0.95 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> age </td>
   <td style="text-align:right;"> 1.42 </td>
   <td style="text-align:right;"> 1.27 </td>
   <td style="text-align:right;"> 1.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SN_tumor_burden </td>
   <td style="text-align:right;"> 1.37 </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.89 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> breslow </td>
   <td style="text-align:right;"> 1.53 </td>
   <td style="text-align:right;"> 1.33 </td>
   <td style="text-align:right;"> 1.76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ulceration - present:absent </td>
   <td style="text-align:right;"> 1.43 </td>
   <td style="text-align:right;"> 1.20 </td>
   <td style="text-align:right;"> 1.70 </td>
  </tr>
</tbody>
</table>

The final model for the prediction of 5-year recurrence from the merged set can be found in Table  \@ref(tab:finalRec5Merged). The c-index of that model  is 0.68 (95 percent c.i. 0.65 to 0.72). Calibration is assessed using a leave-one-center-out cross validation, based on the merged dataset. In this way, information from the 9 development centers of the first objective can be used to impute missing values for `SN_tumor_burden` in the German (validation) dataset. Multiple imputation is again separated between the training and the test set as was done in the original analysis. When leaving out the German dataset, because we cannot use the data from the other centers to impute the missing values, we use the approach of [External validation](#externalVal) (the same holds for the assessment of calibration in the calibrated models for distant metastasis and overall mortality). The results of this analysis can be found in Figure \@ref(fig:leave1outCalMerged) for the 5-year recurrence model, in Figure \@ref(fig:leave1outCalDistantMerged) for the 5-year distant recurrence model and in Figure \@ref(fig:leave1outCalDeathMerged) for the 5-year overall mortality model.


```r
discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3, 3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(datMerged, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(datMerged, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence,
                      datTestCenter$status_RFS_FDA)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  # miDataResultTrain <- mice(data = miDataTrain, m = 20)
  # saveRDS(miDataResultTrain, file.path(".", "miMergedCV",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataResultTrain <- readRDS(file.path(".", "miMergedCV",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  if(any(table(datTestCenter$AJCC_sub_SLNB_8th) == 0))
    miDataTest <- subset(miDataTest, select = -AJCC_sub_SLNB_8th)
  
  # miDataResultTest <- mice(data = miDataTest, m = 20)
  # saveRDS(miDataResultTest, file.path(".", "miMergedCV",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  miDataResultTest <- readRDS(file.path(".", "miMergedCV",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(f.miTrainCenters, TRUE)
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <- cbind(linearPredictor, lpTestCenter)
    
  }
  
  # discrimination[[i]] <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = baselineHazard,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToRecurrence60,
                                status = datTestCenter$status_RFS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalMerged-1.png" alt="\label{fig:leave1outCalMerged}Leave-one-center-out cross validation for prediction of recurrence in the merged dataset"  />
<p class="caption">(\#fig:leave1outCalMerged)\label{fig:leave1outCalMerged}Leave-one-center-out cross validation for prediction of recurrence in the merged dataset</p>
</div>

## Nomogram

The updated nomogram based on the combined dataset can be found in Figure  \@ref(fig:nomogramMerged).


```r
breslow.class <- c(0.1, .15, .25, .4, .6, 1, 1.5, 2.5, 4, 6, 8, 10)
t.class <- c(.05,.1, .2, .4, .8, 1.5, 2.5, 4, 8)
age.class <- c(25, 30, 40, 50, 60, 75)
nom <- nomogram(f.mi,
                lp = FALSE,
                maxscale = 10,
                age = age.class,
                SN_tumor_burden = t.class,
                breslow = breslow.class)
plot(nom,
     total.sep.page = FALSE,
     col.grid =gray(c(.7, .9)))
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/nomogramMerged-1.png" alt="\label{fig:nomogramMerged}Nomogram for 5-year recurrence for the merged dataset"  />
<p class="caption">(\#fig:nomogramMerged)\label{fig:nomogramMerged}Nomogram for 5-year recurrence for the merged dataset</p>
</div>

## Distant metastasis



```r
SDistantMerged <- Surv(datMerged$timeToDistant60,
                       datMerged$status_DMFS_FDA60)
form <- SDistantMerged ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)
f.miDistantFullMerged <- fit.mult.impute(form,
                                         cph,
                                         xtrans = miDataResultMerged,
                                         data = miDataMerged,
                                         n.impute = 20,
                                         pr = FALSE,
                                         fit.reps = TRUE,
                                         y = TRUE,
                                         x = TRUE,
                                         se.fit = TRUE)
SRecurrenceMerged  <- Surv(datMerged$timeToRecurrence60,
                           datMerged$status_RFS_FDA60)
form <- SRecurrenceMerged ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRecMerged <- fit.mult.impute(form,
                                 cph,
                                 xtrans = miDataResultMerged,
                                 data = miDataMerged,
                                 n.impute = 20,
                                 pr = FALSE,
                                 fit.reps = TRUE,
                                 y = TRUE,
                                 x = TRUE,
                                 se.fit = TRUE)
linearPredictor <- f.miRecMerged$linear
dd <- datadist(linearPredictor)
options(datadist = 'dd')
options(digits = 8)
fDistantCalibratedMerged <- cph(SDistantMerged ~ linearPredictor, 
                                y = TRUE,
                                x = TRUE)
cIndexCalibratedMerged <- rcorr.cens(-linearPredictor*coef(fDistantCalibratedMerged), SDistantMerged)
```

For the assessment of distant metastasis we considered a calibrated version of the 5-year recurrence model. The association between distant metastasis and was of the same size (calibration slope of 1.01, 95 percent c.i. 0.87 to 1.16). We compare the performance of the considered model to that of multivariable Cox regression model including all 9 covariates of interest (`sex`, `ulceration`, `no_removed_SNs`, `no_pos_SNs`, `simpleloc`, `histology_simple`, `breslow`, `SN_tumor_burden` and `age`). The full model had a c-index of 0.7 (95 percent c.i. 0.66 to 0.74) while the calibrated model had a c-index of 0.71 (95 percent c.i. 0.69 to 0.73).

Again we perform a leave-one-center-out cross validation to assess the calibration of our final model for distant metastasis. The approach is the same as the one considered for the 5-year recurrence model (Figure \@ref(fig:leave1outCalDistantMerged)).


```r
SRecurrence  <- Surv(datMerged$timeToRecurrence60,
                     datMerged$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(datMerged, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(datMerged, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence60,
                      datTestCenter$status_RFS_FDA60)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  miDataResultTrain <- readRDS(file.path(".", "miMergedCV",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  if(any(table(datTestCenter$AJCC_sub_SLNB_8th) == 0))
    miDataTest <- subset(miDataTest, select = -AJCC_sub_SLNB_8th)
  
  miDataResultTest <- readRDS(file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  SDistantTrainCenter <- Surv(datTrainCenters$timeToDistant60,
                              datTrainCenters$status_DMFS_FDA60)
  lp <- f.miTrainCenters$linear
  fDistant <- cph(SDistantTrainCenter ~ lp, 
                  y = TRUE,
                  x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(fDistant, TRUE) 
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    lpTestCenter <- lpTestCenter*coef(fDistant)
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  # discrimination[[i]] <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  # baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = f.basehazTrainCenters,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToDistant60,
                                status = datTestCenter$status_DMFS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalDistantMerged-1.png" alt="\label{fig:leave1outCalDistantMerged}Leave-one-center-out cross validation for prediction of distant metastasis in the merged dataset"  />
<p class="caption">(\#fig:leave1outCalDistantMerged)\label{fig:leave1outCalDistantMerged}Leave-one-center-out cross validation for prediction of distant metastasis in the merged dataset</p>
</div>

## Overall mortality


```r
SDeathMerged <- Surv(datMerged$timeToFollowup60,
                     datMerged$status_OS_FDA60)
form <- SDeathMerged ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)
f.miDeathFullMerged <- fit.mult.impute(form,
                                       cph,
                                       xtrans = miDataResultMerged,
                                       data = miDataMerged,
                                       n.impute = 20,
                                       pr = FALSE,
                                       fit.reps = TRUE,
                                       y = TRUE,
                                       x = TRUE,
                                       se.fit = TRUE)
SRecurrence  <- Surv(datMerged$timeToRecurrence60,
                     datMerged$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRecMerged <- fit.mult.impute(form,
                                 cph,
                                 xtrans = miDataResultMerged,
                                 data = miDataMerged,
                                 n.impute = 20,
                                 pr = FALSE,
                                 fit.reps = TRUE,
                                 y = TRUE,
                                 x = TRUE,
                                 se.fit = TRUE)
linearPredictor <- f.miRecMerged$linear
dd <- datadist(linearPredictor)
options(datadist = 'dd')
options(digits = 8)
fDeathCalibratedMerged <- cph(SDeathMerged ~ linearPredictor, 
                              y = TRUE,
                              x = TRUE)
cIndexCalibratedMerged <- rcorr.cens(-linearPredictor*coef(fDeathCalibratedMerged), SDeathMerged)
```

For the assessment 5-year overall mortality in the merged dataset we considered a calibrated version of the 5-year recurrence model. We compare the performance of the considered model to that of multivariable Cox regression model including all 9 covariates of interest. The association between recurrence and overall mortality was not different (calibration slope 1.04, 95 percent c.i. 0.88 to 1.20). The full model had a c-index of 0.7 (95 percent c.i. 0.65 to 0.75) while the calibrated model had a c-index of 0.71 (95 percent c.i. 0.69 to 0.74). We assess calibation as previously using leave-one-center-out cross validation.

The results of leave-one-center-out cross validation can be found in Figure \@ref(fig:leave1outCalDeathMerged).


```r
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(datMerged, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(datMerged, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence60,
                      datTestCenter$status_RFS_FDA60)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  
  miDataTrain <- datTrainCenters[miVariables]
  miDataTrain$timeToRecurrence <- log(miDataTrain$timeToRecurrence)
  miDataTrain$timeToDistant <- log(miDataTrain$timeToDistant)
  miDataTrain$timeToFollowup <- log(miDataTrain$timeToFollowup)
  
  miDataResultTrain <- readRDS(file.path(".", "miMergedCV",paste("miDataResultTrain", centers[i], "rds", sep = ".")))
  
  miDataTest <- datTestCenter[miVariables]
  miDataTest$timeToRecurrence <- log(miDataTest$timeToRecurrence)
  miDataTest$timeToDistant <- log(miDataTest$timeToDistant)
  miDataTest$timeToFollowup <- log(miDataTest$timeToFollowup)
  
  if(any(table(datTestCenter$AJCC_sub_SLNB_8th) == 0))
    miDataTest <- subset(miDataTest, select = -AJCC_sub_SLNB_8th)
  
  miDataResultTest <- readRDS(file.path(".", "mi",paste("miDataResultTest", centers[i], "rds", sep = ".")))
  
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResultTrain,
                                      data = miDataTrain,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE)
  
  SDeathTrainCenter <- Surv(datTrainCenters$timeToFollowup60,
                            datTrainCenters$status_OS_FDA60)
  lp <- f.miTrainCenters$linear
  fDeath <- cph(SDeathTrainCenter ~ lp, 
                y = TRUE,
                x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(fDeath, TRUE)
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:20){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResultTest, j))
    lpTestCenter <- lpTestCenter*coef(fDeath)
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  # discrimination[[i]] <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  # baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = f.basehazTrainCenters,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = datTestCenter$timeToFollowup60,
                                status = datTestCenter$status_OS_FDA60)
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       ylab = "Observed proportions",
       xlab = "Predicted probabilities",
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/leave1outCalDeathMerged-1.png" alt="\label{fig:leave1outCalDeathMerged}Leave-one-center-out cross validation for prediction of overall mortality in the merged dataset"  />
<p class="caption">(\#fig:leave1outCalDeathMerged)\label{fig:leave1outCalDeathMerged}Leave-one-center-out cross validation for prediction of overall mortality in the merged dataset</p>
</div>


```r
dat.complete <- complete(miDataResultMerged, 1)
breslowTrunc <- dat.complete$breslow
breslowTrunc[dat.complete$breslow < .2] <- .2
breslowTrunc[dat.complete$breslow > 10] <- 10
SN_tumor_burdenTrunc <- dat.complete$SN_tumor_burden
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden < .2] <- .2
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden > 7] <- 7
ageTrunc <- dat.complete$age
ageTrunc[dat.complete$age < 10] <- 10
ageTrunc[dat.complete$age > 70] <- 70
SRecurrence  <- Surv(datMerged$timeToRecurrence60,
                     datMerged$status_RFS_FDA60)
form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResultMerged,
                           data = miDataMerged,
                           n.impute = 20,
                           pr = FALSE,
                           fit.reps = TRUE,
                           y = TRUE,
                           x = TRUE,
                           se.fit = TRUE)
SDistant <- Surv(datMerged$timeToDistant60,
                 datMerged$status_DMFS_FDA60)
lp <- f.miRec$linear
fDistant <- cph(SDistant ~ lp, 
                y = TRUE,
                x = TRUE)
SDeath <- Surv(datMerged$timeToFollowup60, 
               datMerged$status_OS_FDA60)
fDeath <- cph(SDeath ~ lp,
              y = TRUE, 
              x = TRUE)
baseHazRec <- basehaz(f.miRec)
baseHazDistant <- basehaz(fDistant)
baseHazDeath <- basehaz(fDeath)
baseSurv <- estimtateBaselineSurvival(baseHazRec, 60)
baseSurvDistant <- estimtateBaselineSurvival(baseHazDistant, 60)
baseSurvDeath <- estimtateBaselineSurvival(baseHazDeath, 60)
interceptNom <- attr(nom, "info")$Intercept
pointwiseIncrease <- 1/attr(nom, "info")$sc
points <- 0:21
lp.score<-predict(f.miRec,
                  newdata=data.frame(ulceration=dat.complete$ulceration,
                                     age=ageTrunc,
                                     SN_tumor_burden=SN_tumor_burdenTrunc,
                                     breslow = breslowTrunc))
score <- round((lp.score - interceptNom)/pointwiseIncrease, 0)
h <- hist(score,
          plot = FALSE,
          breaks = 0:22,
          right = FALSE)
recPoints <- interceptNom + pointwiseIncrease*points
distantPoints <- coef(fDistant)*(interceptNom + pointwiseIncrease*points)
deathPoints <- coef(fDeath)*(interceptNom + pointwiseIncrease*points)
recProbs <- 1 - baseSurv^exp(recPoints)
distantProbs <- 1 - baseSurvDistant^exp(distantPoints)
deathProbs <- 1 - baseSurvDeath^exp(deathPoints)
pointsToProbs <- data.frame(points = points,
                            recProbs = recProbs*100,
                            distantProbs = distantProbs*100,
                            deathProbs = deathProbs*100,
                            hist = h$density*100)
pointsToProbs$riskLevel <- c(rep(1, 9),
                             rep(2, 4),
                             rep(3, 3),
                             rep(4, 6))
pointsToProbs$riskLevel <- factor(pointsToProbs$riskLevel, labels = c(" < 25%", " 25-50%", " 50-75%", " > 75%"))
```

The combined results regarding the calibrated models for 5-year distant recurrence and 5-year overall mortality along with the predictions from the 5-year recurrence model can be found in Figure \@ref(fig:combinedPlotMerged).


```r
ggplot(data = pointsToProbs, aes(x = points)) +
  geom_bar(mapping = aes(x = points, y = hist*4, fill = factor(riskLevel)),
           stat = "identity",
           alpha = .7,
           width = 1, 
           position = position_nudge(x = 0.5),
           color = "black") +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(sec.axis = sec_axis(~./4, name = "% of patients", breaks = seq(0, 25, 5))) +
  geom_line(aes(x = points, y = recProbs, col = "Recurrence"), 
            size = 1.2) +
  geom_line(aes(x = points, y = distantProbs, col = "Distant metastasis"), 
            size = 1.2) +
  geom_line(aes(x = points, y = deathProbs, col = "Death"),
            size = 1.2) +
  ylab("5-year probability (%)") + 
  xlab("Risk score") +
  theme_classic() +
  theme(legend.position = c(0.15, 0.75),
        legend.title = element_blank())
```

<div class="figure" style="text-align: center">
<img src="positiveSLNB_files/figure-html/combinedPlotMerged-1.png" alt="\label{fig:combinedPlotMerged}Absolute risks along with risk distribution using the merged dataset"  />
<p class="caption">(\#fig:combinedPlotMerged)\label{fig:combinedPlotMerged}Absolute risks along with risk distribution using the merged dataset</p>
</div>




# Reviewer comments

## Simple risk subgroups

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-21)5-year recurrence-free survival with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.48 </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.56 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.45 </td>
   <td style="text-align:right;"> 0.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.27 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:right;"> 0.33 </td>
  </tr>
</tbody>
</table>

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-21)5-year distant metastasis free survival with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.74 </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.79 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.46 </td>
   <td style="text-align:right;"> 0.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.66 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.37 </td>
  </tr>
</tbody>
</table>

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-21)5-year overall survival with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.48 </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.56 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.45 </td>
   <td style="text-align:right;"> 0.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.27 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:right;"> 0.33 </td>
  </tr>
</tbody>
</table>

## AJCC classifications

```r
eortcData <- read.spss("./Data/EORTC_CLEAN.sav",
                       to.data.frame = TRUE)
```

```
## Warning in read.spss("./Data/EORTC_CLEAN.sav", to.data.frame = TRUE): ./
## Data/EORTC_CLEAN.sav: Very long string record(s) found (record type 7,
## subtype 14), each will be imported in consecutive separate variables
```

```
## re-encoding from UTF-8
```

```
## Warning in read.spss("./Data/EORTC_CLEAN.sav", to.data.frame = TRUE):
## Undeclared level(s) 9 added in variable: RTx
```

```r
eortcData$id <- 1:dim(eortcData)[1] 
## Set dates 1582-10-14 to NA
eortcData$Date_dist_met <- ifelse(eortcData$Date_dist_met == 0, NA, eortcData$Date_dist_met)
eortcData$first_rec <- ifelse(eortcData$first_rec == 0, NA, eortcData$first_rec)
## Transform into dates and create new times to event
eortcData <- eortcData %>%
  mutate(last_FU = spss2Date(last_FU), 
         SLNBdate = spss2Date(SLNBdate),
         first_rec = spss2Date(first_rec),
         Date_dist_met = spss2Date(Date_dist_met),
         recurrence = ifelse(recurrence == "yes", 1, 0)) 
## problematic observations -> typos with dates could be the reason
eortcData[eortcData$id == 344, ]$last_FU <- eortcData[eortcData$id == 344, ]$Date_dist_met
eortcData[eortcData$id == 926, ]$first_rec <- eortcData[eortcData$id == 926, ]$Date_dist_met <- eortcData[eortcData$id == 926, ]$last_FU
eortcData[eortcData$id == 1042, ]$first_rec <- eortcData[eortcData$id == 1042, ]$Date_dist_met <- eortcData[eortcData$id == 1042, ]$last_FU
eortcData <- eortcData %>%
  mutate(timeToFollowup = as.numeric((last_FU - SLNBdate)/30.5),
         timeToRecurrence = as.numeric((first_rec - SLNBdate)/30.5),
         timeToDistant = as.numeric(Date_dist_met - SLNBdate)/30.5)
eortcData$timeToRecurrence <- ifelse(is.na(eortcData$timeToRecurrence),
                                     as.numeric(eortcData$last_FU - eortcData$SLNBdate)/30.5,
                                     eortcData$timeToRecurrence)
eortcData$timeToDistant <- ifelse(is.na(eortcData$timeToDistant),
                                  as.numeric(eortcData$last_FU - eortcData$SLNBdate)/30.5,
                                  eortcData$timeToDistant)
eortcData$clark <- factor(eortcData$clark,
                          levels = c("ii", "iii", "iv", "v")) # dropping declared levels "i" and "unknown" that are never encountered
## exclude weird observations
eortcData <- eortcData %>%
  filter(!(id %in% c(7, 310, 863)))
# glimpse(eortcData)
##################################
# Extract data
##################################
selectedVariables <-c("id", "center", "age", "sex", "simpleloc", "clark", "CLND",
                      "histology_simple", "histology_extended", "breslow",
                      "ulceration", "SLNBdate", "no_removed_SNs", "AJCC_sub_SLNB_8th",
                      "no_pos_SNs", "no_pos_SNs_cat", "SN_tumor_burden",
                      "SN_tumor_burden_extended", "SN_tumor_burden_simple",
                      "no_removed_nonSNs", "no_pos_nonSNs", "no_pos_nonSNs_cat",
                      "EJC_groups", "Risk_classes", "AJCC_7th_extra", "AJCC_8th_extra")
outcomes <- c("status_RFS_FDA", "timeToRecurrence", 
              "status_DMFS_FDA", "timeToDistant",
              "status_OS_FDA", "timeToFollowup",
              "status_RFS_FDA60", "status_DMFS_FDA60",
              "status_OS_FDA60")
statusOutcomesoutcomes <- c("status_RFS_FDA", 
                            "status_DMFS_FDA", 
                            "status_OS_FDA")
dat <- eortcData %>%
  select(c(selectedVariables, outcomes))
dat$AJCC_sub_SLNB_8th <- gsub(" ", "", dat$AJCC_sub_SLNB_8th)
dat$timeToRecurrence60 <- ifelse(eortcData$timeToRecurrence < 60, 
                                 eortcData$timeToRecurrence,
                                 60)
dat$timeToDistant60 <- ifelse(eortcData$timeToDistant < 60,
                              eortcData$timeToDistant,
                              60)
dat$timeToFollowup60 <- ifelse(eortcData$timeToFollowup < 60,
                               eortcData$timeToFollowup,
                               60)

newData <- subset(dat, !is.na(AJCC_7th_extra))
cIndexAJCC_7th_Rec <- cIndexAJCC_7th(data = newData,
                                     timeToOutcome = "timeToRecurrence60",
                                     outcomeStatus = "status_RFS_FDA60")
cIndexAJCC_7th_Distant <- cIndexAJCC_7th(data = newData,
                                     timeToOutcome = "timeToRecurrence60",
                                     outcomeStatus = "status_DMFS_FDA60")
cIndexAJCC_7th_Death <- cIndexAJCC_7th(data = newData,
                                     timeToOutcome = "timeToFollowup",
                                     outcomeStatus = "status_OS_FDA60")

newData <- subset(dat, !is.na(AJCC_8th_extra))
cIndexAJCC_8th_Rec <- cIndexAJCC_8th(data = newData,
                                     timeToOutcome = "timeToRecurrence60",
                                     outcomeStatus = "status_RFS_FDA60")
cIndexAJCC_8th_Distant <- cIndexAJCC_8th(data = newData,
                                     timeToOutcome = "timeToRecurrence60",
                                     outcomeStatus = "status_DMFS_FDA60")
cIndexAJCC_8th_Death <- cIndexAJCC_8th(data = newData,
                                     timeToOutcome = "timeToFollowup",
                                     outcomeStatus = "status_OS_FDA60")
```

Estimating probabilities of AJCC 7th subgroups


```r
newData <- subset(dat, !is.na(AJCC_7th_extra))
SRec <- Surv(newData$timeToRecurrence60,
          newData$status_RFS_FDA60)
form <- SRec ~ AJCC_7th_extra
f.AJCC_7th_Rec <- cph(formula = form,
                      data = newData,
                      x = TRUE,
                      y = TRUE)

runSurvEstimation <- survest(f.AJCC_7th_Rec,
                          newdata = data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra)), times = 60)
resAJCC_7th_Rec <- data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_7th_extra)




# Distant metastasis

SDistant <- Surv(newData$timeToDistant60,
          newData$status_DMFS_FDA60)
form <- SDistant ~ AJCC_7th_extra
f.AJCC_7th_Distant <- cph(formula = form,
                      data = newData,
                      x = TRUE,
                      y = TRUE)

runSurvEstimation <- survest(f.AJCC_7th_Distant,
                          newdata = data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra)), times = 60)
resAJCC_7th_Distant <- data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_7th_extra)





# Overall survival
SDeath <- Surv(newData$timeToFollowup60,
          newData$status_OS_FDA60)
form <- SDeath ~ AJCC_7th_extra
f.AJCC_7th_Death <- cph(formula = form,
                      data = newData,
                      x = TRUE,
                      y = TRUE)

runSurvEstimation <- survest(f.AJCC_7th_Death,
                          newdata = data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra)), times = 60)
resAJCC_7th_Death <- data.frame(AJCC_7th_extra = unique(newData$AJCC_7th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_7th_extra)
```

Estimating probabilities of AJC8 7th subgroups


```r
newData <- subset(dat, !is.na(AJCC_8th_extra))
SRec <- Surv(newData$timeToRecurrence60,
             newData$status_RFS_FDA60)
form <- SRec ~ AJCC_8th_extra
f.AJCC_8th_Rec <- cph(formula = form,
                      data = newData,
                      x = TRUE,
                      y = TRUE)

runSurvEstimation <- survest(f.AJCC_8th_Rec,
                          newdata = data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra)), times = 60)
resAJCC_8th_Rec <- data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_8th_extra)



# Distant metastasis

SDistant <- Surv(newData$timeToDistant60,
                 newData$status_DMFS_FDA60)
form <- SDistant ~ AJCC_8th_extra
f.AJCC_8th_Distant <- cph(formula = form,
                          data = newData,
                          x = TRUE,
                          y = TRUE)

runSurvEstimation <- survest(f.AJCC_8th_Distant,
                          newdata = data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra)), times = 60)
resAJCC_8th_Distant <- data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_8th_extra)




# Overall survival
SDeath <- Surv(newData$timeToFollowup60,
               newData$status_OS_FDA60)
form <- SDeath ~ AJCC_8th_extra
f.AJCC_8th_Death <- cph(formula = form,
                        data = newData,
                        x = TRUE,
                        y = TRUE)

runSurvEstimation <- survest(f.AJCC_8th_Death,
                          newdata = data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra)), times = 60)
resAJCC_8th_Death <- data.frame(AJCC_8th_extra = unique(newData$AJCC_8th_extra), 
                              surv = runSurvEstimation$surv, 
                              lower = runSurvEstimation$lower, 
                              upper = runSurvEstimation$upper) %>%
  arrange(AJCC_8th_extra)

kable(resAJCC_7th_Rec %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year recurrence-free survival with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year recurrence-free survival with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.40 </td>
   <td style="text-align:right;"> 0.16 </td>
   <td style="text-align:right;"> 0.98 </td>
  </tr>
</tbody>
</table>

```r
kable(resAJCC_7th_Distant %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year distant metastasis free survival with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year distant metastasis free survival with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.75 </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.54 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.98 </td>
  </tr>
</tbody>
</table>

```r
kable(resAJCC_7th_Death %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year overall survival with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year overall survival with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.80 </td>
   <td style="text-align:right;"> 0.75 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.54 </td>
   <td style="text-align:right;"> 0.68 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.47 </td>
   <td style="text-align:right;"> 0.56 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;"> 0.14 </td>
   <td style="text-align:right;"> 0.98 </td>
  </tr>
</tbody>
</table>

```r
kable(resAJCC_8th_Rec %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year recurrence-free survival with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year recurrence-free survival with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.73 </td>
   <td style="text-align:right;"> 0.66 </td>
   <td style="text-align:right;"> 0.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.57 </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.36 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.34 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
</tbody>
</table>

```r
kable(resAJCC_8th_Distant %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year distant metastasis free survival with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year distant metastasis free survival with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:right;"> 0.72 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.66 </td>
   <td style="text-align:right;"> 0.54 </td>
   <td style="text-align:right;"> 0.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.65 </td>
   <td style="text-align:right;"> 0.59 </td>
   <td style="text-align:right;"> 0.71 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
</tbody>
</table>

```r
kable(resAJCC_8th_Death %>%
        mutate(surv = round(surv, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year overall survival with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-24)5-year overall survival with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> surv </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.85 </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:right;"> 0.91 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.73 </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.45 </td>
   <td style="text-align:right;"> 0.55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.30 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
</tbody>
</table>

