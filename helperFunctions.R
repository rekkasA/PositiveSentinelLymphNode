# Multiply imputed c index

mi.cindex <- function(f.mi,
                      miDataResult){
  
  cindex <- NULL
  m <- miDataResult$m
  for (i in 1:m){
    
    rc <- Hmisc::rcorr.cens(-predict(f.mi, newdata = complete(miDataResult, i)),
                            f.mi$fits[[i]]$y)
    cindex<-rbind(cindex,
                  c(rc["C Index"], rc["S.D."]/2))
    
  }
  summary(mitools::MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))
  
}

cIndexAJCC_7th <- function(data, timeToOutcome, outcomeStatus){
  
  data$new <- 0
  data[data$AJCC_7th_extra == levels(data$AJCC_7th_extra)[1], ]$new <- 1
  data[data$AJCC_7th_extra == levels(data$AJCC_7th_extra)[2], ]$new <- 2
  data[data$AJCC_7th_extra == levels(data$AJCC_7th_extra)[3], ]$new <- 3
  data[data$AJCC_7th_extra == levels(data$AJCC_7th_extra)[4], ]$new <- 4
  rc <- Hmisc::rcorr.cens(-as.numeric(data$new),
                          Surv(unlist(data[timeToOutcome]), unlist(data[outcomeStatus])))

  cIndex <-  c(rc["C Index"], rc["S.D."]/2)
  
  return(cIndex)
}
cIndexAJCC_8th <- function(data, timeToOutcome, outcomeStatus){
  
  data$new <- 0
  data[data$AJCC_8th_extra == levels(data$AJCC_8th_extra)[1], ]$new <- 1
  data[data$AJCC_8th_extra == levels(data$AJCC_8th_extra)[2], ]$new <- 2
  data[data$AJCC_8th_extra == levels(data$AJCC_8th_extra)[3], ]$new <- 3
  data[data$AJCC_8th_extra == levels(data$AJCC_8th_extra)[4], ]$new <- 4
  data[data$AJCC_8th_extra == levels(data$AJCC_8th_extra)[5], ]$new <- 5
  rc <- Hmisc::rcorr.cens(-as.numeric(data$new),
                          Surv(unlist(data[timeToOutcome]), unlist(data[outcomeStatus])))

  cIndex <-  c(rc["C Index"], rc["S.D."]/2)
  
  return(cIndex)
}

cIndexEJC <- function(miDataResult, data, timeToOutcome, outcomeStatus){
  
  cindex <- NULL
  m <- miDataResult$m
  for (i in 1:m){
    
    completeData <- complete(miDataResult, i)
    completeData$new <- 0
    completeData[completeData$EJC_groups == "No ulceration + low TB", ]$new <- 1
    completeData[completeData$EJC_groups == "No ulceration + high TB", ]$new <- 2
    completeData[completeData$EJC_groups == "Ulceration + low TB", ]$new <- 3
    completeData[completeData$EJC_groups ==  "Ulceration + high TB", ]$new <- 4
    rc <- Hmisc::rcorr.cens(-as.numeric(completeData$new),
                            Surv(unlist(data[timeToOutcome]), unlist(data[outcomeStatus])))
    cindex<-rbind(cindex,
                  c(rc["C Index"], rc["S.D."]/2))
    
  }
  summary(mitools::MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))
  
}



# S(t|x) = S0(t)^exp(lp) ; S0(t) = exp(-Lambda0(t))

estimtateBaselineSurvival <- function(baselineHazard, 
                                      timePoint){
  
  baselineHazard <- subset(baselineHazard, time <= timePoint)
  maxTimeLocation <- dim(baselineHazard)[1]
  baselineSurvival <- exp(-baselineHazard$hazard[maxTimeLocation])
  return(baselineSurvival)
  
}
predictSurvival <- function(baselineHazard, 
                            linearPredictor, 
                            timePoint){
  
  baselineSurvival <- estimtateBaselineSurvival(baselineHazard = baselineHazard,
                                                timePoint = timePoint)
  survivalProbability <- baselineSurvival^exp(linearPredictor)
  
  return(survivalProbability)
  
}


getKaplanMeier <- function(time, status, timePoint){
  S <- Surv(time, status)
  summary(survfit(S ~ 1), times = timePoint)$sur
}

upperKaplanMeier <- function(time, status, timePoint){
  S <- Surv(time, status)
  summary(survfit(S ~ 1), times = timePoint)$upper
}

lowerKaplanMeier <- function(time, status, timePoint){
  S <- Surv(time, status)
  summary(survfit(S ~ 1), times = timePoint)$lower
}



# validate on the new dataset using the model form training set and the multiple imputation result of the test dataset
# the calibration coefficient, estimated in the training set is also required
validateOnDecog <- function(fTraining, fTest, survValidation, miDataResultValidation){
  
  cIndex <- NULL
  linearPredictor <- NULL
  timePoint <- 60
  calibrationCoefficientTraining <- coef(fTest)
  
  for (j in 1:5){
    
    lpValidation <- predict(fTraining, newdata = complete(miDataResultValidation, j))
    rc <- rcorr.cens(-lpValidation, survValidation)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <- cbind(linearPredictor, lpValidation)
    
  }
  
  linearPredictorAverage <- rowMeans(linearPredictor)*calibrationCoefficientTraining
  baselineHazard <- basehaz(fTest, TRUE)
  predictedProbabilities <- predictSurvival(baselineHazard = baselineHazard,
                                            linearPredictor = linearPredictorAverage,
                                            timePoint = timePoint)
  calibrationData <- data.frame(predictedProbabilities, 
                                quintile = cut(predictedProbabilities, 
                                               quantile(predictedProbabilities, seq(0, 1, .2)),
                                               include.lowest = TRUE))
  levels(calibrationData$quintile) <- paste0("q", 1:5)
  
  calibrationData <- data.frame(calibrationData, 
                                time = survValidation[, "time"],
                                status = survValidation[, "status"])
  
  quintileResults <- calibrationData %>%
    group_by(quintile) %>%
    summarise(mean = mean(predictedProbabilities), 
              surv = getKaplanMeier(time,status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  
  calibrationPlot <- ggplot(data = quintileResults, aes(x = mean, y = surv)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper, width = 0)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2, size = 1, alpha = .4) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    xlab("Predicted probability") +
    ylab("Observed frequency") +
    theme_bw()
  
  discrimination <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))
  
  return(list(calibration = calibrationPlot,
              discrimination = discrimination))
  
}


# runSurvEstimation <- survest(f.miRec, newdata = dat, times = 60)
# test <- data.frame(surv = 1 - runSurvEstimation$surv)
# test <- cbind(dat, test) %>%
#   mutate(riskSubgroup = cut(surv, c(0, .25, .5 , .75, 1)))
# 
# kaplanMeier <- survfit(Surv(timeToRecurrence, status_RFS_FDA) ~ riskSubgroup, data = test)
# pp <- summary(kaplanMeier, times = 60)
# 
# data.frame(riskGroup = pp$strata,
#            P = 1 - pp$surv,
#            lower = 1 - pp$upper,
#            upper = 1 - pp$lower)
# 
# kaplanMeierDistant <- survfit(Surv(timeToDistant, status_DMFS_FDA) ~ riskSubgroup, data = test)
# pp <- summary(kaplanMeierDistant, times = 60)
# 
# data.frame(riskGroup = pp$strata,
#            P = 1 - pp$surv,
#            lower = 1 - pp$upper,
#            upper = 1 - pp$lower)
# 
# 
# kaplanMeierDeath <- survfit(Surv(timeToFollowup, status_OS_FDA) ~ riskSubgroup, data = test)
# pp <- summary(kaplanMeierDeath, times = 60)
# 
# data.frame(riskGroup = pp$strata,
#            P = 1 - pp$surv,
#            lower = 1 - pp$upper,
#            upper = 1 - pp$lower)