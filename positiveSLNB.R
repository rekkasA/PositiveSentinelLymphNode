library(foreign)
library(dplyr)
library(naniar)
library(mice)
library(survival)
library(rms)
library(mitools)
source('./helperFunctions.R')

spss2Date <- function(x) as.Date(x/86400, origin = "1582-10-14") # converts spss date to r date


##################################
# Read data
##################################

eortcData <- read.spss("./Data/EORTC.sav",
                       to.data.frame = TRUE)
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

selectedVariables <-c("id", "center", "age", "sex", "simpleloc", "clark",
                      "histology_simple", "histology_extended", "breslow",
                      "ulceration", "SLNBdate", "no_removed_SNs", "AJCC_sub_SLNB_8th",
                      "no_pos_SNs", "no_pos_SNs_cat", "SN_tumor_burden",
                      "SN_tumor_burden_extended",
                      "no_removed_nonSNs", "no_pos_nonSNs", "no_pos_nonSNs_cat",
                      "EJC_groups", "Risk_classes")

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

dat$timeToRecurrence60 <- ifelse(eortcData$timeToRecurrence < 60, 
                                 eortcData$timeToRecurrence,
                                 60)
dat$timeToDistant60 <- ifelse(eortcData$timeToDistant < 60,
                              eortcData$timeToDistant,
                              60)
dat$timeToFollowup60 <- ifelse(eortcData$timeToFollowup < 60,
                               eortcData$timeToFollowup,
                               60)


##################################
# Descriptives
##################################

dd <- datadist(dat)
options(datadist = 'dd')
options(digits = 8)


## Median follow-up time
S <- Surv(dat$timeToFollowup,
          dat$status_OS_FDA == 0)

survfit(S ~ 1,
        data = dat) 

summarySurvFit <- summary(survfit(S ~ 1,
                                  data = dat))
# plot(survfit(S ~ 1,
#              data = dat), conf.int  = F)
# abline(h=.5, col = "red", lty = 2)
# abline(v = 104.07, col = "red", lty = 2)


##################################
# Missingness
##################################

# colnames(dat)[!complete.cases(t(dat))] # length = 15
# gg_miss_var(dat, show_pct = TRUE)
# gg_miss_upset(dat, nsets = 14, nintersects = 500)

miVariables <- c("center", "age", "sex", "simpleloc", "SN_tumor_burden_extended",
                 "histology_simple","breslow", "AJCC_sub_SLNB_8th", 
                 "no_pos_SNs", "no_removed_SNs",
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
# miDataResult <- mice(data = miData,
#                      m = 5)

# saveRDS(miDataResult, file = "miDataResult.rds")

miDataResult <- readRDS("miDataResult.rds")
##################################
# Modeling recurrence-free survival
##################################

S <- Surv(dat$timeToRecurrence,
          dat$status_RFS_FDA)
# plot(survfit(S ~ as.factor(center),
#              data = dat))
# 
# plot(survfit(S ~ 1,
#              data = dat))

## All variables
### No time restrictions
form <- S ~ sex + ulceration + no_removed_SNs + 
  no_pos_SNs + simpleloc + histology_simple  +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
anova(f.mi)

mi.cindex(f.mi, miDataResult)

### Model specification
fastbw(f.mi, rule = "p", sls = .1)

dat.complete <- complete(miDataResult, 1)
describe(dat.complete$breslow)
describe(dat.complete$SN_tumor_burden)

f <- cph(S ~ rcs(breslow, 4),
         data = dat)
f2 <- cph(S ~ log(breslow),
          data = dat)
f4 <- cph(S ~ breslow,
          data = dat)
plot(1, 
     type="n", 
     xlab="breslow", 
     ylab="", 
     xlim=c(0, 15),
     ylim=c(-1, 1))
lines(Predict(f, breslow))
lines(Predict(f2, breslow), col = 2)
lines(Predict(f4, breslow), col = "green")


f <- cph(S ~ rcs(SN_tumor_burden, 4),
         data = dat)
f2 <- cph(S ~ log(SN_tumor_burden),
          data = dat)
plot(1, 
     type="n", 
     xlab="SN_tumor_burden", 
     ylab="", 
     xlim=c(0, 14),
     ylim=c(-2, 2))
lines(Predict(f, SN_tumor_burden))
lines(Predict(f2, SN_tumor_burden), col = 2)

f <- cph(S ~ rcs(age, 4),
         data = dat)
f2 <- cph(S ~ log(age),
          data = dat)
f4 <- cph(S ~ age,
          data = dat)
plot(1, 
     type="n", 
     xlab="SN_tumor_burden", 
     ylab="", 
     xlim=c(10, 80),
     ylim=c(-1, 1))
lines(Predict(f, age))
lines(Predict(f2, age), col = 2)
lines(Predict(f4, age), col = "green")

### Final model 
form <- S ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
f.mi
mi.cindex(f.mi,
          miDataResult)

breslow.class <- c(seq(.2, .6, .2), seq(1, 3, 1), seq(4, 10, 2))
t.class <- c(seq(.2, 1, .2), 1.5, seq(2, 4, .5), 7)
age.class <- c(seq(10, 18, 2), 20, 25, seq(30, 70, 10))
nom <- nomogram(f.mi,
                lp = TRUE,
                maxscale = 10,
                age = age.class,
                SN_tumor_burden = t.class,
                breslow = breslow.class)


plot(nom,
     total.sep.page = TRUE,
     col.grid =gray(c(.7, .9)))

### Results for 5 years
S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)

# plot(survfit(S ~ as.factor(center),
#              data = dat))
# 
# plot(survfit(S ~ 1,
#              data = dat))

form <- S ~ sex + ulceration + no_removed_SNs +
  simpleloc + histology_simple + clark + AJCC_sub_SLNB_8th +
  rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 4) # I get convergence issues otherwise -> why??

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
mi.cindex(f.mi,
          miDataResult)

anova(f.mi)
fastbw(f.mi, rule = "p", sls = .1)

## Final 5-year recurrence free survival model
S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)

form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
# form <- S ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)
summaryFit <- summary(f.mi)
namesSummaryFit <- rownames(summaryFit)
namesSummaryFit <- namesSummaryFit[seq(1, 15, 2)]
summaryFit <- summaryFit[seq(1, 16, 2), c(4, 6, 7)]
rownames(summaryFit) <- namesSummaryFit
colnames(summaryFit)[1] <- "Hazard ratio"

round(mi.cindex(f.mi, miDataResult)[1:4], 2)


## Nomogram for 5-year free survival using final model
breslow.class <- c((2:10)/10, 1 + (1:5)/5, 2.5, 3:5, 6, 8, 10)
t.class <- c((2:10)/10, 1+(1:5)/5, 2.5, 3, 4, 6, 8)
age.class <- c(seq(10, 18, 2), 20, 25, seq(30, 70, 10))
nom <- nomogram(f.mi,
                lp = TRUE,
                maxscale = 10,
                age = age.class,
                SN_tumor_burden = t.class,
                breslow = breslow.class)


plot(nom,
     total.sep.page = TRUE,
     col.grid =gray(c(.7, .9)))

## Bootstrap validation
# f.mi.boot <- f.mi
# opt <- NULL
# slope <- NULL
# intercept <- NULL
# set.seed(0)
# for (i in 1:5){
#   v<-validate(f.mi.boot$fits[[i]], 
#               bw = TRUE, 
#               rule = 'p',
#               sls = 0.1,
#               B = 100,
#               pr = FALSE,
#               type = 'individual')
#   
#   opt <- c(opt, v["Dxy", "optimism"]/2)
#   slope <- c(slope, v["Slope", "test"])
# }
# opt
# mean(opt)
# mean(slope)

## Cross-validation across centers

S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)

form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence,
                      datTestCenter$status_RFS_FDA)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResult,
                                      data = miData,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE,
                                      subset = miData$center != centers[i])
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(f.miTrainCenters, TRUE)
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:5){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResult, j)[miData$center == centers[i], ])
    rc <- rcorr.cens(-lpTestCenter, STestCenter)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  discrimination[[i]] <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))
  
  linearPredictorAverage <- rowMeans(linearPredictor)
  baselineHazard <- basehaz(f.miTrainCenters)
  predictedProbabilities <- predictSurvival(baselineHazard = f.basehazTrainCenters,
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
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}


## Risk distribution
dat.complete <- complete(miDataResult, 1)
breslowTrunc <- dat.complete$breslow
breslowTrunc[dat.complete$breslow < .2] <- .2
breslowTrunc[dat.complete$breslow > 10] <- 10

SN_tumor_burdenTrunc <- dat.complete$SN_tumor_burden
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden < .2] <- .2
SN_tumor_burdenTrunc[dat.complete$SN_tumor_burden > 7] <- 7

ageTrunc <- dat.complete$age
ageTrunc[dat.complete$age < 10] <- 10
ageTrunc[dat.complete$age > 70] <- 70


### Distant metastasis
SDistant <- Surv(dat$timeToDistant60,
                 dat$status_DMFS_FDA60)

# form <- SDistant ~ sex + ulceration +
#   simpleloc + histology_simple + clark + AJCC_sub_SLNB_8th +
#   rcs(breslow, 3) + rcs(SN_tumor_burden, 3) + rcs(age, 3)

form <- SDistant ~ sex + ulceration + log(age) + log(SN_tumor_burden) + log(breslow)

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)

round(mi.cindex(f.mi, miDataResult), 2)

# f.mi.boot <- f.mi
# opt <- NULL
# slope <- NULL
# intercept <- NULL
# set.seed(0)
# for (i in 1:5){
#   v<-validate(f.mi.boot$fits[[i]], 
#               bw = TRUE, 
#               rule = 'p',
#               sls = 0.1,
#               B = 500,
#               pr = FALSE,
#               type = 'individual')
#   
#   opt <- c(opt, v["Dxy", "optimism"]/2)
#   slope <- c(slope, v["Slope", "test"])
# }
# opt
# mean(opt)
# mean(slope)



SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)

form <- SRecurrence ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)

f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResult,
                           data = miData,
                           n.impute = 5,
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


cindex <- NULL
m <- miDataResult$m
for (i in 1:m){
  
  rc <- Hmisc::rcorr.cens(-f.miRec$fits[[i]]$linear.predictors*coef(fDistant),
                          SDistant)
  cindex<-rbind(cindex,
                c(rc["C Index"], rc["S.D."]/2))
  
}
summary(mitools::MIcombine(as.list(cindex[,1]),as.list(cindex[,2]^2)))



## Cross-validation across centers for distant metastasis
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)

form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)

discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  
  formTrainCenters <- update(form, STrainCenters  ~ . )
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResult,
                                      data = miData,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE,
                                      subset = miData$center != centers[i])
  
  SDistantTrainCenter <- Surv(datTrainCenters$timeToDistant60,
                              datTrainCenters$status_DMFS_FDA60)
  
  SDistantTestCenter <- Surv(datTestCenter$timeToDistant60,
                             datTestCenter$status_DMFS_FDA60)
  lp <- f.miTrainCenters$linear
  fDistant <- cph(SDistantTrainCenter ~ lp, 
                  y = TRUE,
                  x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(f.miTrainCenters, TRUE)
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:5){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResult, j)[miData$center == centers[i], ])
    lpTestCenter <- lpTestCenter
    rc <- rcorr.cens(-lpTestCenter, SDistantTestCenter)
    cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
    
    linearPredictor <-cbind(linearPredictor, lpTestCenter)
    
  }
  
  discrimination[[i]] <- summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))

  
  
  linearPredictorAverage <- rowMeans(linearPredictor)*coef(fDistant)
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
              surv = getKaplanMeier(time, status, timePoint),
              lower =  lowerKaplanMeier(time, status, timePoint),
              upper = upperKaplanMeier(time, status, timePoint))
  plot(quintileResults$mean, 
       quintileResults$surv,
       ylim = c(0, 1), 
       xlim = c(0, 1),
       main = as.character(centers[i]), 
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}

discriminationData <- data.frame(do.call(rbind, discrimination))
discriminationData <- discriminationData[, -5]
colMeans(discriminationData)

### Calibration for death
SDeath <- Surv(dat$timeToFollowup60, 
               dat$status_OS_FDA60)
lp <- f.miRec$linear
fDeath <- cph(SDeath ~ lp,
              y = TRUE, 
              x = TRUE)

### Calibration across centers - Death
SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)

form <- SRecurrence ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)

discrimination <- list()
centers <- unique(dat$center)
par(mfrow=c(3,3))
for(i in 1:length(centers)){
  
  datTrainCenters <- subset(dat, center != centers[i])
  STrainCenters <- Surv(datTrainCenters$timeToRecurrence60,
                        datTrainCenters$status_RFS_FDA60)
  datTestCenter <- subset(dat, center == centers[i])
  STestCenter <- Surv(datTestCenter$timeToRecurrence60,
                      datTestCenter$status_RFS_FDA60)
  formTrainCenters <- update(form, STrainCenters  ~ . )
  f.miTrainCenters <- fit.mult.impute(formTrainCenters,
                                      cph,
                                      xtrans = miDataResult,
                                      data = miData,
                                      n.impute = 5,
                                      pr = FALSE,
                                      fit.reps = TRUE,
                                      y = TRUE,
                                      x = TRUE,
                                      se.fit = TRUE,
                                      subset = miData$center != centers[i])
  
  SDeathTrainCenter <- Surv(datTrainCenters$timeToFollowup60,
                            datTrainCenters$status_OS_FDA60)
  lp <- f.miTrainCenters$linear
  fDeath <- cph(SDeathTrainCenter ~ lp, 
                y = TRUE,
                x = TRUE)
  
  lpTrainCenters <- f.miTrainCenters$linear.predictors
  f.basehazTrainCenters <- basehaz(f.miTrainCenters, TRUE)
  timePoint <- 60
  cIndex <- NULL
  linearPredictor <- NULL
  
  for (j in 1:5){
    
    lpTestCenter <- predict(f.miTrainCenters,
                            newdata = complete(miDataResult, j)[miData$center == centers[i], ])
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
       pch = 19)
  arrows(x0 = quintileResults$mean, 
         y0 = quintileResults$lower, 
         x1 = quintileResults$mean, 
         y1 = quintileResults$upper, 
         length=0, angle=90, code=3)
  abline(0, 1, lty = 2)
  
}

### Plot across all outcomes
baseHazRec <- basehaz(f.miRec)
baseHazDistant <- basehaz(fDistant)
baseHazDeath <- basehaz(fDeath)
baseSurv <- estimtateBaselineSurvival(baseHazRec, 60)
baseSurvDistant <- estimtateBaselineSurvival(baseHazDistant, 60)
baseSurvDeath <- estimtateBaselineSurvival(baseHazDeath, 60)


interceptNom <- attr(nom, "info")$Intercept
pointwiseIncrease <- 1/attr(nom, "info")$sc
points <- 0:27

lp.score<-predict(f.miRec,
                  newdata=data.frame(ulceration=dat.complete$ulceration,
                                     age=ageTrunc,
                                     SN_tumor_burden=SN_tumor_burdenTrunc,
                                     breslow = breslowTrunc))
score<-round((lp.score - interceptNom)/pointwiseIncrease, 0)
h<-hist(score,plot=FALSE,breaks=0:28,right=FALSE)

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
pointsToProbs$riskLevel <- c(rep(1, 16),
                             rep(2, 4),
                             rep(3, 8))

pointsToProbs$riskLevel <- factor(pointsToProbs$riskLevel, labels = c("Low risk", "Medium risk", "High risk"))

ggplot(data = pointsToProbs, aes(x = points)) +
  geom_bar(mapping = aes(x = points, y = hist, fill = factor(riskLevel)), stat = "identity") +
  scale_fill_brewer(palette = "Blues") +
  ylab("% of patients")+
  geom_line(aes(x = points, y = recProbs/5, col = "Recurrence"), size = 1.2) +
  geom_line(aes(x = points, y = distantProbs/5, col = "Distant metastasis"), size = 1.2) +
  geom_line(aes(x = points, y = deathProbs/5, col = "Death"), size = 1.2) +
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "5-year probability (%)")) +
  xlab("Risk score") +
  theme_bw() +
  theme(legend.position = c(0.1, 0.8),
        legend.title = element_blank())



##################################
# EXTERNAL VALIDATION
##################################

decogData <- read.spss("./Data/DECOG.sav",
                       to.data.frame = TRUE)
decogData$id <- 1:dim(decogData)[1] 


selectedVariables <-c("id", "age", "sex", "simpleloc",
                      "histology_simple", "breslow",
                      "ulceration", "SLNBdate", "no_removed_SNs", "AJCCsub_8th_SLNB_only_DV",
                      "no_pos_nonSNs", "no_pos_SNs", "no_pos_SNs_cat", "no_pos_nonSNs_cat",
                      "SN_tumor_burden_extended", "SN_tumor_burden_simple",
                      "no_removed_nonSNs", "no_pos_nonSNs", "no_pos_nonSNs_cat",
                      "EJC_groups", "Risk_classes")

outcomes <- c("status_RFS_FDA", "time_RFS", 
              "status_DMFS_FDA", "time_DMFS",
              "status_OS_FDA", "time_OS",
              "status_RFS_FDA60", "time_RFS60", 
              "status_DMFS_FDA60", "time_DMFS60",
              "status_OS_FDA60", "time_OS60")

datValidation <- decogData %>%
  select(c(selectedVariables, outcomes))


# colnames(datValidation)[!complete.cases(t(datValidation))] # length = 15
# gg_miss_var(datValidation, show_pct = TRUE)
# gg_miss_upset(datValidation, nsets = 14, nintersects = 500)


miVariables <- c( "age", "sex", "simpleloc",
                  "histology_simple", "breslow",
                  "ulceration", "SLNBdate",
                  "no_pos_SNs_cat",
                  "SN_tumor_burden_extended",
                  "no_pos_nonSNs_cat", "status_RFS_FDA", "time_RFS", 
                  "status_DMFS_FDA", "time_DMFS",
                  "status_OS_FDA", "time_OS")

miData <- datValidation[miVariables] 

miData$time_RFS <- log(miData$time_RFS)
miData$time_DMFS <- log(miData$time_DMFS)
miData$time_OS <- log(miData$time_O)

# miDataResultVal <- mice(data = miData,
#                         m = 5)

# saveRDS(miDataResultVal, file = "miDataResultValidation.rds")

miDataResultVal <- readRDS("miDataResultValidation.rds")

S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)
form <- S ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)
f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResult,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)

survValidation <- Surv(datValidation$time_RFS60,
                       datValidation$status_RFS_FDA60)

cIndex <- NULL
linearPredictor <- NULL
timePoint <- 60

for (j in 1:5){
  
  lpValidation <- predict(f.mi, newdata = complete(miDataResultVal, j))
  rc <- rcorr.cens(-lpValidation, survValidation)
  cIndex <- rbind(cIndex,c(rc["C Index"],rc["S.D."]/2))
  
  linearPredictor <-cbind(linearPredictor, lpValidation)
  
}

summary(MIcombine(as.list(cIndex[,1]),as.list(cIndex[,2]^2)))


linearPredictorAverage <- rowMeans(linearPredictor)
baselineHazard <- basehaz(f.mi)
predictedProbabilities <- predictSurvival(baselineHazard = baselineHazard,
                                          linearPredictor = linearPredictorAverage,
                                          timePoint = timePoint)
calibrationData <- data.frame(predictedProbabilities, 
                              quintile = cut(predictedProbabilities, 
                                             quantile(predictedProbabilities, seq(0, 1, .2)),
                                             include.lowest = TRUE))
levels(calibrationData$quintile) <- paste0("q", 1:5)

calibrationData <- data.frame(calibrationData, 
                              time = datValidation$time_RFS60,
                              status = datValidation$status_RFS_FDA60)

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
     main = "External validation", 
     pch = 19)
arrows(x0 = quintileResults$mean, 
       y0 = quintileResults$lower, 
       x1 = quintileResults$mean, 
       y1 = quintileResults$upper, 
       length=0, angle=90, code=3)
abline(0, 1, lty = 2)


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


## Distant metastasis

SRecurrence  <- Surv(dat$timeToRecurrence60,
                     dat$status_RFS_FDA60)

form <- SRecurrence ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)

f.miRec <- fit.mult.impute(form,
                           cph,
                           xtrans = miDataResult,
                           data = miData,
                           n.impute = 5,
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

survValidation <- Surv(datValidation$time_DMFS60,
                       datValidation$status_DMFS_FDA60)

pp = validateOnDecog(fTraining = f.miRec, 
                     calibrationCoefficientTraining = coef(fDistant), 
                     survValidation = survValidation, 
                     miDataResultValidation = miDataResultVal)

## Death

SDeath <- Surv(dat$timeToFollowup60, 
               dat$status_OS_FDA60)
fDeath <- cph(SDeath ~ lp,
              y = TRUE, 
              x = TRUE)

survValidation <- Surv(datValidation$time_OS60,
                       datValidation$status_OS_FDA60)

pp = validateOnDecog(fTraining = f.miRec, 
                     calibrationCoefficientTraining = coef(fDeath), 
                     survValidation = survValidation, 
                     miDataResultValidation = miDataResultVal)


##################################
# Tertiary objective
##################################
eortcDataCLND <- subset(eortcData, CLND == "yes")

selectedVariables <-c("id", "center", "age", "sex", "simpleloc", "clark",
                      "histology_simple", "histology_extended", "breslow",
                      "ulceration", "SLNBdate", "no_removed_SNs", "AJCC_sub_SLNB_8th",
                      "no_pos_SNs_cat", "SN_tumor_burden",
                      "SN_tumor_burden_extended",
                      "no_removed_nonSNs", "no_pos_nonSNs", "no_pos_nonSNs_cat",
                      "EJC_groups", "Risk_classes")

outcomes <- c("status_RFS_FDA", "timeToRecurrence", 
              "status_DMFS_FDA", "timeToDistant",
              "status_OS_FDA", "timeToFollowup",
              "status_RFS_FDA60", "status_DMFS_FDA60",
              "status_OS_FDA60")

statusOutcomesoutcomes <- c("status_RFS_FDA", 
                            "status_DMFS_FDA", 
                            "status_OS_FDA")

dat <- eortcDataCLND %>%
  select(c(selectedVariables, outcomes))

dat$timeToRecurrence60 <- ifelse(eortcDataCLND$timeToRecurrence < 60, 
                                 eortcDataCLND$timeToRecurrence,
                                 60)
dat$timeToDistant60 <- ifelse(eortcDataCLND$timeToDistant < 60,
                              eortcDataCLND$timeToDistant,
                              60)
dat$timeToFollowup60 <- ifelse(eortcDataCLND$timeToFollowup < 60,
                               eortcDataCLND$timeToFollowup,
                               60)

miVariables <- c("center", "age", "sex", "simpleloc", "clark", "SN_tumor_burden_extended",
                 "histology_simple","breslow", "AJCC_sub_SLNB_8th",
                 "ulceration", "SLNBdate", "no_removed_SNs",
                 "SN_tumor_burden", "no_removed_nonSNs", "no_pos_nonSNs_cat",
                 "status_RFS_FDA", "timeToRecurrence", 
                 "status_DMFS_FDA", "timeToDistant",
                 "status_OS_FDA", "timeToFollowup")


miData <- dat[miVariables] 

miData$timeToRecurrence <- log(miData$timeToRecurrence)
miData$timeToDistant <- log(miData$timeToDistant)
miData$timeToFollowup <- log(miData$timeToFollowup)

## multiple imputations
# miDataResultCLND <- mice(data = miData,
#                          m = 5)
miDataResultCLND <- readRDS("miDataResultCLND.rds")

S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)

form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow) + no_pos_nonSNs_cat
# form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
# form <- S ~ ulceration + log(age) + SN_tumor_burden_extended + log(breslow)

dd <- datadist(dat)
options(datadist = 'dd')
options(digits = 8)

f.mi <- fit.mult.impute(form,
                        cph,
                        xtrans = miDataResultCLND,
                        data = miData,
                        n.impute = 5,
                        pr = FALSE,
                        fit.reps = TRUE,
                        y = TRUE,
                        x = TRUE,
                        se.fit = TRUE)

round(mi.cindex(f.mi, miDataResultCLND)[1:4], 2)

### Additional value of additional positive lymph modes
lp <- f.mi$linear


no_pos_nonSNs_cat = dat$no_pos_nonSNs_cat
dd <- datadist(no_pos_nonSNs_cat)
options(datadist = 'dd')
f<-cph(S ~ offset(lp) + no_pos_nonSNs_cat)
f

summary(f, no_pos_nonSNs_cat = 0)