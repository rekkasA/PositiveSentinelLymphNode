---
title: "Reviewer comments"
output: 
  bookdown::html_document2:
    keep_md: true
    code_folding: hide
---





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


# Simple risk subgroups




```r
kable(resSimpleRec %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year recurrence probability with EJC subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-4)5-year recurrence probability with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.44 </td>
   <td style="text-align:right;"> 0.58 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> 0.41 </td>
   <td style="text-align:right;"> 0.55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.73 </td>
   <td style="text-align:right;"> 0.67 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
</tbody>
</table>


```r
kable(resSimpleDistant %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year distant metastasis probability with EJC subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-5)5-year distant metastasis probability with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.48 </td>
   <td style="text-align:right;"> 0.40 </td>
   <td style="text-align:right;"> 0.54 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.34 </td>
   <td style="text-align:right;"> 0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.69 </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.74 </td>
  </tr>
</tbody>
</table>


```r
kable(resSimpleDeath %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)),
      caption = "5-year overall mortality probability with EJC subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-6)5-year overall mortality probability with EJC subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> EJC_groups </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> No ulceration + low TB </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.16 </td>
   <td style="text-align:right;"> 0.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> No ulceration + high TB </td>
   <td style="text-align:right;"> 0.41 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + low TB </td>
   <td style="text-align:right;"> 0.35 </td>
   <td style="text-align:right;"> 0.28 </td>
   <td style="text-align:right;"> 0.42 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ulceration + high TB </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.53 </td>
   <td style="text-align:right;"> 0.66 </td>
  </tr>
</tbody>
</table>

# AJCC classifications

```r
eortcData <- read.spss("./Data/EORTC_CLEAN.sav",
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
```
## AJCC 7th subgroups

```r
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

cIndexData <- data.frame(outcome = c("Recurrence", "Distant meatastasis", "Overall mortality"),
                         cIndex = c(cIndexAJCC_7th_Rec[1], cIndexAJCC_7th_Distant[1], cIndexAJCC_7th_Death[1]),
                         lower = c(cIndexAJCC_7th_Rec[1] - qnorm(.975)*cIndexAJCC_7th_Rec[2],
                                   cIndexAJCC_7th_Distant[1] - qnorm(.975)*cIndexAJCC_7th_Distant[2],
                                   cIndexAJCC_7th_Death[1] - qnorm(.975)*cIndexAJCC_7th_Death[2]),
                         upper = c(cIndexAJCC_7th_Rec[1] + qnorm(.975)*cIndexAJCC_7th_Rec[2],
                                   cIndexAJCC_7th_Distant[1] + qnorm(.975)*cIndexAJCC_7th_Distant[2],
                                   cIndexAJCC_7th_Death[1] + qnorm(.975)*cIndexAJCC_7th_Death[2])) %>%
  mutate(cIndex = round(cIndex, 2),
         lower = round(lower, 2),
         upper = round(upper, 2))
```


```r
kable(cIndexData, 
      caption = "c-indices for AJCC 7th classification") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-9)c-indices for AJCC 7th classification</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> outcome </th>
   <th style="text-align:right;"> cIndex </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Recurrence </td>
   <td style="text-align:right;"> 0.61 </td>
   <td style="text-align:right;"> 0.59 </td>
   <td style="text-align:right;"> 0.63 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Distant meatastasis </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Overall mortality </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.59 </td>
   <td style="text-align:right;"> 0.65 </td>
  </tr>
</tbody>
</table>


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
                              probability = 1 - runSurvEstimation$surv, 
                              lower = 1 - runSurvEstimation$upper, 
                              upper = 1 - runSurvEstimation$lower) %>%
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
                                  probability = 1 - runSurvEstimation$surv, 
                                  lower = 1 - runSurvEstimation$upper, 
                                  upper = 1 - runSurvEstimation$lower) %>%
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
                                probability = 1 - runSurvEstimation$surv, 
                                lower = 1 - runSurvEstimation$upper, 
                                upper = 1 - runSurvEstimation$lower) %>%
  arrange(AJCC_7th_extra)
```


```r
kable(resAJCC_7th_Rec %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year recurrence probability with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-11)5-year recurrence probability with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.57 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;"> 0.67 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.84 </td>
  </tr>
</tbody>
</table>


```r
kable(resAJCC_7th_Distant %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year distant metastasis probability with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-12)5-year distant metastasis probability with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.20 </td>
   <td style="text-align:right;"> 0.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.46 </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.57 </td>
   <td style="text-align:right;"> 0.52 </td>
   <td style="text-align:right;"> 0.62 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.85 </td>
  </tr>
</tbody>
</table>


```r
kable(resAJCC_7th_Death %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year overall mortality probability with AJCC-7th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-13)5-year overall mortality probability with AJCC-7th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_7th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0 </td>
   <td style="text-align:right;"> 0.20 </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0 </td>
   <td style="text-align:right;"> 0.40 </td>
   <td style="text-align:right;"> 0.32 </td>
   <td style="text-align:right;"> 0.46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> 0.44 </td>
   <td style="text-align:right;"> 0.53 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.86 </td>
  </tr>
</tbody>
</table>

## AJCC 8th subgroups


```r
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
cIndexData <- data.frame(outcome = c("Recurrence", "Distant meatastasis", "Overall mortality"),
                         cIndex = c(cIndexAJCC_8th_Rec[1], cIndexAJCC_8th_Distant[1], cIndexAJCC_8th_Death[1]),
                         lower = c(cIndexAJCC_8th_Rec[1] - qnorm(.975)*cIndexAJCC_8th_Rec[2],
                                   cIndexAJCC_8th_Distant[1] - qnorm(.975)*cIndexAJCC_8th_Distant[2],
                                   cIndexAJCC_8th_Death[1] - qnorm(.975)*cIndexAJCC_8th_Death[2]),
                         upper = c(cIndexAJCC_8th_Rec[1] + qnorm(.975)*cIndexAJCC_8th_Rec[2],
                                   cIndexAJCC_8th_Distant[1] + qnorm(.975)*cIndexAJCC_8th_Distant[2],
                                   cIndexAJCC_8th_Death[1] + qnorm(.975)*cIndexAJCC_8th_Death[2])) %>%
  mutate(cIndex = round(cIndex, 2),
         lower = round(lower, 2),
         upper = round(upper, 2))
```


```r
kable(cIndexData, 
      caption = "c-indices for AJCC 8th classification") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-15)c-indices for AJCC 8th classification</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> outcome </th>
   <th style="text-align:right;"> cIndex </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Recurrence </td>
   <td style="text-align:right;"> 0.62 </td>
   <td style="text-align:right;"> 0.59 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Distant meatastasis </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.65 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Overall mortality </td>
   <td style="text-align:right;"> 0.63 </td>
   <td style="text-align:right;"> 0.61 </td>
   <td style="text-align:right;"> 0.66 </td>
  </tr>
</tbody>
</table>


```r
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
                              probability = 1 - runSurvEstimation$surv, 
                              lower = 1 - runSurvEstimation$upper, 
                              upper = 1 - runSurvEstimation$lower) %>%
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
                                  probability = 1 - runSurvEstimation$surv, 
                                  lower = 1 - runSurvEstimation$upper, 
                                  upper = 1 - runSurvEstimation$lower) %>%
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
                                probability = 1 - runSurvEstimation$surv, 
                                lower = 1 - runSurvEstimation$upper, 
                                upper = 1 - runSurvEstimation$lower) %>%
  arrange(AJCC_8th_extra)
```


```r
kable(resAJCC_8th_Rec %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year recurrence probability with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-17)5-year recurrence probability with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.27 </td>
   <td style="text-align:right;"> 0.20 </td>
   <td style="text-align:right;"> 0.34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.37 </td>
   <td style="text-align:right;"> 0.23 </td>
   <td style="text-align:right;"> 0.49 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.36 </td>
   <td style="text-align:right;"> 0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.64 </td>
   <td style="text-align:right;"> 0.59 </td>
   <td style="text-align:right;"> 0.68 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.66 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.90 </td>
  </tr>
</tbody>
</table>


```r
kable(resAJCC_8th_Distant %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year distant metastasis probability with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-18)5-year distant metastasis probability with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.34 </td>
   <td style="text-align:right;"> 0.20 </td>
   <td style="text-align:right;"> 0.46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.35 </td>
   <td style="text-align:right;"> 0.29 </td>
   <td style="text-align:right;"> 0.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.58 </td>
   <td style="text-align:right;"> 0.53 </td>
   <td style="text-align:right;"> 0.63 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.68 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.91 </td>
  </tr>
</tbody>
</table>


```r
kable(resAJCC_8th_Death %>%
        mutate(probability = round(probability, 2),
               lower = round(lower, 2),
               upper = round(upper, 2)), 
      caption = "5-year overall mortality probability with AJCC-8th subgroups") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-19)5-year overall mortality probability with AJCC-8th subgroups</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> AJCC_8th_extra </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> IIIA =&lt;1.0mm </td>
   <td style="text-align:right;"> 0.15 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 0.21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIA &gt;1.0mm </td>
   <td style="text-align:right;"> 0.27 </td>
   <td style="text-align:right;"> 0.14 </td>
   <td style="text-align:right;"> 0.38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIB </td>
   <td style="text-align:right;"> 0.30 </td>
   <td style="text-align:right;"> 0.24 </td>
   <td style="text-align:right;"> 0.36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIIC </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 0.45 </td>
   <td style="text-align:right;"> 0.55 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IIID </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.92 </td>
  </tr>
</tbody>
</table>

# Overall probabilities


```r
spss2Date <- function(x) as.Date(x/86400, origin = "1582-10-14") # converts spss date to r date
tableCount <- 0
figureCount <- 0
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
selectedVariables <-c("id", "center", "age", "sex", "simpleloc", "clark", "CLND",
                      "histology_simple", "histology_extended", "breslow",
                      "ulceration", "SLNBdate", "no_removed_SNs", "AJCC_sub_SLNB_8th",
                      "no_pos_SNs", "no_pos_SNs_cat", "SN_tumor_burden",
                      "SN_tumor_burden_extended", "SN_tumor_burden_simple",
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

S <- Surv(dat$timeToRecurrence60,
          dat$status_RFS_FDA60)
form <- S ~ ulceration + log(age) + log(SN_tumor_burden) + log(breslow)
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

dat.complete <- complete(miDataResult, 1)
runSurvEstimation <- survest(f.miRec, newdata = dat.complete, times = 60)
test <- data.frame(surv = 1 - runSurvEstimation$surv)
test <- cbind(dat, test) %>%
  mutate(riskSubgroup = cut(surv, c(0, .25, .5 , .75, 1)))
```


```r
kaplanMeier <- survfit(Surv(timeToRecurrence, status_RFS_FDA) ~ riskSubgroup, data = test)
pp <- summary(kaplanMeier, times = 60)

overallProbabilities <- data.frame(riskGroup = pp$strata,
                                   probability = 1 - pp$surv,
                                   lower = 1 - pp$upper,
                                   upper = 1 - pp$lower) %>%
  mutate(riskGroup = factor(riskGroup, labels = c("(0,0.25]", "(0.25,0.50]", "(0.50,0.75]", "(0.75, 1]")),
         probability = round(probability, 2),
         lower = round(lower, 2),
         upper = round(upper, 2))

kable(overallProbabilities, 
      caption = "Observed 5-year recurrence rates within subgroups risk using our model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-21)Observed 5-year recurrence rates within subgroups risk using our model</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> riskGroup </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (0,0.25] </td>
   <td style="text-align:right;"> 0.13 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> 0.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.25,0.50] </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:right;"> 0.43 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.50,0.75] </td>
   <td style="text-align:right;"> 0.61 </td>
   <td style="text-align:right;"> 0.56 </td>
   <td style="text-align:right;"> 0.66 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.75, 1] </td>
   <td style="text-align:right;"> 0.82 </td>
   <td style="text-align:right;"> 0.73 </td>
   <td style="text-align:right;"> 0.88 </td>
  </tr>
</tbody>
</table>


```r
kaplanMeierDistant <- survfit(Surv(timeToDistant, status_DMFS_FDA) ~ riskSubgroup, data = test)
pp <- summary(kaplanMeierDistant, times = 60)

overallProbabilities <- data.frame(riskGroup = pp$strata,
                                   probability = 1 - pp$surv,
                                   lower = 1 - pp$upper,
                                   upper = 1 - pp$lower) %>%
  mutate(riskGroup = factor(riskGroup, labels = c("(0,0.25]", "(0.25,0.50]", "(0.50,0.75]", "(0.75, 1]")),
         probability = round(probability, 2),
         lower = round(lower, 2),
         upper = round(upper, 2))

kable(overallProbabilities, 
      caption = "Observed 5-year distant metastasis rates within risk subgroups using our model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-22)Observed 5-year distant metastasis rates within risk subgroups using our model</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> riskGroup </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (0,0.25] </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.25,0.50] </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:right;"> 0.26 </td>
   <td style="text-align:right;"> 0.36 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.50,0.75] </td>
   <td style="text-align:right;"> 0.55 </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> 0.60 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.75, 1] </td>
   <td style="text-align:right;"> 0.78 </td>
   <td style="text-align:right;"> 0.69 </td>
   <td style="text-align:right;"> 0.84 </td>
  </tr>
</tbody>
</table>


```r
kaplanMeierDeath <- survfit(Surv(timeToFollowup, status_OS_FDA) ~ riskSubgroup, data = test)
pp <- summary(kaplanMeierDeath, times = 60)

overallProbabilities <- data.frame(riskGroup = pp$strata,
                                   probability = 1 - pp$surv,
                                   lower = 1 - pp$upper,
                                   upper = 1 - pp$lower) %>%
  mutate(riskGroup = factor(riskGroup, labels = c("(0,0.25]", "(0.25,0.50]", "(0.50,0.75]", "(0.75, 1]")),
         probability = round(probability, 2),
         lower = round(lower, 2),
         upper = round(upper, 2))

kable(overallProbabilities, 
      caption = "Observed 5-year overall mortality rates within risk subgroups using our model") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
```

<table class="table table-striped table-hover table-condensed table-responsive" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:unnamed-chunk-23)Observed 5-year overall mortality rates within risk subgroups using our model</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> riskGroup </th>
   <th style="text-align:right;"> probability </th>
   <th style="text-align:right;"> lower </th>
   <th style="text-align:right;"> upper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (0,0.25] </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 0.13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.25,0.50] </td>
   <td style="text-align:right;"> 0.25 </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.30 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.50,0.75] </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.54 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> (0.75, 1] </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.61 </td>
   <td style="text-align:right;"> 0.77 </td>
  </tr>
</tbody>
</table>
