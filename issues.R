##################################
# ISSUES
##################################

## last follow-up before time of recurrence
eortcData %>%
  filter(first_rec > last_FU) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

## Last follow-up before date of distant metastasis
eortcData %>%
  filter(Date_dist_met > last_FU) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

## Last follow-up is the same as date of recurrence
eortcData %>%
  filter(first_rec == last_FU) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

## Last follow-up is the same as date of distant metastasis
eortcData %>%
  filter(Date_dist_met == last_FU) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

## Are these dead or alive?
eortcData %>%
  filter(status == "unknown") %>%
  select(c(id, outcomes, status))

## Has this had a recurrence or not?
eortcData %>%
  filter(status_RFS_FDA == 0 & recurrence == 1) %>%
  select(c(id, outcomes, status))


##################################
# Distant metastasis definitions
##################################

# The difference are the patients with last follow-up same as date of distant metastasis
# plus patient 1042 who had last follow-up before date to distant metastasis

eortcData %>%
  filter(status_DMFS_FDA == 1 & !is.na(Date_dist_met)) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

eortcData %>%
  filter(status_DMFS_FDA == 1 & timeToDistant < timeToFollowup) %>%
  select(id, recurrence, timeToRecurrence, status_RFS_FDA, 
         time_RFS, status_OS_FDA, time_OS, timeToFollowup, 
         first_rec, Date_dist_met, last_FU, status)

# eortcData <- eortcData %>%
#   mutate(distantMetastasis = ifelse(status_DMFS_FDA == 1 & !is.na(Date_dist_met), 1, 0))