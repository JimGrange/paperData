rm(list = ls()) 
setwd("C:/Users/Home/Dropbox/Reliability Manuscript Revision/Exploratory analysis codes/Main Code")

source("functions.R")

library(dplyr)
library(ez)
library(ggplot2)
library(moments)
library(tidyr)
library(Hmisc)
library(ppcor)

# make sure the most recent version of trimr is installed
#devtools::install_github("JimGrange/trimr")
library(trimr)

# import the data
target <- read.csv("raw_target.csv", stringsAsFactors = FALSE)
visual <- read.csv("raw_visual.csv", stringsAsFactors = FALSE)
numeric <- read.csv("raw_numeric.csv", stringsAsFactors = FALSE)
colnames(target) <- c("participant", "trial", "condition", "accuracy", "rt")
colnames(visual) <- c("participant", "trial", "condition", "accuracy", "rt")
colnames(numeric) <- c("participant", "trial", "condition", "accuracy", "rt")

# add accuracy trimming column to each data set & declare each paradigm
target <- mutate(target, paradigm = "target", accTrim = 0)
visual <- mutate(visual, paradigm = "visual", accTrim = 0)
numeric <- mutate(numeric, paradigm = "numeric", accTrim = 0)
#------------------------------------------------------------------------------
# for the 'equal trials' analysis
# visual <- subset(visual, trial < 361)
# numeric <- subset (numeric, trial < 361)

#------------------------------------------------------------------------------
### sort the null trials for each paradigm

## target data first because of the coding error
# trials to remove for participants 1-23
n23 <- c(1, 2, 103, 104, 205, 206, 307, 308)

# trials to remove for participants > 23
n24 <- c(1, 2, 121, 122, 241, 242, 361, 362)

# loop over participants and do the trimming
for(i in 1:nrow(target)){
  if(target$participant[i] <= 23){
    if(target$trial[i] %in% n23) {
      target$condition[i] <- "null"
    }
  }
  if(target$participant[i] > 23){
    if(target$trial[i] %in% n24){
      target$condition[i] <- "null"
    }
  }
}

## visual & numeric null trials
nullTrials <- c(1, 2, 121, 122, 241, 242, 361, 362)

for(i in 1:nrow(visual)){
  if(visual$trial[i] %in% nullTrials){
    visual$condition[i] <- "null"
  }
}

for(i in 1:nrow(numeric)){
  if(numeric$trial[i] %in% nullTrials){
    numeric$condition[i] <- "null"
  }
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### removing the first block of each paradigm to conduct analysis for practice 
### effect

# to be calculated before individual paradigms data are combined 
# (before the null trials are removed)

# target23 <- subset(target, participant < 24)
# target23 <- subset(target23, trial > 102)

# target24 <- subset(target, participant > 23)
# target24 <- subset(target24, trial > 120)

# target <- rbind(target23, target24)
# visual <- subset(visual, trial > 120)
# numeric <- subset(numeric, trial > 120)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### data collation & participant removal checks

# bind all data together
allData <- rbind(target, visual, numeric)

# remove the null trials
allData <- subset(allData, allData$condition != "null")

# check which participants don't have data for all conditions
incompleteRemoval <- completeData(allData)

# check which participants have accuracy too low
accCriterion <- 90
accRemoval <- accuracyRemoval(allData, accCriterion)

# collate all removal, and do the removal
participantsRemoved <- sort(unique(c(incompleteRemoval, accRemoval)))

allData <- allData[!allData$participant %in% participantsRemoved, ]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### overall and per-paradigm proportion of trials removed due to accuracy and 
### RTs trimming

# Part 1 
# 2nd part starts at line 224 (after the trimming is finished)

# take allData, before accuracy trimming (before removing 2 trials after an error),
# and assign the trials length to a vector allTrials 
allTrials <- length(allData$trial)

# additionally subset paradigms from the same allData or calculation of 
# percentage of trials removed per paradigm
allTarget <- subset(allData, paradigm == "target")
allTarget <- length(allTarget$trial)

allVisual <- subset(allData, paradigm == "visual")
allVisual <- length(allVisual$trial)

allNumeric <- subset(allData, paradigm == "numeric")
allNumeric <- length(allNumeric$trial)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# accuracy trimming (remove two trials following an error)
for(i in 3:nrow(allData)){
  allData$accTrim[i] <- allData$accuracy[i - 2] * allData$accuracy[i - 1] 
}

allData <- subset(allData, allData$accTrim == 1)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### main accuracy analysis
accuracy <- allData %>%
  group_by(paradigm, condition, participant) %>%
  summarise(rawAcc = (sum(accuracy) / length(accuracy)) * 100)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### response time analysis

# disable scientific notation
options(scipen = 999)

# trim the slow RTs
rtData <- rtTrimming(allData, sd = 2.5, minRT = 150)

# get the mean RT for each participant
rt <- rtData %>% 
  group_by(paradigm, condition, participant) %>%
  summarise(meanRT = mean(rt))

# how many participants?
nparticipants <- length(unique(rt$participant))

# change paradigm & condition to factor so we can do ANOVA on it
rt$paradigm <- as.factor(rt$paradigm)
rt$condition <- as.factor(rt$condition)

# do the ANOVA
rtANOVA <- ezANOVA(
  data = data.frame(rt), 
  dv = .(meanRT), 
  wid = .(participant), 
  within = .(paradigm, condition), 
  between = NULL, 
  detailed = FALSE
)

rt <- data.frame(rt)

# get the mean RT per cell
meanRT <- rt %>%
  group_by(paradigm, condition) %>%
  summarise(rt = round(mean(meanRT), 0), se = round(sd(meanRT) / 
                                                      sqrt(nparticipants), 0))

# main effect of condition
seq <- rtData %>%
  group_by(condition) %>%
  summarise(meanRT = round(mean(rt), 0), se = round(sd(rt) / 
                                                      sqrt(nparticipants), 0))

# main effect of paradigm
# get the mean RT per cell
paradigm <- rtData %>%
  group_by(paradigm) %>%
  summarise(meanRT = round(mean(rt), 0), se = round(sd(rt) / 
                                                      sqrt(nparticipants), 0))

## t-tests of each paradigm's n-2 repetition cost
targetABA <- subset(rt, rt$paradigm == "target" & rt$condition == "ABA")
targetCBA <- subset(rt, rt$paradigm == "target" & rt$condition == "CBA")
targetTtest <- t.test(targetABA$meanRT, targetCBA$meanRT, paired = TRUE)

visualABA <- subset(rt, rt$paradigm == "visual" & rt$condition == "ABA")
visualCBA <- subset(rt, rt$paradigm == "visual" & rt$condition == "CBA")
visualTtest <- t.test(visualABA$meanRT, visualCBA$meanRT, paired = TRUE)

numericABA <- subset(rt, rt$paradigm == "numeric" & rt$condition == "ABA")
numericCBA <- subset(rt, rt$paradigm == "numeric" & rt$condition == "CBA")
numericTtest <- t.test(numericABA$meanRT, numericCBA$meanRT, paired = TRUE)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Part 2 of calculating proportion of trials removed

# after trimming RTs, assign the length of trimmed data frame to trimmedTrials
trimmedTrials <- length(rtData$trial)

# calculate the overall number of removed trials
removedTrials <- allTrials - trimmedTrials

# subset paradigms from trimmed rtData for calculation of 
# percentage of trials removed per paradigm
# and calculate number of trials removed
trimTarget <- subset(rtData, paradigm == "target")
trimTarget <- length(trimTarget$trial)
removedTarget <- allTarget - trimTarget

trimVisual <- subset(rtData, paradigm =="visual")
trimVisual <- length(trimVisual$trial)
removedVisual <- allVisual - trimVisual

trimNumeric <- subset(rtData, paradigm == "numeric")
trimNumeric <- length(trimNumeric$trial)
removedNumeric <- allNumeric - trimNumeric

# calculate the overall percentage of removed trials
proportionRemoved <- (100*removedTrials)/allTrials
round(proportionRemoved, 1)

# calculate the percetage of trials removed per paradigm
propRemTarget <- (100*removedTarget)/allTarget
round(propRemTarget,1)

propRemVisual <- (100*removedVisual)/ allVisual
round(propRemVisual,1)

propRemNumeric <- (100*removedNumeric)/ allNumeric
round(propRemNumeric,1)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# change paradigm & condition to factor so we can do ANOVA on it
accuracy$paradigm <- as.factor(accuracy$paradigm)
accuracy$condition <- as.factor(accuracy$condition)


# do the ANOVA
accuracyANOVA <- ezANOVA(
  data = data.frame(accuracy), 
  dv = .(rawAcc), 
  wid = .(participant), 
  within = .(paradigm, condition), 
  between = NULL, 
  detailed = FALSE
)

# get the mean accuracy per cell
meanAcc <- accuracy %>%
  group_by(paradigm, condition) %>%
  summarise(meanAcc = round(mean(rawAcc), 2), 
            se = round(sd(rawAcc) / sqrt(nparticipants), 2))

# main effect of condition
seq <- accuracy %>%
  group_by(condition) %>%
  summarise(meanAcc = round(mean(rawAcc), 2), 
            se = round(sd(rawAcc) / sqrt(nparticipants), 2))

# main effect of condition
paradigm <- accuracy %>%
  group_by(paradigm) %>%
  summarise(meanAcc = round(mean(rawAcc), 2), 
            se = round(sd(rawAcc) / sqrt(nparticipants), 2))

## t-tests of each paradigm's n-2 repetition cost
targetABA <- subset(accuracy, accuracy$paradigm == "target" & 
                      accuracy$condition == "ABA")
targetCBA <- subset(accuracy, accuracy$paradigm == "target" & 
                      accuracy$condition == "CBA")
targetTtest <- t.test(targetABA$rawAcc, targetCBA$rawAcc, paired = TRUE)

visualABA <- subset(accuracy, accuracy$paradigm == "visual" & 
                      accuracy$condition == "ABA")
visualCBA <- subset(accuracy, accuracy$paradigm == "visual" & 
                      accuracy$condition == "CBA")
visualTtest <- t.test(visualABA$rawAcc, visualCBA$rawAcc, paired = TRUE)

numericABA <- subset(accuracy, accuracy$paradigm == "numeric" & 
                       accuracy$condition == "ABA")
numericCBA <- subset(accuracy, accuracy$paradigm == "numeric" & 
                       accuracy$condition == "CBA")
numericTtest <- t.test(numericABA$rawAcc, numericCBA$rawAcc, paired = TRUE)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### look at individual differences in the n-2 repetition cost

#---- Response Time
# get a data frame with n-2 repetition cost as the DV. 
wideRt <- spread(rt, condition, meanRT)

n2Cost <- wideRt %>%
  group_by(paradigm,participant) %>%
  summarise(n2Cost=ABA-CBA)

n2Cost$n2Cost <- round(n2Cost$n2Cost, 0)

# load individual differences data
indData <- read.csv("ind_data.csv", stringsAsFactors = FALSE)
wideN2Cost <- spread(n2Cost, paradigm, n2Cost)
corData <- merge(wideN2Cost, indData, by = "participant")

# impute the missing data point for subject 68 (in position 51)
corData$processing[51] <- mean(corData$processing, na.rm = TRUE)

# draw overlapping density functions of n-2 repetition costs
pdf("biDistributions_rt.pdf", width = 8, height = 8)
ggplot(n2Cost, aes(x = n2Cost, colour = paradigm, linetype = paradigm)) + 
  geom_line(stat = "density", size = 1.3) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        panel.background = element_rect(fill = "grey86")) +
  scale_x_continuous(name = "N-2 Repetition Cost (ms)") + 
  scale_y_continuous(name = "Density") 
dev.off()

# get summary of distributions for each paradigm
RtDistributions <- n2Cost %>%
  group_by(paradigm) %>%
  summarise(min = min(n2Cost), 
            max = max(n2Cost), 
            sd = sd(n2Cost), 
            skew = skewness(n2Cost), 
            kurtosis = kurtosis(n2Cost), 
            normality = shapiro.test(n2Cost)$p.value, 
            mean = mean(n2Cost))

#---- Accuracy

# re-calculate mean accuracy per participant/ condition/ paradigm
trimmedAcc <- allData %>%
  group_by(paradigm, condition, participant) %>%
  summarise(rawAcc = (sum(accuracy) / length(accuracy)) * 100)

# change the data frame format to wide
wideTrimmedAcc <- spread(trimmedAcc, condition, rawAcc)

# calculate n-2 repetition cost for accuracy
AccN2Cost <- wideTrimmedAcc %>%
  group_by(paradigm, participant) %>%
  summarise(AccN2Cost = ABA - CBA)

# round to 2 decimal places
AccN2Cost$AccN2Cost <- round(AccN2Cost$AccN2Cost, 2)

# draw overlapping density functions of n-2 repetition costs
pdf("biDistributions_acc.pdf", width = 8, height = 8)
ggplot(AccN2Cost, aes(x = AccN2Cost, colour = paradigm, linetype = paradigm)) + 
  geom_line(stat = "density", size = 1.3) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        panel.background = element_rect(fill = "grey86")) +
  scale_x_continuous(name = "N-2 Repetition Cost (Accuracy)") + 
  scale_y_continuous(name = "Density") 
dev.off()


# get summary of n-2 repetition costs for accuracy distributions for each paradigm
AccDistributions <- AccN2Cost %>%
  group_by(paradigm) %>%
  summarise(min = min(AccN2Cost), 
            max = max(AccN2Cost), 
            sd = sd(AccN2Cost), 
            skew = skewness(AccN2Cost), 
            kurtosis = kurtosis(AccN2Cost), 
            normality = shapiro.test(AccN2Cost)$p.value, 
            mean = mean(AccN2Cost))

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Sequencing effect analysis

# load a .csv file with data of the possible order of paradigms
order <- read.csv("paradigms_order.csv", stringsAsFactors = FALSE)

# change names of columns
colnames(order) <- c("participant", "sixOrders", "threeOrders")

# remove column with data of order of experiment components 
# which include processing speed and RRS order
order <- order[,-2]

# combine the data frame for n-2 repetition cost, the n2Cost and order.csv
orderData <- merge(n2Cost, order, by = "participant")

orderData$threeOrders <- as.numeric(orderData$threeOrders)

# subset orderData based on order 1, 2, and 3
# then assign 1,2,3 depending on which paradigm was conducted first

# order1: target 1st, visual 2nd, numeric 3rd
order1 <- subset(orderData, orderData$threeOrders == 1)

# target is already coded as 1st
# code visual as 2nd
for (i in 1:nrow(order1)){
  if (order1$paradigm[i] == "visual"){
    order1$threeOrders[i] = 2}
}

# code numeric as 3rd
for (i in 1:nrow(order1)){
  if(order1$paradigm[i] == "numeric"){
    order1$threeOrders[i] = 3
  }
}

# order2: visual 1st, numeric 2nd, target 3rd
# subset order2 from orderData
order2 <- subset(orderData, orderData$threeOrders == 2)

# code visual as 1st
for (i in 1:nrow(order2)){
  if(order2$paradigm[i] == "visual"){
    order2$threeOrders[i] = 1
  }
}

# numeric is already coded as 2nd

# code target as 3rd
for (i in 1:nrow(order2)){
  if(order2$paradigm[i] == "target"){
    order2$threeOrders[i] = 3
  }
}

# order3: numeric 1st, target 2nd, visual 3rd
# subset order3 from orderData
order3 <- subset(orderData, orderData$threeOrders == 3)

# code numeric as 1st
for (i in 1:nrow(order3)){
  if(order3$paradigm[i] == "numeric"){
    order3$threeOrders[i] = 1
  }
}

# code target as 2nd
for (i in 1:nrow(order3)){
  if(order3$paradigm[i] == "target"){
    order3$threeOrders[i] = 2
  }
}

# visual is already coded as 3rd

# combine the order1, order2, and order3, wich have correctly coded order
n2CostOrder <- rbind(order1, order2, order3)

# anova
n2CostOrder$paradigm <- as.factor(n2CostOrder$paradigm)
n2CostOrder$threeOrders <- as.factor(n2CostOrder$threeOrders)

# ANOVA for n2cost as DV and threeOrders as IV
orderANOVA <- ezANOVA(
  data = data.frame(n2CostOrder), 
  dv = .(n2Cost), 
  wid = .(participant), 
  within = .(threeOrders), 
  between = NULL, 
  detailed = FALSE
)

meanN2CostOrder <- n2CostOrder %>%
  group_by(threeOrders) %>%
  summarise(meanN2Cost= round(mean(n2Cost), 0), 
            se = round(sd(n2Cost) / sqrt(nparticipants), 0))
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### correlations

#--- Response time

# Overall (mean) response times
rtCor <- rtData %>% 
  group_by(paradigm, participant) %>%
  summarise(meanRT = mean(rt))

wideRtCor <- spread(rtCor, paradigm, meanRT)
wideRtCor <- merge(wideRtCor, indData, by = "participant")

indRtCor <- rcorr(as.matrix(wideRtCor))

# n-2 repetition cost and Ind Diff correlations
indCor <- rcorr(as.matrix(corData))

# partial correlations, controlling for processing speed (as requested
# by reviewer)
partial_target_visual <- pcor.test(corData$target, corData$visual, 
                                   corData$processing)
partial_target_numeric <- pcor.test(corData$target, corData$numeric, 
                                    corData$processing)
partial_visual_numeric <- pcor.test(corData$visual, corData$numeric, 
                                    corData$processing)


#--- Accuracy

# Overall (mean) accuracy
accAve <- allData %>%
  group_by(paradigm, participant) %>%
  summarise(rawAcc = (sum(accuracy) / length(accuracy)) * 100)

wideAcc <- spread(accAve, paradigm, rawAcc)
wideAcc <- merge(wideAcc, indData, by = "participant")
accCor <- rcorr(as.matrix(wideAcc))

# change data frame format to wide
wideAccN2Cost <- spread(AccN2Cost, paradigm, AccN2Cost)
wideAccN2Cor <- merge(wideAccN2Cost, indData, by = "participant")

# impute the missing data point for subject 68 (in position 51)
wideAccN2Cor$processing[51] <- mean(wideAccN2Cor$processing, na.rm = TRUE)

# calulate correlations
indAccN2Cor <- rcorr(as.matrix(wideAccN2Cor))
round(indAccN2Cor$r, 2)
round(indAccN2Cor$P, 3)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# multiple regressions (just for RT)

# normalise Ind Diff scores for regression
corData$rumination <- scale(corData$rumination)
corData$processing <- scale(corData$processing)
nIndCor <- rcorr(as.matrix(corData))
visualReg <- lm(visual ~ rumination + processing, data = corData)
targetReg <- lm(target ~ rumination + processing, data = corData)
numericReg <- lm(numeric ~ rumination + processing, data = corData)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### do the reliability checks

# run the reliability function for response times
set.seed(200)
correlations_rt <- splitHalf(rtData, splitType = "random", nSplits = 500)
colnames(correlations_rt) <- c("Target Detection", "Visual Judgment", 
                               "Numeric Judgment")


# violin plot of the reliability bootstrap
library(vioplot)
pdf("violin Reliability_rt.pdf", width = 8, height = 8)
vioplot(correlations_rt[, 1], correlations_rt[, 2], correlations_rt[, 3], 
        col = "skyblue", names = c("Target Detection", "Visual Judgment", 
                                   "Numeric Judgment"), lwd = 1.5, 
        ylim = c(-0.2, 1))
title(ylab = "Correlation (r)", xlab = "Paradigm")
abline(h = 0.5385, lwd = 2, lty = 2)
dev.off()


# run the reliability function for accuracy
set.seed(200)
correlations_acc <- splitHalf_acc(allData, splitType = "random", nSplits = 500)
colnames(correlations_acc) <- c("Target Detection", "Visual Judgment", 
                                "Numeric Judgment")

# violin plot of the reliability bootstrap
library(vioplot)
pdf("violin Reliability_accuracy.pdf", width = 8, height = 8)
vioplot(correlations_acc[, 1], correlations_acc[, 2], correlations_acc[, 3], 
        col = "skyblue", names = c("Target Detection", "Visual Judgment", 
                                   "Numeric Judgment"), lwd = 1.5, 
        ylim = c(-0.2, 1))
title(ylab = "Correlation (r)", xlab = "Paradigm")
abline(h = 0.5385, lwd = 2, lty = 2)
dev.off()
#------------------------------------------------------------------------------