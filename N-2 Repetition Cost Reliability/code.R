rm(list = ls()) 
setwd("~/Git/paperData/N-2 Repetition Cost Reliability")

source("functions.R")

library(dplyr)
library(ez)
library(ggplot2)
library(moments)
library(tidyr)
library(Hmisc)

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

# get a data frame with n-2 repetition cost as the DV. 
# SURELY there is a more elegant way of doing this...
# ( Agi: see below how I do it, from line 271)
wideRt <- spread(rt, condition, meanRT)

n2Cost <- wideRt %>%
  group_by(paradigm,participant) %>%
  summarise(n2Cost=ABA-CBA)

n2Cost$n2Cost <- round(n2Cost$n2Cost, 0)

#load individual differences data
indData <- read.csv("ind_data.csv", stringsAsFactors = FALSE)
wideN2Cost <- spread(n2Cost, paradigm, n2Cost)
corData <- merge(wideN2Cost, indData, by = "participant")

# draw overlapping density functions of n-2 repetition costs
pdf("biDistributions.pdf", width = 8, height = 8)
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
distributions <- n2Cost %>%
  group_by(paradigm) %>%
  summarise(min = min(n2Cost), 
            max = max(n2Cost), 
            sd = sd(n2Cost), 
            skew = skewness(n2Cost), 
            kurtosis = kurtosis(n2Cost), 
            normality = shapiro.test(n2Cost)$p.value, 
            mean = mean(n2Cost))

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### correlations

# accuracy

accAve <- allData %>%
  group_by(paradigm, participant) %>%
  summarise(rawAcc = (sum(accuracy) / length(accuracy)) * 100)

wideAcc <- spread(accAve, paradigm, rawAcc)
wideAcc <- merge(wideAcc, indData, by = "participant")

# obtain mean average and standard deviations for RRS and processing speed
meanRRS <- mean(wideAcc$rumination)
sdRRS <- sd(wideAcc$rumination)

meanProc <- mean(wideAcc$processing, na.rm = TRUE)
sdProc <- sd(wideAcc$processing, na.rm = TRUE)

accCor <- rcorr(as.matrix(wideAcc))

# normalised Ind Diff scores and accuracy correlations
wideAcc$rumination <- scale(wideAcc$rumination)
wideAcc$processing <- scale(wideAcc$processing)
nAccCor <- rcorr(as.matrix(wideAcc))

# RTs and Ind Diff scores correlations
rtCor <- rtData %>% 
  group_by(paradigm, participant) %>%
  summarise(meanRT = mean(rt))


wideRtCor <- spread(rtCor, paradigm, meanRT)
wideRtCor <- merge(wideRtCor, indData, by = "participant")

indRtCor <- rcorr(as.matrix(wideRtCor))

# normalised Ind Diff scores and RTs correlations
corData$rumination <- scale(corData$rumination)
corData$processing <- scale(corData$processing)
nIndCor <- rcorr(as.matrix(corData))

#------------------------------------------------------------------------------
# n-2 repetition cost and Ind Diff correlations

indCor <- rcorr(as.matrix(corData))

# normalised Ind Diff scores and n-2 rep cost correlations

corData$rumination <- scale(corData$rumination)
corData$processing <- scale(corData$processing)
nIndCor <- rcorr(as.matrix(corData))
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# multiple regressions
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