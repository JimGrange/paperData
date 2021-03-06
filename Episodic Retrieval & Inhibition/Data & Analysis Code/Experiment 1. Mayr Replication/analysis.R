### Analysis script for Agi's experiment (Mayr replication)

#------------------------------------------------------------------------------
### General set-up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Git/lab-book/Agi's PhD/Mayr Replication")

# load necessary functions file & load necessary packages
source("functions.R")
library(dplyr)
library(ggplot2)
library(ez)
library(BayesFactor)
library(BEST)

# load the data file
fullData <- read.csv("allData.csv", header = TRUE)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Error Data Trimming 

# Remove null trials (1st and 2nd trial of each block)

# which trials were null?
nullTrials <- c(1, 2, 121, 122, 241, 242, 361, 362)

# remove these trials
data <- fullData[! fullData$trial %in% nullTrials, ]


# add column to code for removing two trials following an error, initialised 
# to zero
data <- mutate(data, accTrim = 0)

# populate this column
for(i in 3:nrow(data)){
  data$accTrim[i] <- data$accuracy[i - 2] * data$accuracy[i - 1]
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Error Analysis

# First, we need to see the OVERALL accuracy per participant, as we will want
# to remove people who do not meet a certain criteria
accuracy <- data %>%
  group_by(subject) %>%
  summarise(meanAccuracy = (sum(accuracy) / length(accuracy) * 100))

# which subjects fall below our criterion of 90% accuracy?
errorRemove <- accuracy[accuracy$meanAccuracy <= 90, ]$subject

# remove these subjects from the data file
data <- data[! data$subject %in% errorRemove, ]

# remove the 2 trials following an error
data <- subset(data, data$accTrim == 1)

# now produce a data frame with accuracy analysis per condition
accuracy <- data %>%
  group_by(subject, stimRep, sequence) %>%
  summarise(acc = (sum(accuracy) / length(accuracy)))

# now remove all error trials, so we are ready for RT analysis
data <- subset(data, data$accuracy == 1)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### response time analysis

# trim the RTs with fixed standard deviation upper limit
rtData <- getRTs(data = data, minRT = 150, sd = 2.5)

#---------
# calculate the mean BI for each condition
biRT <- rtData %>%
  group_by(subject) %>%
  summarise(rep = abaRep - cbaRep, sw = abaSw - cbaSw)

# get the difference score
biDiff <- biRT$rep - biRT$sw

# get the n--2 repetition cost for response repetitions and export it
# this is used in the mini-meta-analysis in the paper
biRep <- biRT$rep
write.table(biRep, file = "experiment_1_respRep.csv", row.names = FALSE, 
            col.names = FALSE)

#---------
# Bayes Factor Analysis

# what is the Bayes Factor of the difference in BI between stimRep & stimSw?
bfDiff <- ttestBF(x = biDiff, rscale = 0.707)

# plot the progression of the Bayes Factor of the interaction as more 
# subjects are added
bfProg <- plotBF(biDiff, scale = 0.707)
bfProg[, 2] <- log(bfProg[, 2])
bfProg[1, 2] <- 0

pdf("bayesFactor.pdf", width = 8, height = 5)
plot(bfProg[, 1], bfProg[, 2], type = "b", xlab = "Sample Size", 
     ylab = "(Log) Bayes Factor (10)", ylim = c(-2.5, 2.5), pch = 19, 
     col = "skyblue", lwd = 2)

abline(h = 0, lwd = 2)
abline(h = log(6), col = "red", lty = 2, lwd = 2)
abline(h = log(1/6), col = "red", lty = 2, lwd = 2)
text(0, log(3) + 0.4, labels = "Evidence for Alternative", 
     cex = 1.5, col = "red", pos = 4)
text(0, log(1/3) - 0.4, labels = "Evidence for Null", 
     cex = 1.5, col = "red", pos = 4)
dev.off()

#---------
# # BEST analysis
# bestDiff <- BESTmcmc(biDiff)
# bestSummary <- summary(bestDiff)
# 
# # plot the analysis
# pdf("BESTsummary.pdf", width = 8, height = 6)
# plotAll(bestDiff, ROPEeff = c(-0.2, 0.2))
# dev.off()


#---------
### ANOVA analsyis

# disable scientific notation
options(scipen = 999)
options(digits = 4)

#--- accuracy
accuracy <- data.frame(accuracy)
aovACC <- ezANOVA(
  data = accuracy 
  , dv = .(acc)
  , wid = .(subject)
  , within = .(sequence, stimRep)
  , between = NULL 
  , detailed = FALSE)

# summarise the main effects
options(digits = 4)
accuracy %>%
  group_by(stimRep, sequence) %>%
  summarise(accMean = mean(acc))

#--- response time
# transfer to long format for analysis & later plotting, & change necessary
# columns to factors
finalRT <- doLong(rtData)
finalRT$sequence <- factor(finalRT$sequence)
finalRT$stimRep <- factor(finalRT$stimRep)
finalRT$subject <- factor(finalRT$subject)


# summarise the response times
rtSummary <- finalRT %>%
  group_by(stimRep, sequence) %>%
  summarise(meanRT = mean(rt), se = (sd(rt)) / sqrt(length(biDiff)))


# ANOVA using ezANOVA
aovRT <- ezANOVA(
  data = finalRT
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, stimRep)
  , between = NULL
  , detailed = FALSE)


# do the ANOVA the standard way
summary(aov(rt ~ sequence * stimRep + Error(subject / (sequence * stimRep)), 
                                            data = finalRT))
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
### do somple plotting

plotRTs <- rtSummary


levels(plotRTs$stimRep)[levels(plotRTs$stimRep) == "repetition"] <- "Repetition"
levels(plotRTs$stimRep)[levels(plotRTs$stimRep) == "switch"] <- "Switch"
levels(plotRTs$sequence)[levels(plotRTs$sequence) == "aba"] <- "ABA"
levels(plotRTs$sequence)[levels(plotRTs$sequence) == "cba"] <- "CBA"

                       

pd <- position_dodge(0.08)
theme_set(theme_gray(base_size = 18))

pdf("meanRT.pdf", width = 8, height = 8)
plot <- ggplot(plotRTs, (aes(x = sequence, y = meanRT, group = stimRep, 
                               colour = stimRep)))
plot <- plot + geom_errorbar(aes(ymin = meanRT - se, ymax = meanRT + se), 
                             width = 0.05, size = 0.5, position = pd)
plot <- plot + geom_line(aes(linetype = stimRep), position = pd)
plot <- plot + geom_point(aes(shape = stimRep), size = 3.3, position = pd)
plot <- plot + scale_shape_discrete(name = "Response") + 
          scale_linetype_discrete(name = "Response") +
          scale_colour_discrete(name = "Response")
plot <- plot + labs(x = "Task Sequence", y = "Mean Response Time (ms)")
plot + theme(panel.background = element_rect(fill = "grey86")) + 
  annotate("text", x = 1.5, y = 1030, label = "86ms", size = 5) + 
  annotate("text", x = 1.5, y = 950, label = "48ms", size = 5)

dev.off()


#------------------------------------------------------------------------------



#---------
### additional Bayesian analysis

# get the progression of the HDI as sample size increases
# (this takes quite some time!)
hdiProgression <- plotHDI(biDiff)
hdiMax <- max(hdiProgression[, 2:3])
hdiMin <- min(hdiProgression[, 2:3])




pdf("hdiProgression.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# plot the progression of the HDI as sample increases
for(i in 2:nrow(hdiProgression)){
  
  if(i == 2){
    plot(hdiProgression[, 1], hdiProgression[, 2], ylim = c(hdiMin, hdiMax),
         type = "n", xlab = "Sample Size", ylab = "Effect Size (d)")
    
    # add the ROPE
    abline(h = 0.2, lty = 2, col = "red", lwd = 2)
    abline(h = -0.2, lty = 2, col = "red", lwd = 2)
  }
  
  segments(x0 = i, y0 = median(as.numeric(hdiProgression[i, 2:3])), 
           x1 = i, y1 = c(hdiProgression[i, 2], hdiProgression[i, 3]), 
           lwd = 3, col = "skyblue")
  
}

# plot the progression of the HDI width (i.e., precision)
hdiProgression <- data.frame(hdiProgression)
hdiProgression <- mutate(hdiProgression, precision = hdiHigh - hdiLow)

plot(hdiProgression[2:nrow(hdiProgression), 1],  
     hdiProgression[2:nrow(hdiProgression), 4], 
     ylim = c(0, max(hdiProgression$precision)), type = "l", 
     xlab = "Sample Size", ylab = "Precision (0.8 * HDI Width)", 
     lwd = 3, col = "skyblue")

# add the stopping-point at 0.8 * ROPE width
abline(h = 0.8 * 0.4, lty = 2, col = "red", lwd = 2)

dev.off()

