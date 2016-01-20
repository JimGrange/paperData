### Analysis script for Mayr replication

#------------------------------------------------------------------------------
### General set-up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Git/paperData/Episodic Retrieval & Inhibition/Data & Analysis Code")

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
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### response time analysis

# remove all error trials, so we are ready for RT analysis
data <- subset(data, data$accuracy == 1)

# trim the RTs with fixed standard deviation upper limit
rtData <- getRTs(data = data, minRT = 150, sd = 2.5)

#---------
# calculate the mean BI for each condition
biRT <- rtData %>%
  group_by(subject) %>%
  summarise(rep = abaRep - cbaRep, sw = abaSw - cbaSw)

# get the difference score
biDiff <- biRT$rep - biRT$sw


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
     col = "gray48", lwd = 2)

abline(h = 0, lwd = 1)
abline(h = log(6), col = "black", lty = 2, lwd = 2)
abline(h = log(1/6), col = "black", lty = 2, lwd = 2)
text(0, log(3) + 0.4, labels = "Evidence for Alternative", 
     cex = 1.5, col = "black", pos = 4)
text(0, log(1/3) - 0.4, labels = "Evidence for Null", 
     cex = 1.5, col = "black", pos = 4)
dev.off()

#---
## robustness check based on prior

# what is the t-value for the data?
tVal <- as.numeric(t.test(biRT$rep, biRT$sw, paired = TRUE)[['statistic']])

# what are the priors to explore?
priors <- seq(from = 0.01, to = 1.5, length.out = 1000)

# get the Bayes factor for each prior value
robust <- sapply(priors, function(x) 
  exp(ttest.tstat(t = tVal, n1 = nrow(biRT), rscale = x)[['bf']]))

# plot it
pdf("robustPrior.pdf", width = 8, height = 5)
plot(priors, robust, type = "l", lwd = 2, col = "gray48",
     ylim = c(0, max(robust)), xaxt = "n", xlab = "Cauchy Prior Width (r)", 
              ylab = "Bayes Factor (10)")
abline(h = 0, lwd = 1)
abline(h = 6, col = "black", lty = 2, lwd = 2)
axis(1, at = seq(0, 1.5, 0.25))
points(0.707, extractBF(bfDiff, onlybf = TRUE), col = "black", 
       cex = 1.2, pch = 21, bg = "skyblue")
legend(x = 1, y = 13, legend = c("Default Prior", "Stopping Rule"),
       pch = c(21, NA), lty = c(NA, 2), lwd = c(NA, 2), pt.cex = c(1, NA),
       col = c("black", "black"), pt.bg = "skyblue", bty = "n")
dev.off()



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
  summarise(accMean = mean(acc), 
            se = (sd(acc)) / sqrt(length(biDiff)))

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

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Bayesian t-test with Kruschke's t-test
bayesK <- BESTmcmc(y1 = biDiff)

# plot the mean and effect size estimates
pdf("bayesParameter.pdf", width = 8, height = 5)
par(mfrow = c(1, 2))
plot(bayesK, which = c("mean"))
plot(bayesK, which = c("effect"))
dev.off()

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
### plot mean RTs (not in paper)
pdf("meanRT.pdf", width = 8, height = 8)
plot <- ggplot(rtSummary, (aes(x = sequence, y = meanRT, group = stimRep, 
                               colour = stimRep)))
plot <- plot + geom_errorbar(aes(ymin = meanRT, ymax = meanRT + se), 
                             width = 0.05, size = 0.5)
plot <- plot + geom_line(aes(linetype = stimRep))
plot <- plot + geom_point(aes(shape = stimRep), size = 2.3)
plot
dev.off()


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### simulation showing recruitment order does not change stopping rule

# set random seed so user can re-produce plot
set.seed(42)

# how many "experiments" to simulate?
nSims = 50

# get a matrix to store all of the "experiments" in
orderAnalysis <- matrix(0, nrow = length(biDiff), ncol = nSims)

# fill the matrix with random recruitment orders
for(i in 1:nSims){
  
  # simulate new ordering
  orderData <- base::sample(biDiff, size = length(biDiff), 
                            replace = FALSE)
  bf <- plotBF(orderData, scale = 0.707)
  orderAnalysis[, i] <- log(bf[, 2])
  orderAnalysis[1, i] <- 0
}

# plot the progression of the Bayes Factor of the interaction as more 
# subjects are added

pdf("recruitmentOrder.pdf", width = 8, height = 5)

# plot the first "experiment"
plot(seq(1:length(biDiff)), orderAnalysis[, 1], type = "l", 
     xlab = "Sample Size", ylab = "(Log) Bayes Factor (10)", 
     ylim = c(-2.5, 2.5), pch = 19, col = "gray", lwd = 2)

# now do the rest
for(i in 2:nSims){
  lines(seq(1:length(biDiff)), orderAnalysis[, i], type = "l", col = "gray", 
        lwd = 2)
}

# add the plot details
abline(h = 0, lwd = 1)
abline(h = log(6), col = "black", lty = 2, lwd = 2)
abline(h = log(1/6), col = "black", lty = 2, lwd = 2)
text(0, log(3) + 0.4, labels = "Evidence for Alternative", 
     cex = 1.5, col = "black", pos = 4)
text(0, log(1/3) - 0.4, labels = "Evidence for Null", 
     cex = 1.5, col = "black", pos = 4)

dev.off()

#------------------------------------------------------------------------------


