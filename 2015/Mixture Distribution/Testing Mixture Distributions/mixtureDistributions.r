rm(list = ls())
setwd("~/Git/paperData/2015/Mixture Distribution/Testing Mixture Distributions")

source("functions.R")
library(dplyr)

# get the data (& do trimming)
data <- read.csv("data.csv")
data$rt <- round(data$rt / 1000, 3)
data = sdTrim(data, minRT = 0.150, sd = 2.5, perCondition = FALSE, 
              perParticipant = TRUE,  returnType = "raw")

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### do the mixture plots 

# define ratios to include in the analysis
ratios <- c(0.05, 1, 20)
ratioData <- subset(data, data$ratio %in% ratios)

# subset of data to pass to fp package. Then, change names to match 
# fp package function requirements
ratioData <- select(ratioData, rt, ratio, participant)
colnames(ratioData) <- c("rt", "cond", "pp")
# get the fp1 object
res <- fpGet(ratioData[, c("rt", "cond")], 1000, bw = 0.75)  

# call plot
par(las = 1, mfrow = c(1, 3), mar = c(4.2, 4.2, 1, 1))  
plot(res, xlab = "Resonse Time (S)", col = 1) 

# do Bayesian analysis
res <- tapply(1:nrow(ratioData), ratioData$pp, function(X) {
  fpGet(ratioData[X, ], 1000, bw = 0.75)
})  # compute the list
fpAnova(res, stat = "both")

# do box-plot
crosses = fpDensDiff(res)  # first, get the crossing points
boxplot(t(crosses), frame.plot = F, xlab = "Crossing point", 
        ylab = "Condition pair", names = c("3-2", "2-1", "3-1"), 
        horizontal = T)
#------------------------------------------------------------------------------