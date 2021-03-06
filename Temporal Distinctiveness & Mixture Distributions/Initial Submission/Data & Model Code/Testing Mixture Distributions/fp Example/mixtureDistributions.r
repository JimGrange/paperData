rm(list = ls())
setwd("~/Git/paperData/2015/Temporal Distinctiveness & Mixture Distributions/Data & Model Code/Testing Mixture Distributions/fp Example")

source("functions.R")
library(dplyr)

set.seed(1)

#------------------------------------------------------------------------------
### do mixture distribution

# declare the parameters
p <- 0.7
A <- 500
b <- c(A + 120, A + 120)
drift <- c(1.2, 0.6)
s <- 0.3
ter <- 220

# how many subjects?
nSubjects <- 16

# how many trials?
nTrials <- 200

# get the data
z <- NULL
x <- NULL
y <- NULL

# simulate the mixed distribution
for(i in 1:nSubjects){
  dataZ <- simLBA_Mixture(n = nTrials, bs = b, a = A, drift = drift, s = s,
                          ter = ter, p = p)
  dataZ <- mutate(dataZ, cond = 2, pp = i)
  
  z <- rbind(z, dataZ)
}

# simulate the fast distribution x
for(i in 1:nSubjects){
  dataX <- simLBA_Mixture(n = nTrials, bs = b, a = A, drift = drift, s = s, 
                          ter = ter, p = 1)
  dataX <- mutate(dataX, cond = 1, pp = i)
  
  x <- rbind(x, dataX)
}

# simulate the slow distribution y
for(i in 1:nSubjects){
  dataY <- simLBA_Mixture(n = nTrials, bs = b, a = A, drift = drift, s = s, 
                          ter = ter, p = 0)
  dataY <- mutate(dataY, cond = 3, pp = i)
  
  y <- rbind(y, dataY)
}

# collate the data
data <- rbind(z, x, y)
data <- subset(data, data$accuracy == 1)
data <- subset(data, data$rt < 2000)
data <- select(data, rt, cond, pp)
data$rt <- data$rt / 1000


# get the fp1 object
res <- fpGet(data[, c("rt", "cond")], 1000, bw = 0.75)  

# call fixed-point plot
par(las = 1, mfrow = c(1, 3), mar = c(4.2, 4.2, 1, 1))  

# make the plot and save it as PDF
pdf("fpExample.pdf", width = 8, height = 5)
plot(res, xlab = "Resonse Time (S)", col = 1) 
dev.off()
#------------------------------------------------------------------------------