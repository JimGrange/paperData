#------------------------------------------------------------------------------
rm(list = ls())
setwd("~/Git/paperData/2015/Temporal Distinctiveness & Mixture Distributions/Data & Model Code/Model Fit")

# load required libraries & functions
library(dplyr)
source("functions.R")

# get the data (& do trimming & round RTs to nearest integer)
data <- read.csv("data.csv")
data <- subset(data, data$condition == "Memory" & data$errorTrimming == 1 & 
                 data$sequence != "null" & data$sequence == "repeat")

# change the "subject" to participant, just for the SD trimming
colnames(data)[1] <- "participant"

# do the SD trimming
data <- sdTrim(data, minRT = 150, sd = 2.5, perCondition = TRUE, 
               perParticipant = TRUE, omitErrors = FALSE, returnType = "raw")
colnames(data)[1] <- "subject"

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# get the data for QMP
quantiles <- c(.1, .3, .5, .7, .9)

fast <- subset(data, data$ratio == 20)
fast <- prepData(fast)

slow <- subset(data, data$ratio == 0.05)
slow <- prepData(slow)

int <- subset(data, data$ratio == 1)
int <- prepData(int)

data <- list(fast = fast, slow = slow, int = int)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### fit the model

## get the starting parameters separately by condition (A, b, drift, s, ter)
slowParms <- getStartParms(data$slow)
fastParms <- getStartParms(data$fast)
intParms <- getStartParms(data$int)

# sort the starting parameters for overall fit
A <- mean(c(slowParms[1], fastParms[1], intParms[1]))
b <- c(fastParms[2], slowParms[2])
driftFast <- fastParms[3]
driftSlow <- slowParms[3]
driftInt <- intParms[3]
s <- mean(c(slowParms[4], fastParms[4], intParms[4]))
ter <- mean(c(slowParms[5], fastParms[5], intParms[5]))
p <- 0.7

# store the parameters
parameters <- c(log(p), log(A), log(b), driftFast, driftSlow, log(s), 
                log(ter))


# do the fit
fit <- optim(par = parameters, fn = fitFunction, data = data, 
             list(maxit = 10000, parscale = parameters))

x <- fit$par
x[c(1, 2, 3, 4, 7, 8)] <- exp(x[c(1, 2, 3, 4, 7, 8)]) 
x <- round(x, 3)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### assess the fit

# first, extract the parameters
p <- x[1]
A <- x[2]
bFast <- x[3] + A
bSlow <- x[4] + A
driftFast <- x[5]
driftSlow <- x[6]
s <- x[7]
ter <- x[8]

## get model predictions
n = 50000

# fast condition
fastSim <- simLBA_Mixture(n = n, b = c(bFast, bSlow), a = A, 
                          drift = c(driftFast, driftSlow), s = s, ter = ter, 
                          p = 1)
fastSim <- prepData(fastSim)

# slow condition
slowSim <- simLBA_Mixture(n = n, b = c(bFast, bSlow), a = A, 
                          drift = c(driftFast, driftSlow), s = s, ter = ter, 
                          p = 0)
slowSim <- prepData(slowSim)

# mixture distribution
intSim <- simLBA_Mixture(n = n, b = c(bFast, bSlow), a = A, 
                         drift = c(driftFast, driftSlow), s = s, ter = ter, 
                         p = p)
intSim <- prepData(intSim)

## do some plots


# fast
fastHumanCorrectRT <- data$fast$q$correct
fastHumanCorrectCDF <- quantiles * data$fast$p
fastHumanErrorRT <- data$fast$q$error
fastHumanErrorCDF <- quantiles * (1 - data$fast$p)

fastModelCorrectRT <- fastSim$q$correct
fastModelCorrectCDF <- quantiles * fastSim$p
fastModelErrorRT <- fastSim$q$error
fastModelErrorCDF <- quantiles * (1 - fastSim$p)


# intermediate
intHumanCorrectRT <- data$int$q$correct
intHumanCorrectCDF <- quantiles * data$int$p
intHumanErrorRT <- data$int$q$error
intHumanErrorCDF <- quantiles * (1 - data$int$p)

intModelCorrectRT <- intSim$q$correct
intModelCorrectCDF <- quantiles * intSim$p
intModelErrorRT <- intSim$q$error
intModelErrorCDF <- quantiles * (1 - intSim$p)

# slow
slowHumanCorrectRT <- data$slow$q$correct
slowHumanCorrectCDF <- quantiles * data$slow$p
slowHumanErrorRT <- data$slow$q$error
slowHumanErrorCDF <- quantiles * (1 - data$slow$p)

slowModelCorrectRT <- slowSim$q$correct
slowModelCorrectCDF <- quantiles * slowSim$p
slowModelErrorRT <- slowSim$q$error
slowModelErrorCDF <- quantiles * (1 - slowSim$p)

# get min and max RT
minRT <- min(slowHumanCorrectRT, slowHumanErrorRT,
             intHumanCorrectRT, intHumanErrorRT,
             fastHumanCorrectRT, fastHumanErrorRT,
             slowModelCorrectRT, slowModelErrorRT,
             intModelCorrectRT, intModelErrorRT,
             fastModelCorrectRT, fastModelErrorRT)

maxRT <- max(slowHumanCorrectRT, slowHumanErrorRT,
             intHumanCorrectRT, intHumanErrorRT,
             fastHumanCorrectRT, fastHumanErrorRT,
             slowModelCorrectRT, slowModelErrorRT,
             intModelCorrectRT, intModelErrorRT,
             fastModelCorrectRT, fastModelErrorRT)



pdf("qmp2RCI_drift + b.pdf", width = 8, height = 4)
par(mfrow = c(1, 3))
# plot slow
plot(slowHumanCorrectRT, slowHumanCorrectCDF, type = "p", 
     xlim = c(minRT, maxRT), ylim = c(0, 1), main = "0.05", pch = 19, 
     xlab = "Response Time (ms)", ylab = "Defective CDF")
points(slowHumanErrorRT, slowHumanErrorCDF, type = "p")
lines(slowModelCorrectRT, slowModelCorrectCDF, type = "b", lty = 1, lwd = 1, 
      pch = 19, cex = 0.2)
lines(slowModelErrorRT, slowModelErrorCDF, type = "b", lty = 2, lwd = 1,
      pch = 19, cex = 0.2)
legend("topleft", c("Data (Correct)", "Data (Error)",
                    "Model (Correct)", "Model (Error)"), cex = 1, 
       pch = c(19, 1, NA, NA), lty = c(NA, NA, 1, 2),  bty = "n")

# plot intermediate
plot(intHumanCorrectRT, intHumanCorrectCDF, type = "p", 
     xlim = c(minRT, maxRT), ylim = c(0, 1), main = "1", pch = 19, 
     xlab = "Response Time (ms)", ylab = "Defective CDF")
points(intHumanErrorRT, intHumanErrorCDF, type = "p")
lines(intModelCorrectRT, intModelCorrectCDF, type = "b", lty = 1, lwd = 1, 
      pch = 19, cex = 0.2)
lines(intModelErrorRT, intModelErrorCDF, type = "b", lty = 2, lwd = 1,
      pch = 19, cex = 0.2)

# plot fast
plot(fastHumanCorrectRT, fastHumanCorrectCDF, type = "p", 
     xlim = c(minRT, maxRT), ylim = c(0, 1), main = "20", pch = 19, 
     xlab = "Response Time (ms)", ylab = "Defective CDF")
points(fastHumanErrorRT, fastHumanErrorCDF, type = "p")
lines(fastModelCorrectRT, fastModelCorrectCDF, type = "b", lty = 1, lwd = 1, 
      pch = 19, cex = 0.2)
lines(fastModelErrorRT, fastModelErrorCDF, type = "b", lty = 2, lwd = 1, 
      pch = 19, cex = 0.2)
dev.off()
#------------------------------------------------------------------------------
#   p     A       bFast   bSlow     vFast   vSlow   s     ter
round(x, 3)