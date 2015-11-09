#------------------------------------------------------------------------------
# clear workspace:
rm(list=ls(all=TRUE))

setwd("D:/Work/Research/My Papers/In Preparation/Grange et al. (Hangover & RTs)/data & code/Correlation & Descriptives")

data <- read.csv("finalHangoverInfo.csv", header = TRUE)

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### do descriptives

# change negative eBAC to zero
data$eBAC[data$eBAC < 0] <- 0

gender <- table(data$Sex_1F)
age <- mean(data$Age)
  sdAge <- sd(data$Age)
bmi <- mean(data$BMI)
sdBMI <- sd(data$BMI)
ahs <- mean(data$HO_AHS, na.rm = TRUE)
  sdAHS <- sd(data$HO_AHS, na.rm = TRUE)
eBAC <- mean(data$eBAC)
  sdEBAC <- sd(data$eBAC)
units <- mean(data$HO_Units)
  sdUnits <- sd(data$HO_Units)
usual <- mean(data$Usual_units, na.rm = TRUE)
  sdUsual <- sd(data$Usual_units, na.rm = TRUE)
sleep <- mean(data$Epworth)
  sdSleep <- sd(data$Epworth)
#------------------------------------------------------------------------------


------------------------------------------------------------------------------
op <- par(cex.lab = 1.2, font.lab = 2, cex.axis = 1.3)

### do eBAC
par(mfrow = c(3, 3))

# median RT
rtEBAC <- lm(data$diffMedianRT~data$eBAC)
cor <- cor(data$diffMedianRT, data$eBAC)
f <- summary(rtEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p

if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$diffMedianRT, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Median RT", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(rtEBAC, lwd = 4, col = "red", lty = 2)

# sd RT
sdEBAC <- lm(data$sdDiff~data$eBAC)
cor <- cor(data$sdDiff, data$eBAC)
f <- summary(sdEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p

if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$sdDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "SD (RT)", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(sdEBAC, lwd = 4, col = "red", lty = 2)

# accuracy
accEBAC <- lm(data$accDiff~data$eBAC)
cor <- cor(data$accDiff, data$eBAC)
f <- summary(accEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$accDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Accuracy", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(accEBAC, lwd = 4, col = "red", lty = 2)

# mu
muEBAC <- lm(data$muDiff~data$eBAC)
cor <- cor(data$muDiff, data$eBAC)
f <- summary(muEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$muDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Mu", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(muEBAC, lwd = 4, col = "red", lty = 2)

# sigma
sigmaEBAC <- lm(data$sigmaDiff~data$eBAC)
cor <- cor(data$sigmaDiff, data$eBAC)
f <- summary(sigmaEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$sigmaDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Sigma", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(sigmaEBAC, lwd = 4, col = "red", lty = 2)

# tau
tauEBAC <- lm(data$tauDiff~data$eBAC)
cor <- cor(data$tauDiff, data$eBAC)
f <- summary(tauEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$eBAC, data$tauDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Tau", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(tauEBAC, lwd = 4, col = "red", lty = 2)


## EZ model fits

# first get rid of subjects who were not fit by the EZ model
temp <- subset(data, data$boundaryDiff != "NA")

# drift
driftEBAC <- lm(temp$driftDiff~temp$eBAC)
cor <- cor(temp$driftDiff, temp$eBAC)
f <- summary(driftEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$eBAC, temp$driftDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Drift", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(driftEBAC, lwd = 4, col = "red", lty = 2)

# boundary
boundaryEBAC <- lm(temp$boundaryDiff~temp$eBAC)
cor <- cor(temp$boundaryDiff, temp$eBAC)
f <- summary(boundaryEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$eBAC, temp$boundaryDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Boundary", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(boundaryEBAC, lwd = 4, col = "red", lty = 2)

# nonDecision
nonEBAC <- lm(temp$nonDecisionDiff~temp$eBAC)
cor <- cor(temp$nonDecisionDiff, temp$eBAC)
f <- summary(nonEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$eBAC, temp$nonDecisionDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "eBAC", ylab = "Non-Decision", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(nonEBAC, lwd = 4, col = "red", lty = 2)

------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Severity correlations
data <- read.csv("finalHangoverInfo.csv", header = TRUE)
data <- subset(data, data$HO_AHS != "NA")

par(mfrow = c(3, 3))

# median RT
rtEBAC <- lm(data$diffMedianRT~data$HO_AHS)
cor <- cor(data$diffMedianRT, data$HO_AHS)
f <- summary(rtEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$diffMedianRT, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Median RT", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(rtEBAC, lwd = 4, col = "red", lty = 2)

# sd RT
sdEBAC <- lm(data$sdDiff~data$HO_AHS)
cor <- cor(data$sdDiff, data$HO_AHS)
f <- summary(sdEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$sdDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "SD (RT)", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(sdEBAC, lwd = 4, col = "red", lty = 2)

# accuracy
accEBAC <- lm(data$accDiff~data$HO_AHS)
cor <- cor(data$accDiff, data$HO_AHS)
f <- summary(accEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$accDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Accuracy", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(accEBAC, lwd = 4, col = "red", lty = 2)

# mu
muEBAC <- lm(data$muDiff~data$HO_AHS)
cor <- cor(data$muDiff, data$HO_AHS)
f <- summary(muEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$muDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Mu", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(muEBAC, lwd = 4, col = "red", lty = 2)

# sigma
sigmaEBAC <- lm(data$sigmaDiff~data$HO_AHS)
cor <- cor(data$sigmaDiff, data$HO_AHS)
f <- summary(muEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$sigmaDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Sigma", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(sigmaEBAC, lwd = 4, col = "red", lty = 2)

# tau
tauEBAC <- lm(data$tauDiff~data$HO_AHS)
cor <- cor(data$tauDiff, data$HO_AHS)
f <- summary(tauEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(data$HO_AHS, data$tauDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Tau", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(tauEBAC, lwd = 4, col = "red", lty = 2)


## EZ model fits

# first get rid of subjects who were not fit by the EZ model
temp <- subset(data, data$boundaryDiff != "NA")

# drift
driftEBAC <- lm(temp$driftDiff~temp$HO_AHS)
cor <- cor(data$tauDiff, data$HO_AHS)
f <- summary(driftEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$HO_AHS, temp$driftDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Drift", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(driftEBAC, lwd = 4, col = "red", lty = 2)

# boundary
boundaryEBAC <- lm(temp$boundaryDiff~temp$HO_AHS)
cor <- cor(temp$boundaryDiff, temp$HO_AHS)
f <- summary(boundaryEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$HO_AHS, temp$boundaryDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Boundary", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(boundaryEBAC, lwd = 4, col = "red", lty = 2)

# nonDecision
nonEBAC <- lm(temp$nonDecisionDiff~temp$HO_AHS)
cor <- cor(temp$nonDecisionDiff, temp$HO_AHS)
f <- summary(nonEBAC)$fstatistic # get F
p <- unname(pf(f[1],f[2],f[3],lower.tail=F)) # get p
if(p < 0.05){
  sig = "*"
} else {
  sig = ""
}
plot(temp$HO_AHS, temp$nonDecisionDiff, col = "black", pch = 21, bg = "grey", 
     cex = 2, xlab = "AHS", ylab = "Non-Decision", 
     main = paste(round(cor, 3), sig, sep = ""))
abline(nonEBAC, lwd = 4, col = "red", lty = 2)


#------------------------------------------------------------------------------