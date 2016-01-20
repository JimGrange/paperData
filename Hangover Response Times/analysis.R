# Analysis script for Grange, J.A., et al. (under review). The effect of alcohol
#                      hangover onhangover on choice response time.

# Code written by Jim Grange, July 2015. If you have any questions about this
# code, please email grange.jim@gmail.com


#------------------------------------------------------------------------------
# clear workspace:
rm(list=ls(all=TRUE))

#set working directory
setwd("D:/Work/Research/My Papers/In Preparation/Grange et al. (Hangover & RTs)/data & code")

#------------------------------------------------------------------------------
#Load required packages

#source in custom function files
source("functions.r")

#package for ex-Gaussian fitting
if(require(retimes)==FALSE){
  install.packages("retimes", dependencies=TRUE)
}

#package for diffusion modelling
if(require(RWiener)==FALSE){
  install.packages("RWiener", dependencies=TRUE)
}

#package for Bayes Factor analysis
if(require(BayesFactor)==FALSE){
  install.packages("BayesFactor", dependencies=TRUE)
}

#package for changing data from wide format to long format
if(require(reshape2)==FALSE){
  install.packages("reshape2", dependencies=TRUE)
}

#package for summarising data
if(require(dplyr)==FALSE){
  install.packages("dplyr", dependencies=TRUE)
}

#package for plotting
if(require(ggplot2)==FALSE){
  install.packages("ggplot2", dependencies=TRUE)
}

#package for Bayesian estimation
if(require(BEST)==FALSE){
  install.packages("BEST", dependencies=TRUE)
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#get data & set trimming details etc.
untrimmedData <- read.csv("rawData.csv", header=TRUE)
hangoverInfo <- read.csv("rawHangoverInfo.csv", header = TRUE)
conditions <- c("hungover", "control")
accCriterion <- 80 #criterion for accuracy for inclusion in analysis

#set criteria for minimum and maximum acceptable response time
minRT <- 150
maxRT <- 10000

#trim the data
data <- subset(untrimmedData, untrimmedData$RT > minRT &
                 untrimmedData$RT < maxRT)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# This section performs preliminary data sorting, by finding mean RT, median
# RT, standard deviation of RT, and percent accuracy. The data are stored in
# matrices, to be used by later functions for inferential analysis.

# find list of subjects with compelte RT data (some subjects didn't have data
# for both conditions)
completeSubjectsRT <- subjectChecks(data, conditions)

# now only select those who also have hangover data
completeSubjects <- completeSubjectsRT[completeSubjectsRT %in% 
                                         hangoverInfo$HO_Participant_N]

# remove participants from hangover data who don't have RT data
hangoverInfo <- hangoverInfo[which(hangoverInfo$HO_Participant_N %in% completeSubjects), ]

# get matrix of mean correct response times
meanRTs <- getRTs(data, completeSubjects, decimals = 0)

# get matrix of median correct response times
medianRTs <- getMedianRTs(data, completeSubjects)

# get standard deviation of correct RTs
sdRTs <- getSDs(data, completeSubjects)

# get matrix of Errors
meanErrors <- getErrors(data, completeSubjects)

# find position of subjects who don't meet accuracy criterion
accExclusion <- accuracyTrimming(meanErrors, accCriterion)

# find subjects who need to be removed
exclusion <-  accExclusion

# remove subject 66 due to VERY high difference score (position 14)
exclusion <- sort(c(which(meanRTs[, 1] == 66), exclusion))


# subjects who need to be removed due to BMI > 30
bmi <- hangoverInfo$HO_Participant_N[which(hangoverInfo$BMI > 30)]

# this logs the POSITION of the subjects to remove
exclusion <- sort(unique(c(exclusion, which(meanRTs[, 1] %in% bmi))))

# this logs the SUBJECT NUMBER of the subjects to remove
exclusionID <- meanRTs[, 1][exclusion]

# now, remove from the final data matrices those subjects who need to be
# removed, and---if you remove the # before "write.csv" ---save the data files
# to .csv files for storage

# get the final hangover Information database
hangoverExclusion <- which(hangoverInfo$HO_Participant_N %in% exclusionID)
finalHangoverInfo <- hangoverInfo[-hangoverExclusion, ]
write.csv(finalHangoverInfo, "finalHangoverInfo.csv", row.names = FALSE)

# get final mean RT data
finalMeanRT <- meanRTs[-exclusion, ]
finalMeanRT <- data.frame(finalMeanRT)
finalMeanRT$subject <- factor(finalMeanRT$subject)
write.csv(finalMeanRT, "finalMeanRT.csv")

# get final median RT data
finalMedianRT <- medianRTs[-exclusion, ]
finalMedianRT <- data.frame(finalMedianRT)
finalMedianRT$subject <- factor(finalMedianRT$subject)
write.csv(finalMedianRT, "finalMedianRT.csv")

# get final SD data
finalSDs <- sdRTs[-exclusion, ]
finalSDs <- data.frame(finalSDs)
finalSDs$subject <- factor(finalSDs$subject)
write.csv(finalSDs, "finalSD.csv")

# get final Error data
finalError <- meanErrors[-exclusion, ]
finalError <- data.frame(finalError)
finalError$subject <- factor(finalError$subject)
write.csv(finalError, "finalError.csv")

# get final hangover information file
finalInfo <- hangoverInfo[-exclusion, ]
finalInfo <- data.frame(finalInfo)
# finalInfo$HO_Participant_N <- factor(finalInfo$HO_Participant_N)

# who are the complete subjects?
completeSubjects <- finalError[, 1]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# do some summaries (if you want to later plot them)
meanRT <- doSummary(finalMeanRT, id.vars=c("subject"),
                    variable.name = "Condition", value.name = "RT")

medianRT <- doSummary(finalMedianRT, id.vars=c("subject"),
                      variable.name = "Condition", value.name = "RT")

sdRT <- doSummary(finalSDs, id.vars=c("subject"), variable.name = "Condition",
                  value.name = "SD")

accuracy <- doSummary(finalError, id.vars=c("subject"),
                      variable.name = "Condition", value.name = "Accuracy")
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# now perform the inferential statistical tests. This follows two stages:
# 1) Standard paired-samples t-tests (together with Cohen's d);
# 2) Bayes Factors


# 1) Standard t-tests
tMeanRT <- t.test(finalMeanRT[, 2], finalMeanRT[, 3], paired = TRUE)
tMedianRT <- t.test(finalMedianRT[, 2], finalMedianRT[, 3], paired = TRUE)
tSDs <- t.test(finalSDs[, 2], finalSDs[, 3], paired = TRUE)
tError <- t.test(finalError[, 2], finalError[, 3], paired = TRUE)

#1a) Find Cohen's d of the differences
dMeanRT <- mean(finalMeanRT[, 2] - finalMeanRT[, 3]) / sd(finalMeanRT[, 2]
                                                          - finalMeanRT[, 3])

dMedianRT <- mean(finalMedianRT[, 2] - finalMedianRT[, 3]) /
  sd(finalMedianRT[, 2] - finalMedianRT[, 3])

dSDs <- mean(finalSDs[, 2] - finalSDs[, 3]) / sd(finalSDs[, 2] - finalSDs[, 3])

dError <- mean(finalError[, 2] - finalError[, 3]) / sd(finalError[, 2] -
                                                         finalError[, 3])

#2) Calculate Bayes Factors for the differenct dependent variables
bfMeanRT <- ttest.tstat(t = as.numeric(tMeanRT$statistic),
                        n1 = length(completeSubjects))
bfMeanRT <- exp(bfMeanRT[["bf"]]) #take the exponential

bfMedianRT <- ttest.tstat(t = as.numeric(tMedianRT$statistic),
                          n1 = length(completeSubjects))
bfMedianRT <- exp(bfMedianRT[["bf"]])

bfSDs <- ttest.tstat(t = as.numeric(tSDs$statistic),
                     n1 = length(completeSubjects))
bfSDs <- exp(bfSDs[["bf"]])

bfError <- ttest.tstat(t = as.numeric(tError$statistic),
                       n1 = length(completeSubjects))
bfError <- exp(bfError[["bf"]])
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# add correlations with aBIC?

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Conduct the ex-Gaussian analysis
# NOTE: This can take quite a while to run!

# initiate an empty matrix to store all data
exG <- matrix(ncol = 8, nrow = length(completeSubjects))
colnames(exG) = c("Subject", "muHangover", "sigmaHangover", "tauHangover",
                  "Subject", "muControl", "sigmaControl", "tauControl")

# get estimates for hangover data
exgHangover <- subset(data, data$Condition == "hungover")
exG[, 1:4] <- doExG(exgHangover, completeSubjects)

# get estimates for control data
exgControl <- subset(data, data$Condition == "control")
exG[, 5:8] <- doExG(exgControl, completeSubjects)

# delete the extra subject column (just tidying up the data matrix)
exG <- exG[, -5]

# save the estimates to a csv file
write.csv(exG, "ExGaussian Parameters.csv", row.names = FALSE)



### now get some summaries
exG <- read.csv("ExGaussian Parameters.csv", header = TRUE)

mu <- cbind(exG[, 1], exG[, 2], exG[, 5])
colnames(mu) <- c("subject", "hangover", "control")
mu <- data.frame(mu)
mu <- doSummary(mu, id.vars=c("subject"), variable.name = "Condition",
                value.name = "Mu")

sigma <- cbind(exG[, 1], exG[, 3], exG[, 6])
colnames(sigma) <- c("subject", "hangover", "control")
sigma <- data.frame(sigma)
sigma <- doSummary(sigma, id.vars=c("subject"), variable.name = "Condition",
                   value.name = "Sigma")

tau <- cbind(exG[, 1], exG[, 4], exG[, 7])
colnames(tau) <- c("subject", "hangover", "control")
tau <- data.frame(tau)
tau <- doSummary(tau, id.vars=c("subject"), variable.name = "Condition",
                 value.name = "Tau")


#1) do t-tests
muT <- t.test(exG[, 2], exG[, 5], paired = TRUE)
sigmaT <- t.test(exG[, 3], exG[, 6], paired = TRUE)
tauT <- t.test(exG[, 4], exG[, 7], paired = TRUE)

#1a) find effect sizes
dMu <- mean(exG[, 2] - exG[, 5]) / sd(exG[, 2] - exG[, 5])
dSigma <- mean(exG[, 3] - exG[, 6]) / sd(exG[, 3] - exG[, 6])
dTau <- mean(exG[, 4] - exG[, 7]) / sd(exG[, 4] - exG[, 7])

#2) do Bayes factors
bfMu <- ttest.tstat(t = as.numeric(muT$statistic), n1 = nrow(exG))
bfMu <- exp(bfMu[["bf"]]) #take the exponential

bfSigma <- ttest.tstat(t = as.numeric(sigmaT$statistic), n1 = nrow(exG))
bfSigma <- exp(bfSigma[["bf"]]) #take the exponential

bfTau <- ttest.tstat(t = as.numeric(tauT$statistic), n1 = nrow(exG))
bfTau <- exp(bfTau[["bf"]]) #take the exponential
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# add correlation between differences and aBIC?

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#do EZ diffusion modelling!

oldData <- data # transfer data to a new variable name, as we are going to
                # change the data in the variable "data"


# Preliminary analysis (makes use of functions file)
data$RT <- data$RT / 1000 # change the RTs to seconds, rather than milliseconds

#get matrix of RTs, variance, and proportion correct
ezRTs <- getRTs(data, as.numeric((as.character(completeSubjects))), decimals = 3)
ezVariance <- getVariance(data, as.numeric((as.character(completeSubjects))))
ezError <- finalError
ezError[, 2:3] <- round((ezError[, 2:3] / 100), digits=3)


#now we can do the EZ modelling!
#NOTE: This function requires MEAN response time (we have been using median up
#      until now, but it was calculated earlier. It's stored in the variable
#      finalMeanRT)
ez <- doEZ(ezRTs, ezVariance, ezError, completeSubjects)

#now we need to remove those subjects who have negative parameter estimates
ezRemoval <- c(which(ez[, 2]<0), which(ez[, 3]<0), which(ez[, 4]<0),
               which(ez[, 5]<0), which(ez[, 6]<0), which(ez[, 7]<0))
ezRemoval <- unique(ezRemoval)


##HOW MANY PEOPLE REMVOED??
nEZremoval <- length(ezRemoval)

#remove them from the set of parameter estimates
finalEZ <- ez[-ezRemoval, ]

#save the data to a csv file
write.csv(finalEZ, "finalEZ.csv", row.names = FALSE)



#Now produce summaries of the parameter estimates

drift <- cbind(finalEZ[, 1], finalEZ[, 2], finalEZ[, 5])
colnames(drift) <- c("Subject", "hangover", "control")
drift <- data.frame(drift)
drift <- doSummary(drift, id.vars=c("Subject"), variable.name = "Condition",
                   value.name = "Drift")


boundary <- cbind(finalEZ[, 1], finalEZ[, 3], finalEZ[, 6])
colnames(boundary) <- c("Subject", "hangover", "control")
boundary <- data.frame(boundary)
boundary <- doSummary(boundary, id.vars=c("Subject"),
                      variable.name = "Condition", value.name = "Boundary")


nonDecision <- cbind(finalEZ[, 1], finalEZ[, 4], finalEZ[, 7])
colnames(nonDecision) <- c("Subject", "hangover", "control")
nonDecision <- data.frame(nonDecision)
nonDecision <- doSummary(nonDecision, id.vars=c("Subject"),
                         variable.name = "Condition",
                         value.name = "nonDecision")



#now do the analysis

#1) do t-tests
driftT <- t.test(finalEZ[, 2], finalEZ[, 5], paired = TRUE)
boundaryT <- t.test(finalEZ[, 3], finalEZ[, 6], paired = TRUE)
nonDecisionT <- t.test(finalEZ[, 4], finalEZ[, 7], paired = TRUE)

#1a) find effect sizes
dDrift <- mean(finalEZ[, 2] - finalEZ[, 5]) / sd(finalEZ[, 2] - finalEZ[, 5])
dBoundary <- mean(finalEZ[, 3] - finalEZ[, 6]) / sd(finalEZ[, 3]
                                                    - finalEZ[, 6])
dNonDecision <- mean(finalEZ[, 4] - finalEZ[, 7]) / sd(finalEZ[, 4]
                                                       - finalEZ[, 7])


#2) do Bayes factors
bfDrift <- ttest.tstat(t = as.numeric(driftT$statistic), n1 = nrow(finalEZ))
bfDrift <- exp(bfDrift[["bf"]]) #take the exponential

bfBoundary <- ttest.tstat(t = as.numeric(boundaryT$statistic),
                          n1 = nrow(finalEZ))
bfBoundary <- exp(bfBoundary[["bf"]]) #take the exponential

bfNonDecision <- ttest.tstat(t = as.numeric(nonDecisionT$statistic),
                             n1 = nrow(finalEZ), rscale = 0.707)
bfNonDecision <- exp(bfNonDecision[["bf"]]) #take the exponential
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# do correlation between EZ diffusion & eBAC
ezInfo <- finalInfo[-ezRemoval, ]
bac <- ezInfo$eBAC

driftDiff <- finalEZ[, 2] - finalEZ[, 5]
boundDiff <- finalEZ[, 3] - finalEZ[, 6]
nonDiff <- finalEZ[, 4] - finalEZ[, 7]

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Robust EZ analysis
# NOTE: This creates A LOT of text files.

# first generate new files, one for each subject and condition
for(currSub in completeSubjects){

  getData <- subset(data, data$Subject == currSub)

  for(cond in conditions){

    condData <- subset(getData, getData$Condition == cond)
    nTrials <- length(condData$RT)
    temp <- subset(condData, condData$Accuracy == 1)
    finalData <- matrix(0, nrow = length(temp$RT) + 1, ncol = 1)
    finalData[2:nrow(finalData), 1] <- temp$RT
    finalData[1, 1] <- nTrials

    filename <- paste(cond, currSub, ".txt", sep = "")
    write(finalData, file=filename)
  }
}

#DO Robust EZ MODELLING

# The code assumes that the (unique) extension .txt for data files
# The results are written to the file "REZBatch.out"

output    = "REZBatch.out"
datafiles = list.files(pattern="\\.txt$") #selects all files that end with .txt

for (i in 1:length(datafiles)) # loop over data files
{
  # For each file, call the Robust EZ routines and write the results to the
  # output file:
  write(c(datafiles[i], RobustEZ.from.File(datafiles[i])), output,
        ncolumns=5, append = TRUE)
  # Note that append = TRUE, so running this code again (without deleting the
  # output file first) will result in a new output file that contains not only
  # the latest analyses but also the previous analyses.
}

# source data file back in

ezData <- read.table("REZBatch.out")
colnames(ezData) <- c("Condition", "Drift", "Boundary", "NonDecision",
                      "Mixture")
ezData[, 1] <- lapply(ezData[1], as.character)

# save conditions to separate variables
controlEZ <- subset(ezData, grepl("control", ezData[, 1]))
hungoverEZ <- subset(ezData, grepl("hungover", ezData[, 1]))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#Do analysis


#1) Do t-tests
driftT <- t.test(hungoverEZ[, 2], controlEZ[, 2],paired = TRUE)
boundaryT <- t.test(hungoverEZ[, 3], controlEZ[, 3],paired = TRUE)
nonDecisionT <- t.test(hungoverEZ[, 4], controlEZ[, 4],paired = TRUE)


#2) Do Bayes Factors
bfDrift <- ttest.tstat(t = as.numeric(driftT$statistic), n1 = nrow(hungoverEZ))
bfDrift <- exp(bfDrift[["bf"]]) #take the exponential

bfBoundary <- ttest.tstat(t = as.numeric(boundaryT$statistic),
                          n1 = nrow(hungoverEZ))
bfBoundary <- exp(bfBoundary[["bf"]]) #take the exponential

bfNonDecision <- ttest.tstat(t = as.numeric(nonDecisionT$statistic),
                             n1 = nrow(hungoverEZ))
bfNonDecision <- exp(bfNonDecision[["bf"]]) #take the exponential
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# Analysis using RWiener package for diffusion modelling (see Appendix)
newData <- oldData #get the original (trimmed) data file back
  # change column headings. This is required for the RWiener package
  colnames(newData) <- c("Subject", "Condition", "resp", "q")
  newData <- data.frame(newData) #save as data frame

# sort ready for RWiener input
newData$q <- newData$q / 1000 #transpose RTs to seconds
newData$resp[newData$resp == 0] <- 2 #change errors to 2
# change resp to a factor
newData$resp = factor(newData$resp, levels = c(1,2),
                      labels = c("upper", "lower"))

# get hangover data
diffDataHangover <- subset(newData, newData$Condition == "hungover")
  diffDataHangover <- diffDataHangover[, -2] # remove condition column
  diffDataHangover <- diffDataHangover[, c(1,3,2)] # reorder data frame

# get control data
diffDataControl <- subset(newData, newData$Condition == "control")
  diffDataControl <- diffDataControl[, -2] # remove condition column
  diffDataControl <- diffDataControl[, c(1,3,2)] # reorder data frame



# run the diffusion analysis with bias fixed to 0.5
diffHangover_fixedBias <- doWiener_fixedBias(diffDataHangover,
                                             as.numeric(as.character(completeSubjects)))
  # store the best-fitting parameters
  write.csv(diffHangover_fixedBias,
            "Diffusion Parameters fixed bias (Hangover).csv", row.names = FALSE)

diffControl_fixedBias <- doWiener_fixedBias(diffDataControl, 
                                            as.numeric(as.character(completeSubjects)))
  write.csv(diffControl_fixedBias,
            "Diffusion Parmaeters fixed bias (Control).csv", row.names = FALSE)


# get the summaries
drift <- cbind(diffHangover_fixedBias[, 1], diffHangover_fixedBias[, 4],
               diffControl_fixedBias[, 4])
  colnames(drift) <- c("Subject", "hangover", "control")
  drift <- data.frame(drift)
  drift <- doSummary(drift, id.vars=c("Subject"), variable.name = "Condition",
                     value.name = "Drift")


boundary <- cbind(diffHangover_fixedBias[, 1], diffHangover_fixedBias[, 2],
                  diffControl_fixedBias[, 2])
  colnames(boundary) <- c("Subject", "hangover", "control")
  boundary <- data.frame(boundary)
  boundary <- doSummary(boundary, id.vars=c("Subject"),
                        variable.name = "Condition", value.name = "Boundary")


nonDecision <- cbind(diffHangover_fixedBias[, 1], diffHangover_fixedBias[, 3],
                     diffControl_fixedBias[, 3])
  colnames(nonDecision) <- c("Subject", "hangover", "control")
  nonDecision <- data.frame(nonDecision)
  nonDecision <- doSummary(nonDecision, id.vars=c("Subject"),
                           variable.name = "Condition",
                           value.name = "nonDecision")



# 1) Do t-tests (rw stands for "RWiener")
driftT_rw <- t.test(diffHangover_fixedBias[, 4], diffControl_fixedBias[, 4],
                    paired = TRUE)
boundaryT_rw <- t.test(diffHangover_fixedBias[, 2], diffControl_fixedBias[, 2],
                       paired = TRUE)
nonDecisionT_rw <- t.test(diffHangover_fixedBias[, 3], diffControl_fixedBias[, 3],
                          paired = TRUE)

# 1a) Get effect sizes
dDrift_rw <- mean(diffHangover_fixedBias[, 4] - diffControl_fixedBias[, 4]) /
              sd(diffHangover_fixedBias[, 4] - diffControl_fixedBias[, 4])

dBoundary_rw <- mean(diffHangover_fixedBias[, 2] - diffControl_fixedBias[, 2])/
  sd(diffHangover_fixedBias[, 2] - diffControl_fixedBias[, 2])

dNonDecision_rw <- mean(diffHangover_fixedBias[, 3] -
                          diffControl_fixedBias[, 3]) /
  sd(diffHangover_fixedBias[, 3] - diffControl_fixedBias[, 3])

# 2) Do Bayes factors
bfDrift_rw <- ttest.tstat(t = as.numeric(driftT_rw$statistic),
                          n1 = nrow(diffHangover_fixedBias))
bfDrift_rw <- exp(bfDrift_rw[["bf"]]) # take the exponential

bfBoundary_rw <- ttest.tstat(t = as.numeric(boundaryT_rw$statistic),
                             n1 = nrow(diffHangover_fixedBias))
bfBoundary_rw <- exp(bfBoundary_rw[["bf"]]) # take the exponential

bfNonDecision_rw <- ttest.tstat(t = as.numeric(nonDecisionT_rw$statistic),
                                n1 = nrow(diffHangover_fixedBias))
bfNonDecision_rw <- exp(bfNonDecision_rw[["bf"]]) # take the exponential

#------------------------------------------------------------------------------







