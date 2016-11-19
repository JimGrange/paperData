### Analysis script for Agi's experiment (Mayr replication)

#------------------------------------------------------------------------------
### General set-up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Documents/Git/paperData/Episodic Retrieval & Inhibition/Data & Analysis Code/Experiment 1. Mayr Replication/LMM")

# load necessary functions file & load necessary packages
library(dplyr)
library(lme4)
source("functions.R")


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

## format the data frame suitable for trimr processing

# add a new column called condition
data$condition <- paste(data$sequence, "_", data$stimRep, sep = "")

# change subject to participant
names(data)[1] <- "participant"

# trim the RTs with fixed standard deviation upper limit
rtData <- doTrim(data = data, minRT = 150, sd = 2.5)

# do the LMM
lm_out <- lmer(rt ~ sequence * stimRep + (sequence|participant) +
                 (stimRep|participant), data = rtData)
