### Analysis script for episodic retrieval & stimulus repetition

#------------------------------------------------------------------------------
### General set-up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Git/paperData/Episodic Retrieval & Inhibition/Data & Analysis Code/Experiment 3. Stimulus Matching")

# load necessary functions file & load necessary packages
source("functions.R")
library(dplyr)
library(ggplot2)
library(ez)
library(BayesFactor)

# load the data files
match <- read.csv("match_data.csv")
mismatch <- read.csv("mismatch_data.csv")
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### do preparation of data files

## code for n-2 response repetition
for(i in 3:nrow(match)){
  
  if(match$stimulus.CRESP[i] == match$stimulus.CRESP[i - 2]){
    match$respRep[i] <- "yes"
  } else {
    match$respRep[i] <- "no"
  }
}

for(i in 3:nrow(mismatch)){

  if(mismatch$stimulus.CRESP[i] == mismatch$stimulus.CRESP[i - 2]){
    mismatch$respRep[i] <- "yes"
  } else {
    mismatch$respRep[i] <- "no"
  }
}

# remove stimulus repetitions from n-1 and n-2 from mismatch condition
mismatch <- mismatch %>%
  filter(stimN1 == "switch", 
         stimN2 == "switch")

# remove stimN1 and stimN2 from mismatch condition so both data files have
# identical column structure
mismatch <- mismatch %>%
  select(-stimN1, -stimN2)

# remove CRESP column from each data frame
match <- select(match, -stimulus.CRESP)
mismatch <- select(mismatch, -stimulus.CRESP)

# add condition column to each data frame. Then re-order columns
match <- match %>%
  mutate(condition = "match")
match <- match[c(7, 1, 2, 3, 6, 4, 5)]

mismatch <- mismatch %>%
  mutate(condition = "mismatch")
mismatch <- mismatch[c(7, 1, 2, 3, 6, 4, 5)]

# merge data frame
fullData <- rbind(match, mismatch)

# change column names 
colnames(fullData) <- c("condition", "subject", "trial", "sequence", "respRep", 
                        "accuracy", "rt")

# save data frame
write.csv(fullData, "fullData", row.names = FALSE)

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Error Data Trimming 

# Remove null trials (1st and 2nd trial of each block)

# which trials were null?
nullTrials <- c(1, 2, 121, 122, 241, 242)

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
  group_by(subject, condition, respRep, sequence) %>%
  summarise(acc = (sum(accuracy) / length(accuracy)))

# now remove all error trials, so we are ready for RT analysis
data <- subset(data, data$accuracy == 1)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### response time analysis

# trim the RTs with fixed standard deviation upper limit
rtData <- getRTs(data = data, minRT = 150, sd = 2.5)

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
  , within = .(condition, sequence, respRep)
  , between = NULL 
  , detailed = FALSE)


mainSequence <- accuracy %>%
  group_by(sequence) %>%
  summarise(meanAcc = mean(acc))


# how many subjects?
nSubs <- length(unique((finalRT$subject)))

# summarise the main effects
options(digits = 4)
finalAcc <- accuracy %>%
  group_by(condition, respRep, sequence) %>%
  summarise(accMean = mean(acc) * 100, SE = round(sd(acc) / sqrt(nSubs), 4) * 100)

#--- response time
# transfer to long format for analysis & later plotting, & change necessary
# columns to factors
finalRT <- doLong(rtData)
finalRT$condition <- factor(finalRT$condition)
finalRT$sequence <- factor(finalRT$sequence)
finalRT$respRep <- factor(finalRT$respRep)
finalRT$subject <- factor(finalRT$subject)


# summarise the response times
rtSummary <- finalRT %>%
  group_by(condition, respRep, sequence) %>%
  summarise(meanRT = mean(rt))


# ANOVA using ezANOVA
aovRT <- ezANOVA(
  data = finalRT
  , dv = .(rt)
  , wid = .(subject)
  , within = .(condition, sequence, respRep)
  , between = NULL
  , detailed = FALSE)


# summarise the response times

# how many subjects?
nSubs <- length(unique((finalRT$subject)))


mainCondition <- finalRT %>%
  group_by(condition) %>%
  summarise(meanRT = mean(rt))

mainSequence <- finalRT %>%
  group_by(sequence) %>%
  summarise(meanRT = mean(rt))

mainStimRep <- finalRT %>%
  group_by(respRep) %>%
  summarise(meanRT = mean(rt))

conditionSequence <- finalRT %>%
  group_by(condition, sequence) %>%
  summarise(meanRT = mean(rt))

conditionStimRep <- finalRT %>%
  group_by(condition, respRep) %>%
  summarise(meanRT = mean(rt))

sequenceStimRep <- finalRT %>%
  group_by(respRep, sequence) %>%
  summarise(meanRT = mean(rt))

all <- finalRT %>%
  group_by(condition, respRep, sequence) %>%
  summarise(meanRT = mean(rt), SE = round(sd(rt) / sqrt(nSubs)))
all$SE  <- as.numeric(all$SE)


# do plot

pd <- position_dodge(0.08)
theme_set(theme_gray(base_size = 18))

levels(all$respRep)[levels(all$respRep) == "repetition"] <- "Repetition"
levels(all$respRep)[levels(all$respRep) == "switch"] <- "Switch"
levels(all$sequence)[levels(all$sequence) == "aba"] <- "ABA"
levels(all$sequence)[levels(all$sequence) == "cba"] <- "CBA"


pdf("responseTimes.pdf", width = 8, height = 5)
plot <- ggplot(all, aes(x = sequence, y = meanRT, group = respRep, colour = respRep))
plot <- plot + geom_errorbar(aes(ymin = meanRT - SE, ymax = meanRT + SE), width = .15, 
                             size = 0.5, position = pd)
plot <- plot + geom_line(aes(linetype = respRep), position = pd)
plot <- plot + geom_point(aes(shape = respRep), size = 2.3, position = pd)
plot <- plot + scale_x_discrete(name = "Task Sequence") + scale_y_continuous(name = "Response Time (ms)")
plot <- plot + scale_shape_discrete(name = "Response") + 
  scale_linetype_discrete(name = "Response") + 
  scale_colour_discrete(name = "Response")
plot <- plot + theme(panel.background = element_rect(fill = "grey86")) 
plot_1 <- plot + facet_grid(  ~ condition) 
plot_1
dev.off()


#------------------------------------------------------------------------------
# collate data for 2 way ANOVA with n-2 repetition cost as the DV

x <- rtData %>%
  group_by(subject) %>%
  summarise(repMatch = abaRepMatch - cbaRepMatch, 
            swMatch = abaSwMatch - cbaSwMatch, 
            repMismatch = abaRepMismatch - cbaRepMismatch, 
            swMismatch = abaSwMismatch - cbaSwMismatch)


# get the n--2 repetition cost for response repetitions and export it
# this is used in the mini-meta-analysis in the paper
biRepMatch <- x$repMatch
write.table(biRepMatch, file = "experiment_3_respRep_match.csv", 
            row.names = FALSE, col.names = FALSE)

biRepMismatch <- x$repMismatch
write.table(biRepMismatch, file = "experiment_3_respRep_mismatch.csv", 
            row.names = FALSE, col.names = FALSE)



# put into long format
newData <- NULL

for(i in 1:nrow(x)){
  
  # matrix for each subject's data
  subData <- matrix(0, nrow = 4, ncol = 4)
  subData <- data.frame(subData)
  colnames(subData) <- c("subject", "condition", "respRep", "rt")
  
  # populate the data frame
  subData[1, ] <- c(x$subject[i], "match", "repetition", x$repMatch[i])
  subData[2, ] <- c(x$subject[i], "match", "switch", x$swMatch[i])
  subData[3, ] <- c(x$subject[i], "mismatch", "repetition", x$repMismatch[i])
  subData[4, ] <- c(x$subject[i], "mismatch", "switch", x$swMismatch[i])
  
  newData <- rbind(newData, subData)
}

newData$subject <- factor(newData$subject)
newData$condition <- factor(newData$condition)
newData$respRep <- factor(newData$respRep)
newData$rt <- as.numeric(newData$rt)


## calculate the Bayes factor

# first, set the random seed so the BF is reproducible
set.seed(65)
bf <- anovaBF(rt ~ condition * respRep + subject, data = newData,
              whichRandom = "subject")
bf <- recompute(bf, iterations = 100000)
bf
bf[2] / bf[4]

plotData2 <- newData %>%
  group_by(condition, respRep) %>%
  summarise(meanRT = mean(rt), se = sd(rt) / sqrt(nSubs))
colnames(plotData2)[2] <- "Response"
levels(plotData2$Response)[levels(plotData2$Response) == "repetition"] <- "Repetition"
levels(plotData2$Response)[levels(plotData2$Response) == "switch"] <- "Switch"


pdf("n2 cost.pdf", width = 8, height = 8)
p <- ggplot(plotData2, aes(x = condition, y = meanRT, group = Response, 
                           colour = Response))
p<- p + geom_errorbar(aes(ymin = meanRT - se, ymax = meanRT + se), 
                      width = 0.05, size = 0.5, position = pd)
p <- p + geom_line(aes(linetype  = Response), position = pd)
p <- p + geom_point(aes(shape = Response), size = 3, position = pd)
p <- p + labs(x = "Stimulus Condition", y = "N-2 Repetition Cost (ms)")
p <- p + theme(panel.background = element_rect(fill = "grey86"))
p
dev.off()



# Get the response time plot and the n--2 repetition cost plot stitched
library(gridExtra)
pdf("all_rts_Experiment3.pdf", width = 8, height = 8)
grid.arrange(plot_1, p)
dev.off()
#------------------------------------------------------------------------------

