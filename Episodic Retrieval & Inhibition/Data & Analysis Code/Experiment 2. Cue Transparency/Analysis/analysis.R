### Analysis script for Agi's experiment (Mayr replication)

#------------------------------------------------------------------------------
### General set-up

# clear workspace
rm(list = ls())

# set working directory
setwd("~/Git/paperData/Episodic Retrieval & Inhibition/Data & Analysis Code/Experiment 2. Cue Transparency/Analysis")

# load necessary functions file & load necessary packages
source("functions.R")
library(dplyr)
library(ggplot2)
library(ez)
library(BayesFactor)


# load the data files
arrowData <- read.csv("raw_arrows.csv", header = TRUE)
shapeData <- read.csv("raw_shapes.csv", header = TRUE)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### do some processing of the raw files to log for n-2 response repetition

## first, add a "condition" column to log what the current cue type is. Also, 
## add blank column to store response repetitions
arrowData <- arrowData %>% 
  mutate(condition = "arrows", respRep = "")

shapeData <- shapeData %>%
  mutate(condition = "shapes", respRep = "")

## code for n-2 response repetition
for(i in 3:nrow(arrowData)){
  
  if(arrowData$stimulus.CRESP[i] == arrowData$stimulus.CRESP[i - 2]){
    arrowData$respRep[i] <- "yes"
  } else {
    arrowData$respRep[i] <- "no"
  }
  
}

for(i in 3:nrow(shapeData)){
  
  if(shapeData$stimulus.CRESP[i] == shapeData$stimulus.CRESP[i - 2]){
    shapeData$respRep[i] <- "yes"
  } else {
    shapeData$respRep[i] <- "no"
  }
  
}

colnames(arrowData) <- c("subject", "trial", "sequence", "accuracy", "CRESP", 
                         "rt", "condition", "respRep")
arrowData <- arrowData %>%
  select(-CRESP)

colnames(shapeData) <- c("subject", "trial", "sequence", "accuracy", "CRESP", 
                         "rt", "condition", "respRep")
shapeData <- shapeData %>%
  select(-CRESP)


## bind the data together
fullData <- rbind(arrowData, shapeData)
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

# summarise the main effects

mainCondition <- accuracy %>%
  group_by(condition) %>%
  summarise(meanACC = mean(acc))

mainRespRep <- accuracy %>%
  group_by(respRep) %>%
  summarise(meanACC = mean(acc))

intSequenceRespRep <- accuracy %>%
  group_by(respRep, sequence) %>%
  summarise(meanACC = mean(acc))

## get all of the accuracy data for the table
# how many subjects?
nSubs <- length(unique((accuracy$subject)))

options(digits = 4)
allAccuracy <- accuracy %>%
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


# tidy some of the variable names
colnames(all)[2] <- "Response"
levels(all$Response)[levels(all$Response) == "repetition"] <- "Repetition"
levels(all$Response)[levels(all$Response) == "switch"] <- "Switch"
levels(all$sequence)[levels(all$sequence) == "aba"] <- "ABA"
levels(all$sequence)[levels(all$sequence) == "cba"] <- "CBA"


pdf("responseTimes.pdf", width = 8, height = 5)
plot <- ggplot(all, aes(x = sequence, y = meanRT, group = Response, 
                        colour = Response))
plot <- plot + geom_errorbar(aes(ymin = meanRT - SE, ymax = meanRT + SE), 
                             width = .15, size = 0.5, position = pd)
plot <- plot + geom_line(aes(linetype = Response), position = pd)
plot <- plot + geom_point(aes(shape = Response), size = 2.3, position = pd)
plot <- plot + scale_x_discrete(name = "Task Sequence") + 
  scale_y_continuous(name = "Response Time (ms)")
plot <- plot + scale_shape_discrete(name = "Response") + 
  scale_linetype_discrete(name = "Response") + 
  scale_colour_discrete(name = "Response")
plot <- plot + theme(panel.background = element_rect(fill = "grey86")) 
plot + facet_grid(  ~ condition) 
dev.off()

# do the plot again, but store it to a variable
plot <- ggplot(all, aes(x = sequence, y = meanRT, group = Response, 
                        colour = Response))
plot <- plot + geom_errorbar(aes(ymin = meanRT - SE, ymax = meanRT + SE), 
                             width = .15, size = 0.5, position = pd)
plot <- plot + geom_line(aes(linetype = Response), position = pd)
plot <- plot + geom_point(aes(shape = Response), size = 2.3, position = pd)
plot <- plot + scale_x_discrete(name = "Task Sequence") + 
  scale_y_continuous(name = "Response Time (ms)")
plot <- plot + scale_shape_discrete(name = "Response") + 
  scale_linetype_discrete(name = "Response") + 
  scale_colour_discrete(name = "Response")
plot <- plot + theme(panel.background = element_rect(fill = "grey86")) 
plot_1 <- plot + facet_grid(  ~ condition)




#------------------------------------------------------------------------------
# collate data for 2 way ANOVA with n-2 repetition cost as the DV

x <- rtData %>%
  group_by(subject) %>%
  summarise(repArrows = abaRepArrows - cbaRepArrows, 
            swArrows = abaSwArrows - cbaSwArrows, 
            repShapes = abaRepShapes - cbaRepShapes, 
            swShapes = abaSwShapes - cbaSwShapes)


# get the n--2 repetition cost for response repetitions and export it
# this is used in the mini-meta-analysis in the paper
biRepArrows <- x$repArrows
write.table(biRepArrows, file = "experiment_2_respRep_arrows.csv", 
            row.names = FALSE, col.names = FALSE)

biRepShapes <- x$repShapes
write.table(biRepShapes, file = "experiment_2_respRep_shapes.csv", 
            row.names = FALSE, col.names = FALSE)


newData <- NULL

for(i in 1:nrow(x)){
  
  # matrix for each subject's data
  subData <- matrix(0, nrow = 4, ncol = 4)
  subData <- data.frame(subData)
  colnames(subData) <- c("subject", "condition", "respRep", "rt")
  
  # populate the data frame
  subData[1, ] <- c(x$subject[i], "arrows", "repetition", x$repArrows[i])
  subData[2, ] <- c(x$subject[i], "arrows", "switch", x$swArrows[i])
  subData[3, ] <- c(x$subject[i], "shapes", "repetition", x$repShapes[i])
  subData[4, ] <- c(x$subject[i], "shapes", "switch", x$swShapes[i])
  
  newData <- rbind(newData, subData)
  
}

newData$subject <- factor(newData$subject)
newData$condition <- factor(newData$condition)
newData$respRep <- factor(newData$respRep)
newData$rt <- as.numeric(newData$rt)


## Bayes factor
set.seed(65)
bf <- anovaBF(rt ~ condition * respRep + subject, data = newData,
              whichRandom = "subject")
bf <- recompute(bf, iterations = 100000)
bf



plotData2 <- newData %>%
  group_by(condition, respRep) %>%
  summarise(meanRT = mean(rt), se = sd(rt) / sqrt(nSubs))


colnames(plotData2)[2] <- "Response"
levels(plotData2$Response)[levels(plotData2$Response) == "repetition"] <- "Repetition"
levels(plotData2$Response)[levels(plotData2$Response) == "switch"] <- "Switch"


pdf("Experiment2.pdf", width = 8, height = 8)
p <- ggplot(plotData2, aes(x = condition, y = meanRT, group = Response, 
                         colour = Response))
p<- p + geom_errorbar(aes(ymin = meanRT - se, ymax = meanRT + se), 
                      width = 0.05, size = 0.5, position = pd)
p <- p + geom_line(aes(linetype  = Response), position = pd)
p <- p + geom_point(aes(shape = Response), size = 3, position = pd)
p <- p + labs(x = "Cue Type", y = "N-2 Repetition Cost (ms)")
p <- p + theme(panel.background = element_rect(fill = "grey86"))
p
dev.off()

# Get the response time plot and the n--2 repetition cost plot stitched
library(gridExtra)
pdf("all_rts_Experiment2.pdf", width = 8, height = 8)
grid.arrange(plot_1, p)
dev.off()


# log RT analysis--------------------------------------------------------------
new <- finalRT
new$rt <- log(new$rt)

# ANOVA using ezANOVA
aovRTLog <- ezANOVA(
  data = new
  , dv = .(rt)
  , wid = .(subject)
  , within = .(condition, sequence, respRep)
  , between = NULL
  , detailed = FALSE)


newPlot <- all
newPlot$meanRT <- log(newPlot$meanRT)
newPlot$SE <- log(newPlot$SE)

pdf("LOGresponseTimes.pdf", width = 8, height = 5)
plot <- ggplot(newPlot, aes(x = sequence, y = meanRT, group = Response, colour = Response))
plot <- plot + geom_line(aes(linetype = Response), position = pd)
plot <- plot + geom_point(aes(shape = Response), size = 2.3, position = pd)
plot <- plot + scale_x_discrete(name = "Task Sequence") + scale_y_continuous(name = "Response Time (ms)")
plot <- plot + scale_shape_discrete(name = "Response") + 
  scale_linetype_discrete(name = "Response") + 
  scale_colour_discrete(name = "Response")
plot <- plot + theme(panel.background = element_rect(fill = "grey86")) 
plot + facet_grid(  ~ condition) 
dev.off()
