#------------------------------------------------------------------------------
# load packages & source files

rm(list = ls())
setwd("~/Git/paperData/2015/Temporal Distinctiveness & Mixture Distributions/Grange (under review) re-analysis")
source("functions.R")

# package for changing data from wide format to long format
if(require(reshape2)==FALSE){
  install.packages("reshape2", dependencies=TRUE)
}

# package for summarising data
if(require(dplyr)==FALSE){
  install.packages("dplyr", dependencies=TRUE)
}


#package for plotting
if(require(ggplot2)==FALSE){
  install.packages("ggplot2", dependencies=TRUE)
}


#package for ANOVA
if(require(ez)==FALSE){
  install.packages("ez", dependencies=TRUE)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# source the data & declare trimming details etc.
rawData <- read.csv("allData.csv")
colnames(rawData) <- c("subject", "trial", "condition", "task", "sequence",
                       "currentRCI", "previousRCI", "rciChange", "accuracy",
                       "rt")

# add column to store error trimming details (to remove trials afer error)
rawData <- mutate(rawData, accuracyTrim = 0)

# only select the pertinent columns
data <- select(rawData, subject, condition, sequence, currentRCI, rciChange, 
               accuracy, accuracyTrim, rt)

# remove null trials
data <- subset(data, data$sequence != "null")

# add the accuracy trimming column (removing trial after error)
data <- addError(data)

# what is the accuracy threshold for inclusion
accCriterion <- 80

# what is the minimum RT allowed?
minRT <- 150

# how many SDs to trim RTs to?
sdCriterion <- 2.5

# how many trials per participant in the whole condition?
nTrials <- 384
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### do accuracy checks

# perceptual (location) condition
locationData <- subset(data, data$condition == "Perceptual")
locationACC <- data.frame(getAccuracy(locationData))
locationACC <- mutate(locationACC, avAccuracy = 0)
for(i in 1:nrow(locationACC)){
  locationACC$avAccuracy[i] = mean(as.numeric(locationACC[i, 2:9]))
}
#which participants need to be removed?
locationRemoved <- locationACC$subject[which
                                       (locationACC$avAccuracy < accCriterion)]



# Memory (parity) condition
parityData <- subset(data, data$condition == "Memory")
parityACC <- data.frame(getAccuracy(parityData))
parityACC <- mutate(parityACC, avAccuracy = 0)
for(i in 1:nrow(parityACC)){
  parityACC$avAccuracy[i] = mean(as.numeric(parityACC[i, 2:9]))
}

#which participants need to be removed?
parityRemoved <- parityACC$subject[which
                                   (parityACC$avAccuracy < accCriterion)]



#list of all subjects to remove
removedSubjects <- unique(c(locationRemoved, parityRemoved))

# get the positions of the subjects to removed
parityRemoved_position <- which(parityACC$subject %in% removedSubjects)
locationRemoved_position <- which(locationACC$subject %in% removedSubjects)

parityACC_final <- parityACC[-parityRemoved_position, ]
parityACC_final[, 11] <- "Memory"
colnames(parityACC_final)[11] <- "condition"
locationACC_final <- locationACC[-locationRemoved_position, ]
locationACC_final[, 11] <- "Perceptual"
colnames(locationACC_final)[11] <- "condition"

finalACC <- rbind(locationACC_final, parityACC_final)
write.csv(finalACC, "finalAccuracy.csv")
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get mean response times for each condition, after trimming

# only use correct RT preceded by correct RT
errorTrimmed <- subset(data, data$accuracyTrim == 1)

# do location condition
locationData <- subset(errorTrimmed, errorTrimmed$condition == "Perceptual")
locationRT <- getRTs(locationData, minRT, sdCriterion)
locationRT <- data.frame(locationRT)
locationRT[, 10] <- "Perceptual"
colnames(locationRT)[10] <- "condition"
locationRT_final <- locationRT[-locationRemoved_position, ]

# do parity condition
parityData <- subset(errorTrimmed, errorTrimmed$condition == "Memory")
parityRT <- getRTs(parityData, minRT, sdCriterion)
parityRT <- data.frame(parityRT)
parityRT[, 10] <- "Memory"
colnames(parityRT)[10] <- "condition"
parityRT_final <- parityRT[-parityRemoved_position, ]

finalRTs <- rbind(locationRT_final, parityRT_final)
write.csv(finalRTs, "finalRTs.csv")
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Do RT analysis

# first, put the data in long format
rtsLong <- doLong(finalRTs)

# just select the "memory" data
rtsLong <- subset(rtsLong$task == "memory")

# disable scientific notation
options(scipen = 999)
options(digits = 3)

## ANVOA using ezANOVA 
aov.rt <- ezANOVA(
  data = rtsLong
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci, rcichange)
  , between = NULL
  , detailed = F
)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
## examine interaction
repData <- subset(rtsLong, rtsLong$sequence == "repeat")
swData <- subset(rtsLong, rtsLong$sequence == "switch")

repAOV <- ezANOVA(
  data = repData
  , dv = .(rt)
  , wid = .(subject)
  , within = .(rci, rcichange)
  , between = NULL
  ,detailed = F
)

swAOV <- ezANOVA(
  data = swData
  , dv = .(rt)
  , wid = .(subject)
  , within = .(rci, rcichange)
  , between = NULL
  ,detailed = F
)

#------------------------------------------------------------------------------






#------------------------------------------------------------------------------
### do some plotting


# first, put the data in long format (reshape package can do this, but I didn't
# have time to work it out yet!)
rtsLong <- doLong(finalRTs)

# get the means for each condition into a data frame
rts <- summariseRTs(rtsLong)

# re-order factor levels
rts$RCI  <- factor(rts$RCI, levels = c("50", "1000"))
rts$RCIChange <- factor(rts$RCIChange, levels = c("Same", "Different"))

#rename "location" and "parity" tasks to "Perceptual" and "Memory
levels(rts$Task)[levels(rts$Task)=="perceptual"] <- "Perceptual"
levels(rts$Task)[levels(rts$Task)=="memory"] <- "Categorisation"

# just select the "categorisation" data
rts <- subset(rts, rts$Task == "Categorisation")
rts$Task <- as.factor(rts$Task)

# to ensure error bars don't overlap
pd <- position_dodge(.02)


theme_set(theme_gray(base_size = 18))

pdf("grangeData.pdf", width = 8, height = 5)

plot <- ggplot(rts, aes(x = RCI, y = RT, group = Sequence, colour = Sequence)) 
plot <- plot + geom_errorbar(aes(ymin = RT - SE, ymax = RT + SE), width = .05, 
                             size = 0.5)
plot <- plot + geom_line(aes(linetype = Sequence))
plot <- plot + scale_linetype_manual(breaks = c("Repeat", "Switch"), values = c(5,1))
plot <- plot + geom_point(aes(shape = Sequence), size = 2.3)
plot <- plot + scale_x_discrete(name = "Response-Cue Interval (ms)") + 
  scale_y_continuous(name = "Response Time (ms)")



plot + facet_grid(  ~ RCIChange)

dev.off()
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# do separate ANOVAs for each task
memAnova <- subset(rtsLong, rtsLong$task == "memory")
perAnova <- subset(rtsLong, rtsLong$task == "perceptual")


## ANVOA for the peceptual task
aov.per <- ezANOVA(
  data = perAnova
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci, rcichange)
  , between = NULL
  , detailed = F
)

# seqeunce * rci interaction
perAnova %>%
  group_by(rci, sequence) %>%
  summarise(meanRT = mean(rt), SE = (sd(rt) / sqrt(25)))

  # test of switch RT from short to long RCI
  switchRT <- filter(perAnova, sequence == "switch")
  
  switchRT <- ezANOVA(
      data = switchRT
      , dv = .(rt)
      , wid = .(subject)
      , within = .(rci)
      , between = NULL
      , detailed = F
      )

  repRT <- filter(perAnova, sequence == "repeat")
  
  repRT <- ezANOVA(
    data = repRT
    , dv = .(rt)
    , wid = .(subject)
    , within = .(rci)
    , between = NULL
    , detailed = F
  )
  




## ANVOA for the memory task
aov.mem <- ezANOVA(
  data = memAnova
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci, rcichange)
  , between = NULL
  , detailed = F
)

# RCI-same data
same <- filter(memAnova, rcichange == "same")

#2-way interaction of sequence and RCI
sameAnova <- ezANOVA(
  data = same
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci)
  , between = NULL
  , detailed = F
)


# RCI-different data
different <- filter(memAnova, rcichange == "different")

diffAnova <- ezANOVA(
  data = different
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci)
  , between = NULL
  , detailed = F
)

  # RCI effect for repetition trials
  repRT <- filter(different, sequence == "repeat")
  repRT <- ezANOVA(
    data = repRT
    , dv = .(rt)
    , wid = .(subject)
    , within = .(rci)
    , between = NULL
    , detailed = F
  )

  # RCI effect for switch trials
  swRT <- filter(different, sequence == "switch")
  switchRT <- ezANOVA(
    data = swRT
    , dv = .(rt)
    , wid = .(subject)
    , within = .(rci)
    , between = NULL
    , detailed = F
  )
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Do accuracy analysis

# first, put the data in long format
accLong <- doLong(finalACC)

# disable scientific notation
options(scipen = 999)
options(digits = 4)

## ANVOA using ezANOVA 
aov.acc <- ezANOVA(
  data = accLong
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci, rcichange, task)
  , between = NULL
  , detailed = F
)
#------------------------------------------------------------------------------


### log RTs

logRTs <- rtsLong
logRTs$rt <- log(logRTs$rt)

## ANVOA using ezANOVA 
aov.log <- ezANOVA(
  data = logRTs
  , dv = .(rt)
  , wid = .(subject)
  , within = .(sequence, rci, rcichange, task)
  , between = NULL
  , detailed = F
)

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Figure showing difference between decay & temporal distinctiveness account

# get the means for RT (and then we change them to be examplar decay & TD)
decay <- summariseRTs(rtsLong)

# re-order factor levels etc.
decay$RCI  <- factor(rts$RCI, levels = c("50", "1000"))
decay$RCIChange <- factor(rts$RCIChange, levels = c("Same", "Different"))
colnames(decay)[4] <- "Hypothesis"

#rename "location" and "parity" tasks to "Perceptual" and "Memory
levels(decay$Hypothesis)[levels(decay$Hypothesis)=="perceptual"] <- "Decay"
levels(decay$Hypothesis)[levels(decay$Hypothesis)=="memory"] <- "Interference"

# Hypothetical data values
decay$RT <- c(900, 900, 1050, 1050, 1300, 1300, 1150, 1150, 900, 800, 900, 
              1050, 1300, 1300, 1300, 1300)

# to ensure error bars don't overlap
pd <- position_dodge(.02)


theme_set(theme_gray(base_size = 18))

pdf("hypothesisPlot.pdf", width = 9, height = 6)

plot <- ggplot(decay, aes(x = RCI, y = RT, group = Sequence, colour = Sequence)) 
plot <- plot + geom_line(aes(linetype = Sequence))
plot <- plot + scale_linetype_manual(breaks = c("Repeat", "Switch"), values = c(5,1))
plot <- plot + geom_point(aes(shape = Sequence), size = 2.3)
plot <- plot + scale_x_discrete(name = "Response-Cue Interval (ms)") + 
  scale_y_continuous(limits = c(600, 1500), name = "Response Time (ms)")



plot + facet_grid(Hypothesis ~ RCIChange)

dev.off()
#------------------------------------------------------------------------------



