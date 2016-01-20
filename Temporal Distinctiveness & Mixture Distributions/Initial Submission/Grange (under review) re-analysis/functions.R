### functions for analysis

#------------------------------------------------------------------------------
# change data to long format
doLong <- function(data){
  
  locationData <- subset(data, data$condition == "Perceptual")
  parityData <- subset(data, data$condition == "Memory")
  
  # do location data first
  locData <- NULL
  
  for(i in 1:nrow(locationData)){
    
    # matrix for subject data
    # rows = 8 conditions
    # columns = subject, sequence, rci, rciChange, task, rt
    subjectData <- matrix(0, nrow = 8, ncol = 6)
    subjectData <- data.frame(subjectData)
    colnames(subjectData) <- c("subject", "sequence", "rci", "rcichange", 
                               "task", "rt")
    
    # populate the data frame
    subjectData[1, ] <- c(locationData$subject[i], "repeat", "50", "same", 
                          "perceptual", locationData$Rep50Same[i])
    
    subjectData[2, ] <- c(locationData$subject[i], "repeat", "50", "different",
                          "perceptual", locationData$Rep50Different[i])
    
    subjectData[3, ] <- c(locationData$subject[i], "repeat", "1000", "same",
                          "perceptual", locationData$Rep1000Same[i])
    
    subjectData[4, ] <- c(locationData$subject[i], "repeat", "1000", "different", 
                          "perceptual", locationData$Rep1000Different[i])
    
    subjectData[5, ] <- c(locationData$subject[i], "switch", "50", "same", 
                          "perceptual", locationData$Sw50Same[i])
    
    subjectData[6, ] <- c(locationData$subject[i], "switch", "50", "different",
                          "perceptual", locationData$Sw50Different[i])
    
    subjectData[7, ] <- c(locationData$subject[i], "switch", "1000", "same",
                          "perceptual", locationData$Sw1000Same[i])
    
    subjectData[8, ] <- c(locationData$subject[i], "switch", "1000", "different", 
                          "perceptual", locationData$Sw1000Different[i])
    
    #change certain columns to numeric
    subjectData$subject <- as.numeric(subjectData$subject)
    subjectData$rt <- as.numeric(subjectData$rt)
    
    # concatonate with other subjects' data
    locData <- rbind(locData, subjectData)
    
  }
  
  
  # now do parity data
  parData <- NULL
  
  for(i in 1:nrow(parityData)){
    
    # matrix for subject data
    # rows = 8 conditions
    # columns = subject, sequence, rci, rciChange, task, rt
    subjectData <- matrix(0, nrow = 8, ncol = 6)
    subjectData <- data.frame(subjectData)
    colnames(subjectData) <- c("subject", "sequence", "rci", "rcichange", 
                               "task", "rt")
    
    # populate the data frame
    subjectData[1, ] <- c(parityData$subject[i], "repeat", "50", "same", 
                          "memory", parityData$Rep50Same[i])
    
    subjectData[2, ] <- c(parityData$subject[i], "repeat", "50", "different",
                          "memory", parityData$Rep50Different[i])
    
    subjectData[3, ] <- c(parityData$subject[i], "repeat", "1000", "same",
                          "memory", parityData$Rep1000Same[i])
    
    subjectData[4, ] <- c(parityData$subject[i], "repeat", "1000", "different", 
                          "memory", parityData$Rep1000Different[i])
    
    subjectData[5, ] <- c(parityData$subject[i], "switch", "50", "same", 
                          "memory", parityData$Sw50Same[i])
    
    subjectData[6, ] <- c(parityData$subject[i], "switch", "50", "different",
                          "memory", parityData$Sw50Different[i])
    
    subjectData[7, ] <- c(parityData$subject[i], "switch", "1000", "same",
                          "memory", parityData$Sw1000Same[i])
    
    subjectData[8, ] <- c(parityData$subject[i], "switch", "1000", "different", 
                          "memory", parityData$Sw1000Different[i])
    
    #change certain columns to numeric
    subjectData$subject <- as.numeric(subjectData$subject)
    subjectData$rt <- as.numeric(subjectData$rt)
    
    # concatonate with other subjects' data
    parData <- rbind(parData, subjectData)
    
  }
  
  allData <- rbind(locData, parData)
  
  return(allData)
  
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# get means, ses of each condition
# see this: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

summariseRTs <- function(data){
  
  # get final data frame to store everything  
  Sequence <- rep(c("Repeat", "Repeat", "Repeat", "Repeat", "Switch", 
                    "Switch", "Switch", "Switch"), 2)
  RCI <- rep(c("50", "50", "1000", "1000"), 4)
  RCIChange <- rep(c("Same", "Different"), 8)
  Task <- c(rep("perceptual", 8), rep("memory", 8))
  RT <- rep(0, 16)
  SE <- rep(0, 16)
  
  summarised <- data.frame(Sequence, RCI, RCIChange, Task, RT, SE)
  
  
  # do location task first
  tempData <- subset(data, data$task == "perceptual")
  
  repData <- subset(tempData, tempData$sequence == "repeat")
  
  fiftyData <- subset(repData, repData$rci == "50")
  
  sameData <- subset(fiftyData, fiftyData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[1] <- round(mean(rts), 0)
  summarised$SE[1] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(fiftyData, fiftyData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[2] <- round(mean(rts), 0)
  summarised$SE[2] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  thousandData <- subset(repData, repData$rci == "1000")
  
  
  sameData <- subset(thousandData, thousandData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[3] <- round(mean(rts), 0)
  summarised$SE[3] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(thousandData, thousandData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[4] <- round(mean(rts), 0)
  summarised$SE[4] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  swData <- subset(tempData, tempData$sequence == "switch")
  
  fiftyData <- subset(swData, swData$rci == "50")
  
  sameData <- subset(fiftyData, fiftyData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[5] <- round(mean(rts), 0)
  summarised$SE[5] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(fiftyData, fiftyData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[6] <- round(mean(rts), 0)
  summarised$SE[6] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  thousandData <- subset(swData, swData$rci == "1000")
  
  
  sameData <- subset(thousandData, thousandData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[7] <- round(mean(rts), 0)
  summarised$SE[7] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(thousandData, thousandData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[8] <- round(mean(rts), 0)
  summarised$SE[8] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  tempData <- subset(data, data$task == "memory")
  
  repData <- subset(tempData, tempData$sequence == "repeat")
  
  fiftyData <- subset(repData, repData$rci == "50")
  
  sameData <- subset(fiftyData, fiftyData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[9] <- round(mean(rts), 0)
  summarised$SE[9] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(fiftyData, fiftyData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[10] <- round(mean(rts), 0)
  summarised$SE[10] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  thousandData <- subset(repData, repData$rci == "1000")
  
  
  sameData <- subset(thousandData, thousandData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[11] <- round(mean(rts), 0)
  summarised$SE[11] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(thousandData, thousandData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[12] <- round(mean(rts), 0)
  summarised$SE[12] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  swData <- subset(tempData, tempData$sequence == "switch")
  
  fiftyData <- subset(swData, swData$rci == "50")
  
  sameData <- subset(fiftyData, fiftyData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[13] <- round(mean(rts), 0)
  summarised$SE[13] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(fiftyData, fiftyData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[14] <- round(mean(rts), 0)
  summarised$SE[14] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  thousandData <- subset(swData, swData$rci == "1000")
  
  
  sameData <- subset(thousandData, thousandData$rcichange == "same")
  
  rts <- as.numeric(sameData$rt)
  summarised$RT[15] <- round(mean(rts), 0)
  summarised$SE[15] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  diffData <- subset(thousandData, thousandData$rcichange == "different")
  
  rts <- as.numeric(diffData$rt)
  summarised$RT[16] <- round(mean(rts), 0)
  summarised$SE[16] <- round(sd(rts) / sqrt(length(rts)), 0)
  
  
  
  return(summarised)
  
}



#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#create multiple plots with ggplot2
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#----------------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# get matrix of trimmed RTs using SD trimming
getRTs <- function(data, minRT, sdCriterion){
  
  # get rid of fast RTs
  data <- subset(data, data$rt > minRT)
  
  
  # get vector of subject numbers
  subjects <- unique(data$subject)
  
  # initialise matrix to store data in
  rtData <- matrix(0, nrow = length(subjects), ncol = 9)
  colnames(rtData) <- c("subject",
                        "Rep50Same", "Rep50Different", 
                        "Rep1000Same", "Rep1000Different", 
                        "Sw50Same", "Sw50Different", 
                        "Sw1000Same", "Sw1000Different")
  
  
  
  # to keep track of how many subjects have been processed
  i = 1
  
  # loop over subjects
  for(sub in subjects){
    
    # add subject number
    rtData[i, 1] <- sub
    
    #get their data
    subData <- subset(data, data$subject == sub)
    
    # get their "repeat" data
    repData <- subset(subData, subData$sequence == "repeat")
    
    # get "50ms" data
    fiftyData <- subset(repData, repData$currentRCI == 50)
    
    # get "same" data & calculate % correct
    sameData <- subset(fiftyData, fiftyData$rciChange == "same")
    rts <- as.numeric(sameData$rt)
    rtData[i, 2] <- trimRTs(rts, sdCriterion)
    
    # get "different data & calculate % correct
    diffData <- subset(fiftyData, fiftyData$rciChange == "different")
    rts <- as.numeric(diffData$rt)
    rtData[i, 3] <- trimRTs(rts, sdCriterion)
    
    # get "1000ms" data
    thousandData <- subset(repData, repData$currentRCI == 1000)
    
    # get "same" data & calculate % correct
    sameData <- subset(thousandData, thousandData$rciChange == "same")
    rts <- as.numeric(sameData$rt)
    rtData[i, 4] <- trimRTs(rts, sdCriterion)
    
    # get "different data & calculate % correct
    diffData <- subset(thousandData, thousandData$rciChange == "different")
    rts <- as.numeric(diffData$rt)
    rtData[i, 5] <- trimRTs(rts, sdCriterion)
    
    
    
    # get their "switch" data
    swData <- subset(subData, subData$sequence == "switch")
    
    # get "50ms" data
    fiftyData <- subset(swData, swData$currentRCI == 50)
    
    # get "same" data & calculate % correct
    sameData <- subset(fiftyData, fiftyData$rciChange == "same")
    rts <- as.numeric(sameData$rt)
    rtData[i, 6] <- trimRTs(rts, sdCriterion)
    
    # get "different data & calculate % correct
    diffData <- subset(fiftyData, fiftyData$rciChange == "different")
    rts <- as.numeric(diffData$rt)  
    rtData[i, 7] <- trimRTs(rts, sdCriterion)
    
    # get "1000ms" data
    thousandData <- subset(swData, swData$currentRCI == 1000)
    
    # get "same" data & calculate % correct
    sameData <- subset(thousandData, thousandData$rciChange == "same")
    rts <- as.numeric(sameData$rt)
    rtData[i, 8] <- trimRTs(rts, sdCriterion)
    
    # get "different data & calculate % correct
    diffData <- subset(thousandData, thousandData$rciChange == "different")
    rts <- as.numeric(diffData$rt)
    rtData[i, 9] <- trimRTs(rts, sdCriterion)
    
    
    
    # update loop count
    i = i + 1
  }  
  
  return(rtData)
  
}

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# trim RT vector using SD cutoff
trimRTs <- function(rts, sdCriterion){
  
  # get passed vector of RTs as numeric
  meanRT <- mean(rts)
  sdRT <- sd(rts)
  
  upperCutoff <- meanRT + (sdCriterion * sdRT)
  
  finalRTs <- rts[rts < upperCutoff]
  meanRT <- round(mean(finalRTs), digits = 0)
  
  return(meanRT)
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# add 0 after error trials
addError <- function(data){
  
  
  
  # to store final data frame
  addErrorData <- NULL
  
  # get vector of subject numbers
  subjects <- unique(data$subject)
  
  
  # to keep track of how many subjects have been looped over
  j = 1 
  
  
  #loop over all subjects, add the column, and then add to addedColumn df
  for(sub in subjects){
    
    
    #get their data
    subData <- subset(data, data$subject == sub)
    
    # loop over trials and store accuracy trimming
    for(i in 2:nrow(subData)){
      
      # multiply current trial's accuracy by previous trial's
      subData$accuracyTrim[i] = subData$accuracy[i] * subData$accuracy[i - 1]
      
    }
    
    if(j == 1){
      addErrorData = subData
    } else {
      addErrorData = rbind(addErrorData, subData)
    }
    
    # update subject count
    j = j + 1
  }
  
  return(addErrorData)
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# return matrix of accuracy for each subject & condition
getAccuracy <- function(data){
  
  # get vector of subject numbers
  subjects <- unique(data$subject)
  
  # initialise matrix to store data in
  accuracyData <- matrix(0, nrow = length(subjects), ncol = 9)
  colnames(accuracyData) <- c("subject",
                              "Rep50Same", "Rep50Different", 
                              "Rep1000Same", "Rep1000Different", 
                              "Sw50Same", "Sw50Different", 
                              "Sw1000Same", "Sw1000Different")
  
  
  # to keep track of how many subjects have been processed
  i = 1
  
  # loop over subjects
  for(sub in subjects){
    
    # add subject number
    accuracyData[i, 1] <- sub
    
    #get their data
    subData <- subset(data, data$subject == sub)
    
    # get their "repeat" data
    repData <- subset(subData, subData$sequence == "repeat")
    
    # get "50ms" data
    fiftyData <- subset(repData, repData$currentRCI == 50)
    
    # get "same" data & calculate % correct
    sameData <- subset(fiftyData, fiftyData$rciChange == "same")
    accuracyData[i, 2] <- (sum(sameData$accuracy) / nrow(sameData)) * 100
    
    # get "different data & calculate % correct
    diffData <- subset(fiftyData, fiftyData$rciChange == "different")
    accuracyData[i, 3] <- (sum(diffData$accuracy) / nrow(diffData)) * 100
    
    # get "1000ms" data
    thousandData <- subset(repData, repData$currentRCI == 1000)
    
    # get "same" data & calculate % correct
    sameData <- subset(thousandData, thousandData$rciChange == "same")
    accuracyData[i, 4] <- (sum(sameData$accuracy) / nrow(sameData)) * 100
    
    # get "different data & calculate % correct
    diffData <- subset(thousandData, thousandData$rciChange == "different")
    accuracyData[i, 5] <- (sum(diffData$accuracy) / nrow(diffData)) * 100
    
    
    
    # get their "repeat" data
    swData <- subset(subData, subData$sequence == "switch")
    
    # get "50ms" data
    fiftyData <- subset(swData, swData$currentRCI == 50)
    
    # get "same" data & calculate % correct
    sameData <- subset(fiftyData, fiftyData$rciChange == "same")
    accuracyData[i, 6] <- (sum(sameData$accuracy) / nrow(sameData)) * 100
    
    # get "different data & calculate % correct
    diffData <- subset(fiftyData, fiftyData$rciChange == "different")
    accuracyData[i, 7] <- (sum(diffData$accuracy) / nrow(diffData)) * 100
    
    # get "1000ms" data
    thousandData <- subset(swData, swData$currentRCI == 1000)
    
    # get "same" data & calculate % correct
    sameData <- subset(thousandData, thousandData$rciChange == "same")
    accuracyData[i, 8] <- (sum(sameData$accuracy) / nrow(sameData)) * 100
    
    # get "different data & calculate % correct
    diffData <- subset(thousandData, thousandData$rciChange == "different")
    accuracyData[i, 9] <- (sum(diffData$accuracy) / nrow(diffData)) * 100
    
    # update loop count
    i = i + 1
    
    
  }
  
  
  return(accuracyData)
  
}

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------


#------------------------------------------------------------------------------