#------------------------------------------------------------------------------
### check which participants don't have data for all conditions
completeData <- function(data){
  
  removedparticipants <- NULL
  
  # loop over all participants
  for(i in 1:max(data$participant)){
    
    # get their data
    currData <- subset(data, data$participant == i)
    
    # get the data for each condition
    targetData <- subset(currData, currData$paradigm == "target")
    visualData <- subset(currData, currData$paradigm == "visual")
    numericData <- subset(currData, currData$paradigm == "numeric")
    
    # if any of these data sets are empty, then add the current participant to
    # the removal list
    if(nrow(targetData) == 0 | nrow(visualData) == 0 | nrow(numericData) == 0){
      removedparticipants <- c(removedparticipants, i)
    }
    
  }
  
  return(removedparticipants)
  
}
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
### find participants to remove due to excess error
accuracyRemoval <- function(data, accCriterion){
  
  removedparticipants <- NULL
  
  
  # get the accuracy per condition
  targetAcc <- subset(data, data$paradigm == "target")
  targetAcc <- targetAcc %>%
    group_by(participant) %>%
    summarise(meanAcc = (sum(accuracy) / length(accuracy)) * 100)
  removedparticipants <- c(removedparticipants, 
                       targetAcc$participant[targetAcc$meanAcc < accCriterion])
  
  
  visualAcc <- subset(data, data$paradigm == "visual")
  visualAcc <- visualAcc %>%
    group_by(participant) %>%
    summarise(meanAcc = (sum(accuracy) / length(accuracy)) * 100)
  removedparticipants <- c(removedparticipants, 
                       visualAcc$participant[visualAcc$meanAcc < accCriterion])
  
  
  
  numericAcc <- subset(data, data$paradigm == "numeric")
  numericAcc <- numericAcc %>%
    group_by(participant) %>%
    summarise(meanAcc = (sum(accuracy) / length(accuracy)) * 100)
  removedparticipants <- c(removedparticipants, 
                       numericAcc$participant[numericAcc$meanAcc < accCriterion])
  
  
  removedparticipants <- sort(unique(removedparticipants))
  
  a <- 61
  removedparticipants <- c(removedparticipants, a)
  
  return(removedparticipants)
  
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
### do the response time trimming
rtTrimming <- function(data, minRT, sd){
  
  # remove all errors from the data file
  data <- subset(allData,accuracy == 1)
  
  # do each paradigm's trimming
  targetData <- subset(allData, allData$paradigm == "target")
  targetData <- sdTrim(targetData, minRT = 150, perCondition = TRUE, 
                       perParticipant = TRUE, sd = 2.5, returnType = "raw")
  
  visualData <- subset(data, data$paradigm == "visual")
  visualData <- sdTrim(visualData, minRT = 150, perCondition = TRUE, 
                     perParticipant = TRUE, sd = 2.5, returnType = "raw")
  
  numericData <- subset(data, data$paradigm == "numeric")  
  numericData <- sdTrim(numericData, minRT = 150, perCondition = TRUE, 
                       perParticipant = TRUE, sd = 2.5, returnType = "raw")
  
  trimmedData <- rbind(targetData, visualData, numericData)
  
  return(trimmedData)
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### split-half reliability for response times
splitHalf <- function(data, splitType = "oddEven", nSplits = 1){
  
  # add column to data that will determine which "half" the data is in
  data <- mutate(data, half = 0)
  
  # data frame to store correlation data
  targetCor <- numeric(nSplits)
  visualCor <- numeric(nSplits)
  numericCor <- numeric(nSplits)
  corData <- data.frame(targetCor, visualCor, numericCor)
  
  ## do odd-even split
  if(splitType == "oddEven"){
    
    # assign the "half" (1 for odd trials, 0 for even)
    data$half <- data$trial %% 2
    
    # do the different paradigms (this is a hack; there must be an elegant
    # way of doing this)
    
    #--- target data
    curData <- subset(data, data$paradigm == "target")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanRT - cba1$meanRT
    bi0 <- aba0$meanRT - cba0$meanRT
    
    # store the correlation
    corData$targetCor <- round(cor(bi1, bi0), 3)
    
    
    #--- visual data
    curData <- subset(data, data$paradigm == "visual")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanRT - cba1$meanRT
    bi0 <- aba0$meanRT - cba0$meanRT
    
    # store the correlation
    corData$visualCor <- round(cor(bi1, bi0), 3)
    
    
    #--- numeric data
    curData <- subset(data, data$paradigm == "numeric")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanRT = mean(rt))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanRT - cba1$meanRT
    bi0 <- aba0$meanRT - cba0$meanRT
    
    # store the correlation
    corData$numericCor <- round(cor(bi1, bi0), 3)
    
    
    return(corData)
    
  } # end of odd/even split
  
  ## do multiple splits
  if(splitType == "random"){
    
    for(i in 1:nSplits){
      
      # assign a 1 or zero randomly to the "half" column
      data$half <- base::sample(c(1, 0), size = nrow(data), replace = TRUE,
                                prob = c(0.5, 0.5))
      
      #--- target data
      curData <- subset(data, data$paradigm == "target")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanRT - cba1$meanRT
      bi0 <- aba0$meanRT - cba0$meanRT
      
      # store the correlation
      corData$targetCor[i] <- round(cor(bi1, bi0), 3)
      
      
      #--- visual data
      curData <- subset(data, data$paradigm == "visual")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanRT - cba1$meanRT
      bi0 <- aba0$meanRT - cba0$meanRT
      
      # store the correlation
      corData$visualCor[i] <- round(cor(bi1, bi0), 3)
      
      
      #--- numeric data
      curData <- subset(data, data$paradigm == "numeric")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanRT = mean(rt))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanRT - cba1$meanRT
      bi0 <- aba0$meanRT - cba0$meanRT
      
      # store the correlation
      corData$numericCor[i] <- round(cor(bi1, bi0), 3)
      
    } # end of nSplits loop
    
    return(corData) 
    
  } # end of random split
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### split-half reliability for accuracy
splitHalf_acc <- function(data, splitType = "oddEven", nSplits = 1){
  
  # add column to data that will determine which "half" the data is in
  data <- mutate(data, half = 0)
  
  # data frame to store correlation data
  targetCor <- numeric(nSplits)
  visualCor <- numeric(nSplits)
  numericCor <- numeric(nSplits)
  corData <- data.frame(targetCor, visualCor, numericCor)
  
  ## do odd-even split
  if(splitType == "oddEven"){
    
    # assign the "half" (1 for odd trials, 0 for even)
    data$half <- data$trial %% 2
    
    # do the different paradigms (this is a hack; there must be an elegant
    # way of doing this)
    
    #--- target data
    curData <- subset(data, data$paradigm == "target")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanAcc - cba1$meanAcc
    bi0 <- aba0$meanAcc - cba0$meanAcc
    
    # store the correlation
    corData$targetCor <- round(cor(bi1, bi0), 3)
    
    
    #--- visual data
    curData <- subset(data, data$paradigm == "visual")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanAcc - cba1$meanAcc
    bi0 <- aba0$meanAcc - cba0$meanAcc
    
    # store the correlation
    corData$visualCor <- round(cor(bi1, bi0), 3)
    
    
    #--- numeric data
    curData <- subset(data, data$paradigm == "numeric")
    
    aba1 <- curData %>%
      filter(condition == "ABA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    aba0 <- curData %>%
      filter(condition == "ABA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba1 <- curData %>%
      filter(condition == "CBA" & half == 1) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    cba0 <- curData %>%
      filter(condition == "CBA" & half == 0) %>%
      group_by(condition, half, participant) %>%
      summarise(meanAcc = mean(accuracy))
    
    # calculate the n-2 repetition cost for each half
    bi1 <- aba1$meanAcc - cba1$meanAcc
    bi0 <- aba0$meanAcc - cba0$meanAcc
    
    # store the correlation
    corData$numericCor <- round(cor(bi1, bi0), 3)
    
    
    return(corData)
    
  } # end of odd/even split
  
  ## do multiple splits
  if(splitType == "random"){
    
    for(i in 1:nSplits){
      
      # assign a 1 or zero randomly to the "half" column
      data$half <- base::sample(c(1, 0), size = nrow(data), replace = TRUE,
                                prob = c(0.5, 0.5))
      
      #--- target data
      curData <- subset(data, data$paradigm == "target")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanAcc - cba1$meanAcc
      bi0 <- aba0$meanAcc - cba0$meanAcc
      
      # store the correlation
      corData$targetCor[i] <- round(cor(bi1, bi0), 3)
      
      
      #--- visual data
      curData <- subset(data, data$paradigm == "visual")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanAcc - cba1$meanAcc
      bi0 <- aba0$meanAcc - cba0$meanAcc
      
      # store the correlation
      corData$visualCor[i] <- round(cor(bi1, bi0), 3)
      
      
      #--- numeric data
      curData <- subset(data, data$paradigm == "numeric")
      
      aba1 <- curData %>%
        filter(condition == "ABA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      aba0 <- curData %>%
        filter(condition == "ABA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba1 <- curData %>%
        filter(condition == "CBA" & half == 1) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      cba0 <- curData %>%
        filter(condition == "CBA" & half == 0) %>%
        group_by(condition, half, participant) %>%
        summarise(meanAcc = mean(accuracy))
      
      # calculate the n-2 repetition cost for each half
      bi1 <- aba1$meanAcc - cba1$meanAcc
      bi0 <- aba0$meanAcc - cba0$meanAcc
      
      # store the correlation
      corData$numericCor[i] <- round(cor(bi1, bi0), 3)
      
    } # end of nSplits loop
    
    return(corData) 
    
  } # end of random split
  
}
#------------------------------------------------------------------------------
