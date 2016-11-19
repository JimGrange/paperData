doTrim <- function(data, minRT, sd, digits = 3){
  
  trimmedData <- data
  
  ###-------------
    # get the list of participant numbers
    participant <- sort(unique(trimmedData$participant))
    
    # get the list of experimental conditions
    conditionList <- unique(trimmedData$condition)
    
    # trim the data to remove trials below minRT
    trimmedData <- subset(trimmedData, trimmedData$rt > minRT)
    
    # change the variable name for sd (as this is an R function)
    stDev <- sd
      
      # initialise variable to keep trimmed data in
      finalData <- NULL
      
      # loop over all participants
      for(currSub in participant){
        
        # loop over all conditions
        for(currCond in conditionList){
          
          # get the relevant data
          tempData <- subset(trimmedData, trimmedData$condition == currCond &
                               trimmedData$participant == currSub)
          
          # find the cutoff
          curMean <- mean(tempData$rt)
          curSD <- sd(tempData$rt)
          curCutoff <- curMean + (stDev * curSD)
          
          # perform the trim
          curData <- subset(tempData, tempData$rt < curCutoff)
          
          # store the data
          finalData <- rbind(finalData, curData)
        }
      }
      
      return(finalData)
  
} # end of function
