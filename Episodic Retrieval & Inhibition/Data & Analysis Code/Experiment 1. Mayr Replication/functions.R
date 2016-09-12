### Functions script

#------------------------------------------------------------------------------
### Trim RTs using SD criterion
getRTs <- function(data, minRT, sd){
  
  
  # get rid of fast RTs
  data <- subset(data, data$rt > minRT)
  
  # get vector of unique subject numbers
  subjects <- unique(data$subject)
  
  # initialise matrix to store data in 
  rtData <- matrix(0, nrow = length(subjects), ncol = 5)
  colnames(rtData) <- c("subject", 
                        "abaRep", "abaSw", 
                        "cbaRep", "cbaSw")
  
  # to keep track of how many subjects have been processed
  i = 1
  
  # loop over each subject
  for(sub in subjects){
    
    # add subject number to data matrix
    rtData[i, 1] <- sub
    
    # get this subject's data
    subData <- subset(data, data$subject == sub)
    
    # get their ABA data first
    abaData <- subset(subData, subData$sequence == "ABA")
    
      # get their stimRepetition == yes data
      repData <- subset(abaData, abaData$stimRep == "yes")
    
      # get the trimmed mean RT
      rts <- as.numeric(repData$rt)
      rtData[i, 2] <- trimRTs(rts, sd)
    
      # get their stimRepetition == no data
      swData <- subset(abaData, abaData$stimRep == "no")
      
      # get the trimmed mean RT
      rts <- as.numeric(swData$rt)
      rtData[i, 3] <- trimRTs(rts, sd)
    
    
    # now get their CBA data 
    cbaData <- subset(subData, subData$sequence == "CBA")
    
      # get their stimRepetition == yes data
      repData <- subset(cbaData, cbaData$stimRep == "yes")
      
      # get the trimmed mean RT
      rts <- as.numeric(repData$rt)
      rtData[i, 4] <- trimRTs(rts, sd)
      
      # get their stimRepetition == no data
      swData <- subset(cbaData, cbaData$stimRep == "no")
      
      # get the trimmed mean RT
      rts <- as.numeric(swData$rt)
      rtData[i, 5] <- trimRTs(rts, sd)
    
    # update looping count
    i <- i + 1
        
  } # end of subject loop
  
  # change to data frame
  rtData <- as.data.frame(rtData)
  
  return(rtData)
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### sub-function of "getRTs" that actually does the trimming
trimRTs <- function(rts, sd){
  
  # get passed vector of RTs as numeric
  meanRT <- mean(rts)
  sdRT <- sd(rts)
  
  # what is the upper cutoff?
  upperCutoff <- meanRT + (sd * sdRT)
  
  # now trim these RTs
  finalRTs <- rts[rts < upperCutoff]
  meanRT <- round(mean(finalRTs), digits = 0)
  
  return(meanRT)
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Change data to long format

doLong <- function(data){
  
  # to store all data
  finalData <- NULL 
  
  
  # loop over each subject (indexed by row)
  for(i in 1:nrow(data)){

    # matrix for each subject's data
    subjectData <- matrix(0, nrow = 4, ncol = 4)
    subjectData <- data.frame(subjectData)
    colnames(subjectData) <- c("subject", "sequence", "stimRep", "rt")
    
    # populate the data frame
    subjectData[1, ] <- c(data$subject[i], "aba", "repetition", data$abaRep[i])
    subjectData[2, ] <- c(data$subject[i], "aba", "switch", data$abaSw[i])
    subjectData[3, ] <- c(data$subject[i], "cba", "repetition", data$cbaRep[i])
    subjectData[4, ] <- c(data$subject[i], "cba", "switch", data$cbaSw[i])
    
    # change certain columns to numeric
    subjectData$subject = as.numeric(subjectData$subject)
    subjectData$rt = as.numeric(subjectData$rt)
    
    # concatonate with other subjects' data
    finalData <- rbind(finalData, subjectData)
  
  }
  
  return(finalData)
  
}

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Plot progression of Bayes factor of interaction as subjects are added
plotBF <- function(data, scale = 1){
  
  # blank matrix to store data in
  bfProg <- matrix(0, ncol = 2, nrow = length(data))
  colnames(bfProg) <- c("N", "BF")
  
  for(i in 1:length(data)){
    
    if(i == 1){
      bfProg[i, 1] <- 1
      bfProg[i, 2] <- 0
    }
    
    if(i > 1){
      tempData <- data[1:i]
      bfProg[i, 1] <- i
      bf <- ttestBF(x = tempData, rscale = scale)
      bfProg[i, 2] <- exp(bf@bayesFactor$bf)
    }
  }
  return(bfProg)
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### plot the BEST HDI as sample size increases
plotHDI <- function(data){
  
  # blank matrix to store data in
  hdiProg <- matrix(0, ncol = 3, nrow = length(data))
  colnames(hdiProg) <- c("N", "hdiLow", "hdiHigh")
  
  for(i in 1:length(data)){
    
    if(i == 1){
      hdiProg[i, 1] <- 1
      hdiProg[i, 2] <- 0
      hdiProg[i, 3] <- 0
    }
    
    if(i > 1){
      tempData <- data[1:i]
      hdiProg[i, 1] <- i
      best <- BESTmcmc(tempData)
      best <- summary(best)
      hdiProg[i, 2] <- best[25]
      hdiProg[i, 3] <- best[30]
    }
  }
  return(hdiProg)
}
#------------------------------------------------------------------------------
