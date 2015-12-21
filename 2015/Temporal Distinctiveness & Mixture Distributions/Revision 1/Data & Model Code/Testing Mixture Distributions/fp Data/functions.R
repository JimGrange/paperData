#------------------------------------------------------------------------------
### fixed-point code adapted from the "fp" package
# http://www.leendertvanmaanen.com/resources/fp/example.html

fp1 <- setClass("fp1", representation(dens="array", diff="data.frame", 
                                      dat="data.frame"))

plot.fp1 <- function(x, ..., ylab=c("Density","Density difference"), 
                     xlim=NULL) {
  options(warn=-1)
  if (is.null(xlim)) xlim <- range(x@dens[[1]]$x)
  plot(1,type="n",
       ylim=c(0,max(unlist(lapply(x@dens, function(X) {max(X$y)})))),
       xlim=xlim,
       bty="n", ylab=ylab[1], ...)
  for (i in 1:length(x@dens)) {
    lines(x@dens[[i]], col=i)
  }
  legend("top", c("1: 0.05", "2: 1.00", "3: 20.00"), lty = c(1, 1, 1), 
         col = c(1, 2, 3), bty = "n")
  
  matplot(x@dens[[1]]$x,x@diff, type='l', pch=16,
          lwd=3, bty="n", ylab=ylab[2], xlim=xlim, ...); abline(h=0, lty=1)
  legend("top", c("3-1", "2-1", "3-2"), bty = "n", 
         lty = c(1, 2, 3), 
         lwd = c(2, 2, 2))
}
setMethod("plot","fp1", function(x, ...) plot.fp1(x, ...))

fpPlot <- function(...) {
  plot.fp1(...) # just a wrapper
}

fpGet <- function(dat, n=512, bw='nrd0') {
  # dat: nx2 dataframe or matrix with in col 1: RT; col 2: condition
  if (is.matrix(dat)) dat <- as.data.frame(dat)
  
  rng <- range(dat[[1]])
  dens <- tapply(dat[[1]], dat[[2]], density, from=rng[1], to=rng[2], n=n, 
                 bw=bw)
  diff <- NULL
  for (i in 2:length(dens)) {
    for (j in 1:(i-1)) {
      
      # difference (density-based method)
      tmp <- data.frame(dens[[i]]$y - dens[[j]]$y)
      names(tmp) <- paste(i,j,sep='-')
      if (is.null(diff)) diff <- tmp else diff <- cbind(diff, tmp)
    }
  }
  fp1(dens=dens, diff=diff, dat=dat)
}

fpDensDiff <- function(object) {
  if (!is.list(object)) {
    object <- list(object)
  }
  
  .get.diffs <- function(X) {
    lwr <- min(unlist(lapply(X@dens, function(Y) {which.max(diff(Y$y))})))
    #min(unlist(lapply(X@dens, function(Y) {which.max(Y$y)})))
    upr <- max(unlist(lapply(X@dens, function(Y) {which.min(diff(Y$y))})))
    
    sapply(1:ncol(X@diff), function(i) {
      y <- X@diff[,i][lwr:upr]
      x <- X@dens[[1]]$x[lwr:upr]
      y1 <- y[!is.na(y)&y!=Inf&y!=0]
      x1 <- x[!is.na(y)&y!=Inf&y!=0]
      index <- which.min(abs(y1))
      x1[index]
    })
  }
  
  roots <- lapply(object, .get.diffs)
  root <- array(dim=c(dim(object[[1]]@diff)[2],length(object)))
  for (i in 1:length(roots)) {
    root[,i] <- unlist(roots[[i]])
  }
  root
}

fpAnova <- function(object, stat="BF", na.rm=TRUE) {
  require(BayesFactor)
  bf <- p <- NULL
  tmp <- fpDensDiff(object)
  tmp <- data.frame(x=c(tmp), cross=factor(1:nrow(tmp)),
                    pp=factor(rep(1:ncol(tmp),each=nrow(tmp))))
  # because tmp is a pp x cross array, we need to test across pp whether the ratios 
  # differ.
  if (na.rm) tmp <- tmp[!is.na(tmp$cross),]
  if (stat=="BF"|stat=="both") {
    bf <- anovaBF(x~cross+pp, whichRandom="pp", data=tmp, progress=F)
  }
  if (stat=="p"|stat=="both") {
    p <- summary(aov(x~cross+Error(pp/cross), data=tmp))
  }
  list(BF=bf, p=p)
}

dnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: quantiles
  #mean/sd: vector of 2
  #p: mixture prop
  p*dnorm(x, mean[1], sd[1]) + (1-p)*dnorm(x, mean[2], sd[2])
}

pnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: probabilities
  #mean/sd: vector of 2
  #p: mixture prop
  p*pnorm(x, mean[1], sd[1]) + (1-p)*pnorm(x, mean[2], sd[2])
}

qnormMix <- function(x, mean=c(0,1), sd=c(1,1), p=1) {
  #x: quantiles
  #mean/sd: vector of 2
  #p: mixture prop
  p*qnorm(x, mean[1], sd[1]) + (1-p)*qnorm(x, mean[2], sd[2])
}

rnormMix <- function(n, mean=c(0,1), sd=c(1,1), p=1) {
  #n number of obs
  #mean/sd: vector of 2
  #p: mixture prop
  ifelse(sample(0:1,n, replace=T, prob=c(p,1-p)),rnorm(n,mean[1], sd[1]), 
         rnorm(n,mean[2], sd[2]))
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### do standard deviation response time trimming
# adapted from my "trimr" R package:
# https://cran.r-project.org/web/packages/trimr/index.html

sdTrim <- function(data, minRT, sd, perCondition = TRUE, perParticipant = TRUE,
                   omitErrors = TRUE, returnType = "mean", digits = 3){
  
  ###-------------
  if(perCondition == FALSE & perParticipant == FALSE){
    # change the variable name for sd (as this is an R function)
    stDev <- sd
    
    # remove errors if the user has asked for it
    if(omitErrors == TRUE){
      trimmedData <- subset(data, data$accuracy == 1)
    } else {
      trimmedData <- data
    }
    
    # get the list of participant numbers
    participant <- sort(unique(trimmedData$participant))
    
    # get the list of experimental conditions
    conditionList <- unique(trimmedData$condition)
    
    # trim the data to remove trials below minRT
    trimmedData <- subset(trimmedData, trimmedData$rt > minRT)
    
    # what is the mean & SD of the whole group's data?
    meanRT <- mean(trimmedData$rt)
    sdRT <- sd(trimmedData$rt)
    
    # what is the cut-off value?
    cutoff <- meanRT + (stDev * sdRT)
    
    # remove these rts
    trimmedData <- subset(trimmedData, trimmedData$rt < cutoff)
    
    
    # if the user asked for trial-level data, return immediately to user
    if(returnType == "raw"){
      return(trimmedData)
    }
    
    # if the user has asked for means, then split the data into separate
    # conditions, and display the means per condition.
    if(returnType == "mean"){
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      
      # loop over all conditions, and over all subjects, and find mean RT
      
      j <- 2 # to keep track of conditions looped over. Starts at 2 as this is
      # where the first condition's column is.
      
      for(currCondition in conditionList){
        
        # get the current condition's data
        tempData <- subset(trimmedData, trimmedData$condition == currCondition)
        
        
        
        #now loop over all participants
        i <- 1
        
        for(currParticipant in participant){
          
          # get that participant's data
          participantData <- subset(tempData,
                                    tempData$participant == currParticipant)
          
          # calculate & store their mean response time
          finalData[i, j] <- round(mean(participantData$rt), digits = digits)
          
          # update participant counter
          i <- i + 1
        }
        
        # update nCondition counter
        j <- j + 1
        
      } # end of condition loop
      
      return(finalData)
      
    } ## end MEAN sub-function
    
    
    # if the user has asked for medians, then split the data into separate
    # conditions, and display the medians per condition.
    if(returnType == "median"){
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      
      # loop over all conditions, and over all subjects, and find mean RT
      
      j <- 2 # to keep track of conditions looped over. Starts at 2 as this is
      # where the first condition's column is.
      
      for(currCondition in conditionList){
        
        # get the current condition's data
        tempData <- subset(trimmedData, trimmedData$condition == currCondition)
        
        
        #now loop over all participants
        i <- 1
        
        for(currParticipant in participant){
          
          # get that participant's data
          participantData <- subset(tempData,
                                    tempData$participant == currParticipant)
          
          # calculate & store their mean response time
          finalData[i, j] <- round(median(participantData$rt), digits = digits)
          
          # update participant counter
          i <- i + 1
        }
        
        
        # update nCondition counter
        j <- j + 1
        
      } # end of condition loop
      
      return(finalData)
    }
    
  } # end of perCell == FALSE & perParticipant == FALSE
  
  
  ###-------------
  if(perCondition == TRUE & perParticipant == FALSE){
    
    # change the variable name for sd (as this is an R function)
    stDev <- sd
    
    # remove errors if the user has asked for it
    if(omitErrors == TRUE){
      trimmedData <- subset(data, data$accuracy == 1)
    } else {
      trimmedData <- data
    }
    
    # get the list of participant numbers
    participant <- sort(unique(trimmedData$participant))
    
    # get the list of experimental conditions
    conditionList <- unique(trimmedData$condition)
    
    # trim the data to remove trials below minRT
    trimmedData <- subset(trimmedData, trimmedData$rt > minRT)
    
    ### do "raw"
    if(returnType == "raw"){
      
      # initialise variable to keep trimmed data in
      finalData <- NULL
      
      # loop over each condition
      for(cond in conditionList){
        
        # get the data, & find cutoff
        curData <- subset(trimmedData, trimmedData$condition == cond)
        curMean <- mean(curData$rt)
        curSD <- sd(curData$rt)
        curCutoff <- curMean + (stDev * curSD)
        
        # trim the data
        curData <- subset(curData, curData$rt < curCutoff)
        
        # bind the data
        finalData <- rbind(finalData, curData)
      }
      
      return(finalData)
    }
    
    ### do "mean"
    if(returnType == "mean"){
      
      ## first, find the cutoff for each condition, and remove the necessary
      ## trials
      
      # initialise variable to keep trimmed data in
      tempData <- NULL
      
      for(cond in conditionList){
        # get the data, & find cutoff
        curData <- subset(trimmedData, trimmedData$condition == cond)
        curMean <- mean(curData$rt)
        curSD <- sd(curData$rt)
        curCutoff <- curMean + (stDev * curSD)
        
        # trim the data
        curData <- subset(curData, curData$rt < curCutoff)
        
        # bind the data
        tempData <- rbind(tempData, curData)
      }
      
      # change variable names
      trimmedData <- tempData
      tempData <- NULL
      
      ## now loop over each subject and calculate their average
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      # loop over conditions & subjects and calculate their average
      
      # to index over conditions. It starts at 2 because this is the first
      # column in the data frame containing condition information
      j <- 2
      
      for(curCondition in conditionList){
        
        # get the current condition's data
        tempData <- subset(trimmedData, trimmedData$condition == curCondition)
        
        #now loop over all participants
        i <- 1
        
        for(currParticipant in participant){
          
          # get that participant's data
          participantData <- subset(tempData,
                                    tempData$participant == currParticipant)
          
          # calculate & store their mean response time
          finalData[i, j] <- round(mean(participantData$rt), digits = digits)
          
          # update participant counter
          i <- i + 1
        }
        
        # update nCondition counter
        j <- j + 1
      }
      
      return(finalData)
    }
    
  } # end of perCell == TRUE & perParticipant == FALSE
  
  
  ###-------------
  if(perCondition == FALSE & perParticipant == TRUE){
    
    # change the variable name for sd (as this is an R function)
    stDev <- sd
    
    # remove errors if the user has asked for it
    if(omitErrors == TRUE){
      trimmedData <- subset(data, data$accuracy == 1)
    } else {
      trimmedData <- data
    }
    
    # get the list of participant numbers
    participant <- sort(unique(trimmedData$participant))
    
    # get the list of experimental conditions
    conditionList <- unique(trimmedData$condition)
    
    # trim the data to remove trials below minRT
    trimmedData <- subset(trimmedData, trimmedData$rt > minRT)
    
    
    ### do "raw"
    if(returnType == "raw"){
      
      # initialise variable to keep trimmed data in
      finalData <- NULL
      
      # loop over each subject
      for(currSub in participant){
        
        # get the current subject's data
        curData <- subset(trimmedData, trimmedData$participant == currSub)
        
        # find their mean, sd, & cutoff
        curMean <- mean(curData$rt)
        curSD <- sd(curData$rt)
        curCutoff <- curMean + (stDev * curSD)
        
        # trim the data
        curData <- subset(curData, curData$rt < curCutoff)
        
        # bind the data
        finalData <- rbind(finalData, curData)
      }
      
      return(finalData)
    }
    
    
    ### do "mean"
    if(returnType == "mean"){
      
      # initialise variable to keep trimmed data in
      tempData <- NULL
      
      # loop over each subject
      for(currSub in participant){
        
        # get the current subject's data
        curData <- subset(trimmedData, trimmedData$participant == currSub)
        
        # find their mean, sd, & cutoff
        curMean <- mean(curData$rt)
        curSD <- sd(curData$rt)
        curCutoff <- curMean + (stDev * curSD)
        
        # trim the data
        curData <- subset(curData, curData$rt < curCutoff)
        
        # bind the data
        tempData <- rbind(tempData, curData)
      }
      
      # change variable names
      trimmedData <- tempData
      tempData <- NULL
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      # loop over conditions & subjects and calculate their average
      
      # to index over conditions. It starts at 2 because this is the first
      # column in the data frame containing condition information
      j <- 2
      
      for(curCondition in conditionList){
        
        # get the current condition's data
        tempData <- subset(trimmedData, trimmedData$condition == curCondition)
        
        #now loop over all participants
        i <- 1
        
        for(currParticipant in participant){
          
          # get that participant's data
          participantData <- subset(tempData,
                                    tempData$participant == currParticipant)
          
          # calculate & store their mean response time
          finalData[i, j] <- round(mean(participantData$rt), digits = digits)
          
          # update participant counter
          i <- i + 1
        }
        
        # update nCondition counter
        j <- j + 1
      }
      
      return(finalData)
      
    }
    
    
    ### do "median"
    if(returnType == "median"){
      
      # initialise variable to keep trimmed data in
      tempData <- NULL
      
      # loop over each subject
      for(currSub in participant){
        
        # get the current subject's data
        curData <- subset(trimmedData, trimmedData$participant == currSub)
        
        # find their mean, sd, & cutoff
        curMean <- mean(curData$rt)
        curSD <- sd(curData$rt)
        curCutoff <- curMean + (stDev * curSD)
        
        # trim the data
        curData <- subset(curData, curData$rt < curCutoff)
        
        # bind the data
        tempData <- rbind(tempData, curData)
      }
      
      # change variable names
      trimmedData <- tempData
      tempData <- NULL
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      # loop over conditions & subjects and calculate their average
      
      # to index over conditions. It starts at 2 because this is the first
      # column in the data frame containing condition information
      j <- 2
      
      for(curCondition in conditionList){
        
        # get the current condition's data
        tempData <- subset(trimmedData, trimmedData$condition == curCondition)
        
        #now loop over all participants
        i <- 1
        
        for(currParticipant in participant){
          
          # get that participant's data
          participantData <- subset(tempData,
                                    tempData$participant == currParticipant)
          
          # calculate & store their median response time
          finalData[i, j] <- round(median(participantData$rt), digits = digits)
          
          # update participant counter
          i <- i + 1
        }
        
        # update nCondition counter
        j <- j + 1
      }
      
      return(finalData)
    }
    
    
  } # end of perCell == FALSE & perParticipant == TRUE
  
  
  
  ###-------------
  if(perCondition == TRUE & perParticipant == TRUE){
    # change the variable name for sd (as this is an R function)
    stDev <- sd
    
    # remove errors if the user has asked for it
    if(omitErrors == TRUE){
      trimmedData <- subset(data, data$accuracy == 1)
    } else {
      trimmedData <- data
    }
    
    # get the list of participant numbers
    participant <- sort(unique(trimmedData$participant))
    
    # get the list of experimental conditions
    conditionList <- unique(trimmedData$condition)
    
    # trim the data to remove trials below minRT
    trimmedData <- subset(trimmedData, trimmedData$rt > minRT)
    
    ### do "raw"
    if(returnType == "raw"){
      
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
    }
    
    
    
    ### do "mean"
    if(returnType == "mean"){
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      # intialise looping variable for subjects
      i <- 1
      
      # loop over all subjects
      for(currSub in participant){
        
        # intialise looping variable for conditions. It starts at 2 because the
        # first column in the data file containing condition information is the
        # second one.
        j <- 2
        
        # loop over all conditions
        for(currCond in conditionList){
          
          # get the relevant data
          tempData <- subset(trimmedData, trimmedData$participant == currSub &
                               trimmedData$condition == currCond)
          
          # find the cutoff
          curMean <- mean(tempData$rt)
          curSD <- sd(tempData$rt)
          curCutoff <- curMean + (stDev * curSD)
          
          # trim the data
          curData <- subset(tempData, tempData$rt < curCutoff)
          
          # find the average, and add to the data frame
          finalData[i, j] <- round(mean(curData$rt), digits = digits)
          
          # update condition loop counter
          j <- j + 1
        }
        
        # update participant loop counter
        i <- i + 1
      }
      
      return(finalData)
      
    }
    
    
    ### do "median"
    if(returnType == "median"){
      
      # ready the final data set
      finalData <- matrix(0, nrow = length(participant),
                          ncol = length(conditionList))
      
      # give the columns the condition names
      colnames(finalData) <- conditionList
      
      # add the participant column
      finalData <- cbind(participant, finalData)
      
      # convert to data frame
      finalData <- data.frame(finalData)
      
      # intialise looping variable for subjects
      i <- 1
      
      # loop over all subjects
      for(currSub in participant){
        
        # intialise looping variable for conditions. It starts at 2 because the
        # first column in the data file containing condition information is the
        # second one.
        j <- 2
        
        # loop over all conditions
        for(currCond in conditionList){
          
          # get the relevant data
          tempData <- subset(trimmedData, trimmedData$participant == currSub &
                               trimmedData$condition == currCond)
          
          # find the cutoff
          curMean <- mean(tempData$rt)
          curSD <- sd(tempData$rt)
          curCutoff <- curMean + (stDev * curSD)
          
          # trim the data
          curData <- subset(tempData, tempData$rt < curCutoff)
          
          # find the average, and add to the data frame
          finalData[i, j] <- round(median(curData$rt), digits = digits)
          
          # update condition loop counter
          j <- j + 1
        }
        
        # update participant loop counter
        i <- i + 1
      }
      
      return(finalData)
    }
    
    
  } # end of perCell == TRUE & perParticipant == TRUE
  
  
} # end of function

#------------------------------------------------------------------------------