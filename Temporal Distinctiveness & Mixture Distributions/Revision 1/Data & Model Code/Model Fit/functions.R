#------------------------------------------------------------------------------
### optimisation function
fitFunction <- function(parameters, data){
  
  # extract the parameters
  p <- exp(parameters[1])
  A <- exp(parameters[2])
  bFast <- exp(parameters[3]) + A
  bSlow <- exp(parameters[4]) + A
  drift1Fast <- parameters[5]
  drift2Fast <- 1 - drift1Fast
  drift1Slow <- parameters[6]
  drift2Slow <- 1 - driftSlow
  s <- exp(parameters[7])
  ter <- exp(parameters[8])

  
  # check against minimum and maximum parameter values, and return large
  # discrepancy value if boundaries are breached
  if(p < 0 | p > 1){
    return(.Machine$double.xmax)
  }
  
  ### Do fit for ratio == 20 ####--------------------
  
  # call standard defective CDF with driftFast
  
  # evaluate defective CDF at infinity for each response accumulator
  # to get proportion of these responses
  fastOut <- list(p = numeric(2))
  fastOut$pfail <- prod(pnorm(-c(drift1Fast, drift2Fast) / s))
  
  fastOut$p[1] <- defectiveCDF(t = Inf, A = A, b = bFast, 
                               drift = c(drift1Fast, drift2Fast), s = s)
  fastOut$p[2] <- defectiveCDF(t = Inf, A = A, b = bFast, 
                               drift = c(drift2Fast, drift1Fast), s = s)
  
  
  ## do correct trials
  # evaluate defective CDF for each element of the quantiles
  fastCorrect <- defectiveCDF(t = data$fast$q[, 1] - ter, A = A, b = bFast, 
                              drift = c(drift1Fast, drift2Fast), s = s)
  # take the difference between each successive quantile value
  fastCorrect <- diff(c(0, fastCorrect, fastOut$p[1]))
  
  ## do error trials
  # evaluate defective CDF for each element of the quantiles
  fastError <- defectiveCDF(t = data$fast$q[, 2] - ter, A = A, b = bFast, 
                            drift = c(drift2Fast, drift1Fast), s = s)
  # take the difference between each successive quantile value
  fastError <- diff(c(0, fastError, fastOut$p[2]))
  
  # store model data in new data frame
  # first, replace all negative values with 1e-10
  fastCorrect[which(fastCorrect <= 1e-10)] <- 1e-10
  fastError[which(fastError <= 1e-10)] <- 1e-10
  
  fastModel <- data.frame(correct = fastCorrect, error = fastError)
  
  # multiply the number of trials in each bin by the CDF difference
  # between each quantile value, and then take the negative sum of its log
  fastFit <- -sum(data$fast$pb * log(fastModel))
  
  #--------------------------------------------------
  
  ### Do fit for ratio == 0.05 ####------------------
  
  # call standard defective CDF with driftSlow
  
  # evaluate defective CDF at infinity for each response accumulator
  # to get proportion of these responses
  slowOut <- list(p = numeric(2))
  slowOut$pfail <- prod(pnorm(-c(drift1Slow, drift2Slow) / s))
  
  slowOut$p[1] <- defectiveCDF(t = Inf, A = A, b = bSlow, 
                               drift = c(drift1Slow, drift2Slow), s = s)
  slowOut$p[2] <- defectiveCDF(t = Inf, A = A, b = bSlow, 
                               drift = c(drift2Slow, drift1Slow), s = s)
  
  
  ## do correct trials
  # evaluate defective CDF for each element of the quantiles
  slowCorrect <- defectiveCDF(t = data$slow$q[, 1] - ter, A = A, b = bSlow, 
                              drift = c(drift1Slow, drift2Slow), s = s)
  # take the difference between each successive quantile value
  slowCorrect <- diff(c(0, slowCorrect, slowOut$p[1]))
  
  ## do error trials
  # evaluate defective CDF for each element of the quantiles
  slowError <- defectiveCDF(t = data$fast$q[, 2] - ter, A = A, b = bSlow, 
                            drift = c(drift2Slow, drift1Slow), s = s)
  # take the difference between each successive quantile value
  slowError <- diff(c(0, slowError, slowOut$p[2]))
  
  # store model data in new data frame
  # first, replace all negative values with 1e-10
  slowCorrect[which(slowCorrect <= 1e-10)] <- 1e-10
  slowError[which(slowError <= 1e-10)] <- 1e-10
  
  slowModel <- data.frame(correct = slowCorrect, error = slowError)
  
  # multiply the number of trials in each bin by the CDF difference
  # between each quantile value, and then take the negative sum of its log
  slowFit <- -sum(data$slow$pb * log(slowModel))
  
  #--------------------------------------------------
  
  
  ### Do fit for ratio == 1 ####---------------------
  
  # create a new defective CDF function that takes 4 drift rates:
  # drift1Fast, drift2Fast, drift1Slow, drift2Slow
  # eventually the other fits can use this, too
  
  # evaluate defective CDF at infinity for each response accumulator
  # to get proportion of these responses
  intOut <- list(p = numeric(2))
  intOut$pfail <- prod(pnorm(-c(drift1Fast, drift2Slow, 
                                drift1Slow, drift2Slow) / s))
  
  intOut$p[1] <- defectiveCDF_Mixed(t = Inf, A = A, b = c(bFast, bSlow), 
                                    drift = c(drift1Fast, drift2Fast, 
                                              drift1Slow, drift2Slow), s = s, 
                                    p = p)
  intOut$p[2] <- defectiveCDF_Mixed(t = Inf, A = A, b = c(bFast, bSlow), 
                                    drift = c(drift2Fast, drift1Fast, 
                                              drift2Slow, drift1Slow), s = s, 
                                    p = p)
  
  
  ## do correct trials
  # evaluate defective CDF for each element of the quantiles
  intCorrect <- defectiveCDF_Mixed(t = data$int$q[, 1] - ter, A = A, 
                                   b = c(bFast, bSlow), 
                                   drift = c(drift1Fast, drift2Fast, 
                                             drift1Slow, drift2Slow), s = s, 
                                   p = p)
  # take the difference between each successive quantile value
  intCorrect <- diff(c(0, intCorrect, intOut$p[1]))
  
  ## do error trials
  # evaluate defective CDF for each element of the quantiles
  intError <- defectiveCDF_Mixed(t = data$int$q[, 2] - ter, A = A, 
                                 b = c(bFast, bSlow), 
                                 drift = c(drift2Fast, drift1Fast, 
                                           drift2Slow, drift1Slow), s = s, 
                                 p = p)
  # take the difference between each successive quantile value
  intError <- diff(c(0, intError, intOut$p[2]))
  
  # store model data in new data frame
  # first, replace all negative values with 1e-10
  intCorrect[which(intCorrect <= 1e-10)] <- 1e-10
  intError[which(intError <= 1e-10)] <- 1e-10
  
  intModel <- data.frame(correct = intCorrect, error = intError)
  
  # multiply the number of trials in each bin by the CDF difference
  # between each quantile value, and then take the negative sum of its log
  intFit <- -sum(data$int$pb * log(intModel))
  
  #--------------------------------------------------
    
  
  fit <- fastFit + slowFit + intFit

  
  return(fit)
  
} # end of function
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get starting parameters for the model 
getStartParms <- function(data){
  
  # drift rate 
  drift <- qnorm((1 - data$p), sd = 0.3 * sqrt(2)) * 2
  drift <- 0.5 + (0.5 * drift)
  drift <- (1 - drift)
  
  # ter 
  ter <- min(data$q) * .9
  
  # A
  A <- mean(as.numeric(data$q[4, ]) - as.numeric(data$q[2, ])) * 2
  
  # b
  b <- ter / 4
  
  # s
  s <- 0.3
  
  # return parameters
  parms <- c(A, b, drift, s, ter)
  
  return(parms)
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### get the defective CDF of the PDF (Equation 3) for MIXED distributions
# drift is a vector with 4 elements: 
# drift1Fast, drift2Fast, drift1Slow, drift2Slow
# b is a vector with 2 elements:
# bFast, bSlow
# p is an additional parameter controls the mixture probability
defectiveCDF_Mixed <- function(t, A, b, drift, s, p){
  
  # vector to store output. t can be a vector (i.e., of quantile values)
  output <- numeric(length(t))
  
  # what's the lower and upper limit of the integration?
  bounds <- c(0, t)
  
  # loop over each interval of t and try integration. There are error catches
  # throughout in those cases where the bounds are breached
  for(i in 1:length(t)){
    
    # set tmp to NULL. This will be replaced with a numeric value if
    # numerical integration succeeds.
    tmp <- NULL
    
    repeat{
      
      # catches if current bound is greater than next
      # (Why would this occur?)
      if(bounds[i] >= bounds[i + 1]){
        output[i] <- 0
        break
      }
      
      # try numerical integration
      tmp <- try(integrate(f = defectivePDF_Mixed, lower = bounds[i], 
                           upper = bounds[i + 1], A = A, b = b, drift = drift, 
                           s = s, p = p)$value, silent = TRUE)
      
      # check whether it worked
      if(is.numeric(tmp)){
        output[i] <- tmp
        break
      }
      
      # if not, try "smart" lower bound
      if(bounds[i] <= 0){
        bounds[i] <- max(c((mean(b) - 0.98 * A) / 
                             (max(mean(drift), drift[1]) + 2 * s), 0))
        next
      }
      
      # try "smart" upper bound
      if(bounds[i + 1] == Inf){
        bounds[i + 1] <- 0.02 * mean(b) / (mean(drift) - 2 * s)
        next
      }
      
      # if nothing works, return an error.
      stop("Error in defective CDF")
      
    } # end of repeat function
    
  } # end of loop over t
  
  # calculate cumulative value
  return(cumsum(output))
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### defective PDF for MIXED distributions
# drift is a vector with 4 elements: 
# drift1Fast, drift2Fast, drift1Slow, drift2Slow
# b is a vector with 2 elements:
# bFast, bSlow
# p is an additional parameter controls the mixture probability
defectivePDF_Mixed <- function(t, A, b, drift, s, p){
  
  # do fast distribution
  cdfFast2 <- cdf(t = t, A = A, b = b[1], drift = drift[2], s = s)
  pdfFast1 <- PDF(t = t, A = A, b = b[1], drift = drift[1], s = s)
  
  # do slow distribution
  cdfSlow2 <- cdf(t = t, A = A, b = b[2], drift = drift[4], s = s)
  pdfSlow1 <- PDF(t = t, A = A, b = b[2], drift = drift[3], s = s)
  
  # now do the likelihoods
  fast <- p * pdfFast1 * (1 - cdfFast2)
  slow <- ((1 - p) * pdfSlow1) * (1 - cdfSlow2)
  
  return(fast + slow)
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### get the defective CDF of the PDF (Equation 3)
# drift is again a vector with the drifts of all accumulators in the race
defectiveCDF <- function(t, A, b, drift, s){
  
  # vector to store output. t can be a vector (i.e., of quantile values)
  output <- numeric(length(t))
  
  # what's the lower and upper limit of the integration?
  bounds <- c(0, t)
  
  # loop over each interval of t and try integration. There are error catches
  # throughout in those cases where the bounds are breached
  for(i in 1:length(t)){
    
    # set tmp to NULL. This will be replaced with a numeric value if
    # numerical integration succeeds.
    tmp <- NULL
    
    repeat{
      
      # catches if current bound is greater than next
      # (Why would this occur?)
      if(bounds[i] >= bounds[i + 1]){
        output[i] <- 0
        break
      }
      
      # try numerical integration
      tmp <- try(integrate(f = defectivePDF, lower = bounds[i], 
                           upper = bounds[i + 1], A = A, b = b, drift = drift, 
                           s = s)$value, silent = TRUE)
      
      # check whether it worked
      if(is.numeric(tmp)){
        output[i] <- tmp
        break
      }
      
      # if not, try "smart" lower bound
      if(bounds[i] <= 0){
        bounds[i] <- max(c((b - 0.98 * A) / 
                             (max(mean(drift), drift[1]) + 2 * s), 0))
        next
      }
      
      # try "smart" upper bound
      if(bounds[i + 1] == Inf){
        bounds[i + 1] <- 0.02 * b / (mean(drift) - 2 * s)
        next
      }
      
      # if nothing works, return an error.
      stop("Error in defective CDF")
      
    } # end of repeat function
    
  } # end of loop over t
  
  # calculate cumulative value
  return(cumsum(output))
  
} # end of function
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get the defective PDF
# Equation 3 in Brown & Heathcote

# here, drift is a vector with v1 and v2
defectivePDF <- function(t, A, b, drift, s){
  
  # get the recioprocal of the CDF for the second accumulator
  tmp <- 1 - cdf(t = t, A = A, b = b, drift = drift[2], s = s)
  
  # multiply this by the PDF for the first accumulator
  tmp <- tmp * PDF(t = t, A = A, b = b, drift = drift[1], s = s)
  
  return(tmp)
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get the PDF for first-passage time 
# (related to Equation 2 Brown & Heathcote)
PDF = function(t, A, b, drift, s) {
  
  # ts
  zs <- t * s
  
  # tv
  zu <- t * drift
  
  # b - tv
  chiminuszu <- b - zu
  
  # b - tv / ts
  chizu <- chiminuszu / zs
  
  # b - tv - A / ts
  chizumax <- (chiminuszu - A) / zs
  
  return((drift * (pnorm(chizu) - pnorm(chizumax)) + 
            s * (dnorm(chizumax) - dnorm(chizu))) / A)
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get the CDF for first-passage time 
# (related to Equation 1 Brown & Heathcote)
cdf <- function(t, A, b, drift, s){
  
  # ts
  zs <-  t * s 
  
  # tv
  zu = t * drift
  
  # b - tv
  chiminuszu <- b - zu
  
  # b - A - tv
  xx <- chiminuszu - A
  
  # b - tv / ts
  chizu <- chiminuszu / zs
  
  # b - A - tv / ts
  chizumax = xx / zs
  
  tmp1 = zs *(dnorm(chizumax) - dnorm(chizu))
  tmp2 = xx * pnorm(chizumax) - chiminuszu * pnorm(chizu)
  
  return(1 + (tmp1 + tmp2) / A)
  
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### prepare the data for QMP
prepData <- function(data, quantiles = c(.1, .3, .5, .7, .9)){
  
  # get average quantile values
  q <- getQuantiles(data, quantiles)
  
  # get number of trials for each accuracy
  n <- getN(data)
  
  # get the proportion of correct responses
  p <- n$correct / (n$correct + n$error)
  
  # get the number of RTs in each bin
  bins <- c(.1, .2, .2, .2, .2, .1)
  nCorrect <- rep(n$correct, length(quantiles) + 1) * bins
  nError <- rep(n$error, length(quantiles) + 1) * bins
  pb <- data.frame(correct = nCorrect, error = nError)
  
  return(list(quantiles = quantiles, p = p, q = q, n = n, pb = pb))
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### get total number of trials for correct & error trials
getN <- function(data){
  
  # get a list of subjects
  subjects <- sort(unique(data$subject))
  
  # get data vectors
  correct <- numeric(length(subjects))
  error <- numeric(length(subjects))
  
  # loop over subjects
  for(i in 1:length(subjects)){
    
    # who is the current subject?
    curSub <- subjects[i]
    
    # get their data
    correctSub <- subset(data, subject == curSub & data$accuracy == 1)
    errorSub <- subset(data, subject == curSub & data$accuracy == 0)
    
    correct[i] <- nrow(correctSub)
    error[i] <- nrow(errorSub)
    
  }
  
  return(data.frame(correct = mean(correct), error = mean(error)))
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### get the average quantile values for each accuracy level
getQuantiles <- function(data, quantiles){
  
  # get a list of subjects
  subjects <- sort(unique(data$subject))
  
  # get a matrix for the CDFs for each subject for each accuracy level
  corCDF <- matrix(0, nrow = length(subjects), ncol = length(quantiles))
  errCDF <- matrix(0, nrow = length(subjects), ncol = length(quantiles))
  
  # loop over each subject and populate the matrices
  for(i in 1:length(subjects)){
    
    # get the current subject
    curSub <- subjects[i]
    
    # get their data
    correctData <- subset(data, data$subject == curSub & data$accuracy == 1)
    errorData <- subset(data, data$subject == curSub & data$accuracy == 0)
    
    # populate the matrices
    corCDF[i, ] <- as.numeric(quantile(correctData$rt, probs = quantiles))
    
    if(nrow(errorData) >= 1){
      errCDF[i, ] <- as.numeric(quantile(errorData$rt, probs = quantiles))
    }
    
  }
  
  corCDF <- apply(corCDF, 2, mean)
  errCDF <- apply(errCDF, 2, mean)
  
  return(data.frame(correct = corCDF, error = errCDF))
  
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### generate data from the LBA (for an example)
simLBA_Mixture <- function(n, bs, a = 1, drift, s, ter, p){
  
  # initialise a matrix to store the data
  data <- matrix(0, nrow = n, ncol = 2)
  
  # loop over every trial and generate an RT and accuracy
  for(i in 1:n){
    
    m <<- i
    
    # get the distribution to sample from
    # 1 = Fast, 2 = Slow
    position <- sample(c(1, 2), 1, prob = c(p, 1 - p), replace = TRUE)
    
    # get the drift rate for the current trial, as a function of mixture prob p
    v <- drift[position]
    
    # get b for the current trial, as a function of mixture prob p
    b <- bs[position]
    
    # accumulator for correct response
    a1 <- (b - runif(1, min = 0, max = a)) / rnorm(1, mean = v, sd = s) 
    # accumulator for error response
    a2 <- (b - runif(1, min = 0, max = a)) / rnorm(1, mean = 1 - v, sd = s)
    
    # if both accumulators are negative, the following code re-generates
    # finishing times until one is positive
    while(max(c(a1, a2)) < 0){
      # accumulator for correct response
      a1 <- (b - runif(1, min = 0, max = a)) / rnorm(1, mean = v, sd = s) 
      # accumulator for error response
      a2 <- (b - runif(1, min = 0, max = a)) / rnorm(1, mean = (1 - v), sd = s)
    }
    
    # if the losing accumulator is still negative, set it to large 
    # value so it won't be selected
    if(a1 < 0){
      a1 <- 1e10
    }
    if(a2 < 0){
      a2 <- 1e10
    }
    
    # which accumulator is the winner?
    response <- which.min(c(a1, a2))
    
    # add ter to the winning accumulator's finishing time, and store it
    data[i, 1] <- round(min(c(a1, a2)), 3) + ter
    
    # add the accuracy
    if(response == 1){
      data[i, 2] <- 1
    } else {
      data[i, 2] <- 0
    }
    
  } # end of trial loops
  
  colnames(data) <- c("rt", "accuracy")
  data <- data.frame(data)
  
  data <- mutate(data, subject = 1)
  
  return(data)
  
} # end of function
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