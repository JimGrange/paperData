#Custom functions for various anlytical steps in hangover paper. 
#Written by Jim Grange, September 2014. 
#Any questions, please email j.a.grange@keele.ac.uk or grange.jim@gmail.com

#----------------------------------------------------------------#
#check which subjects have complete data sets
subjectChecks <- function(data, conditions){
  
  #first, check I have the same subject numbers in each condition
  hangoverData <- subset(data, data$Condition==conditions[1])
  controlData <- subset(data, data$Condition==conditions[2])
  
  #find the unique subject numbers
  hangoverSubjects <- unique(hangoverData$Subject)
  controlSubjects <- unique(controlData$Subject)
  
  
  #only select subject numbers that are present in BOTH data sets. 
  completeSubjects <-NULL
  for(i in 1:max(data$Subject)){ #for all possible subject numbers....
    
    if(i %in% hangoverSubjects & i %in% controlSubjects) #if the current subject number is present in BOTH data sets
    {completeSubjects = c(completeSubjects, i)} #store this subject number in completeSubjects vector
    
  }  
  #return a vector with complete subject numbers
  return(completeSubjects)
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#get mean RTs for each condition
getRTs <- function(data, completeSubjects, decimals){
  
  #how many are there?
  nSubjects <- length(completeSubjects)
  
  #empty matrix for mean response times
  rts <- matrix(0, nrow = nSubjects, ncol = 3) #cols = subject number, hungover rts, control rts
  colnames(rts) <- c("subject", "hungover", "control")
  
  
  ##############
  #loop over every subject, get the mean rt for each condition
  
  #reset counter for subjects
  i = 1
  
  #for each subject...
  for(sub in completeSubjects){
    
    #update Matrix with subject number
    rts[i, 1] <- sub
    
    #get that subject's data
    subData <- subset(data, data$Subject==sub)
    
    
    #get hangover data
    hungoverData <- subset(subData, subData$Condition=="hungover")
    
    
    #what was their mean response time?
    rtData <- subset(hungoverData, hungoverData$Accuracy==1) #only select correct trials
    rts[i, 2] <- round(mean(rtData$RT), digits = decimals)
    
    #get control data
    controlData <- subset(subData, subData$Condition=="control")
    
    #what was their mean response time?
    rtData <- subset(controlData, controlData$Accuracy==1)
    rts[i, 3] <- round(mean(rtData$RT), digits = decimals)
    
    #update counter
    i = i + 1                  
  }
  
  return(rts)
  
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#get median RTs for each condition
getMedianRTs <- function(data, completeSubjects){
  
  #how many are there?
  nSubjects <- length(completeSubjects)
  
  #empty matrix for mean response times
  rts <- matrix(0, nrow = nSubjects, ncol = 3) #cols = subject number, hungover rts, control rts
  colnames(rts) <- c("subject", "hungover", "control")
  
  
  ##############
  #loop over every subject, get the mean rt for each condition
  
  #reset counter for subjects
  i = 1
  
  #for each subject...
  for(sub in completeSubjects){
    
    #update Matrix with subject number
    rts[i, 1] <- sub
    
    #get that subject's data
    subData <- subset(data, data$Subject==sub)
    
    
    #get hangover data
    hungoverData <- subset(subData, subData$Condition=="hungover")
    
    
    #what was their median response time?
    rtData <- subset(hungoverData, hungoverData$Accuracy==1) #only select correct trials
    rts[i, 2] <- round(median(rtData$RT), digits = 0)
    
    #get control data
    controlData <- subset(subData, subData$Condition=="control")
    
    #what was their median response time?
    rtData <- subset(controlData, controlData$Accuracy==1)
    rts[i, 3] <- round(median(rtData$RT), digits = 0)
    
    #update counter
    i = i + 1                  
  }
  
  return(rts)
  
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#get standard deviation of RTs
getSDs <- function(data, completeSubjects){
 
  #how many are there?
  nSubjects <- length(completeSubjects)
  
  #empty matrix for mean response times
  sds <- matrix(0, nrow = nSubjects, ncol = 3) #cols = subject number, hungover rts, control rts
  colnames(sds) <- c("subject", "hungover", "control")
  
  ##############
  #loop over every subject, get the sd rt for each condition
  
  #reset counter for subjects
  i = 1
  
  #for each subject...
  for(sub in completeSubjects){
    
    #update Matrix with subject number
    sds[i, 1] <- sub
    
    #get that subject's data
    subData <- subset(data, data$Subject==sub)
    
    
    #get hangover data
    hungoverData <- subset(subData, subData$Condition=="hungover")
    
    
    #what was their sd response time?
    sdData <- subset(hungoverData, hungoverData$Accuracy==1) #only select correct trials
    sds[i, 2] <- round(sd(sdData$RT), digits = 2)
    
    #get control data
    controlData <- subset(subData, subData$Condition=="control")
    
    #what was their mean response time?
    sdData <- subset(controlData, controlData$Accuracy==1)
    sds[i, 3] <- round(sd(sdData$RT), digits = 2)
    
    #update counter
    i = i + 1                  
  }
  
  return(sds)
  
}

#----------------------------------------------------------------#
#get the variance of RTs
getVariance <- function(data, completeSubjects){
  
  #how many are there?
  nSubjects <- length(completeSubjects)
  
  #empty matrix for mean response times
  variance <- matrix(0, nrow = nSubjects, ncol = 3) #cols = subject number, hungover rts, control rts
  colnames(variance) <- c("subject", "hungover", "control")
  
  ##############
  #loop over every subject, get the sd rt for each condition
  
  #reset counter for subjects
  i = 1
  
  #for each subject...
  for(sub in completeSubjects){
    
    #update Matrix with subject number
    variance[i, 1] <- sub
    
    #get that subject's data
    subData <- subset(data, data$Subject==sub)
    
    
    #get hangover data
    hungoverData <- subset(subData, subData$Condition=="hungover")
    
    
    #what was their sd response time?
    varData <- subset(hungoverData, hungoverData$Accuracy==1) #only select correct trials
    variance[i, 2] <- round(var(varData$RT), digits = 4)
    
    #get control data
    controlData <- subset(subData, subData$Condition=="control")
    
    #what was their mean response time?
    varData <- subset(controlData, controlData$Accuracy==1)
    variance[i, 3] <- round(var(varData$RT), digits = 4)
    
    #update counter
    i = i + 1                  
  }
  
  return(variance)
  
}
#----------------------------------------------------------------#



#----------------------------------------------------------------#
#get mean Errors for each condition
getErrors <- function(data, completeSubjects){
  
  #how many are there?
  nSubjects <- length(completeSubjects)
  
  #empty matrix for mean response times
  errors <- matrix(0, nrow = nSubjects, ncol = 3) #cols = subject number, hungover rts, control rts
  colnames(errors) <- c("subject", "hungover", "control")
  
  
  ##############
  #loop over every subject, get the mean rt for each condition
  
  #reset counter for subjects
  i = 1
  
  #for each subject...
  for(sub in completeSubjects){
    
    #update Matrix with subject number
    errors[i, 1] <- sub
    
    #get that subject's data
    subData <- subset(data, data$Subject==sub)
    
    
    #get hangover data
    hungoverData <- subset(subData, subData$Condition=="hungover")
    
    #what was their accuracy? Store it in the matrix
    errors[i, 2] <- round(((sum(hungoverData$Accuracy) / nrow(hungoverData)) * 100), digits = 2)
    
    #get control data
    controlData <- subset(subData, subData$Condition=="control")
    
    #what was their accuracy? Store it in the matrix
    errors[i, 3] <- round(((sum(controlData$Accuracy) / nrow(controlData)) * 100), digits = 2)
    
    #update counter
    i = i + 1                  
  }
  
  return(errors)
  
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#find subject lines who don't meet accuracy criterion
accuracyTrimming <- function(errors, criterion){
    
  removed <- which(meanErrors[, 3] < criterion | meanErrors[, 2] < criterion)
  
  return(removed)
  
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#Run BEST Bayesian analysis
doBEST <- function(y1, y2, dv){
  
  filename <- paste(dv, "BEST.Rdata")
  
  mcmcChain = BESTmcmc(y1 , y2) 
  
  # Plot the results of the Bayesian analysis:
  postInfo = BESTplot( y1 , y2 , mcmcChain , pairsPlot=FALSE)
  
  # Show detailed summary info on console:
  show(postInfo) 
  # You can save the plot(s) using the pull-down menu in the R graphics window,
  # or by using the following:
  # saveGraph( file="BESTexample" , type="eps" )
   saveGraph( file=dv, type="eps" )
    
  # Save the data and results for future use:
  save( y1, y2, mcmcChain, postInfo, file=filename )
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#do ex-Gaussian analysis
doExG <- function(data, subjects){
  
  #get empty matrix to store data
  exgData <- matrix(nrow = length(subjects), ncol = 4)
    colnames(exgData) <- c("subject", "mu", "sigma", "tau")
  
  #reset looping counter
  i = 1
  
  #loop over subjects, get correct RTs, and find exG parameters
  for(currSub in subjects){
    
    #store subject number
    exgData[i, 1] <- currSub
    
    #get current subject's data
    currData <- subset(data, data$Subject == currSub)
    
    #only include correct trials
    currData <- subset(currData, currData$Accuracy == 1)
    
    #get the RTs
    currRTs <- currData$RT
    
    #do the ex gaussian fitting
    currParameters <- timefit(currRTs, iter = 1000)
    
    #store the results (rounded to nearest millisecond)
    exgData[i, 2:4] <- round(currParameters@par, 0)
    
    #update looping counter
    i = i + 1
  }
  
  return(exgData)
  
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#do Wiender diffusion analysis
doWiener <- function(data, subjects){
  
  pastDiff <- c(0,0,0,0)
  
  #create empty matrix
  finalData <- matrix(nrow = length(subjects), ncol = 5)
  colnames(finalData) <- c("subject", "boundary", "nonDecision", "starting", "drift")
  
  #initiate indexing
  i = 1
  
  #do each subject
  for(currSub in subjects){
    

    
    finalData[i, 1] <- currSub
    
    #get current subject's data
    subData <- subset(data, data$Subject==currSub)
    
    #remove subject column
    subData <- subData[, -1]
    
    #estimate parameters
    try(diff <- optim(c(1, .1, .1, 1), wiener_deviance, dat=subData, method="Nelder-Mead"))
  
    
    comparison <- sum(pastDiff - diff$par)
    
    
    #store best fitting parameters
    finalData[i, 2:5] <- diff$par
    
      if(comparison==0){finalData[i, 2:5] <- c(0,0,0,0)}
    
    i = i + 1
    
    #pass parameters to new parameter, for comparison
    pastDiff <- diff$par
    
  }
  
  return(finalData)
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#do RWiener with bias fixed to 0.5 (probably sensible)
doWiener_fixedBias <- function(data, subjects){
  
  
  #first, here's the new wiener function
  wiener_deviance2 <- function(x, dat) {
    beta <- 0.5
    wiener_deviance(c(x[1],x[2],beta,x[3]), dat)
  }
  
  pastDiff <- c(0,0,0)
  
  #create empty matrix
  finalData <- matrix(nrow = length(subjects), ncol = 4)
  colnames(finalData) <- c("subject", "boundary", "nonDecision", "drift")
  
  #initiate indexing
  i = 1
  
  #do each subject
  for(currSub in subjects){
    
    
    
    finalData[i, 1] <- currSub
    
    #get current subject's data
    subData <- subset(data, data$Subject==currSub)
    
    #remove subject column
    subData <- subData[, -1]
    
    #estimate parameters
    try(diff <- optim(c(1, .1, 1), wiener_deviance2, dat=subData, method="Nelder-Mead"))
    try(diff <- optim(diff$par, wiener_deviance2, dat=subData, method = "BFGS", hessian=TRUE))
    
    comparison <- sum(pastDiff - diff$par)
    
    
    #store best fitting parameters
    finalData[i, 2:4] <- diff$par
    
    if(comparison==0){finalData[i, 2:4] <- c(0,0,0)}
    if(diff$convergence>0){finalData[i, 2:4] <- c(0,0,0)}
    
    i = i + 1
    
    #pass parameters to new parameter, for comparison
    pastDiff <- diff$par
    
  }
  
  return(finalData) 
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#take mean data, return melted data frame ready for ggplotting
doSummary <- function(data, id.vars, variable.name, value.name){
  
  meltedData <- melt(data, id.vars = id.vars, variable.name = variable.name, 
                     value.name = value.name)
  
  #how many subjects?
  N <- nrow(meltedData) / 2
  
  if(value.name=="RT"){
  
    mean <- select(meltedData, Condition, RT) %>%
    group_by(Condition) %>%
    summarise(mean(RT))
  
    se <- select(meltedData, Condition, RT) %>%
      group_by(Condition) %>%
      summarise(sd(RT) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
  
  }
 
  if(value.name=="Accuracy"){
    
    mean <- select(meltedData, Condition, Accuracy) %>%
      group_by(Condition) %>%
      summarise(mean(Accuracy))
    
    se <- select(meltedData, Condition, Accuracy) %>%
      group_by(Condition) %>%
      summarise(sd(Accuracy) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  
  if(value.name=="SD"){
    
    mean <- select(meltedData, Condition, SD) %>%
      group_by(Condition) %>%
      summarise(mean(SD))
    
    se <- select(meltedData, Condition, SD) %>%
      group_by(Condition) %>%
      summarise(sd(SD) / sqrt(N))
    
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
  }
  
  if(value.name=="Mu"){
    
    mean <- select(meltedData, Condition, Mu) %>%
      group_by(Condition) %>%
      summarise(mean(Mu))
    
    se <- select(meltedData, Condition, Mu) %>%
      group_by(Condition) %>%
      summarise(sd(Mu) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  if(value.name=="Sigma"){
    
    mean <- select(meltedData, Condition, Sigma) %>%
      group_by(Condition) %>%
      summarise(mean(Sigma))
    
    se <- select(meltedData, Condition, Sigma) %>%
      group_by(Condition) %>%
      summarise(sd(Sigma) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  if(value.name=="Tau"){
    
    mean <- select(meltedData, Condition, Tau) %>%
      group_by(Condition) %>%
      summarise(mean(Tau))
    
    se <- select(meltedData, Condition, Tau) %>%
      group_by(Condition) %>%
      summarise(sd(Tau) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  
  if(value.name=="Drift"){
    
    mean <- select(meltedData, Condition, Drift) %>%
      group_by(Condition) %>%
      summarise(mean(Drift))
    
    se <- select(meltedData, Condition, Drift) %>%
      group_by(Condition) %>%
      summarise(sd(Drift) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  if(value.name=="Boundary"){
    
    mean <- select(meltedData, Condition, Boundary) %>%
      group_by(Condition) %>%
      summarise(mean(Boundary))
    
    se <- select(meltedData, Condition, Boundary) %>%
      group_by(Condition) %>%
      summarise(sd(Boundary) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  if(value.name=="nonDecision"){
    
    mean <- select(meltedData, Condition, nonDecision) %>%
      group_by(Condition) %>%
      summarise(mean(nonDecision))
    
    se <- select(meltedData, Condition, nonDecision) %>%
      group_by(Condition) %>%
      summarise(sd(nonDecision) / sqrt(N))
    
    Condition <- levels(mean$Condition)
    mean <- mean$mean
    se <- se$sd
    
    finalData <- data.frame(Condition, mean, se)
    
    return(finalData)
    
  }
  
  

}
#----------------------------------------------------------------#



#----------------------------------------------------------------#
#do EZ modelling
doEZ <- function(meanRT, varRT, pCorrect, completeSubjects){
  
  #initialise empty matrix to store parameters
  ezParameters <- matrix(nrow = length(completeSubjects), ncol = 7)
  colnames(ezParameters) <- c("Subject", "Drift_H", "Boundary_H", "Non-Decision_H", 
                              "Drift_C", "Boundary_C", "Non-Decision_C")
  

  
  
  #start loop for model fitting
  for(i in 1:length(completeSubjects)){
    
    #store subject number
    ezParameters[i, 1] <- meanRT[i, 1]
    
    #Do hangover condition first
    mRT <- as.numeric(meanRT[i, 2])
    vRT <- as.numeric(varRT[i, 2])
    pC <- as.numeric(pCorrect[i, 2])
    
    #run ez model
    ez <- ezIndividual(pC, vRT, mRT) 
    ez <- as.numeric(ez)
    
    #store parameters
    ezParameters[i, 2:4] <- ez
    
    
    #now do control condition
    mRT <- meanRT[i, 3]
    vRT <- varRT[i, 3]
    pC <- pCorrect[i, 3]
    
    #run ez model
    ez <- ezIndividual(pC, vRT, mRT) 
    ez <- as.numeric(ez)
    
    #store parameters
    ezParameters[i, 5:7] <- ez
    
  }
  
  return(ezParameters)
}
#----------------------------------------------------------------#


#----------------------------------------------------------------#
#provide EZ estimates 
#Written by E-J Wagenmakers
ezIndividual = function(Pc, VRT, MRT, s=.1)
{
  s2 = s^2
  # The default value for the scaling parameter s equals .1
  
  if (Pc == 0)
    cat("Oops, Pc == 0!\n")
  if (Pc == 0.5)
    cat("Oops, Pc == .5!\n")
  if (Pc == 1)
    cat("Oops, Pc == 1!\n")
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge-correction is required.
  
  L = qlogis(Pc)
  # The function "qlogis" calculates the logit.
  x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
  v = sign(Pc-0.5)*s*x^(1/4)
  # This gives drift rate.
  
  a = s2*qlogis(Pc)/v
  # This gives boundary separation.
  
  y   = -v*a/s2
  MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
  Ter = MRT-MDT
  # This gives nondecision time.
  
  return(list(v, a, Ter))
}
#----------------------------------------------------------------#


#------------------------------------------------------------------------------
get.vaTer = function(Pc, VRT, MRT, s=.1)
{
  s2 = s^2
  # The default value for the scaling parameter s equals .1
  
  if (Pc == 0)
    cat("Oops, Pc == 0!\n")
  if (Pc == 0.5)
    cat("Oops, Pc == .5!\n")
  if (Pc == 1)
    cat("Oops, Pc == 1!\n")
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge-correction is required.
  
  L = qlogis(Pc)
  # The function "qlogis" calculates the logit.
  x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
  v = sign(Pc-0.5)*s*x^(1/4)
  # This gives drift rate.
  
  a = s2*qlogis(Pc)/v
  # This gives boundary separation.
  
  y   = -v*a/s2
  MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
  Ter = MRT-MDT
  # This gives nondecision time.
  
  return(c(v, a, Ter))
}
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

# revised includes parameter p: mix proportion, mixture of Exg & Fixed Unif.

erf <- function(x) { 2 * pnorm(x * sqrt(2)) - 1 }

# exg
exgf <- function(tau,xm,s,t) {
  t1 = 1.0/tau
  t4 = s^2
  t5 = tau^2
  tmp=(xm-t)*t1+t4/t5/2.0
  t10 = exp(tmp)
  t17 = sqrt(2.0)
  t21 = erf((tau*t-tau*xm-t4)/s*t1*t17/2.0)
  t1*t10*(t21+1.0)/2.0 
  
}
# deriv m
exgfdm <-function(tau,xm,s,t) {
  t1 = tau*t
  t3 = tau*xm
  t5 = s^2
  t7 = tau^2
  t8 = 1.0/t7
  tmp=-(2.0*t1-2.0*t3-t5)*t8/2.0
  t11 = exp(tmp)
  t12 = sqrt(pi)
  t13 = t12*s
  t14 = t1-t3-t5
  t15 = 1.0/s
  t18 = sqrt(2.0)
  t22 = erf(t14*t15/tau*t18/2.0)
  t24 = t14^2
  tmp=(-t24/t5*t8/2.0)
  t29 = exp(tmp)
  -t11*(-t13*t22-t13+t29*t18*tau)*t8/t12*t15/2.0
}
# deriv s
exgfds<-function(tau,xm,s,t) {
  t1 = tau*t
  t3 = tau*xm
  t5 = s^2
  t7 = tau^2
  t8 = 1.0/t7
  tmp=(-(2.0*t1-2.0*t3-t5)*t8/2.0)
  t11 = exp(tmp)
  t13 = sqrt(pi)
  t14 = t5*s*t13
  t15 = t1-t3-t5
  t19 = sqrt(2.0)
  t23 = erf(t15/s/tau*t19/2.0)
  t25 = t15^2
  t26 = 1.0/t5
  tmp=(-t25*t26*t8/2.0)
  t30 = exp(tmp)
  t31 = t30*t19
  -t11*(-t14*t23-t14+t31*tau*t5+t31*t7*t-t31*t7*xm)/t7/tau/t13*t26/2.0
}
# deriv t
exgfdt<-function(tau,xm,s,t) {
  t1 = tau*t
  t3 = tau*xm
  t5 = s^2
  t7 = tau^2
  t8 = 1.0/t7
  tmp = (-(2.0*t1-2.0*t3-t5)*t8/2.0)
  t11 = exp(tmp)
  t12 = sqrt(pi)
  t13 = t7*t12
  t14 = t1-t3-t5
  t18 = sqrt(2.0)
  t22 = erf(t14/s/tau*t18/2.0)
  t24 = t12*tau
  t31 = t12*t5
  t33 = t14^2
  tmp = (-t33/t5*t8/2.0)
  t38 = exp(tmp)
  t44 = t7^2
  t11*(-t13*t22-t13+t24*t*t22+t24*t-t24*xm*t22-t24*xm-t31*t22-t31+t38*t18*s*tau)/t44/t12/2.0
}

# exgauss dens
exgauss<-function(xm,s,tau,t) {exgf(tau,xm,s,t)}

# uniform density
unif<-function(a,b,t) {(t-t+1)/(b-a)} # t-t+1 = 1 nonsense. but useful when t is a vector

# mixture
# log likelihood
loglexgm<-function(par,x,y)
{
  m=par[1]
  s=par[2]
  t=par[3]
  p1=par[4]
  p2=1-p1
  a=y[1]
  b=y[2]
  L=-sum(log(exgauss(m,s,t,x)*p1+unif(a,b,x)*p2))
  print(L)
  L
}


postp<-function(par,x,y)
{
  m=par[1]
  s=par[2]
  t=par[3]
  a=y[1]
  b=y[2]
  p1=par[4]
  p2=1-p1
  xtmp=exgauss(m,s,t,x)*p1+unif(a,b,x)*p2
  xtmp1=exgauss(m,s,t,x)
  xtmp2=unif(a,b,x)
  cbind(p1*xtmp1/xtmp,p2*xtmp2/xtmp)
}



# mixture
# derivatives
loglexggm<-function(par,x,y)
{
  g=c(0,0,0,0)
  m=par[1]
  s=par[2]
  t=par[3]
  a=y[1]
  b=y[2]
  p1=par[4]
  p2=1-p1
  xtmp=exgauss(m,s,t,x)*p1+unif(a,b,x)*p2
  xtmp1=exgauss(m,s,t,x)
  xtmp2=unif(a,b,x)
  g[1]=-sum(p1*exgfdm(t,m,s,x)/xtmp)       
  g[2]=-sum(p1*exgfds(t,m,s,x)/xtmp) 
  g[3]=-sum(p1*exgfdt(t,m,s,x)/xtmp)
  g[4]=-sum((xtmp1-xtmp2)/xtmp)
  g
}

# mixture
loglexg2m<-function(par,x,y)
{
  m=par[1]
  s=par[2]
  t=par[3]
  a=y[1]
  b=y[2]
  p1=par[4]
  p2=1-p1
  L=-sum( log(exgauss(m,s,t,x)*p1+unif(a,b,x)*p2))
  g=loglexgg(par,x)
  attr(L,"gradient")=g
  L
}

# only exg
# log likelihood
loglexg<-function(par,x)
{
  m=par[1]
  s=par[2]
  t=par[3]
  -sum(log(exgauss(m,s,t,x)))
}
# only exg
# derivatives
loglexgg<-function(par,x)
{
  g=c(0,0,0)
  m=par[1]
  s=par[2]
  t=par[3]
  xtmp=exgf(t,m,s,x)
  g[1]=-sum(exgfdm(t,m,s,x)/xtmp)       
  g[2]=-sum(exgfds(t,m,s,x)/xtmp) 
  g[3]=-sum(exgfdt(t,m,s,x)/xtmp)
  g
}

# only exg
loglexg2<-function(par,x)
{
  m=par[1]
  s=par[2]
  t=par[3]
  L=-sum(log(exgauss(m,s,t,x)))
  g=loglexgg(par,x)
  attr(L,"gradient")=g
  L
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
source("fdexgun.R") # the disfit routines; these should be in your working directory

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
RobustEZ.from.File = function(file.name)
{
  #========================================================================== 
  # Read in data from file.name
  #==========================================================================
  
  d    = scan(file = file.name)
  N    = d[1]
  rt   = d[2:length(d)]
  ac   = length(rt)/N
  
  pars = Get.Robust.vaTer(rt, ac, min_U=min(rt), max_U=max(rt),start_p_EG=.95)
  
  return(pars)
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Get.Robust.vaTer = function(dat, Pc, min_U=min(dat), max_U=max(dat), start_p_EG=.95, s=.1)
{
  #========================================================================== 
  # Get cleaned-up values for EZ pars
  #==========================================================================
  meanvarp = Get.Robust.MV(dat, min_U, max_U, start_p_EG)
  MRT      = meanvarp[1]
  VRT      = meanvarp[2]
  p_EG	   = meanvarp[3]
  # MRT & VRT cleaned up, now plug these in to EZ equations
  s2 = s^2
  # The default value for the scaling parameter s equals .1
  if (Pc == 0)
    cat("Oops, Pc == 0!\n")
  if (Pc == 0.5)
    cat("Oops, Pc == .5!\n")
  if (Pc == 1)
    cat("Oops, Pc == 1!\n")
  # If Pc equals 0, .5, or 1, the method will not work, and
  # an edge-correction is required.
  L = qlogis(Pc)
  # The function "qlogis" calculates the logit.
  x = L*(L*Pc^2 - L*Pc + Pc - 0.5)/VRT
  v = sign(Pc-0.5)*s*x^(1/4)
  # This gives drift rate.
  a = s2*qlogis(Pc)/v
  # This gives boundary separation.
  y = -v*a/s2
  MDT = (a/(2*v))*(1-exp(y))/(1+exp(y))
  Ter = MRT-MDT
  # This gives nondecision time.
  return(c(v, a, Ter, p_EG))
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Get.Robust.MV = function(dat, min_U, max_U, start_p_EG)
{
  #========================================================================== 
  # Get cleaned-up values for RT mean and RT variance
  #==========================================================================
  pars     = Est.ExGaussUnif(dat, min_U, max_U, start_p_EG)
  EG_mean  = pars[1] + pars[3]
  EG_var   = pars[2]^2 + pars[3]^2
  meanvarp = c(EG_mean, EG_var, pars[4])
  return(meanvarp)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=0.1; max_U=1.2
# dat  = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# Est.ExGaussUnif(dat,.1,1.2,.9)
# Get.Robust.MV(dat, .1,1.2,.9)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Est.ExGaussUnif = function(dat, min_U, max_U, start_p_EG)
{
  #========================================================================== 
  # Get Estimates for Ex-Gauss/Unif RT distribution
  #==========================================================================
  start.pars = c(MMest.ExGauss(dat), start_p_EG)
  
  lo   = c(0.1, .010, .010, 0.60)      # lower bounds; note bound on mix proportion 
  up   = c(1,   .300, .300, 0.9999999) # upper bounds
  y    = c(min_U,max_U) # fixed unif pars
  
  # quasi newton
  res = try(optim(start.pars, loglexgm, loglexggm, method="L-BFGS-B",
                  lower=lo, upper=up, hessian=T, x=dat, y=y), TRUE)
  
  if (class(res) == "try-error") #optimization failed
  {
    pars = c(NA, NA, NA, NA)   
  }
  
  if (class(res) != "try-error") #optimization worked
  {
    pars = res$par[1:4]
  }
  
  return(pars)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=.1; max_U=1.2
# dat = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# Est.ExGaussUnif(dat,.1,1.2,.9)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
MMest.ExGauss = function(dat)
{
  #========================================================================== 
  # Get Method of Moments Estimates for ExGauss distribution
  # cf. Heathcote (1996)
  #==========================================================================
  M1    = mean(dat)
  M2    = var(dat) #division by n-1
  M3    = (length(dat)/(length(dat)-1)) * mean ( (dat-M1)^3 ) #effectively division by n-1
  
  if (M3 >= 0)
    tau  = (0.5*M3)^(1/3)
  if (M3 < 0)
    tau  = 0
  mu    = M1 - tau
  if ( (M2 - tau^2)>0 )
  {
    sigma = sqrt(M2 - tau^2) 
  }
  if ( ((M2 - tau^2)<=0) || (sigma<=0) || (tau<=0) ) #in case of weird values
  {
    tau   = 0.8 * sd(dat)
    mu    = M1 - tau
    #sigma = sqrt(M2 - tau^2) # this is the same as:
    sigma = 0.6 * sd(dat)	    # (correcting a small mistake in Heathcote 1996)
  }
  pars = c(mu, sigma, tau)
  return(pars)
}
# to check:
# mu=.4; sigma=.035; tau=.065
# dat = Gen.ExGauss(1000, mu, sigma, tau)
# MMest.ExGauss(dat)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Gen.ExGaussUnif = function(N, mu, sigma, tau, p_EG, min_U, max_U)
{
  #========================================================================== 
  # Generate data from ExGauss distribution Contaminated with a Uniform 
  #==========================================================================
  EG = rnorm(round(p_EG*N), mean=mu, sd=sigma) +	rexp(round(p_EG*N), rate=1/tau)	
  Un = runif(N-round(p_EG*N), min=min_U, max=max_U)
  dat = c(EG,Un)
  
  return(dat)
}
# to check:
# mu=.4; sigma=.035; tau=.065; p_EG=.8; min_U=.1; max_U=1.2
# dat = Gen.ExGaussUnif(1000, mu, sigma, tau, p_EG, min_U, max_U)
# plot(dat)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
Gen.ExGauss = function(N, mu, sigma, tau)
{
  #========================================================================== 
  # Generate data from ExGauss distribution
  #==========================================================================
  dat = rnorm(N, mean=mu, sd=sigma) +	rexp(N, rate=1/tau)	
  return(dat)
}
# to check:
# mu=.4; sigma=.035; tau=.065
# dat = Gen.ExGauss(1000, mu, sigma, tau)
# mu+tau # theoretical mean
# mean(dat)
# sigma^2 + tau^2 # theoretical variance
# var(dat)
# 2*tau^3 # theoretical 3rd moment
# (length(dat)/(length(dat)-1)) * mean ( (dat-mean(dat))^3 )
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#------------------------------------------------------------------------------