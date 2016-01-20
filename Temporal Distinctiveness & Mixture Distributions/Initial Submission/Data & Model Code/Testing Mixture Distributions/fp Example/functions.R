#------------------------------------------------------------------------------
### fixed-point code adapted from the "fp" package
# http://www.leendertvanmaanen.com/resources/fp/example.html

fp1 <- setClass("fp1", representation(dens="array", diff="data.frame", 
                                      dat="data.frame"))

plot.fp1 <- function(x, ..., ylab=c("Density","Density difference"), 
                     xlim=NULL) {
  
  par(mfrow = c(1, 2))
  options(warn=-1)
  if (is.null(xlim)) xlim <- range(x@dens[[1]]$x)
  plot(1,type="n",
       ylim=c(0,max(unlist(lapply(x@dens, function(X) {max(X$y)})))),
       xlim=xlim,
       bty="n", ylab=ylab[1], ...)
  for (i in 1:length(x@dens)) {
    lines(x@dens[[i]], col=i)
  }
  legend("topright", c("x", "z", "y"), lty = c(1, 1, 1), 
         col = c(1, 2, 3), bty = "n")
  
  matplot(x@dens[[1]]$x,x@diff, type='l', pch=16,
          lwd=3, bty="n", ylab=ylab[2], xlim=xlim, ...); abline(h=0, lty=1)
  legend("topleft", c("z-x", "y-z", "y-x"), bty = "n", 
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