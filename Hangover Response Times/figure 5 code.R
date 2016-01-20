# clear workspace:
rm(list=ls(all=TRUE))

par(mfrow=c(1,1))

library(retimes)

set.seed(34) #this ensures you get the same plot as in the paper. Remove this line
#if you want a new random estimate!
x = rexgauss(100000, mu = 397, sigma = 42, tau = 127) #control
y = rexgauss(100000, mu = 401, sigma = 47, tau = 176) #hangover


#plot the density function
plot(density(x), xaxt="n", yaxt="n", xlim=c(0, 2000), main="", xlab=("Response Time (ms)"), lwd=4)
axis(1, xaxp=c(0, 2000, 10))
#add the hangover density function
lines(density(y), lwd=4, lty=2, col="black")

#add the legend
legend("topleft", c("Control", "Hangover"), lty=c(1, 2), lwd=c(3,3),
       col=c("black", "black"), bty="n")