##Code to compute overlap of two quasi-gaussian phenology curves.
##Author:Ian Breckheimer 
##code modified from http://stats.stackexchange.com/questions/12209/percentage-of-overlapping-regions-of-two-normal-distributions
##8 June 2014

## Returns proportion in flower from regression coefficients
## Two functions to allow functional forms to differ
f1 <- function(x,a1,b1,c1) { exp(a1 + b1 * x + c1 * x^2)/(1 + exp(a1 + b1 * x + c1 * x^2)) }
f2 <- function(x,a2,b2,c2) { exp(a2 + b2 * log(x) + c2 * log(x)^2)/(1 + exp(a2 + b2 * log(x) + c2 * log(x)^2)) }

## Returns curves as densities
f1_dens <- function(x,a1,b1,c1) { 
  y <- f1(x,a1,b1,c1)
  yi <- integrate(f1, -Inf, +Inf, a1=a1, b1=b1, c1=c1)
  return(y/yi[[1]])
}

f2_dens <- function(x,a2,b2,c2) { 
  y <- f2(x,a2,b2,c2)
  yi <- integrate(f2, 1e-10, +Inf, a2=a2, b2=b2, c2=c2)
  return(y/yi[[1]])
}

## Computes the intersection between them
min_f1f2 <- function(x, a1, b1, c1, a2, b2, c2) {
  f1 <- f1(x,a1,b1,c1) 
  f2 <- f2(x,a2,b2,c2)
  pmin(f1, f2)
}

min_f1f2_dens <- function(x, a1, b1, c1, a2, b2, c2) {
  f1 <- f1_dens(x,a1,b1,c1) 
  f2 <- f2_dens(x,a2,b2,c2)
  pmin(f1, f2)
}

## Fake coefficients
a1 <- -65
b1 <- 40
c1 <- -6

a2 <- -2
b2 <- 7
c2 <- -6
xs <- seq(0.0001,5,by=0.01)

y1s <- f1(xs,a1,b1,c1)
y2s <- f2(xs,a2,b2,c2)
y1d <- f1_dens(xs,a1,b1,c1)
y2d <- f2_dens(xs,a2,b2,c2)

yi <- min_f1f2(xs, a1, b1, c1, a2, b2, c2)
xp <- c(xs, xs[1])
yp <- c(yi, yi[1])

yid <- min_f1f2_dens(xs, a1, b1, c1, a2, b2, c2)
xpd <- c(xs, xs[1])
ypd <- c(yid, yid[1])

##Computes overlap
over <- integrate(min_f1f2, 0, Inf, a1=a1, b1=b1, c1=c1, a2=a2, b2=b2, c2=c2)
over_dens <- integrate(min_f1f2_dens, 0, Inf, a1=a1, b1=b1, c1=c1, a2=a2, b2=b2, c2=c2)

##Plots the two distributions and overlap in their proportion or density.
par(mfrow=c(1,2))
plot(xs, y1s, type="n", ylim=c(0, 1),xlab="x", ylab="Proportion")
polygon(xp, yp, col="gray")
lines(xs, y1s, lty="solid")
lines(xs, y2s, lty="dotted")
title(main=paste("Overlap: ",round(over[[1]],2)))

plot(xs, y1s, type="n", ylim=c(0, max(y1d,y2d)),xlab="x", ylab="Density")
polygon(xpd, ypd, col="gray")
lines(xs, y1d, lty="solid")
lines(xs, y2d, lty="dotted")
title(main=paste("Overlap: ",round(over_dens[[1]],2)))
