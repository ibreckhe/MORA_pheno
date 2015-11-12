##Script to estimate overlap coefficients for visitor and flower phenology
##Ian Breckheimer
##November 13th, 2015

####Sets up workspace####
library(dplyr)
library(ggmcmc)

##Brings in fit model objects and converts them to data frames.
load("./scratch/visitors_jags_output_2011_2014.Rdata")
visit_model <- out2
load("./scratch/jags_output_2011_2014.Rdata")
flower_model <- out2

##Converts them to data frames
visit_samples <- filter(ggs(visit_model),Parameter %in% c("height_int","height_slope",
                                         "opt_int","opt_slope",
                                         "width_int","width_slope"))
visit_samples$Parameter <- factor(visit_samples$Parameter)
visit_params <- tidyr::spread(visit_samples,Parameter,value)
visit_means <- colMeans(visit_params)[3:8]

flower_samples <- filter(ggs(flower_model),Parameter %in% c("height_int","height_slope",
                                                          "opt_int","opt_slope",
                                                          "width_int","width_slope"))
flower_samples$Parameter <- factor(flower_samples$Parameter)
flower_params <- tidyr::spread(flower_samples,Parameter,value)
flower_means <- colMeans(flower_params)[3:8]

params <- cbind(visit_params[3:8],flower_params[,3:8])

####Workhorse Functions ####
inv.logit <- function(x){exp(x)/(1+exp(x))}

##Functional forms of the fit relationships for flowers and visitors
f1 <- function(x1,x2,height_int,height_slope,
                       opt_int,opt_slope,
                       width_int,width_slope) {
  exp(exp(width_int + width_slope * x2) * -1 * 
  (x1 -(opt_int + opt_slope * x2))^2 + (height_int + height_slope * x2))
  }
f2 <- function(x1,x2,height_int2,height_slope2,
                       opt_int2,opt_slope2,
                       width_int2,width_slope2) {
               inv.logit(exp(width_int2 + width_slope2 * x2) * -1 * 
               (x1 -(opt_int2 + opt_slope2 * x2))^2 + (height_int2 + height_slope2 * x2))
}

##Density functions
f1_dens <- function(x1,x2,height_int,height_slope,
                            opt_int,opt_slope,
                            width_int,width_slope) { 
  y <- f1(x1,x2,height_int,height_slope,
                opt_int,opt_slope,
                width_int,width_slope)
  yi <- integrate(f1, -Inf, +Inf, x2,height_int=height_int,height_slope=height_slope,
                                  opt_int=opt_int, opt_slope=opt_slope,
                                  width_int=width_int,width_slope=width_slope)
                                       
                                      
  return(y/yi[[1]])
}


f2_dens <- function(x1,x2,height_int2,height_slope2,
                            opt_int2,opt_slope2,
                            width_int2,width_slope2) { 
  y <- f2(x1,x2,height_int2,height_slope2,
                opt_int2,opt_slope2,
                width_int2,width_slope2)
  yi <- integrate(f2, -Inf, +Inf, x2,height_int2=height_int2,height_slope2=height_slope2,
                                       opt_int2=opt_int2, opt_slope2=opt_slope2,
                                       width_int2=width_int2,width_slope2=width_slope2)
  return(y/yi[[1]])
}

##Function to find the minimum of these two functions.
min_f1f2_dens <- function(x1,x2, height_int,height_slope,
                                 opt_int,opt_slope,
                                 width_int, width_slope,
                                 height_int2,height_slope2,
                                 opt_int2,opt_slope2,
                                 width_int2,width_slope2) {
  f1 <- f1_dens(x1,x2,height_int,height_slope,
                      opt_int,opt_slope,
                      width_int,width_slope) 
  f2 <- f2_dens(x1,x2,height_int2,height_slope2,
                      opt_int2,opt_slope2,
                      width_int2,width_slope2)
  pmin(f1, f2)
}

par(mfrow=c(1,1))
xtest <- seq(-5,8,by=0.1)
ytest1 <- f1(xtest,-2,visit_means[1],visit_means[2],
             visit_means[3],visit_means[4],
             visit_means[5],visit_means[6])
ytest2 <- f2(xtest,-2,flower_means[1],flower_means[2],
             flower_means[3],flower_means[4],
             flower_means[5],flower_means[6])
plot(xtest,ytest1,type="l",col="black",ylim=c(0,2))
points(xtest,ytest2,type="l",col="red")


ytest3 <- f1_dens(xtest,-2,visit_means[1],visit_means[2],
                  visit_means[3],visit_means[4],
                  visit_means[5],visit_means[6])
ytest4 <- f2_dens(xtest,-2,flower_means[1],flower_means[2],
                  flower_means[3],flower_means[4],
                  flower_means[5],flower_means[6])
points(xtest,ytest3,type="l",col="black",lty=2)
points(xtest,ytest4,type="l",col="red",lty=2)

ytest5 <- min_f1f2_dens(xtest,-2,visit_means[1],visit_means[2],
                                 visit_means[3],visit_means[4],
                                 visit_means[5],visit_means[6],
                                 flower_means[1],flower_means[2],
                                 flower_means[3],flower_means[4],
                                 flower_means[5],flower_means[6])
points(xtest,ytest5,type="l",col="purple")

##Function to calculate overlap coefficients.
measure_overlap <- function(params,x2val){
  integrate(min_f1f2_dens, -Inf, Inf, x2=x2val,height_int=params[1],height_slope=params[2],
                                      opt_int=params[3], opt_slope=params[4],
                                      width_int=params[5],width_slope=params[6],
                                      height_int2=params[7],height_slope2=params[8],
                                      opt_int2=params[9], opt_slope2=params[10],
                                      width_int2=params[11],width_slope2=params[12])$value
}

##Calculates overlap for each MCMC sample.
x2seq <- seq(-2,2,by=0.5)
overmat <- matrix(NA,nrow=nrow(params),ncol=length(x2seq))
for (i in 1:length(x2seq)){
  overmat[,i] <- apply(params,FUN=measure_overlap,MARGIN=1,x2val=x2seq[i])
}

