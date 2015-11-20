##Script to estimate overlap coefficients for visitor and flower phenology
##Ian Breckheimer
##November 13th, 2015

####Sets up workspace####
library(dplyr)
library(rjags)
library(ggmcmc)
library(doParallel)
library(foreach)
source("./scripts/overlap_functions.R")

##Brings in fit model objects and converts them to data frames.
load("./scratch/visitors_jags_output_2011_2014.Rdata")
visit_model <- out2
load("./scratch/jags_output_2011_2014.Rdata")
flower_model <- out2

##Processes model output to measure phenological mismatch for each site.

##Names and days for measurement
site_names <- c("Mowich","Paradise","Sunrise","Chinook Pass")
measure_sdd <- seq(100,260,by=20)

##Function that does the heavy-lifting.
jags_mismatch <- function(visit_model,flower_model,group_index,measure_sdd,
                          site_names){
  params <- prep_mcmc_vars(visit_model,flower_model,group_index=group_index,
                           mod1_param_names = c("height_int_site","height_slope",
                                                "opt_int","opt_slope",
                                                "width_int","width_slope"),
                           mod2_param_names = c("height_int","height_slope",
                                                "opt_int_site","opt_slope",
                                                "width_int","width_slope"))
  mismatch <- measure_mismatch(params=params,x2sdd=measure_sdd)
  mismatch$Site <- site_names[i]
  return(mismatch)
}

##Registers a parallel backend.
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

##Runs the computation in parallel.
mismatch <- foreach(i=1:length(site_names),.combine='rbind') %dopar% {
                  jags_mismatch(visit_model,flower_model,group_index=i,measure_sdd=measure_sdd,
                                site_names=site_names)
}
stopCluster(cl)

##Graphs the results.
month_breaks <- c(121,152,182,213,244)
month_labels <- c("May","June","July","Aug","Sept")

#pdf("./figs/overlap_sdd_paradise_2015.pdf",width=4,height=4)
ggplot(data=mismatch)+
  geom_ribbon(aes(x=SDD,ymax=overlap_upr,ymin=overlap_lwr),fill="grey20",alpha=0.2)+
#  geom_ribbon(aes(x=SDD,ymax=overlap_q10,ymin=overlap_q90),fill="grey20",alpha=0.4)+
  geom_ribbon(aes(x=SDD,ymax=overlap_q25,ymin=overlap_q75),fill="grey20",alpha=0.6)+
  geom_line(aes(x=SDD,y=overlap_q50))+
  scale_x_continuous(breaks = month_breaks,labels = month_labels )+
  xlab("Last Snow Melt")+
  ylab("Phenological Match")+
  facet_wrap(facets=~Site)+
  ylim(c(0,1))+
  theme_bw()
#dev.off()

##Graphs functions and overlap at 2011 and 2015 Snow Disappearance Dates
paramelt2015 <- 139
paramelt2011 <- 211
xseq <- seq(1,365,by=1) / 100 - 2
x2_2015 <- paramelt2015 / 100 - 2
x2_2011 <- paramelt2011 / 100 - 2

params <- prep_mcmc_vars(visit_model,flower_model,group_index=2,
                         mod1_param_names = c("height_int_site","height_slope",
                                              "opt_int","opt_slope",
                                              "width_int","width_slope"),
                         mod2_param_names = c("height_int","height_slope",
                                              "opt_int_site","opt_slope",
                                              "width_int","width_slope"))
visit_params <- params[,1:6]
flower_params <- params[,7:12]

visit_fun_2011 <- function(x) { f1_dens(x1=xseq, x2=x2_2011, 
                              height_int=x[1],height_slope=x[2],
                              opt_int=x[3],opt_slope=x[4],
                              width_int=x[5],width_slope=x[6])}
visits_2011 <- apply(visit_params, FUN=visit_fun_2011,MARGIN = 1)

visit_fun_2015 <- function(x) { f1_dens(x1=xseq, x2=x2_2015, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
visits_2015 <- apply(visit_params, FUN=visit_fun_2015,MARGIN = 1)

flower_fun_2011 <- function(x) { f2_dens(x1=xseq, x2=x2_2011, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
flowers_2011 <- apply(flower_params, FUN=flower_fun_2011,MARGIN = 1)

flower_fun_2015 <- function(x) { f2_dens(x1=xseq, x2=x2_2015, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
flowers_2015 <- apply(flower_params, FUN=flower_fun_2015,MARGIN = 1)

visits_quant_2011 <- data.frame(t(apply(visits_2011,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
visits_quant_2011$group <- "Visitors"
visits_quant_2011$Year <- 2011


visits_quant_2015 <- data.frame(t(apply(visits_2015,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
visits_quant_2015$group <- "Visitors"
visits_quant_2015$Year <- 2015

flowers_quant_2011 <- data.frame(t(apply(flowers_2011,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
flowers_quant_2011$group <- "Flowers"
flowers_quant_2011$Year <- 2011

flowers_quant_2015 <- data.frame(t(apply(flowers_2015,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
flowers_quant_2015$group <- "Flowers"
flowers_quant_2015$Year <- 2015

flowers_visitors <- rbind(visits_quant_2015,visits_quant_2011,
                          flowers_quant_2015,flowers_quant_2011)
flowers_visitors$DOY <- rep((xseq + 2) * 100,4)
colnames(flowers_visitors) <- c("lwr","lwr25","median","upr75","upr","group","year","doy")

pdf("./figs/curves_2011_2015_paradise.pdf",width=8,height=4)
ggplot(data=flowers_visitors)+
  geom_line(aes(x=doy,y=median,color=group))+
  geom_ribbon(aes(x=doy,ymin=lwr,ymax=upr,fill=group),alpha=0.2)+
  facet_wrap(facets=~year)+
  xlab("Day of Year")+
  ylab("Density")+
  xlim(c(0,366))+
  theme_bw()
dev.off()
