##Script to model plant phenology in Flickr data##
##Ian Breckheimer
##10 November 2015

####Sets up workspace####
library(dplyr)
library(ggplot2)
library(mgcv)
library(rjags)

####Brings in data####
f_all <- read.csv("./data/MORA_flickr_classified_2011_2014.csv")

###Removes photos with more than 20 photos per unique coordinate####
f_all$locfact <- as.factor(paste(f_all$long,f_all$lat,sep="by"))
f_all_grp <- group_by(f_all,locfact)
nphotos <- summarize(f_all_grp,nphotos=n())
small_locs <- nphotos$locfact[which(nphotos$nphotos < 20)]


###Creates a new dataset indicating whether there are any flowers at all.
f_all <- filter(f_all,Phenophase!="np" &
                      Species != "other" &
                      days_since_snow > -100 &
#                      sdd_pred > 150 &
#                     sdd_pred < 210 &
#                      year != 2012 &
#                      year != 2014 &
                      locfact %in% small_locs)

f_grp <- group_by(f_all,id,owner,year,datetaken_DOY,
                  nearest_center,days_since_snow, relev_30m,
                  canopy_pct,srad_noc,sdd_pred,X_UTM,Y_UTM)
f_flwr <- summarise(f_grp,FlwrYN = any(Phenophase=="flowering"))
f_flwr$Flwr <- as.numeric(f_flwr$FlwrYN)
f_flwr$yearfact <- as.factor(f_flwr$year)


####Plots data.####
ggplot(f_flwr,groups=yearfact)+
  geom_point(aes(x=days_since_snow,y=Flwr,color=nearest_center),
             position=position_jitter(height=0.1),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr,color=nearest_center),method="glm",formula=y~poly(x,2),se=TRUE,
              method.args=list(family="binomial"))+
  facet_wrap(facets=~nearest_center)+
  theme_bw()

####Generalized additive model####
gam2011 <-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=8)+nearest_center,
                               family=binomial("logit"), data = subset(f_flwr,year==2011))
summary(gam2011)
plot(gam2011,pages=1,residuals=TRUE)

gam_all<-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=7)+yearfact+
                 nearest_center*datetaken_DOY,
                family=binomial("logit"), data = f_flwr)
summary(gam_all)
plot(gam_all,pages=1,residuals=TRUE)


####Plots GAM Predictions####
gam_pred_data <- expand.grid(sdd_pred=seq(110,230,by=1),datetaken_DOY=seq(100,350,by=1),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"),yearfact=2014)
gam_pred_data$pFlwr <- predict(gam_all,newdata=gam_pred_data,type="response")
gam_pred_data$SE <- predict(gam_all,newdata=gam_pred_data,type="response",se.fit=TRUE)$se.fit

##Finds the maximum week for each prediction.

# Function to find the values that define 68.4% of the area under the curve.####
gam_width <- function(week,pred,threshold=0.1586553){
  
  # Converts to numeric vectors
  week <- as.numeric(week)
  pred <- as.numeric(pred)
  
  # Gets the bin width from the first pred interval.
  bin_width <- week[2] - week[1]
  
  # Make sure that all predictions are positive and none are missing.
  stopifnot(anyNA(pred)==FALSE,
            anyNA(week)==FALSE,
            any(pred<0)==FALSE)
  
  # Total area under the curve.
  total_area <- sum(pred*bin_width,na.rm=TRUE)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(pred)*bin_width/total_area
  cumprop_down <- rev(cumsum(rev(pred))*bin_width)/total_area
  
  # Finds the indices of the first and last values greater than 0.158
  lwr_index <- min(which(cumprop_up >= threshold))
  upr_index <- max(which(cumprop_down >= threshold))
  
  # Finds the corresponding values of dss.
  lwr_bound <- week[lwr_index]
  upr_bound <- week[upr_index]
  bounds <- c(lwr_bound,upr_bound)
  names(bounds) <- c("lwr_bound","upr_bound")
  
  # Output
  return(bounds)
}

##Computes optimum and width across the model fit.
gam_pred_grp <- group_by(gam_pred_data,sdd_pred,nearest_center)
gam_pred_max <- summarise(gam_pred_grp,DOYmax=datetaken_DOY[which.max(pFlwr)],
                          daystart=gam_width(datetaken_DOY,pFlwr,threshold=0.25)[1],
                          dayend=gam_width(datetaken_DOY,pFlwr,threshold=0.25)[2])

pdf("./figs/Flickr_gam_flwrprob_2011_2014_2d_fits.pdf",width=5,height=4)
ggplot(gam_pred_data)+
  geom_raster(aes(y=datetaken_DOY,x=sdd_pred,fill=pFlwr),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
#  geom_point(aes(y=datetaken_DOY,x=sdd_pred,shape=as.factor(Flwr)),
#             data=f_flwr,alpha=0.1)+
  scale_shape(solid = FALSE)+
  geom_line(aes(y=DOYmax,x=sdd_pred),size=0.5,linetype=3,data=gam_pred_max)+
  geom_line(aes(y=daystart,x=sdd_pred),size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(y=dayend,x=sdd_pred),size=0.5,linetype=2,data=gam_pred_max)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(facets=~nearest_center)+
  coord_cartesian()+
  xlab("Day of Snow Melt")+
  ylab("Day of Year")+
  theme_bw()
dev.off()

pdf("./figs/Flickr_gam_2d_se.pdf",width=5,height=4)
ggplot(gam_pred_data)+
  geom_raster(aes(y=week,x=melt_week,fill=SE),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  #  geom_point(aes(y=week,x=melt_week,size=nusers),
  #             data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
  #  scale_shape(solid = FALSE)+
  geom_line(aes(x=melt_week,y=weekmax),size=0.5,linetype=3,data=gam_pred_max)+
  geom_line(aes(x=melt_week,y=weekstart),size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(x=melt_week,y=weekend),size=0.5,linetype=2,data=gam_pred_max)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(facets=~nearest_center)+
  coord_cartesian()+
  xlab("Week of Snow Melt")+
  ylab("Week of Year")+
  theme_bw()
dev.off()

####Fits a parametric model using JAGS####

## Subsamples data for model development
rand <- runif(dim(f_flwr)[1],0,1)
f_samp <- f_flwr[rand <= 1,]

## Prep data.
obsy <- f_samp$Flwr
x1 <- as.numeric(scale(f_samp$datetaken_DOY/100,scale=FALSE))
obsx2 <- as.numeric(scale(f_samp$sdd_pred/100,scale=FALSE))
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
groups <- as.numeric(factor(f_samp$owner))
ngroups <- length(unique(groups))
sites <- as.numeric(factor(f_samp$nearest_center))
nsites <- length(unique(sites))

# simple plot of data..
obsy_t <- obsy[obsy==1]
x1_t <- x1[obsy==1]
obsx2_t <- obsx2[obsy==1]
obsy_n <- obsy[obsy==0]
x1_n <- x1[obsy==0]
obsx2_n <- obsx2[obsy==0]

par(mfrow=c(1,1))
plot(obsx2_n,x1_n,col="grey80",pch=20,cex=0.2)
points(obsx2_t,x1_t,col="red",pch=20,cex=0.3)

## specify model
cat("
    model {
    ## Priors
    height_int ~ dunif(-20,20)
    height_slope ~ dunif(-5,5)
    opt_int ~ dunif(-20,20)
    opt_slope ~ dunif(-5,5)
    width_int ~ dunif(-50,20)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)
#   site_opt_int_sd ~ dgamma(1,0.1)
#   site_opt_int_tau <- 1 / (site_opt_int_sd * site_opt_int_sd)
    
    ## Random effects priors
    for (j in 1:ngroups){
    height_grp[j] ~ dnorm(0,group_tau)
    }
    for (k in 1:nsites){
     opt_int_site[k] ~ dnorm(0,site_opt_int_tau)
    }

    
    ## Likelihood
    for (i in 1:n){
      x2true[i] ~ dnorm(0,pop_taux)
      x2[i] ~ dnorm(x2true[i],tau_obs)
      height[i] <- height_int + height_slope * x2true[i] + height_grp[group[i]] 
      opt[i] <- opt_int + opt_slope * x2true[i] + opt_int_site[site[i]]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror.txt")

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = groups, site = sites,
               ngroups = ngroups, n = length(x1))

# initiate model
mod2 <- jags.model("./scratch/xerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)

# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_int","opt_slope", "opt_int_site",
                                      "width_int","width_slope",
                                      "site_opt_int_sd","group_sd"))
gelman.diag(out2)
save(out2,file="./scratch/jags_output_2011_2014.Rdata",compress=TRUE)
