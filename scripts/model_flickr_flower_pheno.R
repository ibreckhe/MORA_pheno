##Script to model plant phenology in Flickr data##
##Ian Breckheimer
##10 November 2015

####Sets up workspace####
library(Hmisc)
library(dplyr)
library(ggplot2)
library(mgcv)
library(rjags)

####Brings in data####
f_all <- read.csv("./data/MORA_flickr_classified_2011_2015.csv")

###Removes photos with more than 20 photos per unique coordinate####
f_all$locfact <- as.factor(paste(f_all$long,f_all$lat,sep="by"))
f_all_grp <- group_by(f_all,locfact)
nphotos <- summarize(f_all_grp,nphotos=n())
small_locs <- nphotos$locfact[which(nphotos$nphotos < 20)]


###Creates a new dataset indicating whether there are any flowers at all.
f_all <- filter(f_all,Phenophase!="np"&
                      Species != "other" &
                      days_since_snow > -100 &
                       datetaken_DOY > 90 &
                       datetaken_DOY < 330 &
#                      year != 2012 &
#                      year != 2014 &
                      locfact %in% small_locs)

###Relative abundance of Erythronium and Anemone in each year.

f_grp <- group_by(f_all,id,owner,year,datetaken_DOY,
                  nearest_center,days_since_snow, relev_30m,
                  canopy_pct,srad_noc,sdd_pred,X_UTM,Y_UTM)
f_flwr <- summarise(f_grp,FlwrYN = any(Phenophase=="flowering"),
                          EarlyYN = any(Species %in% c("erythronium montanum","anemone occidentalis") & 
                                          Phenophase == "flowering"),
                          LateYN = any(Species %in% c("ligusticum grayi","polygonum bistortoides") & 
                                    Phenophase == "flowering"),
                          NoEarly = any(Species %nin% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering"),
                          NoLate = any(Species %nin% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering"),
                    
                          nSpp = sum(as.numeric(Phenophase=="flowering")))
f_flwr$Flwr <- as.numeric(f_flwr$FlwrYN)
f_flwr$Early <- as.numeric(f_flwr$EarlyYN)
f_flwr$Late <- as.numeric(f_flwr$LateYN)
f_flwr$notEarly <- as.numeric(f_flwr$NoEarly)
f_flwr$notLate <- as.numeric(f_flwr$NoLate)
f_flwr$yearfact <- as.factor(f_flwr$year)

##Computes the proportion of photos with early and late flowers in every year.
f_yr <- group_by(f_flwr,year)
f_prop <- summarise(f_yr,pFlwr=mean(Flwr),
                         nFlwr=sum(Flwr),
                         nEarly=sum(Early),
                         nLate=sum(Late),
                         nNotEarly=sum(notEarly),
                         nNotLate=sum(notLate),
                         pEarly=nEarly/nFlwr,
                         pLate=nLate/nFlwr)
f_prop_2012_2014 <- f_prop[f_prop$year %in% 2012:2014,c("nEarly","nLate","nNotEarly","nNotLate")]
sums_2012_2014 <- colSums(f_prop_2012_2014)

f_prop_2015 <- f_prop[f_prop$year == 2015,c("nEarly","nLate","nNotEarly","nNotLate")]
sums_2015 <- colSums(f_prop_2015)

##Chi-squared test of the proportion of photos with early and late flowering species.
early_ctable <- matrix(NA,ncol=2,nrow=2)
early_ctable[1,1] <- sums_2012_2014["nEarly"]
early_ctable[1,2] <- sums_2015["nEarly"]
early_ctable[2,1] <- sums_2012_2014["nNotEarly"]
early_ctable[2,2] <- sums_2015["nNotEarly"]
colnames(early_ctable) <- c("Typical Melt","Early Melt")
rownames(early_ctable) <- c("Early Flowering","Not Early Flowering")

chisq.test(early_ctable,simulate.p.value=TRUE,B=10000)

late_ctable <- matrix(NA,ncol=2,nrow=2)
late_ctable[1,1] <- sums_2012_2014["nLate"]
late_ctable[1,2] <- sums_2015["nLate"]
late_ctable[2,1] <- sums_2012_2014["nNotLate"]
late_ctable[2,2] <- sums_2015["nNotLate"]
colnames(late_ctable) <- c("Typical Melt","Early Melt")
rownames(late_ctable) <- c("Late Flowering","Not Late Flowering")

chisq.test(late_ctable,simulate.p.value=TRUE,B=10000)

##Computes the proportion of unique users that have flower photos in every year.
f_owner <- group_by(f_all,owner,year)
f_owner_stats <- summarise(f_owner,FlwrYN = any(Phenophase=="flowering"),
                    EarlyYN = any(Species %in% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering"),
                    LateYN = any(Species %in% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering"),
                    NoFlwr = all(Phenophase != "flowering" & Species != "other"),
                    NoEarly = any(Species %nin% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering"),
                    NoLate = any(Species %nin% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering"),
                    nSpp = sum(as.numeric(Phenophase=="flowering")))
f_owner_stats$Flwr <- as.numeric(f_owner_stats$FlwrYN)
f_owner_stats$noFlwr <- as.numeric(f_owner_stats$NoFlwr)
f_owner_stats$Early <- as.numeric(f_owner_stats$EarlyYN)
f_owner_stats$Late <- as.numeric(f_owner_stats$LateYN)
f_owner_stats$notEarly <- as.numeric(f_owner_stats$NoEarly)
f_owner_stats$notLate <- as.numeric(f_owner_stats$NoLate)
f_owner_stats$yearfact <- as.factor(f_owner_stats$year)

##Computes the proportion of photos with early and late flowers in every year.
f_owner_yr <- group_by(f_owner_stats,year)
f_owner_prop <- summarise(f_owner_yr,pFlwr=mean(Flwr),
                    nFlwr=sum(Flwr),
                    nNoFlwr=sum(noFlwr),
                    nEarly=sum(Early),
                    nLate=sum(Late),
                    nNotEarly=sum(notEarly),
                    nNotLate=sum(notLate),
                    pEarly=nEarly/nFlwr,
                    pLate=nLate/nFlwr)
f_owner_prop_2012_2014 <- f_owner_prop[f_owner_prop$year %in% 2012:2014,c("nEarly","nLate","nNotEarly","nNotLate")]
sums_owner_2012_2014 <- colSums(f_owner_prop_2012_2014)

f_owner_prop_2015 <- f_owner_prop[f_owner_prop$year == 2015,c("nEarly","nLate","nNotEarly","nNotLate")]
sums_owner_2015 <- colSums(f_owner_prop_2015)

##Chi-squared test of the proportion of photos with early and late flowering species.
early_owner_ctable <- matrix(NA,ncol=2,nrow=2)
early_owner_ctable[1,1] <- sums_owner_2012_2014["nEarly"]
early_owner_ctable[1,2] <- sums_owner_2015["nEarly"]
early_owner_ctable[2,1] <- sums_owner_2012_2014["nNotEarly"]
early_owner_ctable[2,2] <- sums_owner_2015["nNotEarly"]
colnames(early_owner_ctable) <- c("Typical Melt","Early Melt")
rownames(early_owner_ctable) <- c("Early Flowering","Not Early Flowering")

chisq.test(early_owner_ctable,simulate.p.value=TRUE,B=10000)

late_owner_ctable <- matrix(NA,ncol=2,nrow=2)
late_owner_ctable[1,1] <- sums_owner_2012_2014["nLate"]
late_owner_ctable[1,2] <- sums_owner_2015["nLate"]
late_owner_ctable[2,1] <- sums_owner_2012_2014["nNotLate"]
late_owner_ctable[2,2] <- sums_owner_2015["nNotLate"]
colnames(late_owner_ctable) <- c("Typical Melt","Early Melt")
rownames(late_owner_ctable) <- c("Late Flowering","Not Late Flowering")

chisq.test(late_owner_ctable,simulate.p.value=TRUE,B=10000)

####Plots data.####

pdf("./figs/alldat_classified_2011_2015_yearly.pdf",width=6,height=10)
ggplot(f_flwr)+
  geom_point(aes(x=sdd_pred,y=datetaken_DOY,color=as.factor(FlwrYN)),size=0.1)+
  geom_abline(slope=1,intercept=0,linetype=1)+
  scale_color_grey(start=0.8,end=0.2,guide=FALSE)+
  geom_smooth(aes(x=sdd_pred,y=datetaken_DOY,color=as.factor(FlwrYN)),
              method="lm",formula=y~x,se=TRUE)+
  xlab("Predicted Snow Disappearance Day")+
  ylab("Day of Year")+
  facet_grid(facets=year~.)+
  theme_bw()
dev.off()

ggplot(f_flwr,groups=yearfact)+
  geom_point(aes(x=days_since_snow,y=Flwr,color=nearest_center),
             position=position_jitter(height=0.1),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr,color=nearest_center),method="glm",formula=y~poly(x,2),se=TRUE,
              method.args=list(family="binomial"))+
#  facet_grid(facets=year~.)+
  theme_bw()

ggplot(f_flwr,groups=yearfact)+
  geom_point(aes(x=days_since_snow,y=nSpp),
             position=position_jitter(height=0.1),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr),method="gam",formula=y~s(x),se=TRUE,
              method.args=list(family="nb",optimizer="perf"))+
#  facet_grid(facets=nearest_center~.)+
  ylim(c(0,7))+
  theme_bw()


####Generalized additive model####
gam2011 <-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=8)+nearest_center,
                               family=binomial("logit"), data = subset(f_flwr,year==2011))
summary(gam2011)
plot(gam2011,pages=1,residuals=TRUE)

gam_all<-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=10)+yearfact+
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

pdf("./figs/Flickr_gam_flwrprob_2011_2015_2d_fits.pdf",width=5,height=4)
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
x1 <- f_samp$datetaken_DOY/100 - 2
obsx2 <- f_samp$sdd_pred/100 - 2
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
groups <- as.numeric(factor(f_samp$owner))
ngroups <- length(unique(groups))
year <- as.numeric(factor(f_samp$year))
sites <- as.numeric(factor(f_samp$nearest_center))
nsites <- length(unique(sites))

ggdat <- data.frame(FlwrYN=obsy,DOY=x1,SDD=obsx2,year=factor(f_samp$year),sites=factor(f_samp$nearest_center))
ggplot(ggdat)+
  geom_point(aes(x=SDD,y=DOY,color=as.factor(FlwrYN)),size=0.1)+
  geom_abline(slope=1,intercept=0,linetype=1)+
  scale_color_grey(start=0.8,end=0.2)+
  geom_smooth(aes(x=SDD,y=DOY,color=as.factor(FlwrYN)),
                  method="lm",formula=y~x,se=TRUE)+
#  facet_grid(facets=.~sites)+
  theme_bw()

# simpler plot of data..
obsy_t <- obsy[obsy==1]
x1_t <- x1[obsy==1]
obsx2_t <- obsx2[obsy==1]
obsy_n <- obsy[obsy==0]
x1_n <- x1[obsy==0]
obsx2_n <- obsx2[obsy==0]

par(mfrow=c(1,1))
plot(obsx2_n,x1_n,col="black",pch=20,cex=0.2)
points(obsx2_t,x1_t,col="red",pch=20,cex=0.3)
abline(0,1,lty=2)

## specify model with common parameters across all sites.
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
   
    ## Random effects priors
    for (j in 1:ngroups){
    height_grp[j] ~ dnorm(0,group_tau)
    }
    
    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_int + height_slope * x2true[i] + height_grp[group[i]] 
    opt[i] <- opt_int + opt_slope * x2true[i]
    width[i] <- exp(width_int + width_slope * x2true[i]) * -1
    logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_common.txt")


## specify model with a separate opt intercept for each site.
cat("
    model {
    ## Priors
    height_int ~ dunif(-20,20)
    height_slope ~ dunif(-5,5)
    opt_slope ~ dunif(-5,5)
    width_int ~ dunif(-50,20)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)
    site_opt_int_sd ~ dgamma(1,0.1)
    site_opt_int_tau <- 1 / (site_opt_int_sd * site_opt_int_sd)
    site_opt_slope_sd ~ dgamma(1,0.1)
    site_opt_slope_tau <- 1 / (site_opt_slope_sd * site_opt_slope_sd)
    
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
      opt[i] <- opt_slope * x2true[i] + opt_int_site[site[i]]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_siteopt.txt")

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = groups, site = sites,
               ngroups = ngroups,nsites=nsites, n = length(x1))

# initiate model
mod2 <- jags.model("./scratch/xerror_common.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)

# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int",
                                      "width_int","width_slope",
                                      "group_sd"))
gelman.diag(out2)
save(out2,file="./scratch/jags_flower_output_common_2011_2015.Rdata",compress=TRUE)
load("./scratch/jags_flower_output_common_2011_2015.Rdata")


# initiate model
mod3 <- jags.model("./scratch/xerror_siteopt.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)

# simulate posterior
out3 <- coda.samples(mod3, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int_site",
                                      "width_int","width_slope",
                                      "site_opt_int_sd","group_sd"))
gelman.diag(out3)
save(out3,file="./scratch/jags_flower_output_siteopt_2011_2015.Rdata",compress=TRUE)

####Plots JAGS Predictions####
jags_pred_data <- expand.grid(sdd_pred=seq(110,230,by=1),datetaken_DOY=seq(100,350,by=1),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"),yearfact=2014)

jags_fit_fun <- function(height_int,height_slope,
                         opt_int,opt_slope,
                         width_int,width_slope,
                         x1vec,x2vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  
  #Calculates the value of the response.
  height <- height_int + height_slope * x2vec
  opt <- opt_int + opt_slope * x2vec
  width<- exp(width_int + width_slope * x2vec) * -1
  mu <- exp(width * (x1vec - opt)^2 + height)
  return(mu)
}

jags_meds <- summary(out2)$quantiles[,3]
jags_pred_data$predmu <- jags_fit_fun(height_int=jags_meds["height_int"],
                                      height_slope=jags_meds["height_slope"],
                                      width_int=jags_meds["width_int"],
                                      width_slope=jags_meds["width_slope"],
                                      opt_int=jags_meds["opt_int"],
                                      opt_slope=jags_meds["opt_slope"],
                                      x1vec=jags_pred_data$datetaken_DOY/100 - 2,
                                      x2vec=jags_pred_data$sdd_pred/100 - 2)

pdf("./figs/jags_fit_flower_common.pdf",width=6,height=4)
ggplot(jags_pred_data)+
  geom_raster(aes(y=datetaken_DOY,x=sdd_pred,fill=predmu),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  geom_point(aes(y=datetaken_DOY,x=sdd_pred,
                 shape=factor(Flwr),color=factor(Flwr)),
             data=f_samp,alpha=0.5)+
  scale_color_grey(start=0,end=1)+
  xlim(c(110,230))+
  ylim(c(100,350))+
  xlab("Day of Snow Melt")+
  ylab("Day of Year")+
  guides(fill=guide_colorbar("Predicted Prob."),
         size=guide_legend("Flowering"))+
  theme_bw()
dev.off()
