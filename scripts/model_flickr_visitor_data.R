##Script to model the relationship between visitors and snow on Mt. Rainier.
##Ian Breckheimer
##27 October 2015

####Sets up workspace####
library(ggplot2)
library(mgcv)
library(dplyr)

####Loads and munges data####
d_f <- read.csv("./data/MORA_flickr_metadata_all_2009_2015_cleaned.csv")

##Calculates the number of unique users per week for every year.##
d_f$datetaken <- as.POSIXct(d_f$datetaken)
d_f$week <- as.numeric(format(d_f$datetaken,format="%U"))
d_grp <- group_by(d_f,week,year,nearest_center)
d_sum <- summarize(d_grp,
                   nphotos=length(id),
                   nusers=length(unique(owner)),
                   mean_dss=mean(days_since_snow),
                   mean_doy=mean(datetaken_DOY),
                   mean_travel_t=mean(acc_times))
##Calculates the number of unique visitors for every day since snow melt.
d_grp2 <- group_by(d_f,days_since_snow,year,nearest_center)
d_sum_snow <- summarize(d_grp2,
                        nphotos=length(id),
                        nusers=length(unique(owner)),
                        mean_doy=mean(datetaken_DOY),
                        mean_sdd=mean(sdd_pred),
                        mean_travel_t=mean(acc_times))

##Does the same to get yearly totals.
d_grp3 <- group_by(d_f,nearest_center,year)
d_sum_year<- summarize(d_grp3,
                        nphotos=length(id),
                        nusers=length(unique(owner)),
                        mean_doy=mean(datetaken_DOY),
                        mean_travel_t=mean(acc_times))

##Adds back groups with zero records.
weeks <- expand.grid(week = 1:53, year = 2009:2015,nearest_center=levels(d_grp$nearest_center))
d_sum_complete <- left_join(weeks,d_sum, by=c("week","year","nearest_center"))
d_sum_complete <- d_sum_complete %>% mutate(nphotos = ifelse(is.na(nphotos),0,nphotos),
                                            nusers = ifelse(is.na(nusers),0,nusers),
                                            study = "Flickr")
days <- expand.grid(days_since_snow = -100:200, year = 2009:2015,nearest_center=levels(d_grp2$nearest_center))
d_sum_snow_complete <- left_join(days,d_sum_snow, by=c("days_since_snow","year","nearest_center"))
d_sum_snow_complete <- d_sum_snow_complete %>% mutate(nphotos = ifelse(is.na(nphotos),0,nphotos),
                                                 nusers = ifelse(is.na(nusers),0,nusers),
                                                 study = "Flickr")

##Melt dates for paradise SNOTEL 2009 - 2015
#para_dates <- c(202,205,240,209,199,204,151)
#para_elevs <- rep(1563,7)

#Adds predicted melt dates at parking lots from snow regression model.
para_preds <- c(177,191,210,191,180,181,139)
sunr_preds <- c(172,199,211,192,181,183,131)
tipsoo_preds <- c(175,173,207,187,171,171,107)
mow_preds <- c(177,187,203,185,178,176,149)
snow_preds <- c(para_preds,sunr_preds,tipsoo_preds,mow_preds)
snow_years <- rep(2009:2015,4)
snow_wks <- (snow_preds / 365) * 52
snow_locs <- factor(rep(c("Paradise","Sunrise","Tipsoo","Mowich"),each=7))
snow_elevs <- rep(c(1563,1950,1626,1506),4)

##Appends those values to data frame.
d_snow <- data.frame(d_sum_complete[1:28,])
d_snow[,] <- rep(NA,nrow(d_snow)*ncol(d_snow))
d_snow$week <- snow_wks
d_snow$year <- snow_years
d_snow$nearest_center <- snow_locs
d_snow$study <- rep("Snow",28)

d_sum$study <- "Flickr"
d_sum_snotel <- rbind(d_sum_complete,d_snow)

##Removes tipsoo lake in 2009, where there are no photos.
d_sum_snotel <- filter(d_sum_snotel,nearest_center!="Tipsoo" | year != 2009)
d_sum_snow_complete <-  filter(d_sum_snow_complete,nearest_center!="Tipsoo" | year != 2009)

####Model Fitting####

##Fits generalized additive models of unique Flickr users per week and location.
d_data_week <- group_by(filter(d_sum_snotel,study=="Flickr"),year,nearest_center)
d_week_fits <- d_data_week %>% do(mod=gam(nusers~s(week,k=5),family=nb(),optimizer="perf",data = . ))
d_week_grp <- group_by(d_week_fits,year,nearest_center)
d_week_preds <- d_week_grp %>% do(as.data.frame(predict(.$mod[[1]],newdata=data.frame(week=seq(0,53,by=0.01)),
                                                         type="response",se=TRUE)))
d_week_preds$week <- rep(seq(0,53,by=0.01),times=dim(d_week_fits)[2])

week_fit_max <- d_week_preds %>% 
                group_by(year,nearest_center) %>% 
                summarise(max_week=week[which.max(fit)])
week_fit_snow <- left_join(week_fit_max,d_snow,by=c("year","nearest_center"))[,1:4]
week_fit_snow <- left_join(week_fit_snow,d_sum_year,by=c("year","nearest_center"))


##Generalized additive models of unique Flickr users by DSS.
d_data_dss <- group_by(filter(d_sum_snow_complete,study=="Flickr"),year,nearest_center)
d_dss_fits <- d_data_dss %>% do(mod=gam(nusers~s(days_since_snow,k=5),family=nb(),optimizer="perf",data = . ))
d_dss_grp <- group_by(d_dss_fits,year,nearest_center)
d_dss_preds <- d_dss_grp %>% do(as.data.frame(predict(.$mod[[1]],newdata=data.frame(days_since_snow=seq(-100,200,by=0.1)),
                                                        type="response",se=TRUE)))
d_dss_preds$dss <- rep(seq(-100,200,by=0.1),times=dim(d_dss_fits)[2])

dss_fit_max <- d_dss_preds %>% 
  group_by(year,nearest_center) %>% 
  summarise(max_dss=dss[which.max(fit)])
dss_fit_snow <- left_join(dss_fit_max,d_snow,by=c("year","nearest_center"))[,1:4]
dss_fit_snow <- left_join(dss_fit_snow,d_sum_year,by=c("year","nearest_center"))

##Tries to fit all years simultaneously
weeks_snowmelt <- left_join(d_data_week[,1:5],d_snow[,1:3],by=c("year","nearest_center"))
colnames(weeks_snowmelt) <- c("week","year","nearest_center","nphotos","nusers","melt_week")

##Estimates day of year for each week and melt week.
weeks_snowmelt$DOY <- NA
weeks_snowmelt$DOY[weeks_snowmelt$year==2012] <- (weeks_snowmelt$week[weeks_snowmelt$year==2012] / (366/7)) * 366 - 3.5
weeks_snowmelt$DOY[weeks_snowmelt$year!=2012] <- (weeks_snowmelt$week[weeks_snowmelt$year!=2012] / (365/7)) * 365 - 3.5

weeks_snowmelt$melt_DOY <- NA
weeks_snowmelt$melt_DOY[weeks_snowmelt$year==2012] <- (weeks_snowmelt$melt_week[weeks_snowmelt$year==2012] / (366/7)) * 366 - 3.5
weeks_snowmelt$melt_DOY[weeks_snowmelt$year!=2012] <- (weeks_snowmelt$melt_week[weeks_snowmelt$year!=2012] / (365/7)) * 365 - 3.5


mod_week_all <- gam(nusers~s(DOY,melt_DOY,bs="tp",k=12)+
                           s(DOY,nearest_center,bs="fs",k=20),
                    family=nb(),optimizer="perf",data = weeks_snowmelt )
summary(mod_week_all)
plot(mod_week_all,pages=1,residuals=TRUE)

gam_pred_data <- expand.grid(DOY=seq(1,365,by=1),melt_DOY=seq(90,240,by=1),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"))
gam_pred_data$Users <- predict(mod_week_all,newdata=gam_pred_data,type="response")
gam_pred_data$SE <- predict(mod_week_all,newdata=gam_pred_data,type="response",se.fit=TRUE)$se.fit

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
gam_pred_grp <- group_by(gam_pred_data,melt_DOY,nearest_center)
gam_pred_max <- summarise(gam_pred_grp,weekmax=DOY[which.max(Users)],
                                       weekstart=gam_width(DOY,Users,threshold=0.25)[1],
                                       weekend=gam_width(DOY,Users,threshold=0.25)[2])

##Adds flowering window
flwr_start <- 23
flwr_end <- 49
flwr_max <- 36

gam_pred_max$flwrstart <- gam_pred_max$melt_DOY + flwr_start
gam_pred_max$flwrend <- gam_pred_max$melt_DOY + flwr_end
gam_pred_max$flwrmax <- gam_pred_max$melt_DOY + flwr_max


pdf("./figs/Flickr_gam_2d_fits.pdf",width=5,height=4)
ggplot(gam_pred_data)+
  geom_raster(aes(y=DOY,x=melt_DOY,fill=Users),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  #geom_point(aes(y=DOY,x=melt_DOY,size=nphotos),
  #           data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
  scale_shape(solid = FALSE)+
  geom_line(aes(x=melt_DOY,y=weekmax),size=0.5,linetype=3,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=weekstart),size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=weekend),size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=flwrmax),color="white",size=0.5,linetype=3,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=flwrstart),color="white",size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=flwrend),color="white",size=0.5,linetype=2,data=gam_pred_max)+
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
  geom_raster(aes(y=DOY,x=melt_DOY,fill=SE),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  #  geom_point(aes(y=DOY,x=melt_DOY,size=nusers),
  #             data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
  #  scale_shape(solid = FALSE)+
  geom_line(aes(x=melt_DOY,y=weekmax),size=0.5,linetype=3,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=weekstart),size=0.5,linetype=2,data=gam_pred_max)+
  geom_line(aes(x=melt_DOY,y=weekend),size=0.5,linetype=2,data=gam_pred_max)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_wrap(facets=~nearest_center)+
  coord_cartesian()+
  xlab("Week of Snow Melt")+
  ylab("Week of Year")+
  theme_bw()
dev.off()

##Determines sensitivity across all years.

##Week Model
week_lm <- lm(max_week~week,weights=nusers,data=week_fit_snow)
summary(week_lm)
confint(week_lm,parm=c("week"),level=0.95)

##DSS Model
sens_lm <- lm(peak_day~melt_day,weights=nusers,data=dss_fit_snow)
summary(sens_lm)
confint(sens_lm,parm=c("melt_day"),level=0.95)

####Makes plots of fit models.####

##Preps color palette.
library(RColorBrewer)
pal <- brewer.pal(7,name="Dark2")

##Function to add alpha
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pal[c(1,2,4,5,6)] <- add.alpha(pal[c(1,2,4,5,6)],alpha=0.4)


pdf("./figs/visitors_weekly.pdf",width=5.25,height=4)
ggplot(data=d_sum_snotel,aes(x=week,y=nusers,color=factor(year),
                             size=nphotos))+
  stat_smooth(method="gam",formula=y~s(x,k=7),method.args=list(family=nb(),optimizer="perf"),se=TRUE)+
  #  geom_vline(aes(xintercept=week,color=factor(year)),linetype=1,
  #             data=subset(d_sum_snotel,study=="Snow"))+
  scale_color_manual(values=pal)+
  geom_point(alpha=0)+
  ylab("Number of Unique Users")+
  xlab("Week of Year")+
  xlim(c(15,45))+
  #facet_grid(facets=year~nearest_center)+
  theme_bw()
dev.off()

pdf("./figs/visitors_dss.pdf",width=5.25,height=4)
ggplot(data=d_sum_snow_complete,aes(x=days_since_snow,y=nusers,color=factor(year),
                                    size=nphotos))+
  geom_point(alpha=0.05)+
  stat_smooth(method="gam",formula=y~s(x,k=5),method.args=list(family=nb(),optimizer="perf"),se=FALSE)+
  geom_vline(aes(xintercept=0,color=factor(year)),linetype=1)+
  scale_color_manual(values=pal)+
  ylab("Number of Unique Users")+
  xlab("Days Since Snow Melt")+
  ylim(c(0,25))+
  facet_grid(facets=year~nearest_center)+
  theme_bw()
dev.off()

####Plots snow melt against peak visitation.
pdf("./figs/visitor_sensitivity_weighted.pdf",width=5.25,height=4)
ggplot(data=week_fit_snow)+
  geom_point(aes(x=week,y=max_week,color=factor(year),
                 size=nusers))+
  geom_abline(slope=1,intercept=0)+
  #  stat_smooth(aes(x=week,y=max_week,color=nearest_center),method="lm",
  #              formula=y~x,se=FALSE)+
  stat_smooth(aes(x=week,y=max_week,weight=nusers),method="lm",
              formula=y~x,se=TRUE,fullrange=FALSE)+
  scale_color_manual(values=pal)+
  xlab("Week of Snow Melt")+
  ylab("Week of Peak Visitation")+
  ylim(c(15,45))+
  xlim(c(15,45))+
  theme_bw()
dev.off()

####Plots snow melt against peak visitation.
ggplot(data=dss_fit_snow) +
  geom_point(aes(x=melt_day,y=peak_day,color=factor(year),size=nusers),alpha=0.4)+
  geom_abline(slope=1,intercept=0)+
  #  stat_smooth(aes(x=melt_day,y=peak_day,color=nearest_center),method="lm",
  #              formula=y~x,se=FALSE)+
  stat_smooth(aes(x=melt_day,y=peak_day,weight=nusers),method="lm",
              formula=y~x,se=TRUE,fullrange=TRUE)+
  xlab("Snow Melt Day")+
  ylab("Peak Visit Day")+
  ylim(c(105,315))+
  xlim(c(105,315))+
  theme_bw()

####Models a parametric approach using JAGS ####
library(rjags)

##Preps data.

##Sample of data for model development.
sample_ind <- runif(length(weeks_snowmelt$nusers),0,1)
weeks_samp <- weeks_snowmelt[sample_ind <= 1,]

obsy <- weeks_samp$nusers
x1 <- as.numeric(scale(weeks_samp$DOY/100,scale=FALSE))
x1scale <- attributes(scale(weeks_samp$DOY/100,scale=FALSE))$`scaled:center`
obsx2 <- as.numeric(scale(weeks_samp$melt_DOY/100,scale=FALSE))
x2scale <- attributes(scale(weeks_samp$DOY/100,scale=FALSE))$`scaled:center`
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
year <- as.numeric(as.factor(weeks_samp$year))
site <- as.numeric(weeks_samp$nearest_center)
n <- length(obsy)
nyears <- length(unique(year))
nsites <- length(unique(site))

# write model
cat("
    model{
    ## Priors
    height_slope ~ dnorm(0,1)T(-5,5)
    opt_int ~ dunif(-20,20)
    opt_slope ~ dunif(-5,5)
    width_int ~ dunif(-20,20)
    width_slope ~ dnorm(0,1)T(-5,5)
    year_sd ~ dgamma(1,0.1)
    year_tau <- pow(year_sd,-2)
    site_height_sd ~ dgamma(1,0.1)
    site_height_tau <- pow(site_height_sd,-2)
    err_sd ~ dunif((sd_obs - 0.1),(sd_obs + 0.1))
    tau_obs <- 1 / (err_sd * err_sd)
    r ~ dgamma(1,0.1)
    
    ## Likelihood
    for (j in 1:nyr){
     height_yr[j] ~ dnorm(0,year_tau)
    }
    for (k in 1:nsites){
     height_site[k] ~ dnorm(0,site_height_tau)
    }
    for (i in 1:n){
      x2true[i] ~ dnorm(0,pop_taux)
      x2[i] ~ dnorm(x2true[i],tau_obs)
      height[i] <- height_slope * x2true[i] + height_yr[year[i]] + height_site[site[i]]
      opt[i] <- opt_int + opt_slope * x2true[i]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      log(mu[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      p[i] <- r/(r+mu[i])
      y[i] ~ dnegbin(p[i],r)
    }
    }
    ",
    fill=TRUE, file="./scratch/xyerror_nb.txt")

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,year = year,site=site,
               nyr=nyears,nsites=nsites, n = n)

# initiate model
mod2 <- jags.model("./scratch/xyerror_nb.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)

# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_slope","height_yr",
                                      "opt_int","opt_slope","opt_slope_site", 
                                      "width_int","width_slope", "height_site",
                                      "year_sd","site_height_sd"))
gelman.diag(out2)
save(out2,file="./scratch/visitors_jags_output_2011_2014.Rdata",compress=TRUE)


####Visualizes the JAGS fit ####
jags_pred_data <- expand.grid(DOY=seq(1,365,length.out=100),
                              melt_DOY=seq(120,220,length.out=100))

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
                                     x1vec=jags_pred_data$DOY/100-x1scale,
                                     x2vec=jags_pred_data$melt_DOY/100-x2scale)
ggplot(jags_pred_data)+
  geom_raster(aes(y=DOY,x=melt_DOY,fill=predmu),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  geom_point(aes(y=DOY,x=melt_DOY,size=nusers),
             data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
  theme_bw()