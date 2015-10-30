##Script to model the relationship between visitors and snow on Mt. Rainier.
##Ian Breckheimer
##27 October 2015

####Sets up workspace####
library(ggplot2)
library(mgcv)
library(dplyr)

####Loads and munges data data####
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

dss_fit_snow$melt_day <- dss_fit_snow$week * 7
dss_fit_snow$peak_day <- dss_fit_snow$melt_day + dss_fit_snow$max_dss

##Tries to fit all years simultaneously
weeks_snowmelt <- left_join(d_data_week[,1:5],d_snow[,1:3],by=c("year","nearest_center"))
colnames(weeks_snowmelt) <- c("week","year","nearest_center","nphotos","nusers","melt_week")

mod_week_all <- gam(nusers~s(week,melt_week,bs="tp",k=10)+
                           s(week,nearest_center,bs="fs",k=7),
                    family=nb(),optimizer="perf",data = weeks_snowmelt )
summary(mod_week_all)

gam_pred_data <- expand.grid(week=seq(1,53,by=0.5),melt_week=seq(10,35,by=0.5),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"))
gam_pred_data$Users <- predict(mod_week_all,newdata=gam_pred_data,type="response")
gam_pred_data$SE <- predict(mod_week_all,newdata=gam_pred_data,type="response",se.fit=TRUE)$se.fit

ggplot(gam_pred_data)+
  geom_raster(aes(y=week,x=melt_week,fill=Users),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
#  geom_point(aes(y=week,x=melt_week,size=nusers),
#             data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
#  scale_shape(solid = FALSE)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_grid(facets=.~nearest_center)+
  coord_cartesian()+
  xlab("Week of Snow Melt")+
  ylab("Week of Year")+
  theme_bw()

ggplot(gam_pred_data)+
  geom_raster(aes(y=week,x=melt_week,fill=SE),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  #  geom_point(aes(y=week,x=melt_week,size=nusers),
  #             data=subset(weeks_snowmelt,nusers>0),alpha=0.1)+
  #  scale_shape(solid = FALSE)+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  facet_grid(facets=.~nearest_center)+
  coord_cartesian()+
  xlab("Week of Snow Melt")+
  ylab("Week of Year")+
  theme_bw()
  

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


