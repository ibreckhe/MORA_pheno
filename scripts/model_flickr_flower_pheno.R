##Script to model plant phenology in Flickr data using generalized aditive models.

####Sets up workspace####
library(dplyr)
library(ggplot2)
library(mgcv)

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
#                      year != 2012 &
#                      year != 2014 &
                      locfact %in% small_locs)

f_grp <- group_by(f_all,id,year,datetaken_DOY,
                  nearest_center,days_since_snow,
                  sdd_pred)
f_flwr <- summarise(f_grp,FlwrYN = any(Phenophase=="flowering"))
f_flwr$Flwr <- as.numeric(f_flwr$FlwrYN)


####Plots data.####
ggplot(f_flwr)+
  geom_point(aes(x=days_since_snow,y=Flwr),
             position=position_jitter(height=0.2),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr),method="glm",formula=y~poly(x,2),se=TRUE,
              method.args=list(family="binomial"))+
  theme_bw()

####Generalized additive model####
gam2011 <-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=8)+nearest_center,
                               family=binomial("log"), data = subset(f_flwr,year==2011))
summary(gam2011)
plot(gam2011,pages=1,residuals=TRUE)

gam_all<-  gam(Flwr~s(sdd_pred,datetaken_DOY,bs="tp",k=8)+year+nearest_center,
                family=binomial("log"), data = f_flwr)
summary(gam_all)
plot(gam_all,pages=1,residuals=TRUE)


####Plots GAM Predictions####
gam_pred_data <- expand.grid(sdd_pred=seq(110,230,by=1),datetaken_DOY=seq(100,350,by=1),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"),year=2014)
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
