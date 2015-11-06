##Script to explore and visualize unclassified flickr data.
##Author:Ian Breckheimer
##Date: 10-9-2015

####Sets up workspace####

##Loads required packages.
library(ggmap)
library(leaflet)
library(rgdal)
library(maptools)
library(raster)
library(MBA)
library(dplyr)
library(mgcv)

setwd("/Users/ian/code/MORA_pheno")

###Loads data.
d <- read.csv("MORA_flickr_metadata_all_2009_2015_covar_sdd.csv")

d$yearfact <- as.factor(d$year)

d$datetaken <- strptime(d$datetaken,format="%Y-%m-%d")
d$datetaken <- as.POSIXct(d$datetaken) # convert to POSIXct to use in data frames / ddply
d$datetaken_DOY <- as.numeric(strftime(d$datetaken,format="%j"))
d$days_since_snow <- round((d$datetaken_DOY - d$sdd_pred),digits=0)

##Measures distance from major centers of activity.
paradise <- c(46.785919, -121.736449)
sunrise <- c(46.914453, -121.643671)
tipsoo <- c(46.869605, -121.519758)
mowich <- c(46.932694, -121.863214)
points_ll <- data.frame(rbind(paradise,sunrise,tipsoo,mowich))
colnames(points_ll) <- c("lat","long")
coordinates(points_ll) <- ~long+lat
points_ll$longitude <- coordinates(points_ll)[,1]
points_ll$latitude <- coordinates(points_ll)[,2] 
proj4string(points_ll) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
points_utm <- spTransform(points_ll,
                          CRSobj=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
points_utm$X_UTM <- coordinates(points_utm)[,1]
points_utm$Y_UTM <- coordinates(points_utm)[,2]

##Computes distances
dist_fun <- function(x1,y1,x2,y2){
  sqrt((x2-x1)^2 + (y2-y1)^2)
}

dist_para <- dist_fun(x1=d$X_UTM,y1=d$Y_UTM,
                        x2=points_utm$X_UTM[1],y2=points_utm$Y_UTM[1])
dist_sun <- dist_fun(x1=d$X_UTM,y1=d$Y_UTM,
                        x2=points_utm$X_UTM[2],y2=points_utm$Y_UTM[2])
dist_tips <- dist_fun(x1=d$X_UTM,y1=d$Y_UTM,
                        x2=points_utm$X_UTM[3],y2=points_utm$Y_UTM[3])
dist_mow <- dist_fun(x1=d$X_UTM,y1=d$Y_UTM,
                        x2=points_utm$X_UTM[4],y2=points_utm$Y_UTM[4])

##Finds minimum distance for each point.
dist <- data.frame(d_para=dist_para,d_sun=dist_sun,d_tips=dist_tips,d_mow=dist_mow)
min_dist <- apply(dist,FUN=min,MARGIN=1)
dist_which <- unlist(apply(dist,FUN=function(x){ ifelse(all(is.na(x)),NA,which.min(x))},MARGIN=1))
dist_fact <- factor(dist_which,labels=c("Paradise","Sunrise","Tipsoo","Mowich"))

##Adds to data
d$min_dist_center <- min_dist
d$nearest_center <- dist_fact

##Filters out photos outside of the scope of interest and those where uncertainty is high.
d_f <- filter(d,days_since_snow > -200 & 
                days_since_snow < 250 & 
                sdd_range < 20 &
                acc_times < 5 &
                elevation > 1200 &
                elevation < 2400)

d2009 <- d_f[d_f$year==2009,]
d2010 <- d_f[d_f$year==2010,]
d2011 <- d_f[d_f$year==2011,]
d2012 <- d_f[d_f$year==2012,]
d2013 <- d_f[d_f$year==2013,]
d2014 <- d_f[d_f$year==2014,]
d2015 <- d_f[d_f$year==2015,]

write.csv(d2009,"./data/flickr_subalp_2009.csv")
write.csv(d2010,"./data/flickr_subalp_2010.csv")
write.csv(d2011,"./data/flickr_subalp_2011.csv")
write.csv(d2012,"./data/flickr_subalp_2012.csv")
write.csv(d2013,"./data/flickr_subalp_2013.csv")
write.csv(d2014,"./data/flickr_subalp_2014.csv")
write.csv(d2015,"./data/flickr_subalp_2015.csv")

##Writes filtered data to disk.
write.csv(d_f,"./data/MORA_flickr_metadata_all_2009_2015_cleaned.csv")

####Compares density distributions of photos across the year.####
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

Year <- d_f$yearfact
svg(file="./figs/photo_density_2009_2015.svg",width=4.5,height=3)
ggplot(data=d_f)+
  geom_density(aes(x=datetaken_DOY,color=Year,linetype=Year),adjust=2.5)+
  scale_color_manual(values=pal)+
  xlab("Day of Year")+
  theme_bw()
dev.off()

svg(file="./figs/photo_density_2009_2015_dss.svg",width=4.5,height=3)
ggplot(data=d_f)+
  geom_density(aes(x=days_since_snow,color=Year,linetype=Year),adjust=2.5)+
  scale_color_manual(values=pal)+
  xlab("Days Since Snow Melt")+
  theme_bw()
dev.off()

####Calculates the number of unique users per week for every year.####
d_f$week <- as.numeric(format(d_f$datetaken,format="%U"))
d_grp <- group_by(d_f,week,year)
d_sum <- summarize(d_grp,
                   mean_travel_t=mean(acc_times),
                   nphotos=length(id),
                   nusers=length(unique(owner)))
d_grp2 <- group_by(d_f,days_since_snow,year)
d_sum_snow <- summarize(d_grp2,
                        nphotos=length(id),
                        nusers=length(unique(owner)))
d_grp3 <- group_by(d_f,grid1000,year)
d_sum_grid <- summarize(d_grp3,
                        nphotos=length(id),
                        nusers=length(unique(owner)))

####Measures linear trend in each grid.
d_grp4 <- group_by(d_sum_grid,grid1000)
grid_glm <- do(d_grp4,mod_lin=lm(nusers~year,data=.),
                      mod_null=lm(nusers~1,data=.))
grid_compare <- grid_glm %>% do(grid = unlist(as.character(.$grid1000)),
                                aov = unlist(anova(.$mod_null,.$mod_lin)$`Pr(>F)`[2]))
grid_coefs <- grid_glm %>% summarise(int=coef(mod_lin)[1],
                                     slope=coef(mod_lin)[2])
grid_stats <- data.frame(grid=unlist(grid_compare$grid),
                         pval=unlist(grid_compare$aov),
                         intercept=grid_coefs$int,
                         slope=grid_coefs$slope)

####Maps coefficients.
gridxy <- read.csv("./data/MORA_1000m_grid_2009cent.csv")
grid_stats_xy <- merge(grid_stats,gridxy,by.x="grid",by.y="ID",all.x=TRUE)
coordinates(grid_stats_xy) <- ~X+Y
proj4string(grid_stats_xy) <- "+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"
grid_stats_ll <- spTransform(grid_stats_xy,CRSobj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
grid_stats_ll@data$long <- coordinates(grid_stats_ll)[,1]
grid_stats_ll@data$lat <- coordinates(grid_stats_ll)[,2]


MBaccessToken <- "pk.eyJ1IjoiaWJyZWNraGUiLCJhIjoidVNHX1VpRSJ9.9fPQ1A3rdxyCAzPkeYSYEQ"
MBurlTemplate <- "https://a.tiles.mapbox.com/v4/ibreckhe.map-z05003mi/{z}/{x}/{y}.png?access_token="
MBTemplate <- paste(MBurlTemplate,MBaccessToken,sep="")

parkbound <- readShapePoly("~/Dropbox/Research/MORA_community_turnover/data/MORA_boundary_UTM10.shp")
proj4string(parkbound) <- "+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"
parkbound_latlng <- spTransform(parkbound,CRSobj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

##Transforms values to the range 0-1, maintaining zero as 0.5
slopes <- grid_stats_ll$slope
slopes[is.na(slopes) | is.nan(slopes) ] <- 0
trans_fun_low <- function(x, ...){
  scaled <- (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  return(scaled/2)
}
trans_fun_high <- function(x, ...){
  scaled <- (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  return((scaled/2)+0.5)
}
slopes_trans <- slopes
slopes_trans[slopes<0] <- trans_fun_low(slopes[slopes<0])
slopes_trans[slopes>=0] <- trans_fun_high(slopes[slopes>=0])
plot(density(slopes_trans))

##Makes color ramp.
cramp <- colorRamp(c("blue","white","red"))
ramp <- ramp(slopes_trans)
cols <- rgb(ramp,maxColorValue=256)

pthresh <- 0.1
grid_filter <- subset(grid_stats_ll,pval < pthresh)
cols_filter <- cols[which(grid_stats_ll$pval < pthresh)]
slopes_filter <- slopes[which(grid_stats_ll$pval < pthresh)]

m <- leaflet() %>% addTiles(MBTemplate) %>% setView(lat=46.854039, lng=-121.760366, zoom = 11)
m <- m %>% addPolygons(data=parkbound_latlng,
                       color="white",fillColor="white",
                       opacity=0.7,fillOpacity=0.1)


m <- m %>%  addCircleMarkers(lng=grid_stats_ll@data$long,
                          lat=grid_stats_ll@data$lat,
                          color=cols,
                          radius=slopes * 5,
                          fillOpacity=0.2,
                          weight=0.1)
m <- m %>%  addCircleMarkers(lng=grid_filter@data$long,
                             lat=grid_filter@data$lat,
                             color=cols_filter,
                             radius=slopes_filter * 5,
                             fillOpacity=0.8,
                             weight=0.5)
m <- m %>%  addCircleMarkers(lng=grid_filter@data$long,
                             lat=grid_filter@data$lat,
                             stroke=TRUE,
                             color="black",
                             fill=FALSE,
                             radius=slopes_filter * 5,
                             weight=4)
m

##Melt dates for paradise SNOTEL 2009 - 2015
para_dates <- c(202,205,240,209,199,204,151)
para_elevs <- rep(1563,7)

##Appends those values to data frame.
#Predicted melt date from regression model.
para_preds <- para_preds <- c(177,191,210,191,180,181,139)
para_wks <- (para_preds / 365) * 52
para_elevs <- rep(1563,length(para_preds))

##Appends those values to data frame.
d_para <- data.frame(d_sum[1:7,])
d_para[,] <- rep(NA,nrow(d_para)*ncol(d_para))
d_para$week <- para_wks
d_para$year <- 2009:2015
d_para$study <- rep("SNOTEL",7)

d_sum$study <- "Flickr"
d_sum_snotel <- rbind(d_sum,d_para)

ggplot(data=d_sum_snotel,aes(x=week,y=nusers,color=factor(year),
                      size=nphotos))+
  geom_point()+
  stat_smooth(method="gam",formula=y~s(x),se=FALSE)+
  geom_vline(aes(xintercept=week,color=factor(year)),linetype=1,
             data=subset(d_sum_snotel,study=="SNOTEL"))+
  scale_color_manual(values=pal)+
  ylab("Number of Unique Users")+
  xlab("Week of Year")+
  #facet_grid(facets=year~.)+
  theme_bw()

ggplot(data=d_sum_snow,aes(x=days_since_snow,y=nusers,color=factor(year),
                             size=nphotos))+
  geom_point(alpha=0.05)+
  stat_smooth(method="gam",formula=y~s(x),method.args=list(family="poisson"),se=TRUE)+
  geom_vline(aes(xintercept=0,color=factor(year)),linetype=1)+
  scale_color_manual(values=pal)+
  ylab("Number of Unique Users")+
  xlab("Days Since Snow Melt")+
  #facet_grid(facets=year~.)+
  theme_bw()

####Simple model of unique Flickr Users.
week_mod <- gam(nusers~s(week,factor(year),bs="fs"),family="poisson",data=d_sum_snotel)
plot(week_mod)

dss_mod <- gam(nusers~s(days_since_snow,factor(year),bs="fs"),family="poisson",data=d_sum_snow)
plot(dss_mod)

##Predicts peak visitation day based on the model.
year <- factor(2009:2015)
week <- seq(0,52,by=0.01)
week_year <- rep(week,times=length(year))
years <- rep(year,length.out=length(week_year))
pred_dat <- data.frame(week=week_year,year=years)

nuser_pred <- exp(predict(week_mod,newdata=pred_dat))
pred_dat$nusers <- nuser_pred
pred_dat_grp <- group_by(pred_dat,year)
pred_dat_max <- summarise(pred_dat_grp,max_wk=week[which.max(nusers)])

##Plots peak vs. paradise melt.
para_preds <- c(177,191,210,191,180,181,139)
para_wks <- (para_preds / 365) * 52
plot(pred_dat_max$max_wk~para_wks,xlim=c(18,40),ylim=c(18,40),
     xlab="Snow Melt (Week)",ylab="Peak Visitation (Week)")
visit_mod <- lm(pred_dat_max$max_wk~para_wks)
abline(0,1,lty=2)
abline(coef=coef(visit_mod),lty=1)

####

####Makes an interactive map####
MBaccessToken <- "pk.eyJ1IjoiaWJyZWNraGUiLCJhIjoidVNHX1VpRSJ9.9fPQ1A3rdxyCAzPkeYSYEQ"
MBurlTemplate <- "https://a.tiles.mapbox.com/v4/ibreckhe.map-z05003mi/{z}/{x}/{y}.png?access_token="
MBTemplate <- paste(MBurlTemplate,MBaccessToken,sep="")
m <- leaflet() %>% addTiles(MBTemplate) %>% setView(lat=46.854039, lng=-121.760366, zoom = 11)

parkbound <- readShapePoly("~/Dropbox/Research/MORA_community_turnover/data/MORA_boundary_UTM10.shp")
proj4string(parkbound) <- "+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"
parkbound_latlng <- spTransform(parkbound,CRSobj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

m <- m %>% addPolygons(data=parkbound_latlng,
                       color="white",fillColor="white",
                       opacity=0.7,fillOpacity=0.1)

m <- m %>% addCircleMarkers(lng=d2015$long,
                            lat=d2015$lat,
                            color=rgb(1,0.5,0.5,1),
                            radius=1,
                            fillOpacity=0.9,
                            weight=0.5)

m <- m %>% addCircleMarkers(lng=d2011$long,
                            lat=d2011$lat,
                            color=rgb(0.5,1,0.5,1),
                            radius=1,
                            fillOpacity=0.9,
                            weight=0.5)
m

paradise <- c(46.785919, -121.736449)
sunrise <- c(46.914453, -121.643671)
tipsoo <- c(46.869605, -121.519758)
mowich <- c(46.932694, -121.863214)
points_ll <- data.frame(rbind(paradise,sunrise,tipsoo,mowich))
colnames(points_ll) <- c("lat","long")
coordinates(points_ll) <- ~long+lat
points_ll$longitude <- coordinates(points_ll)[,1]
points_ll$latitude <- coordinates(points_ll)[,2] 
proj4string(points_ll) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
points_utm <- spTransform(points_ll,
                          CRSobj=CRS("+proj=utm +zone=10 +ellps=GRS80 +units=m +no_defs"))
points_utm$X_UTM <- coordinates(points_utm)[,1]
points_utm$Y_UTM <- coordinates(points_utm)[,2] 