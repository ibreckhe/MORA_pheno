##Script to predict snow disapearance dates of meadowatch sites.
##Author:Ian Breckheimer
##Date: 14 February 2014

##Loads required packages
library(raster)
library(ggmap)
library(dplyr)

##Brings in the 3m environmental data
setwd("/Users/ian/GIS/")
env3 <- brick("MORA_snow_pred.grd",
              crs=CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

names(env3) <- c("canopy_pct","elevation","slope","srad","stream_dist","asp",
                     "cold_air_ind","srad_noc","relev_30m","X_UTM","Y_UTM")
access <- raster("/Volumes/ib_working/GIS/access_time_hrs_wgs84.tif")

##Brings in the spatial grid.
setwd("/Volumes/ib_working/GIS/")
grid100 <- raster("MORA_grid_100m.tif")
grid500 <- raster("MORA_grid_500m.tif")
grid1000 <- raster("MORA_grid_1000m.tif")
grid5000 <- raster("MORA_grid_500m.tif")

##Brings in the flickr metadata.
msites <- read.csv("~/code/MORA_pheno/MORA_flickr_metadata_all_2009_2015.csv")

##Converts locations to a spatial points data frame.
coordinates(msites) <- ~UTME+UTMN
proj4string(msites) <- "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

##Samples the covariate data at the points.
msites_covar <- extract(env3,msites,method="simple",sp=TRUE)

##Samples the access time raster as well.
msites_covar@data$acc_times <- extract(access,cbind(msites_covar@data$long,msites_covar@data$lat),
                                       method='simple')

msites_covar@data$grid100 <- extract(grid100,coordinates(msites_covar),
                                       method='simple')
msites_covar@data$grid500 <- extract(grid500,coordinates(msites_covar),
                                     method='simple')
msites_covar@data$grid1000 <- extract(grid1000,coordinates(msites_covar),
                                     method='simple')
msites_covar@data$grid5000 <- extract(grid5000,coordinates(msites_covar),
                                     method='simple')

##Writes to csv file
setwd("~/code/MORA_pheno")
write.csv(msites_covar@data,"MORA_flickr_metadata_all_2009_2015_covar.csv",row.names=FALSE)
