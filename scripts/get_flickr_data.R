####Script to format REST JSON-formatted data from flickr into a .csv that can be used in Microsoft Excel.
####Author: Ian Breckheimer.
####Created: 27 January 2013
####Updated: 10-7-2015

###Flickr API call: https://api.flickr.com/services/rest/?method=flickr.photos.search&api_key=ef98b511cc80c5ad0e0e9b54d9dc0050&min_taken_date=2014-01-01&max_taken_date=2015-01-01&bbox=-121.993775%2C46.695609%2C-121.517243%2C+47.014907&accuracy=13&extras=geo%2C+date_upload%2C+date_taken%2C+url_o%2C+url_z%2C+tags&format=json&nojsoncallback=1&auth_token=72157651259053427-8b7ae1abd4ab78b4&api_sig=161c57669f96a10d085565396ba4a329
##loads required packages.
library(jsonlite)
library(RCurl)
library(maptools)
library(rgdal)
library(dplyr)
library(ggmap)

##Sets the working directory
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/Flickr/api_processing")

##Function to scrape data from the Flickr API.
get_flickr_data <- function(api_key,min_taken_date,max_taken_date,accuracy){
          ##Builds the API call
          api_list <- list(start="https://api.flickr.com/services/rest/?method=flickr.photos.search&api_key=",
                           key=api_key,
                           taken_sep="&min_taken_date=",
                           min_taken_date=min_taken_date,
                           taken_max_sep="&max_taken_date=",
                           max_taken_date=max_taken_date,
                           bbox="&bbox=-121.993775%2C46.695609%2C-121.423867%2C+47.014907",
                           acc_sep="&accuracy=",
                           accuracy=accuracy,
                           extras="&extras=geo%2C+date_upload%2C+date_taken%2C+url_o%2C+url_z%2C+tags",
                           per_page="&per_page=250&page=",
                           page=1,
                           format="&format=json&nojsoncallback=1")
          
          api_string <- paste(api_list,collapse="")
          
          ##calls the API:
          result <- getURLContent(api_string)
          
          ##loads the data, converting it from JSON format.
          result_conv <- fromJSON(result)
          result_pages <- result_conv$photos$pages
          result_photos <- as.numeric(result_conv$photos$total)
          
          ##Loops through the rest of the pages
          print(paste("Found ",result_photos," photos on ",result_pages," pages, downloading metadata...",sep=""))
          data <- list()
          for (i in 1:result_pages){
            print(paste("Page ",i,"..."))
            api_list2 <- api_list
            api_list2$page <- i
            api_string2 <- paste(api_list2,collapse="")
            result2 <- fromJSON(getURLContent(api_string2))
            data2 <- result2$photos$photo
            data <- rbind(data,data2)
          }
          data <- unique(data)
          return(data)
}

##Default arguments
api_key <- "5feeec895890f7c0a6d604e25207a4de"
api_secret <- "00c14c0da44e5352"
accuracy <- "16"

##Gets data for 2009.
min_taken_date <- "2009-01-01"
max_taken_date <- "2009-8-01"
data_2009_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2009_1$year <- 2009

min_taken_date <- "2009-8-1"
max_taken_date <- "2009-12-31"
data_2009_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2009_2$year <- 2009

##Gets data for 2010.
min_taken_date <- "2010-01-01"
max_taken_date <- "2010-8-01"
data_2010_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2010_1$year <- 2010

min_taken_date <- "2010-8-2"
max_taken_date <- "2010-12-31"
data_2010_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2010_2$year <- 2010

##Gets data for 2011.
min_taken_date <- "2011-01-01"
max_taken_date <- "2011-8-01"
data_2011_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2011_1$year <- 2011

min_taken_date <- "2011-08-02"
max_taken_date <- "2011-12-31"
data_2011_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2011_2$year <- 2011

##Gets data for 2012.
min_taken_date <- "2012-01-01"
max_taken_date <- "2012-08-01"
data_2012_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2012_1$year <- 2012

min_taken_date <- "2012-08-02"
max_taken_date <- "2012-12-31"
data_2012_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2012_2$year <- 2012

##Gets data for 2013.
min_taken_date <- "2013-01-01"
max_taken_date <- "2013-08-01"
data_2013_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2013_1$year <- 2013

min_taken_date <- "2013-08-02"
max_taken_date <- "2013-12-31"
data_2013_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2013_2$year <- 2013

##Gets data for 2014.
min_taken_date <- "2014-01-01"
max_taken_date <- "2014-08-01"
data_2014_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2014_1$year <- 2014

min_taken_date <- "2014-08-02"
max_taken_date <- "2014-12-31"
data_2014_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2014_2$year <- 2014

min_taken_date <- "2015-01-01"
max_taken_date <- "2015-08-01"
data_2015_1 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2015_1$year <- 2015
min_taken_date <- "2015-08-02"
max_taken_date <- "2015-12-31"
data_2015_2 <- get_flickr_data(api_key,min_taken_date,max_taken_date,accuracy)
data_2015_2$year <- 2015

data_2011 <- rbind(data_2011_1,data_2011_2)
data_2015 <- rbind(data_2015_1,data_2015_2)

data_all <- rbind(data_2009_1,
                  data_2009_2,
                  data_2010_1,
                  data_2010_2,
                  data_2011_1,
                  data_2011_2,
                  data_2012_1,
                  data_2012_2,
                  data_2013_1,
                  data_2013_2,
                  data_2014_1,
                  data_2014_2,
                  data_2015_1,
                  data_2015_2)

##Gets rid of duplicate records.
data_all <- unique(data_all)
data_2011 <- unique(data_2011)
data_2015 <- unique(data_2015)

##creates a stable URL.
data_all$URL <- paste("http://flickr.com/photos/",data_all$owner,"/",data_all$id,sep="")
data_2015$URL <- paste("http://flickr.com/photos/",data_2015$owner,"/",data_2015$id,sep="")
data_2011$URL <- paste("http://flickr.com/photos/",data_2011$owner,"/",data_2011$id,sep="")


##Adds UTM coordinates.
data_all$latitude <- as.numeric(data_all$latitude)
data_all$longitude <- as.numeric(data_all$longitude)

##Filters nonsense coordinates.
data_all <- subset(data_all,latitude > 44 & latitude < 48)
data_all <- subset(data_all, longitude > -124 & longitude < -120)

coordinates(data_all) <- ~longitude+latitude
proj4string(data_all) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
data_all_UTM <- spTransform(data_all,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_all_UTM$UTME <- coordinates(data_all_UTM)[,1]
data_all_UTM$UTMN <- coordinates(data_all_UTM)[,2]
data_all_UTM$lat <- coordinates(data_all)[,1]
data_all_UTM$long <- coordinates(data_all)[,2]

data_2015$latitude <- as.numeric(data_2015$latitude)
data_2015$longitude <- as.numeric(data_2015$longitude)
coordinates(data_2015) <- ~longitude+latitude
proj4string(data_2015) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
data_2015_UTM <- spTransform(data_2015,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_2015_UTM$UTME <- coordinates(data_2015_UTM)[,1]
data_2015_UTM$UTMN <- coordinates(data_2015_UTM)[,2]
data_2015_UTM$lat <- coordinates(data_2015)[,1]
data_2015_UTM$long <- coordinates(data_2015)[,2]

data_2011$latitude <- as.numeric(data_2011$latitude)
data_2011$longitude <- as.numeric(data_2011$longitude)
coordinates(data_2011) <- ~longitude+latitude
proj4string(data_2011) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
data_2011_UTM <- spTransform(data_2011,CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
data_2011_UTM$UTME <- coordinates(data_2011_UTM)[,1]
data_2011_UTM$UTMN <- coordinates(data_2011_UTM)[,2]
data_2011_UTM$lat <- coordinates(data_2011)[,1]
data_2011_UTM$long <- coordinates(data_2011)[,2]

##Writes the converted data to disk
write.csv(data_2015_UTM@data,"MORA_flickr_metadata_all_2015.csv")
write.csv(data_2011_UTM@data,"MORA_flickr_metadata_all_2011.csv")
write.csv(data_all_UTM@data,"MORA_flickr_metadata_all_2009_2015.csv")


