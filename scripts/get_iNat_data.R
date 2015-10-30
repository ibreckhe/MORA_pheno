##Script to grab iNaturalist data for Mt. Rainier National Park

####Sets up worspace####
setwd("~/code/MORA_pheno/")
library(rinat)
library(dplyr)

####Grabs all geotagged taxon records for Mt. Rainier National Park
bbox <- c(46.695609,-121.993775, 47.014907,-121.423867)
plantID <- 47126
animalID <- 1
fungiID <- 47170

plants <- get_inat_obs(taxon_id=plantID,bounds=bbox,maxresults=5000)
animals <- get_inat_obs(taxon_id=animalID,bounds=bbox,maxresults=5000)
fungi <- get_inat_obs(taxon_id=fungiID,bounds=bbox,maxresults=5000)

all <- rbind(plants,animals,fungi)
all$date <- as.POSIXct(all$datetime,format="%Y-%m-%d %H:%M:%S")

####Filters records####
all_f <- filter(all,date > as.POSIXct("2009-01-01"))

####Writes to disk####
write.csv(all_f,"./data/iNat_all_2009_2015.csv",row.names=FALSE)
