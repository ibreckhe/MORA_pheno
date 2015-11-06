##Script to join classified flickr data for 2011 to 2015.
##Ian Breckheimer
##November 6th, 2015

####Sets up workspace####
library(dplyr)

####Read and munge data####

##Brings in all of the subalpine records from 2009 to 2015##
f <- as.tbl(read.csv("./data/MORA_flickr_metadata_all_2009_2015_cleaned.csv"))
f$datePOSIX <- as.POSIXct(f$datetaken)

##Brings in all of the classified data
c2011 <- as.tbl(read.csv("./data/MORA_flickr_subalp_2011_classified_cleaned.csv"))
c2012 <- as.tbl(read.csv("./data/MORA_flickr_classified_2012_PresOnly.csv"))
c2013 <- as.tbl(read.csv("./data/MORA_flickr_classified_all_2013_covars.csv"))
c2014 <- as.tbl(read.csv("./data/MORA_flickr_2014_classified_cleaned.csv"))

####Munges data for 2011####
f2011 <- filter(f,year==2011)
f2011$UniqueID <- paste(f2011$owner,f2011$id,sep=".")

c2011$UniqueID <- paste(c2011$owner,c2011$id,sep=".")
d2011_cols <- select(c2011,UniqueID,Collector=Collector,Species,Phenophase,Notes)

d2011 <- left_join(f2011,d2011_cols,by="UniqueID")

####Munges data for 2012####
##Brings in the 4-letter species codes
scodes <- as.tbl(read.csv("./data/2009_2013_species_codes.csv"))

##Merge species codes with 2012 data.
c2012 <- left_join(c2012,scodes,by=c("Species" = "Code"))

##Creates a date field
c2012$datePOSIX <- as.POSIXct(c2012$DateTaken,format="%m/%d/%y")

##Creates a unique factor for location title and owner.
c2012$UniqueID <- paste(c2012$Title,c2012$datePOSIX,c2012$Latitude,c2012$Longitude,sep=".")

f2012 <- filter(f,year==2012)
f2012$UniqueID <- paste(f2012$title,f2012$datePOSIX,f2012$lat,f2012$long,sep=".")

##Joins observations to new data.
c2012_cols <- select(c2012,UniqueID,Collector,Species = Name,Phenophase,Notes)
d2012 <- left_join(f2012,c2012_cols,by="UniqueID")

####Munges data for 2013####
f2013 <- filter(f,year==2013)
f2013$UniqueID <- paste(f2013$owner,f2013$id,sep=".")

c2013$UniqueID <- paste(c2013$owner,c2013$id,sep=".")
d2013_cols <- select(c2013,UniqueID,Species=SPECIES,Phenophase=PHEN_PHASE)
d2013 <- left_join(f2013,d2013_cols,by="UniqueID")
d2013$Collector <- "Sam"
d2013$Notes <- NA

####Munges data for 2014####
f2014 <- filter(f,year==2014)
f2014$UniqueID <- paste(f2014$owner,f2014$id,sep=".")

c2014$UniqueID <- paste(c2014$owner,c2014$id,sep=".")
d2014_cols <- select(c2014,UniqueID,Collector=CollectorName,Species,Phenophase,Notes)

d2014 <- left_join(f2014,d2014_cols,by="UniqueID")

####Merges data from different years.####
d_classed <- rbind(d2011,d2012,d2013,d2014)

###Converts species names to lowercase
d_classed$Phenophase <- as.factor(tolower(d_classed$Phenophase))
d_classed$Species <- as.factor(tolower(d_classed$Species))

###Corrects some misspelled species names.
d_classed$Species[d_classed$Species=="eucephalus ledophyllus"] <- "aster ledophyllus"
d_classed$Species[d_classed$Species=="potentilla flavbellifolia"] <- "potentilla flabellifolia"
d_classed$Species[d_classed$Species=="eucephalus ledophyllus"] <- "aster ledophyllus"
d_classed$Species[d_classed$Species=="nothocalais alpstris"] <- "microseris alpestris"

focal <- c("lupinus arcticus",
           "polygonum bistortoides",
           "castilleja parviflora",
           "anemone occidentalis",
           "valeriana sitchensis",
           "erigeron peregrinus",
           "ligusticum grayi",
           "pedicularis bracteosa",
           "erythronium montanum",
           "microseris alpestris",
           "other",
           "nf",
           "np",
           NA)

##Filters output.
`%nin%` <- Negate(`%in%`)
d_classed$Species <- as.character(d_classed$Species)
d_classed$Species[which(d_classed$Species %nin% focal)] <- "other"

d_classed$Species <- as.factor(d_classed$Species)
d_classed <- distinct(d_classed)

####Writes output to disk####
write.csv(d_classed,"./data/MORA_flickr_classified_2011_2014.csv",row.names=FALSE) 


