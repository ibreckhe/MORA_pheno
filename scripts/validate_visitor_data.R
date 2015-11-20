##Script to compare monthly flickr users to monthly visitor data.
##Ian Breckheimer
##14 November 2015

####Sets up workspace####
library(ggplot2)
library(mgcv)
library(dplyr)

####Loads and munges data####
d_f <- read.csv("./data/MORA_flickr_metadata_all_2009_2015_cleaned.csv")
v_p <- read.csv("./data/MORA_Visitor_Use_Monthly_2009_2014.csv")


##Calculates the number of unique users per month for every year.##
v_p$Year <- factor(v_p$Year)
d_f$datetaken <- as.POSIXct(d_f$datetaken)
d_f$month <- as.numeric(format(d_f$datetaken,format="%m"))
d_grp <- group_by(d_f,year,month)
d_sum <- summarize(d_grp,
                   nphotos=length(id),
                   nusers=length(unique(owner)),
                   mean_dss=mean(days_since_snow),
                   mean_doy=mean(datetaken_DOY),
                   mean_travel_t=mean(acc_times))
d_sum$Year<- factor(d_sum$year)

####Merges visitor info with flickr data.####
d_sum <- left_join(d_sum,v_p,by=c("month"="Numeric.Month","Year"="Year"))

d_sum <- filter(d_sum,year != 2015)

###Plots relationship for each year####
pdf("./figs/flickr_visitor_validation_monthly.pdf",width=7,height=4)
ggplot(data=d_sum)+
  geom_point(aes(x=nusers,y=Recreation.Visitors/1000,color=Year))+
  geom_smooth(aes(x=nusers,y=Recreation.Visitors/1000,color=Year),
              method="lm",formula=y~x-1,se=TRUE)+
#  scale_x_continuous(breaks=seq(0,150,by=50),trans="log")+
#  scale_y_continuous(trans="log")+
  xlab("Unique Flickr Users")+
  xlim(c(0,124))+
  ylab("Recreation Visitors (thousands)")+
  facet_wrap(facets=~Year)+
  theme_bw()
dev.off()

###Linear model####
visit_lm <- lm(Recreation.Visitors~nusers*yearfact - 1,data=d_sum)
  
  
