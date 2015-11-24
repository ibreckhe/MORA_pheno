##Script to make non-parametric graphs of visitors and flowers.
##Ian Brecheimer
##20 November 2016

##Sets up Workspace
library(Hmisc)
library(dplyr)
library(ggplot2)

##Loads data
flowers <- read.csv("./data/MORA_flickr_classified_2011_2015.csv")
visitors <-read.csv("./data/MORA_flickr_metadata_all_2009_2015_cleaned.csv")

##Filters flower data.
flowers$locfact <- as.factor(paste(flowers$long,flowers$lat,sep="by"))
flowers_grp <- group_by(flowers,locfact)
nphotos <- summarize(flowers_grp,nphotos=n())
small_locs <- nphotos$locfact[which(nphotos$nphotos < 20)]


###Creates a new dataset indicating whether there are any flowers at all.
flowers <- filter(flowers,Phenophase!="np"&
                  Species != "other" &
#                  days_since_snow > -100 &
#                  datetaken_DOY > 90 &
#                  datetaken_DOY < 305 &
                  #                      year != 2012 &
                  #                      year != 2014 &
                  locfact %in% small_locs)

##Reduces to one record per flower.
f_grp <- group_by(flowers,id,owner,year,datetaken_DOY,elevation,
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
                    ANOC = any(Species== "anemone occidentalis" & 
                                   Phenophase == "flowering"),
                    ERMO = any(Species == "erythronium montanum" & 
                                 Phenophase == "flowering"),
                    CAPA = any(Species == "castilleja parviflora" & 
                                 Phenophase == "flowering"),
                    VASI = any(Species == "valeriana sitchensis" & 
                                 Phenophase == "flowering"),
                    ERPE = any(Species == "erigeron peregrinus" & 
                                 Phenophase == "flowering"),
                    LUAR= any(Species == "lupinus arcticus" & 
                                 Phenophase == "flowering"),
                    POBI= any(Species == "polygonum bistortoides" & 
                                Phenophase == "flowering"),
                    MIAL= any(Species == "microseris alpestris" & 
                                Phenophase == "flowering"),
                    PEBR= any(Species == "pedicularis bracteosa" & 
                                Phenophase == "flowering"),
                    LIGR= any(Species == "ligusticum grayi" & 
                                Phenophase == "flowering"),
                    nSpp = sum(as.numeric(Phenophase=="flowering")))
f_flwr$Flwr <- as.numeric(f_flwr$FlwrYN)
f_flwr$Early <- as.numeric(f_flwr$EarlyYN)
f_flwr$Late <- as.numeric(f_flwr$LateYN)
f_flwr$notEarly <- as.numeric(f_flwr$NoEarly)
f_flwr$notLate <- as.numeric(f_flwr$NoLate)
f_flwr$yearfact <- as.factor(f_flwr$year)
f_flwr$elevfact <- as.factor(f_flwr$elevation >= 1650)

##Creates the plot.
month_breaks <- c(1,32,60,91,121,152,182,213,244,274,305,335)
month_labels <- c("Jan","","March","","May","",
                  "July","","Sept","","Nov","")
pdf("./figs/photo_density_2011_2015.pdf",width=4,height=10)
ggplot(f_flwr)+
  geom_density(aes(x=datetaken_DOY),fill="grey40",adjust=3)+
#  geom_point(aes(x=datetaken_DOY,y=Flwr,color=yearfact),alpha=0.1,
#             position=position_jitter(height=0.1))+
#  geom_smooth(aes(x=datetaken_DOY,y=Flwr,color=yearfact),
#              method="glm",formula=y~poly(x,2),se=TRUE,
#              method.args=list(family="binomial"))+
  scale_x_continuous(breaks = month_breaks,labels = month_labels)+
  xlab("")+
  facet_grid(facets=yearfact~.)+
  theme_bw()
dev.off()

pdf("./figs/flower_div_2011_2015.pdf",width=4,height=10)
ggplot(f_flwr)+
   geom_point(aes(x=datetaken_DOY,y=Flwr,group=yearfact),color="black",alpha=0.1,
              position=position_jitter(height=0.1))+
   geom_smooth(aes(x=datetaken_DOY,y=Flwr,group=yearfact),color="black",
               method="glm",formula=y~poly(x,2),se=TRUE,
               method.args=list(family="binomial"),data=f_flwr)+
  scale_x_continuous(breaks = month_breaks,labels = month_labels)+
  scale_y_continuous(limits=c(0,1),labels=function(x) {format(x,nsmall=0)}) +
  xlab("")+
  ylab("Diversity of Focal Species")+
  facet_grid(facets=yearfact~.)+
  theme_bw()
dev.off()

##Plot by DSS.
pdf("./figs/photo_DSS_2011_2015.pdf",width=4,height=10)
ggplot(f_flwr)+
  geom_density(aes(x=days_since_snow),fill="grey40",adjust=2.2)+
  # geom_point(aes(x=days_since_snow,y=as.numeric(LUAR)),color="purple",alpha=0,
  #            position=position_jitter(height=0.1))+
  # geom_smooth(aes(x=days_since_snow,y=as.numeric(LUAR)),color="purple",
  #             method="glm",formula=y~poly(x,2),se=FALSE,
  #             method.args=list(family="binomial"))+
  # geom_point(aes(x=days_since_snow,y=as.numeric(LIGR)),color="green",alpha=0,
  #            position=position_jitter(height=0.1))+
  # geom_smooth(aes(x=days_since_snow,y=as.numeric(LIGR)),color="green",
  #             method="glm",formula=y~poly(x,2),se=FALSE,
  #             method.args=list(family="binomial"))+
  # geom_point(aes(x=days_since_snow,y=as.numeric(CAPA)),color="red",alpha=0,
  #            position=position_jitter(height=0.1))+
  # geom_smooth(aes(x=days_since_snow,y=as.numeric(CAPA)),color="red",
  #             method="gam",formula=y~poly(x,2),se=FALSE,
  #             method.args=list(family="binomial"))+
  # geom_point(aes(x=days_since_snow,y=as.numeric(ERMO)),color="orange",alpha=0,
  #            position=position_jitter(height=0.1))+
  # geom_smooth(aes(x=days_since_snow,y=as.numeric(ERMO)),color="orange",
  #             method="gam",formula=y~poly(x,2),se=FALSE,
  #             method.args=list(family="binomial"))+
  scale_x_continuous(limits=c(-30,140))+
  scale_y_continuous(labels=function(x) {format(x,nsmall=3)}) +
  xlab("Days Since Snow")+
  ylab("Photo Density")+
  facet_grid(facets=year~.)+
  theme_bw()
dev.off()

##Plot by DSS of individual species.
pdf("./figs/flower_DSS_2011_2015.pdf",width=4,height=10)
ggplot(f_flwr)+
#  geom_density(aes(x=days_since_snow),fill="grey40",adjust=2.2)+
geom_point(aes(x=days_since_snow,y=as.numeric(LUAR)),color="purple",alpha=0,
           position=position_jitter(height=0.1))+
geom_smooth(aes(x=days_since_snow,y=as.numeric(LUAR)),color="purple",
            method="glm",formula=y~poly(x,2),se=FALSE,
            method.args=list(family="binomial"))+
geom_point(aes(x=days_since_snow,y=as.numeric(LIGR)),color="green",alpha=0,
           position=position_jitter(height=0.1))+
geom_smooth(aes(x=days_since_snow,y=as.numeric(LIGR)),color="green",
            method="glm",formula=y~poly(x,2),se=FALSE,
            method.args=list(family="binomial"))+
geom_point(aes(x=days_since_snow,y=as.numeric(CAPA)),color="red",alpha=0,
           position=position_jitter(height=0.1))+
geom_smooth(aes(x=days_since_snow,y=as.numeric(CAPA)),color="red",
            method="gam",formula=y~poly(x,2),se=FALSE,
            method.args=list(family="binomial"))+
geom_point(aes(x=days_since_snow,y=as.numeric(ERMO)),color="orange",alpha=0,
           position=position_jitter(height=0.1))+
geom_smooth(aes(x=days_since_snow,y=as.numeric(ERMO)),color="orange",
            method="gam",formula=y~poly(x,2),se=FALSE,
            method.args=list(family="binomial"))+
scale_x_continuous(limits=c(-30,140))+
  scale_y_continuous(labels=function(x) {format(x,nsmall=3)}) +
  xlab("Days Since Snow")+
  ylab("Photo Density")+
  facet_grid(facets=year~.)+
  theme_bw()
dev.off()
