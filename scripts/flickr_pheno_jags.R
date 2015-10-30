######################
####Script to fit distributions and gaussian curves to Flickr data#####
####Author: Ian Breckheimer
####Date: 14 April 2014
######################

####Sets up the workspace####
library(plyr)
library(ggplot2)
library(rjags)
library(ggmcmc)
library(foreach)
library(doParallel)
setwd("~/Dropbox/Lab/EcoForecasting_SDD_Phenology (1)/Data&Analysis/Phenology/Flickr&iNaturalist/cleaned")

####Sources my functions.
source("~/code/MORA_microclimate/phenology_functions.R")

####Reads in the data####
fsites <- read.csv("flickr_inat_all_classed_2001_2013.csv")

####Filters the data.
fsites <- fsites[-c(which(fsites$Species != "NF" & fsites$dss < -30)),]
fsites <- fsites[-c(which(fsites$Species != "NF" & fsites$can_pct > 0.8)),]
fsites <- fsites[-c(which(is.na(fsites$dss))),]
#fsites <- fsites[-c(which(fsites$Owner=="88227046@N00")),]

####Standardizes phenology stages
fsites$Phenological.Stage[which(fsites$Phenological.Stage %in% c("Flowering",
                                                                 "Flowering ",
                                                                 "FLOWERING"))] <- "Flowering"
fsites$Phenological.Stage[which(fsites$Phenological.Stage %in% c("Budding",
                                                                 "Budding ",
                                                                 "BUDDING"))] <- "Budding"
fsites$Phenological.Stage[which(fsites$Phenological.Stage %in% c("Ripening Fruit",
                                                                 "Ripening Fruit ",
                                                                 "Ripening Fruit  ",
                                                                 "RIPENING FRUIT"))] <- "Ripening Fruit"
fsites$Phenological.Stage[which(fsites$Phenological.Stage %in% c("Releasing Seed",
                                                                 "Releasing Seed ",
                                                                 "RELEASING SEED"))] <- "Releasing Seed"

##Converts date taken to julian day.
fsites$datetaken <- strptime(fsites$datetaken,format="%m/%d/%y")
fsites$datetaken <- as.POSIXct(fsites$datetaken) # convert to POSIXct to use in data frames / ddply
fsites$datetaken_DOY <- as.numeric(strftime(fsites$datetaken,format="%j"))

##Creates a unique factor for photos taken on a particular 
##date by a particular person.
fsites$ownerdate <- factor(paste(fsites$Owner,fsites$datetaken_DOY,sep="_"))

##Creates a unique factor for each photo, since there are multiple records per photo.
fsites$photo <- factor(paste(fsites$ownerdate,fsites$UTM_E,
                             fsites$UTM_N,fsites$ID,fsites$Title,sep="_"))

# Focal species with at least 100 records
summary(fsites$Species)[order(summary(fsites$Species))]

spp <- c("Anemone occidentalis",
         "Castilleja parviflora",
         "Erigeron peregrinus",
         "Erythronium montanum",
         "Lupinus arcticus",
         "Pedicularis bracteosa",
         "Polygonum bistortoides",
         "Valeriana sitchensis") 

##Removes duplicate records of photos.
fsites_d <- subset(fsites, select=-c(ID.Person))
drop_d <- which(duplicated(fsites_d))
fsites_nd <- fsites_d[-drop_d,]

##Drops records of non-focal species
spp <- c("Anemone occidentalis",
         "Castilleja parviflora",
         "Erigeron peregrinus",
         "Erythronium montanum",
         "Ligusticum grayi",
         "Lupinus arcticus",
         "Microseris alpestris",
         "Pedicularis bracteosa",
         "Polygonum bistortoides",
         "Valeriana sitchensis",
         "NF")
nspp <- unique(fsites_nd$Species)[!unique(fsites_nd$Species) %in% spp]
drop_spp <- which(fsites_nd$Species %in% nspp)
drop_fl <- which(fsites_nd$Species %in% spp[1:10] & fsites_nd$Flowering == 0 )
fsites_sp <- fsites_nd[-c(drop_spp,drop_fl),]

##Plots flowering and not flowering photos by year.
ggplot(data=subset(fsites_sp,Year<2013))+
  geom_point(aes(x=DOY,y=Flowering,color=Year),
             position=position_jitter(height=0.05),
             size=1)+
  scale_color_distiller(type="qual")

##Drops photos out of the DSS range of interest.
drop_dss <- which(fsites_sp$dss < -30 | fsites_sp$dss > 130)
fsites_ds <- fsites_sp[-c(drop_dss),]

##Creates flowering presence-absence records for each species.
fsites_m <- fsites_ds[,c("Species","photo","ownerdate","Flowering")]
fsites_d <- dcast(fsites_m,formula=photo~Species,)

##Assumes species is absent from photo if not recorded
fsites_d[is.na(fsites_d)] <- 0

##Puts data back in long format
fsites_m2 <- melt(fsites_d)
colnames(fsites_m2) <- c("photo","Species","FlwrYN")
fsites_fspp <- merge(fsites_m2,fsites_sp[,c("photo","ownerdate","dss",
                                           "datetaken_DOY","Year","elev",
                                           "UTM_E","UTM_N")],all.x=TRUE)
fsites_fspp <- fsites_fspp[-which(fsites_fspp$Species == "NF"),]
fsites_comp <- fsites_fspp[complete.cases(fsites_fspp),]

##Plots data.
par(mfrow=c(1,1))
plot(fsites_comp$dss,jitter(fsites_comp$FlwrYN,amount=0.1),col=fsites_comp$Species,pch=".")

##Subsets data.
set.seed(38)
fsites_sub <- fsites_comp[sample(1:dim(fsites_comp)[1],size=8000,replace=FALSE),]

##Preps sample of data for jags.
dss <- scale(fsites_sub$dss)
dss_scale <- attr(dss,"scaled:scale")
dss_center <- attr(dss,"scaled:center")
species <- factor(fsites_sub$Species)
flowering <- fsites_sub$FlwrYN
photo <- factor(fsites_sub$photo)
visit <- factor(fsites_sub$ownerdate)
year <- fsites_sub$Year

##Preps full data for jags.
dss <- scale(fsites_comp$dss)
dss_scale <- attr(dss,"scaled:scale")
dss_center <- attr(dss,"scaled:center")
species <- factor(fsites_comp$Species)
flowering <- fsites_comp$FlwrYN
photo <- factor(fsites_comp$photo)
visit <- factor(fsites_comp$ownerdate)
year <- fsites_comp$Year

###GLM fit for comparison.
fl_glm_allspp <- glm(flowering~dss+I(dss^2)*species,family="binomial")

fl_jags_allspp <- fit.jags.mixed.allspp.fl(x=dss,y=flowering,groups=photo,observers=visit,species=species,
                         years=year,nsamples=1000)
fl_jags_up <- update.jags.mixed(fl_jags_allspp,n.update=5000,n.iter=5000,thin=5,
                  params=c("height.s","opt.s","width.s","width.s.mu",
                           "width.s.sigma","height.s.mu",
                           "height.s.sigma","opt.s.mu","opt.s.sigma",
                           "height.g.sigma","opt.g.sigma","width.g.sigma"))
fl_jags_up <- update.jags.mixed(fl_jags_up,n.update=30000,n.iter=30000,thin=30,
                                params=c("height.s","opt.s","width.s","width.s.mu",
                                         "width.s.sigma","height.s.mu",
                                         "height.s.sigma","opt.s.mu","opt.s.sigma",
                                         "height.g.sigma","opt.g.sigma","width.g.sigma"))
save(fl_jags_up,file="fl_jagsoutput.Rdata")
load(file="fl_jagsoutput.Rdata")

## Makes a function matrix for each species.
fl_fun_mat <- make.fun.matrix(fl_jags_up$out,par.names=c("width.s","opt.s","height.s"),
                              n.samples=1000)
fl_fun_mat_med <- make.med.fun.matrix(fl_jags_up$out,par.names=c("width.s","opt.s","height.s"))

## Makes community mean function matrix.
fl_fun_mat_com <- make.community.fun.matrix(fl_jags_up$out,n.samples=1000)
fl_fun_com_med <- make.med.fun.matrix(fl_jags_up$out,par.names=c("width.s.mu","opt.s.mu","height.s.mu"),
                                      n.samples=1,spp=FALSE)[1,1][[1]]

##Creates predictions based on median parameter estimates.
jagsestimates.allspp <- foreach(i=1:length(levels(species)),.combine="rbind") %do% {
  newx <- seq(-5,5,length.out=10000)
  newx_unscaled <- newx * dss_scale + dss_center
  pred.fun <- fl_fun_mat_med[,i][[1]]
  pred_y <- pred.fun(newx)
  pred_y_unscaled <- pred_y * dss_scale + dss_center
  out_preds <- data.frame(species=levels(species)[i],dss=newx_unscaled,y=pred_y)
  colnames(out_preds) <- c("Species","dss","y")
  (out_preds)
}

##Creates credible intervals for each curve.
quantfun <- function(x){quantile(x,probs=c(0.025,0.25,0.50,0.75,0.975))}
jagscurves.allspp <- foreach(i=1:length(levels(species)),.combine="rbind") %do% {
  newx <- seq(-5,5,length.out=10000)
  newx_unscaled <- newx * dss_scale + dss_center
  pred.y <- matrix(NA,ncol=length(newx),nrow=dim(fl_fun_mat)[1])
  for(j in 1:dim(fl_fun_mat)[1]){
    pred.fun <- fl_fun_mat[j,i][[1]]
    pred.y[j,] <- pred.fun(newx)
  }
  out_curves <- t(apply(pred.y,FUN=quantfun,MARGIN=2))
  out_frame <- data.frame(Species=levels(species)[i],
                          x=newx_unscaled,
                          pred_y=out_curves[,3],
                          y_lwr50=out_curves[,2],
                          y_upr50=out_curves[,4],
                          y_lwr95=out_curves[,1],
                          y_upr95=out_curves[,5])
  (out_frame)
}

##Creates credible intervals for height,optimum, and tolerance.
quantfun <- function(x){quantile(x,probs=c(0.025,0.25,0.50,0.75,0.975))}
jagsintervals.allspp <- foreach(i=1:length(levels(species)),.combine="rbind") %do% {
  newx <- seq(-5,5,length.out=10000)
  newx_unscaled <- newx * dss_scale + dss_center
  intervals_mat <- matrix(NA,nrow=dim(fl_fun_mat)[1],ncol=5)
  for(j in 1:dim(fl_fun_mat)[1]){
    pred.fun <- fl_fun_mat[j,i][[1]]
    preds_df <- data.frame(dss=newx_unscaled,pred=pred.fun(newx))
    intervals <- obs_intervals(preds_df)
    intervals_mat[j,1] <- intervals[1]
    intervals_mat[j,2] <- intervals[2]
    intervals_mat[j,3] <- intervals[2] - intervals[1]
    intervals_mat[j,4] <- preds_df$dss[which.max(preds_df$pred)]
    intervals_mat[j,5] <- max(preds_df$pred)
  }
  out_intervals <- t(apply(intervals_mat,FUN=quantfun,MARGIN=2))
  out_frame <- data.frame(Species=levels(species)[i],
                          Stat=c("lwr","lwr25","est","upr75","upr"),
                          lwr_bound = out_intervals[1,],
                          upr_bound = out_intervals[2,],
                          width = out_intervals[3,],
                          opt = out_intervals[4,],
                          height = out_intervals[5,])
  (out_frame)
}
jagsintervals_m.allspp <- melt(jagsintervals.allspp,id.vars=c("Species","Stat"))
jagsintervals_c.allspp <- dcast(jagsintervals_m.allspp,formula=Species~Stat+variable)
write.csv(jagsintervals_c.allspp,file="fl_jagsparams_allspp.csv",row.names=FALSE)

##Creates median predictions for the community estimate.
newx <- seq(-5,5,length.out=10000)
newx_unscaled <- newx * dss_scale + dss_center
pred.fun <- fl_fun_com_med
pred_y <- pred.fun(newx)
pred_y_unscaled <- pred_y * dss_scale + dss_center
out_preds <- data.frame(species="Community",dss=newx_unscaled,y=pred_y)
colnames(out_preds) <- c("Species","dss","y")
jagsestimates.com <- out_preds
jagsestimates.com.spp <- rbind(jagsestimates.allspp,jagsestimates.com)

##Creates credible intervals for the community estimate.
pred.y <- matrix(NA,ncol=length(newx),nrow=dim(fl_fun_mat_com)[1])
for(j in 1:dim(fl_fun_mat_com)[1]){
  pred.fun <- fl_fun_mat_com[j,1][[1]]
  pred.y[j,] <- pred.fun(newx)
}
out_curves <- t(apply(pred.y,FUN=quantfun,MARGIN=2))
out_frame <- data.frame(Species="Community",
                        x=newx_unscaled,
                        pred_y=out_curves[,3],
                        y_lwr50=out_curves[,2],
                        y_upr50=out_curves[,4],
                        y_lwr95=out_curves[,1],
                        y_upr95=out_curves[,5])
jagscurves.com <- out_frame
jagscurves.com.spp <- rbind(jagscurves.allspp,jagscurves.com)

##Community estimates of optimum and duration.
intervals_mat <- matrix(NA,nrow=dim(fl_fun_mat_com)[1],ncol=5)
for(j in 1:dim(fl_fun_mat_com)[1]){
  pred.fun <- fl_fun_mat_com[j,1][[1]]
  preds_df <- data.frame(dss=newx_unscaled,pred=pred.fun(newx))
  intervals <- obs_intervals(preds_df)
  intervals_mat[j,1] <- intervals[1]
  intervals_mat[j,2] <- intervals[2]
  intervals_mat[j,3] <- intervals[2] - intervals[1]
  intervals_mat[j,4] <- preds_df$dss[which.max(preds_df$pred)]
  intervals_mat[j,5] <- max(preds_df$pred)
}
out_intervals <- t(apply(intervals_mat,FUN=quantfun,MARGIN=2))
jagsintervals.com <- data.frame(Species="Community",
                                Stat=c("lwr","lwr25","est","upr75","upr"),
                                lwr_bound = out_intervals[1,],
                                upr_bound = out_intervals[2,],
                                width = out_intervals[3,],
                                opt = out_intervals[4,],
                                height = out_intervals[5,])
jagsintervals_m.com <- melt(jagsintervals.com,id.vars=c("Species","Stat"))
jagsintervals_c.com <- dcast(jagsintervals_m.com,formula=Species~Stat+variable)
jagsintervals_c.com.spp <- rbind(jagsintervals_c.allspp,jagsintervals_c.com)
write.csv(jagsintervals_c.com.spp,file="fl_jagsparams_allspp.csv",row.names=FALSE)

##Data for plotting.
fsites_com <- fsites_comp
fsites_com$Species <- "Community"
fsites_com <- rbind(fsites_comp,fsites_com)

##Plots the fit curves and credible intervals.
pdf("fl_curves_jags.pdf",width=10,height=8)
ggplot(jagscurves.com.spp)+
  facet_wrap(~Species)+
  geom_point(aes(x=dss,y=FlwrYN,color=Species),data=fsites_com,size=0.8,
             position=position_jitter(w=0.01,h=0.01))+
  geom_ribbon(aes(x=x,ymin=y_lwr95,ymax=y_upr95,bg=Species),alpha=0.5,data=jagscurves.com.spp)+
  geom_line(aes(x=dss,y=y,color=Species),data=jagsestimates.com.spp)+
  geom_errorbarh(aes(x=est_lwr_bound,xmin=lwr_lwr_bound,xmax=upr_lwr_bound,y=-0.11,yend=-0.11,col=Species),
                 data=jagsintervals_c.com.spp,lwd=0.5,height=0.01)+
  geom_errorbarh(aes(x=est_upr_bound,xmin=lwr_upr_bound,xmax=upr_upr_bound,y=-0.09,yend=-0.09,col=Species),
                 data=jagsintervals_c.com.spp,lwd=0.5,height=0.01)+
  geom_errorbarh(aes(x=est_opt,xmin=lwr_opt,xmax=upr_opt,y=-0.04,yend=-0.04,col=Species),
                 data=jagsintervals_c.com.spp,lwd=0.5,height=0.01)+
  geom_segment(aes(x=est_lwr_bound,xend=est_upr_bound,y=-0.06,yend=-0.06,col=Species),
               data=jagsintervals_c.com.spp,alpha=0.5,lwd=2)+
  geom_point(aes(x=est_opt,y=-0.06,color=Species),pch=21,fill="white",size=3,
             lwd=0.5,data=jagsintervals_c.com.spp)+
  coord_cartesian(xlim=c(-5,99),ylim=c(-0.2,1.1))+
  theme_bw()+
  theme(legend.key=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_blank())
dev.off()


