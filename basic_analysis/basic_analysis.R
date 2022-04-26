#------------------------------------------------------------------------
#script for basic analysis of the data I use
#------------------------------------------------------------------------

#load packages and sources
library(animation)
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/colormaps.R")
# ------------------------------------------------------------------------------------------ 
# 00. Read model data
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")

for (region in regions){
  
print(region)
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
figpath=paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/",region)

setwd(path)
X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST

#read ncdf files into RB objects
X.anom.RB=brick(paste0(X, "_",region,"_anom.nc"))
X.data.RB = brick(paste0(X, "_",region,".nc")) 
Y.anom.RB =brick(paste0(Y,"_",region,"_anom.nc"))
Y.data.RB = brick(paste0(Y,"_",region,".nc"))



# ------------------------------------------------------------------------------------------ 
# 01. Preprocessing: Define regions, crop data, area weighting and mean
# ------------------------------------------------------------------------------------------

#crop data to smaller region
if (region=="Atlantic"){
  Y.extent = extent(c(-70, 0, 20, 60))
  X.extent = Y.extent + 40
} else if (region=="Pacific") {
  Y.extent = extent(c(140,240, 20, 60))
  X.extent = Y.extent + 40
}

X.data.region.RB = crop(X.data.RB, X.extent)
Y.data.region.RB = crop(Y.data.RB, X.extent)
X.anom.region.RB = crop(X.anom.RB, X.extent)
Y.anom.region.RB = crop(Y.anom.RB, X.extent)


#SVD
#X.svd=svd(values(area.weight.RB(X.anom.region.RB,svd=TRUE)))

#area weighting
Y.anom.mean.aw=mean.area.weighted.RB(Y.anom.region.RB)
Y.data.mean.aw=mean.area.weighted.RB(Y.data.region.RB)


# ------------------------------------------------------------------------------------------ 
# 02. analysis
# ------------------------------------------------------------------------------------------
#set up for figures
setwd(figpath)

if(region=="Pacific") {map="world2";fig.height=350;fig.width=470}
if(region=="Atlantic") {map="world";fig.height=350;fig.width=380}

#spatial raster for figures
X.spatial.raster=raster(X.data.region.RB,1)
Y.spatial.raster=raster(Y.data.region.RB,1)

for (month in c(1:12)){
  #calculate climatologies
  X.climatology.RB=mean(subset(X.data.region.RB,seq(month,24000,12)))
  Y.climatology.RB=mean(subset(Y.data.region.RB,seq(month,24000,12)))
  
  png(paste0("PSL_climatology_",month,"_",region,".png"),height=fig.height,width=fig.width)
  par(mar=c(2,2,1,1),oma=c(0,0,0,0))
  image(X.climatology.RB/100,zlim=c(990,1040),col=color.warm.cold,main="")
  maps::map(map,add=T,interior=F)
  dev.off()
  
  png(paste0("SST_climatology_",month,"_",region,".png"),height=fig.height,width=fig.width)
  par(mar=c(2,2,1,1),oma=c(0,0,0,0))
  image(Y.climatology.RB,zlim=c(-2,35),col=color.warm.cold,main="")
  maps::map(map,add=T,interior=F,fill=T)
  dev.off()
  
  #plot climatologies together
  png(paste0("SST_PSL_climatology_",month,"_",region,".png"),height=fig.height,width=fig.width)
  X.plot.RB=X.climatology.RB/100
  Y.plot.RB=Y.climatology.RB
  levels=pretty(values(X.plot.RB),n=6);levels=levels[which(levels!=0)]
  lty=rep("solid",length(levels));lty[which(levels<1013)]="dashed"
  par(mar=c(2,2,1,1),oma=c(0,0,0,0))
  image(Y.plot.RB,col=color.warm.cold,zlim=c(-2,35))
  contour(X.plot.RB,levels=levels,lty="solid",add=T,lwd=2,labcex=1)
  maps::map(map,add=T,interior=F)
  dev.off()
}

for (month in c(1:12)){
  #calculate sd
  X.sd.RB=X.spatial.raster
  Y.sd.RB=Y.spatial.raster
  
  values(X.sd.RB)=rowSds(values(subset(X.data.region.RB,seq(month,24000,12))))
  values(Y.sd.RB)=rowSds(values(subset(Y.data.region.RB,seq(month,24000,12))))
  
  png(paste0("PSL_sd_",month,"_",region,".png"),height=fig.height,width=fig.width)
  par(mar=c(2,2,1,1),oma=c(0,0,0,0))
  image(X.sd.RB/100,zlim=c(0,10),col=color.high.low,main="")
  maps::map(map,add=T,interior=F)
  dev.off()
  
  png(paste0("SST_sd_",month,"_",region,".png"),height=fig.height,width=fig.width)
  par(mar=c(2,2,1,1),oma=c(0,0,0,0))
  image(Y.sd.RB,zlim=c(0,2),col=color.high.low,main="")
  maps::map(map,add=T,interior=F,fill=T)
  dev.off()

}
}




#plot examples
{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/examples"))
year=1830;month=1
X.plot.RB=X.anom.region.RB/100
Y.plot.RB=Y.anom.region.RB
levels=pretty(c(-20,20),n=20);levels=levels[which(levels!=0)]
lty=rep("solid",length(levels));lty[which(levels<0)]="dashed"
date.seq=as.Date(substring(text = names(X.plot.RB), first = 2, last = 11), "%Y.%m.%d")
years=as.numeric(format(date.seq, "%Y"))
months=as.numeric(format(date.seq, "%m"))
plot.idx=which(years==year & months==month)

png(paste0("case_study_",year,"-",month,".png"),height=4*fig.height,width=3*fig.width)
par(mfrow=c(4,3),mar=c(0,2,4,0),oma=c(3,2,0,1),cex.main=3,cex.axis=2)
for (i in 11:0) {
  cur.idx=plot.idx-i
  image(Y.plot.RB,cur.idx,col=color.warm.cold,main=NULL,zlim=c(-4,4),add=F,axes=T,labels=F)
  title(substring(date.seq[cur.idx],first=1,last=7),line=0.2)
  axis(1,labels=(if(i<3) labels=TRUE else labels=FALSE))
  axis(2,labels=(if(i%%3==2) labels=TRUE else labels=FALSE))
  contour(subset(X.plot.RB,cur.idx),levels=levels,lty=lty,add=T,lwd=2,labcex=1.5)
  maps::map(map,add=T,interior=F)}
dev.off()}




#plot case studies
setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/",region,"/animations"))
X.plot.RB=X.anom.region.RB/100
Y.plot.RB=Y.anom.region.RB
levels=pretty(c(-20,20),n=20);levels=levels[which(levels!=0)]
lty=rep("solid",length(levels));lty[which(levels<0)]="dashed"
  
for (i in 6000:12000){
  cur.idx=i
  png(paste0("SST_PSL_anom_",names(Y.plot.RB)[cur.idx],".png"),height=fig.height,width=fig.width)
  par(mfrow=c(1,1),mar=c(2,2,2,0.1),oma=c(0,0,0,0))
  plot(Y.plot.RB,cur.idx,col=color.warm.cold,zlim=c(-4,4),legend = T)
  contour(subset(X.plot.RB,cur.idx),levels=levels,lty=lty,add=T,lwd=2,labcex=1.5,)
  maps::map(map,add=T,interior=F)
  dev.off()
}



#plot basemaps
setwd(path)
X.data.Atlantic.RB = brick(paste0(X, "_Atlantic.nc")) 
X.data.Pacific.RB = brick(paste0(X, "_Pacific.nc")) 

Atlantic.extent = define.region.extent("Atlantic",NULL)
Pacific.extent = define.region.extent("Pacific",NULL)

X.Atlantic.region.RB = crop(X.data.Atlantic.RB, Atlantic.extent+40)
X.Pacific.region.RB = crop(X.data.Pacific.RB, Pacific.extent+40)

X.spatial.raster.atlantic=raster(X.Atlantic.region.RB,1)
X.spatial.raster.pacific=raster(X.Pacific.region.RB,1)

{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
png(paste0("basemap_atlantic.png"),height=350)
values(X.spatial.raster.atlantic)=NA
par(mar=c(2,2,2,2),oma=c(0,0,0,0))
plot(X.spatial.raster.atlantic,useRaster=T)
maps::map("world",add=T)
plot(Atlantic.extent+40,add=T,col=4,lwd=1.5)
plot(Atlantic.extent,add=T,col=2,lwd=1.5)
dev.off()}

{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
png(paste0("basemap_pacific.png"),height=350)
values(X.spatial.raster.pacific)=NA
par(mar=c(2,2,2,2),oma=c(0,0,0,0))
plot(X.spatial.raster.pacific)
maps::map("world2",add=T)
plot(Pacific.extent+40,add=T,col=4,lwd=1.5)
plot(Pacific.extent,add=T,col=2,lwd=1.5)
dev.off()}


{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
  png(paste0("basemap_atlantic_subregions.png"),height=350,width=400)
  values(X.spatial.raster.atlantic)=NA
  par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  plot(crop(X.spatial.raster.atlantic,define.region.extent("Atlantic")+20),useRaster=T)
  maps::map("world",add=T,interior=F)
  plot(define.region.extent("Atlantic","Sargasso"),col=2,add=T,lwd=1.5)
  plot(define.region.extent("Atlantic","Gulfstream"),col=3,add=T,lwd=1.5)
  plot(define.region.extent("Atlantic","North"),col=4,add=T,lwd=1.5)
  plot(define.region.extent("Atlantic","Spain"),col="orange",add=T,lwd=1.5)
dev.off()}

{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
  png(paste0("basemap_pacific_subregions.png"),height=350,width=500)
  values(X.spatial.raster.pacific)=NA
  par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  plot(crop(X.spatial.raster.pacific,define.region.extent("Pacific")+20),useRaster=T)
  maps::map("world2",add=T,interior=F)
  plot(define.region.extent("Pacific","West"),col=2,add=T,lwd=1.5)
  plot(define.region.extent("Pacific","East"),col=4,add=T,lwd=1.5)
  plot(define.region.extent("Pacific","NorthNP"),col=3,add=T,lwd=1.5)
  plot(define.region.extent("Pacific","Central"),col="orange",add=T,lwd=1.5)
  
dev.off()}

{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
  png(paste0("basemap_atlantic.png"),height=350,width=350)
  values(X.spatial.raster.atlantic)=NA
  par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  plot(define.region.extent("Atlantic")+40,col=2,lwd=1.5)
  plot(define.region.extent("Atlantic"),col=4,add=T,lwd=1.5)
  maps::map("world",add=T,interior=F)
  dev.off()}

{setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/"))
  png(paste0("basemap_pacific.png"),height=350,width=450)
  values(X.spatial.raster.pacific)=NA
  par(mar=c(2,2,0,0),oma=c(0,0,0,0))
  plot(define.region.extent("Pacific")+40,col=2,lwd=1.5)
  plot(define.region.extent("Pacific"),col=4,add=T,lwd=1.5)  
  maps::map("world2",add=T,interior=F)
  dev.off()}
