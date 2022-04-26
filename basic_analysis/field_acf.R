#load packages and sources
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)

source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/additional_functions/filled_contour3.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/colormaps.R")

# ------------------------------------------------------------------------------------------ 
# 00. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")

for (region in regions){
  
subregion=NULL

if(region=="Pacific") {map="world2";fig.height=350;fig.width=510}
if(region=="Atlantic") {map="world";fig.height=350;fig.width=450}


#set directories
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/",region)
setwd(path)

# ------------------------------------------------------------------------------------------ 
# 00. Read processed data and SVDs
# ------------------------------------------------------------------------------------------

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))


Y.extent=define.region.extent(region,NULL)
X.extent=Y.extent+40

X.region.RB = crop(X.data.RB, Y.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)


# ------------------------------------------------------------------------------------------ 
# 01. Calculate SVD
# ------------------------------------------------------------------------------------------

X.region.svd=values(area.weight.RB(X.region.RB,svd=TRUE))


#calculate total SVDs and store in separate objects
#X.svd=svd(X.region.svd)
#X.svd.u=X.svd$u
#X.svd.v=X.svd$v
#X.svd.d=X.svd$d
#X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 

#calculate SVDs for each month separately
{X.svd.mon=list()
  X.svd.mon.u=list()
  X.svd.mon.v=list()
  X.svd.mon.d=list()
  for (mon in 1:12){
    month=month.abb[mon]
    X.svd.mon.temp=svd(X.region.svd[,seq(mon,length(X.region.svd[1,]),12)])
    X.svd.mon.temp$d=X.svd.mon.temp$d^2/sum(X.svd.mon.temp$d^2)
    X.svd.mon[[month]] <- X.svd.mon.temp
    X.svd.mon.v[[month]]<-X.svd.mon.temp$v
    X.svd.mon.u[[month]]<-X.svd.mon.temp$u
    X.svd.mon.d[[month]]<-X.svd.mon.temp$d
    
    #to account for the arbitrary sign change in certain months in EOF 1 and 2:
    if(region=="Atlantic" & mon %in% c(3,12)) {X.svd.mon.v[[month]][,1]<- -X.svd.mon.temp$v[,1] ; X.svd.mon.u[[month]][,1]<- -X.svd.mon.temp$u[,1]}
    if(region=="Atlantic" & mon %in% c(1,2,3,7,12))  {X.svd.mon.v[[month]][,2]<- -X.svd.mon.temp$v[,2] ; X.svd.mon.u[[month]][,2]<- -X.svd.mon.temp$u[,2]}
    if(region=="Pacific" & !(mon %in% c(2,7,8))) {X.svd.mon.v[[month]][,1]<- -X.svd.mon.temp$v[,1] ;X.svd.mon.u[[month]][,1]<- -X.svd.mon.temp$u[,1]}    }
    if(region=="Pacific" & mon %in% c(11,3,8,11,12)){X.svd.mon.v[[month]][,2]<- -X.svd.mon.temp$v[,2] ; X.svd.mon.u[[month]][,2]<- -X.svd.mon.temp$u[,2]}
  
    
  }


#store all monthly PC ts in one array
X.svd.mon.all.v=matrix(NA,nrow=dim(X.svd.v)[1],ncol=dim(X.svd.mon.v[[1]])[2])
for (i in 1:12){
  idx=seq(i,length(X.svd$v[,1]),12)
  X.svd.mon.all.v[idx,]=X.svd.mon.v[[i]]
}

# ------------------------------------------------------------------------------------------ 
# 02. Calculate ACF and CCF
# ------------------------------------------------------------------------------------------

gridcell.acf.jan=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=Y.region.RB,month=1,nr.lags=24)
gridcell.acf.jul=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=Y.region.RB,month=7,nr.lags=24)

#gridcell cross-correlation with EOF 1
X.PC1=X.svd.mon.all.v[,1]
gridcell.ccf.jan=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=X.PC1,month=1,nr.lags=24)
gridcell.ccf.jul=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=X.PC1,month=7,nr.lags=24)

X.PC2=X.svd.mon.all.v[,2]
gridcell.ccf2.jan=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=X.PC2,month=1,nr.lags=24)
gridcell.ccf2.jul=ccf.monthly.gridcell(Y.RB=Y.region.RB,X.RB=X.PC2,month=7,nr.lags=24)

setwd(figpath)

for (i in 1:25){
  png(paste0("gridcell_acf_",region,"_Jan_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.acf.jan,i,col=color.high.low.negative,zlim=c(-1,1),main="")
  maps::map(map,add=T,interior=F,fill=T)
  box("inner")
  dev.off()
}

for (i in 1:25){
  png(paste0("gridcell_acf_",region,"_Jul_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.acf.jul,i,col=color.high.low.negative,zlim=c(-1,1),main="")
  maps::map(map,add=T,interior=F,fill=T)
  box("inner")
  dev.off()
}

for (i in 1:25){
  png(paste0("gridcell_ccf_pc1_",region,"_Jan_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.ccf.jan,i,col=color.high.low.negative,zlim=c(-1,1),main="")
  maps::map(map,add=T,interior=F,fill=T)
  box("inner")
  dev.off()
}

for (i in 1:25){
  png(paste0("gridcell_ccf_pc1_",region,"_Jul_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.ccf.jul,i,col=color.high.low.negative,zlim=c(-1,1),main="")
  maps::map(map,add=T,interior=F,fill=T)
  box("inner")
  dev.off()
}

for (i in 1:25){
  png(paste0("gridcell_ccf_pc2_",region,"_Jan_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.ccf2.jan,i,col=color.high.low.negative,zlim=c(-1,1),main="",axes=F)
  box("inner")
  maps::map(map,add=T,interior=F,fill=T)
  dev.off()
}

for (i in 1:25){
  png(paste0("gridcell_ccf_pc2_",region,"_Jul_lag",i-1,".png"),height=fig.height,width=fig.width)
  par(mar=c(0,0,0,0),oma=c(0,0,0,0))
  image(gridcell.ccf2.jul,i,col=color.high.low.negative,zlim=c(-1,1),main="")
  maps::map(map,add=T,interior=F,fill=T)
  box("inner")
  dev.off()
}
}
  
