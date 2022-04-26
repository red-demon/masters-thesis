#script for calculating and plotting monthly EOFs


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
  
  
  #set directories
  path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/basic_analysis/",region)
  outpath="/net/h2o/climphys1/rdaenzer/output/"
  setwd(path)
  
  # ------------------------------------------------------------------------------------------ 
  # 00. Read processed data
  # ------------------------------------------------------------------------------------------
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  
  #read ncdf files into RB objects
  X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
  Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
  
    
  Y.extent=define.region.extent(region,NULL)
  X.extent=Y.extent+40
  
  X.region.RB = crop(X.data.RB, X.extent)
  Y.region.RB = crop(Y.data.RB, Y.extent)
  
  #area weighting and centering
  
  X.region.aw.RB=area.weight.RB(X.region.RB)
  Y.region.aw.RB=area.weight.RB(Y.region.RB)
  
  Y.mean.aw=mean.area.weighted.RB(Y.region.RB)
  
  
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
      
      #to account for the arbitrary sign change in certain months:
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
  

  
  #---------------------------------------------------------------------------------
  # 03. Plotting leading EOFs and corresponding PCs
  #---------------------------------------------------------------------------------
  setwd(figpath)
  

  
  #worldmaps
  if(region=="Pacific") {map="world2";fig.height=350;fig.width=450}
  if(region=="Atlantic") {map="world";fig.height=350;fig.width=360}
  
  
  #plot total EOFs 1-5:
  X.spatial.raster=raster(X.region.RB,1)
  for (i in 1:5){
    figname=paste0("EOF_",paste(i,X,region,sep="_"),"_Jan.png")
    png(figname,height=fig.height,width=fig.width)
    par(mar=c(2,2,2,0.5),oma=c(0,0,0,0))
    values(X.spatial.raster)=X.svd.mon.u[[1]][,i]
    image(X.spatial.raster,col=color.warm.cold,zlim=c(-0.1,0.1),main=paste0(i,". EOF"),cex.main=1.1)
    maps::map(map,add=T,interior=F)
    dev.off()}
  
  for (i in 1:5){
    figname=paste0("EOF_",paste(i,X,region,sep="_"),"_Jul.png")
    png(figname,height=fig.height,width=fig.width)
    par(mar=c(2,2,2,0.5),oma=c(0,0,0,0))
    values(X.spatial.raster)=X.svd.mon.u[[7]][,i]
    image(X.spatial.raster,col=color.warm.cold,zlim=c(-0.1,0.1),main=paste0(i,". EOF"),cex.main=1.1)
    maps::map(map,add=T,interior=F)
    dev.off()}
  
  #save variance explained
  var.expl=matrix(NA,nrow=2,ncol=5,dimnames=list(c("Jan","Jul"),c("EOF1","EOF2","EOF3","EOF4","EOF5")))
  var.expl[1,]=X.svd.mon.d[[1]][1:5]
  var.expl[2,]=X.svd.mon.d[[7]][1:5]
  
  write.csv(var.expl,paste0("EOF_var_explained_",region,".csv"))
}
  
  