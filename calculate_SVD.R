#load packages and sources
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")

# ------------------------------------------------------------------------------------------ 
# 00. Read model data
# ------------------------------------------------------------------------------------------
for (region in c("Atlantic","Pacific")){

print(region)
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
outpath="/net/h2o/climphys1/rdaenzer/output"
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))

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

X.region.RB = crop(X.data.RB, X.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)

#area weighting and centering

X.region.aw.RB=area.weight.RB(X.region.RB)
Y.region.aw.RB=area.weight.RB(Y.region.RB)

Y.mean.aw=mean.area.weighted.RB(Y.region.RB)

#---------------------------------------------------------------------------------
# 02. SVD total and monthly
#---------------------------------------------------------------------------------

#calculate total SVDs and store in separate objects
X.region.svd=values(area.weight.RB(X.region.RB,svd=TRUE))
X.svd=svd(X.region.svd)
X.svd.u=X.svd$u
X.svd.v=X.svd$v
X.svd.d=X.svd$d
X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 

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
  }
}

#store all monthly PC ts in one array
X.svd.mon.all.v=matrix(NA,nrow=dim(X.svd.v)[1],ncol=dim(X.svd.mon.v[[1]])[2])
for (i in 1:12){
  idx=seq(i,length(X.svd$v[,1]),12)
  X.svd.mon.all.v[idx,]=X.svd.mon.v[[i]]
}





#------------------------------------------------------------------------------------------
#03. save SVDs as R files
#------------------------------------------------------------------------------------------
setwd(outpath)

saveRDS(X.svd, file = paste0("X_svd_",region,".rds"))
saveRDS(X.svd.mon, file = paste0("X_svd_mon_",region,".rds"))
saveRDS(X.svd.mon.all.v,file=paste0("X_svd_mon_all_v_",region,".rds"))

saveRDS(X.region.RB,file=paste0("X_region_RB_",region,".rds"))
saveRDS(Y.region.RB,file=paste0("Y_region_RB_",region,".rds"))

saveRDS(X.region.aw.RB,file=paste0("X_region_aw_RB_",region,".rds"))
saveRDS(Y.region.aw.RB,file=paste0("Y_region_aw_RB_",region,".rds"))

saveRDS(Y.mean.aw,file=paste0("Y_mean_aw_",region,".rds"))
}

#------------------------------------------------------------------------------------------
#03. read SVDs and data as R files
#------------------------------------------------------------------------------------------
setwd(outpath)
X.svd=readRDS(file = paste0("X_svd_",region,".rds"))
X.svd.mon=readRDS(file = paste0("X_svd_mon_",region,".rds"))
X.svd.mon.all.v=readRDS(file=paste0("X_svd_mon_all_v_",region,".rds"))

X.region.RB=readRDS(file=paste0("X_region_RB_",region,".rds"))
Y.region.RB=readRDS(file=paste0("Y_region_RB_",region,".rds"))

X.region.aw.RB=readRDS(file=paste0("X_region_aw_RB_",region,".rds"))
Y.region.aw.RB=readRDS(file=paste0("Y_region_aw_RB_",region,".rds"))

X.mean.aw=readRDS(file=paste0("X_mean_aw_",region,".rds"))
Y.mean.aw=readRDS(file=paste0("Y_mean_aw_",region,".rds"))
