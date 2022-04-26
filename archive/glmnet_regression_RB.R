# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment PC Regression simple script:
# Dynamically adjust regional SST or TREFHT mean temperatures using PSL as predictor.
# Data are read in using the Raster Package
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 22.7.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(fields)
library(glmnet)
# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
region="Atlantic"
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/",region)
figpath="/net/h2o/climphys1/rdaenzer/figures/"
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST
season=c("DJF","MAM","JJA","SON")

for (i in season){
  #read ncdf files into RB objects
  X.data.RB = brick(paste(X,region,i,"anom.nc",sep="_")) 
  Y.data.RB = brick(paste(Y,region,i,"anom.nc",sep="_"))
  
  # ------------------------------------------------------------------------------------------ 
  # 2. Define parameters and regions, crop data
  # ------------------------------------------------------------------------------------------

  #crop data to smaller region
  if (region=="Atlantic"){
    Y.extent = extent(c(-70, 0, 0, 70))
    X.extent = Y.extent + 40
  } else if (region=="Pacific") {
    Y.extent = extent(c(140,240, 0, 70))
    X.extent = Y.extent + 40
  }
  
  X.region.RB = crop(X.data.RB, X.extent)
  Y.region.RB = crop(Y.data.RB, Y.extent)
  X.region=values(X.region.RB)
  Y.region=values(Y.region.RB)
  
  #extract coordinatess
  Y.coord=coordinates(Y.region.RB)
  Y.lon=seq(min(Y.coord[,1]),max(Y.coord[,1]),2)
  Y.lat=seq(min(Y.coord[,2]),max(Y.coord[,2]),2)
  X.coord=coordinates(X.region.RB)
  X.lon=seq(min(X.coord[,1]),max(X.coord[,1]),2)
  X.lat=seq(min(X.coord[,2]),max(X.coord[,2]),2)
  
  #dimension of regional domain incl. time
  X.region.dim=dim(X.region.RB)
  Y.region.dim=dim(Y.region.RB)

  
  #calculate region mean
  X.mean=colMeans(X.region,na.rm=TRUE)
  Y.mean=colMeans(Y.region,na.rm=TRUE)
  
  # ------------------------------------------------------------------------------------------ 
  # 3.Prediction using glmnet
  # ------------------------------------------------------------------------------------------
  
  #choose period for the prediction and training
  train.idx=seq(1,1000) #indices for training step
  pred.idx=seq(1001,1999) #indices for prediction
  alpha=0.1
  
  model.cv.glmnet=cv.glmnet(t(X.region[,train.idx]),Y.mean[train.idx], alpha=alpha)
  Y.hat=predict(model.cv.glmnet,t(X.region[,pred.idx]))
  
  #plotting
  figname=paste("glmnet_regression",Y,region,i, sep="_")
  jpeg(paste0(figpath,figname,".jpg"), width = 600, height = 350)
  plot(Y.mean[pred.idx],Y.hat,main=paste("Predicted (cv.glmnet) vs. observed mean",Y,"anomalies,",region,i,sep=" "),
       xlab="Mean SST anomalies, observed [K]", ylab="Mean SST anomalies, predicted [K]")
  correlation=cor(Y.hat,Y.mean[pred.idx])
  text(0.2, -0.3, labels = paste("cor=",correlation),col="red")
  dev.off()
}
  
