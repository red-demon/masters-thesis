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
  X.data.RB = brick(paste0(X, "_",region,"_", i, "_anom.nc")) 
  Y.data.RB = brick(paste0(Y,"_",region,"_", i, "_anom.nc"))
  
  # ------------------------------------------------------------------------------------------ 
  # 2. Define parameters and regions, crop data
  # ------------------------------------------------------------------------------------------
  #choose period for the prediction and training
  train.idx=seq(1,1000) #indices for training step
  pred.idx=seq(1001,2000) #indices for prediction
  
  #crop data to smaller region
  if (region=="Atlantic"){
    Y.extent = extent(c(-70, 0, 0, 70))
    X.extent = Y.extent + 40
  } else if (region=="Pacific") {
    Y.extent = extent(c(-140,-110, 0, 70))
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
  
  #-------------------------------------------------------------------------------------
  # 2. SVD
  #-----------------------------------------------------------------------------
  #calculate SVD
  #adjust area size
  area.vec=values(raster::area(X.region.RB))
  X.region.area.weighted = values(X.region.RB) * area.vec^2
  
  #calculate SVD
  X.svd=svd(X.region.area.weighted) #finding SVD of the matrix
  X.svd.u=X.svd$u #spatial EOFs of the initial matrix
  X.svd.v=X.svd$v # PCs or time series of the EOFs
  X.svd.d=X.svd$d # variance of each eof
  
  # standardize variance
  X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 
  
  
  # ------------------------------------------------------------------------------------------ 
  # 3. Regression 
  # ------------------------------------------------------------------------------------------
  #regress temp on SVD PCs
  predictor=cbind(X.svd.v[train.idx,1],X.svd.v[train.idx,2],X.svd.v[train.idx,3],X.svd.v[train.idx,4],X.svd.v[train.idx,5])
  model.eof=lm(Y.mean[train.idx]~predictor)
  
  #test predicted data on subset of data
  Y.pred= X.svd.v[pred.idx,1:5] %*% model.eof$coefficients[2:6]
  Y.pred= model.eof$coefficients[1]+Y.pred
  
  #plotting
  figname=paste("pc_regression",Y,region,i, sep="_")
  jpeg(paste0(figpath,figname,".jpg"), width = 600, height = 350)
  plot(Y.pred,Y.mean[pred.idx],main="Observed vs. predicted temperature anomalies from PC regression",
       xlab="Predicted values from PC regression [K]", ylab="Observed Mean Temp. anomalies [K]")
  model.fit=lm(Y.mean[pred.idx]~Y.pred)
  abline(model.fit$coefficients,col="red")
  correlation=cor(Y.pred,Y.mean[pred.idx])
  text(mean(Y.pred), min(Y.mean), labels = paste("cor=",correlation),col="red")
  dev.off()
  
  # ------------------------------------------------------------------------------------------ 
  # 4. How does correlation change with increasing EOF as predictors? 
  # ------------------------------------------------------------------------------------------
  #how many eofs to test?
  rm(predictor)
  rm(model.eof)
  eofs=1000
  correlation.numbereofs=rep(NA,eofs)
  predictor=X.svd.v[train.idx,1]
  
  
  #first iteration regression
  model.eof=lm(Y.mean[train.idx]~predictor)
  #test predicted data on subset of data
  Y.pred= X.svd.v[pred.idx,1] * model.eof$coefficients[2]
  Y.pred= model.eof$coefficients[1]+Y.pred
  correlation.numbereofs[1]=cor(Y.pred,Y.mean[pred.idx])
  
  #regress temp on SVD PCs
  for (j in seq(2,eofs)){
    print("eofs:") 
    print(j)
    predictor=cbind(predictor,X.svd.v[train.idx,j]) #add new eof to predictor-array
    model.eof=lm(Y.mean[train.idx]~predictor) #regression based on predictor-array
    
    #test predicted data on subset of data
    Y.pred= X.svd.v[pred.idx,1:j] %*% model.eof$coefficients[2:(j+1)] # multiply subset of eofs with coefficients to obtain predictions
    Y.pred= model.eof$coefficients[1]+Y.pred #add intercept to predictions
    correlation.numbereofs[j]=cor(Y.pred,Y.mean[pred.idx]) # calculate correlation and store in correlation-array
  }
  
  #plotting
  #plot correlation in dependence of EOFs
  figname=paste("corr_eof",Y,region,i,sep="_")
  jpeg(paste0(figpath,figname,".jpg"), width = 600, height = 350)
  plot(seq(1,eofs),correlation.numbereofs[1:eofs], type='l', main="Correlation between predicted and observed Temperature for different #EOFS",
       xlab="# EOFs", ylab="Correlation")  
  dev.off()
  
}
  