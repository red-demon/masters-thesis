# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment GLMNET Regression simple script:
# Dynamically adjust regional SST or TREFHT mean temperatures using lagged PSL data as predictor.
# Data are read in using the Raster Package
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 12.8.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(glmnet)


# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/time_lag.R")
# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
regions=c("Pacific","Atlantic")
period.type="seasonal" #seasonal or monthly

for (region in regions) {
  print(region)
  path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/",region)
  figpath="/net/h2o/climphys1/rdaenzer/figures"
  outpath="/net/h2o/climphys1/rdaenzer/output"
  setwd(path)
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  
  #read ncdf files into RB objects
  X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
  Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
  
  # ------------------------------------------------------------------------------------------ 
  # 2. Define regions, crop data
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
  
  # -----------------------------------------------------------------------------------------
  # 3. Extract lagged monthly or seasonal data
  #------------------------------------------------------------------------------------------
  if (period.type=="monthly") period.names=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  if (period.type=="seasonal") period.names=c("djf","mam","jja","son")
  if (period.type=="monthly")  lag.names=seq(1,12)
  if (period.type=="seasonal") lag.names=seq(1,4)
  
  correlation=matrix(NA,nrow=length(lag.names),ncol=length(period.names),dimnames=list(lag.names,period.names))
  
  for (i in seq(1,length(period.names))){
    
    if (period.type=="monthly") Y.period=i
    if (period.type=="seasonal") Y.period=period.names[i]
    
    print(Y.period)
    
    for (j in lag.names){
      print(j)
      #iterate trough lag from month 1 to 12
      lag=j
      
      X.lag.RB=extract.Xlag.period.RB(X.region.RB,period=period.type,Y.period=Y.period,lag=lag)
      Y.lag.RB=extract.Ylag.period.RB(Y.region.RB,period=period.type,Y.period=Y.period,lag=lag)
      
      Y.lag=values(Y.lag.RB)
      
      #calculate region mean
      Y.mean=colMeans(Y.lag,na.rm=TRUE)
      
      #-------------------------------------------------------------------------------------
      # 2. SVD
      #-----------------------------------------------------------------------------
      #choose period for the prediction and training
      train.idx=seq(1,1000) #indices for training step
      pred.idx=seq(1001,1998) #indices for prediction
      alpha=1
      
      #adjust area size
      area.vec=values(raster::area(X.lag.RB))
      X.lag.areaweighted = values(X.lag.RB) * area.vec^2
      
      #svd regression
      model.cv.glmnet=cv.glmnet(t(X.lag.areaweighted[,train.idx]),Y.mean[train.idx], alpha=alpha)
      Y.hat=predict(model.cv.glmnet,t(X.lag.areaweighted[,pred.idx]))
      
      model.linear=lm(Y.hat~Y.mean[pred.idx])
      correlation[j,i]=cor(Y.mean[pred.idx],Y.pred)
    }
    
  }
  
  # save correlation as RDS file
  filename=paste("glmnet_regr_lagcor",Y,region,period.type,"anom.RDS",sep="_")
  file=file.path(outpath,filename)
  saveRDS(correlation,file=file)
}

