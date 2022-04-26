# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment PC Regression simple script:
# Dynamically adjust regional SST or TREFHT mean temperatures using PSL as predictor.
# Data are read in using the Raster Package
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 29.7.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")

# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
regions=c("Pacific","Atlantic")
period.type="monthly" #seasonal or monthly

for (region in regions) {
  print(region)
  path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
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
    Y.extent = extent(c(-90, 20,20, 70))
    X.extent = Y.extent + 40
  } else if (region=="Pacific") {
    Y.extent = extent(c(140,240, 20, 70))
    X.extent = Y.extent + 40
  }
  
  X.region.RB = crop(X.data.RB, X.extent)
  Y.region.RB = crop(Y.data.RB, Y.extent)
  
  # -----------------------------------------------------------------------------------------
  # 3. Extract lagged monthly or seasonal data
  #------------------------------------------------------------------------------------------
  if (period.type=="monthly") period.names=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  if (period.type=="seasonal") period.names=c("djf","mam","jja","son")
  if (period.type=="monthly")  lags=seq(0,12)
  if (period.type=="seasonal") lags=seq(0,4)
  
  correlation=matrix(NA,nrow=length(lags),ncol=length(period.names),dimnames=list(lags,period.names))
  
  for (i in 1:length(period.names)){

    if (period.type=="monthly") Y.period=i
    if (period.type=="seasonal") Y.period=period.names[i]
    
    print(Y.period)
    
    for (lag in lags){
      print(lag)
      #iterate trough lag from month 1 to 12
  
      X.lag.RB=extract.Xlag.period.RB(X.region.RB,period=period.type,Y.period=Y.period,lag=lag)
      Y.lag.RB=extract.Ylag.period.RB(Y.region.RB,period=period.type,Y.period=Y.period,lag=lag)
      
      #adjust area size
      X.area.vec=values(raster::area(X.lag.RB))
      X.lag.areaweighted = values(X.lag.RB) * X.area.vec^2
      
      Y.area.vec=values(raster::area(Y.lag.RB))
      Y.lag.areaweighted=values(Y.lag.RB) * Y.area.vec^2
      
      #calculate region mean
      Y.mean=colMeans(Y.lag.areaweighted,na.rm=TRUE)
      #-------------------------------------------------------------------------------------
      # 2. SVD
      #-----------------------------------------------------------------------------
      #choose period for the prediction and training
      train.idx=seq(2,1000) #indices for training step
      pred.idx=seq(1001,length(Y.mean)) #indices for prediction
      
      #svd
      X.svd=svd(X.lag.areaweighted)
      X.svd.v=X.svd$v
      
      #svd regression
      Y.pred=svd.regression(Y.mean,X.svd.v,predictors=5,train.years=train.idx,pred.years=pred.idx)
      
      #calculate correlation betw. Y.pred and Y.mean for prediction Period
      correlation[lag+1,i]=cor(Y.mean[pred.idx],Y.pred)
    }
    
  }
  
  #-------------------------------------------------------------------------------------
  # 4. Plotting and save data
  #-----------------------------------------------------------------------------
  setwd(outpath)
  outname=paste0("pc_regr_lag_corr",Y,"_",region,"_",period.type,".RDS")
  saveRDS(correlation,outname)
  
  setwd(figpath)
  {figname=paste0("02_pc_regr_lag_corr_",paste(Y,region,period.type,sep="_"),".pdf")
    pdf(figname, width = 6*4, height = 5.5*3)
    par(mfrow=c(4, 3))
    for (i in 1:length(period.names)){
      plot(lags,correlation1[,i],main=paste("PC regression lag correlations (pred vs. obs) mean SST anomalies,",period.type,period.names[i],sep=" "), 
           ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag", ylab="linear correlation",type="l")}
    dev.off()
  }
  
}

