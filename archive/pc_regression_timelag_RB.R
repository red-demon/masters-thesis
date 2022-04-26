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
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")

# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
region="Atlantic"
nr.lags=24

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

#adjust area size
X.area.vec=values(raster::area(X.region.RB))
X.area.weighted = values(X.region.RB) * X.area.vec^2

Y.area.vec=values(raster::area(Y.region.RB))
Y.area.weighted=values(Y.region.RB) * Y.area.vec^2

#calculate Y regional mean
Y.mean=colMeans(Y.area.weighted,na.rm=TRUE)
# -----------------------------------------------------------------------------------------
# 3. Lagged SVD Regression
#------------------------------------------------------------------------------------------
period.names=month.abb
lags=seq(0,nr.lags)

#create empty matrix for storing lagged correlations with dimensions lags x periods
correlation=matrix(NA,nrow=length(lags),ncol=length(period.names),dimnames=list(lags,period.names))

for (i in 1:length(period.names)){

  print(period.names[i])
  
  start.idx=i+nr.lags
  Y.idx=seq(start.idx,length(Y.mean),12)
  Y.mean.mon=Y.mean[Y.idx]
  
  for (lag in seq(0,nr.lags)){
    print(lag)
    #iterate trough lag from month 1 to 12
    
    X.idx=Y.idx-lag
    X.lag.area.weighted=X.area.weighted[,X.idx]
    
    #choose period for the prediction and training
    train.years=seq(1,1000) #indices for training step
    pred.years=seq(1001,length(Y.mean.mon)) #indices for prediction
    
    #svd regression
    X.svd=svd(X.lag.area.weighted)
    X.svd.v=X.svd$v
    Y.pred=svd.regression(Y.mean.mon,X.svd.v,predictors=5,train.years=train.years,pred.years=pred.years)
    
    #calculate correlation betw. Y.pred and Y.mean for prediction Period
    correlation[lag+1,i]=cor(Y.mean.mon[pred.years],Y.pred)
    
  }
  
}






#-------------------------------------------------------------------------------------
# 4. Cumulative Lagged SVD regression
#-----------------------------------------------------------------------------
lags=seq(0,nr.lags)
months=c(1,2,3,4,5,6)
train.years=1:1000
pred.years=1001:1998
predictors=5

correlation.lag=array(NA,12)
for (i in 1:12){
  res=svd.regression.lag.cumulative(Y.ts=Y.mean,X.ts=X.area.weighted,predictors=predictors,
        train.years=train.years,pred.years=pred.years,month.idx=i,months=lags,type="lag")
  Y.obs=res$Y.obs
  Y.hat=res$Y.hat
  correlation.lag[i]=cor(Y.obs[pred.years],Y.hat[pred.years])
}


correlation.months=array(NA,12)
for (i in 1:12){
  res=svd.regression.lag.cumulative(Y.ts=Y.mean,X.ts=X.area.weighted,predictors=predictors,
                                    train.years=train.years,pred.years=pred.years,month.idx=i,months=months,type="months")
  Y.obs=res$Y.obs
  Y.hat=res$Y.hat
  correlation.months[i]=cor(Y.obs[pred.years],Y.hat[pred.years])
}



#-------------------------------------------------------------------------------------
# 5. Plotting and save data
#-------------------------------------------------------------------------------------
setwd(outpath)
outname=paste0("pc_regr_lag_corr",Y,"_",region,"_",period.type,".RDS")
saveRDS(correlation,outname)

#plot monthly pc regression lag correlation
setwd(figpath)
{figname=paste0("02_pc_regr_lag_corr_",paste(Y,region,period.type,sep="_"),".pdf")
  pdf(figname, width = 6*4, height = 5.5*3)
  par(mfrow=c(4, 3))
  for (i in 1:length(period.names)){
    plot(lags,correlation[,i],main=paste("PC regression lag correlations (pred vs. obs) mean SST anomalies,",period.type,period.names[i],sep=" "), 
         ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag", ylab="linear correlation",type="l")}
  dev.off()
}


#plot cumulative lag correlation
{pdf(paste0("02_pc_regr_cumulative_lag_corr_",region,".pdf"), width = 6*4, height = 5.5*3)
  par(mfrow=c(2, 1))
  plot(correlation.lag,type="l",ylim=c(0,1),xaxt="n", main="PC Regr. with cumulative lagged (24m) time series from 5 leading PSL EOFs",
       ylab="correlation of obs. and pred. mean SST",xlab=NA)
  axis(side=1,at=c(1:12),labels=month.abb)
  
  plot(correlation.months,type="l",ylim=c(0,1),xaxt="n", main="PC Regr. with fixed (m:1-6) time series from 5 leading PSL EOFs",
       ylab="correlation of obs. and pred. mean SST",xlab=NA)
  axis(side=1,at=c(1:12),labels=month.abb)
  dev.off()}
  