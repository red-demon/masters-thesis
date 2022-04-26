#-------------------------------------------------------------------------------
# Calculate cross correlations and auto correlations of mean SST and PSL EOFs
#-------------------------------------------------------------------------------

#load packages and sources
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)

source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
region="Atlantic"
type="ccf" #acf or ccf

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
  Y.extent = extent(c(-90, 20, 20, 70))
  X.extent = Y.extent + 40
} else if (region=="Pacific") {
  Y.extent = extent(c(140,240, 20, 70))
  X.extent = Y.extent + 40
}

X.region.RB = crop(X.data.RB, X.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)
#----------------------------------------------------------------------------------------
# 3. calculate SVDs
# ------------------------------------------------------------------------------------
#adjust area size
X.centered.RB=center.data.RB(X.region.RB)
Y.centered.RB=center.data.RB(Y.region.RB)
X.area.weighted.RB=area.weighting.RB(X.centered.RB)
Y.area.weighted.RB=area.weighting.RB(Y.centered.RB)

Y.mean.area.weighted=colMeans(values(Y.area.weighted.RB),na.rm=TRUE)

X.svd.v=svd(values(X.region.area.weighted.RB))$v
#-----------------------------------------------------------------------------------
# 4. calculate ccf and acf
#----------------------------------------------------------------------------------

#calculate lag correlation

X.ts=X.region.area.weighted
X.ts2=X.svd.v[,1]

#calculate lag correlation from SVDs calculated monthly
lag.cor = sapply(X = 1:12, FUN = function(mon.ix) ccf.svd.monthly(Y.ts = Y.mean,X.ts=X.ts, mon.ix = mon.ix, nr.lags = 48,eof=1))

#calculate lag correlation from SVDs calculated over whole dataset
lag.cor2 = sapply(X = 1:12, FUN = function(mon.ix) ccf.monthly(Y.ts = Y.mean.area.weighted,X.ts=X.ts2, mon.ix = mon.ix, nr.lags = 48))



#plotting
setwd(figpath)

{pdf(paste0("01_ccf_",region,".pdf"), width = 6*4, height = 5.5*3)
par(mfrow=c(4, 3))
for (i in 1:12) {
  plot(x = 0:48, y = lag.cor[,i], 
       ylim = c(-1, 1), xlim = c(0, 50), type="l",xlab = "time lag [months]", ylab = "Cross Correlation", 
       main=paste0("CCF of Atlantic mean SST and 1. PSL EOF (EOFs calculated per month), ",month.abb[i]))
}
dev.off()}

{pdf(paste0("01_ccf_",region,"2.pdf"), width = 6*4, height = 5.5*3)
  par(mfrow=c(4, 3))
  for (i in 1:12) {
    plot(x = 0:48, y = lag.cor2[,i], 
         ylim = c(-1, 1), xlim = c(0, 50), type="l",xlab = "time lag [months]", ylab = "Cross Correlation", 
         main=paste0("CCF of Atlantic mean SST and 1. PSL EOF (EOFs calculated over all months), ",month.abb[i]))
  }
dev.off()}

{jpeg(paste0("01_ccf_",region,"_allmonths.pdf"), width = 6*4, height = 5.5*3)
plot(x = 0:48, y = lag.cor[,1], 
     ylim = c(-1, 1), xlim = c(0, 50), type="l",xlab = "time lag [months]", ylab = "Cross Correlation", 
     main=paste0("CCF of Atlantic mean SST and 1. PSL EOF (EOFs calculated over all months), ",month.abb[i]))
for (i in 2:12) {
  lines(x = 0:48, y = lag.cor[,i],col=i)} 
dev.off()}
