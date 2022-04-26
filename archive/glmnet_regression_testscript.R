
# ----------------------------------------------------------------------------------------------
# Testscript for testing glmnet_regression on model control data by using different predictors
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 29.9.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)
library(hydroGOF)
library(glmnet)
library(SpatialEpi)

library(foreach)
library(doParallel)
library(bigmemory)
registerDoParallel(cores = 4)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_00_preprocessing_4DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_00_preprocessing_SPACETIME_4DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_01_GRIDCELL_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_01_SPACE_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ_RB.R")

# other functions to read:
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/convert.to.eurocentric.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/frenchcolormap.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/project_raster.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/gridcorts.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/02_elasticnet_regression.R")

# ------------------------------------------------------------------------------------------ 
# 00. Read processed data and SVDs
# ------------------------------------------------------------------------------------------
region="Atlantic"
subregion=NULL
print(region)

#set directories
figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/elasticnet_regression/",region)
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
setwd(figpath)

#define subregion, crop data to regional extent
Y.extent=define.region.extent(region,subregion)
X.extent=Y.extent+40

if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)

X.region.RB = crop(X.data.RB, X.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)

#area weighting and centering
X.region.aw.RB=area.weight.RB(X.region.RB)
Y.region.aw.RB=area.weight.RB(Y.region.RB)

Y.mean.aw=mean.area.weighted.RB(Y.region.RB)


# ------------------------------------------------------------------------------------------ 
# 01. elasticnet regression for Y.mean
# ------------------------------------------------------------------------------------------


#elasticnet regression for Y.mean using first 100 EOFs as predictors
month=1
lags=0:12
train.years=1:1000
pred.years=1001:1998
s="lambda.min"
alpha=0.1

#extract and define Y.obs and predictor
Y.obs.cur.month=extract.lag.period.idx(Y.mean.aw,month,0)
predictor=t(values(X.region.aw.RB)) #define predictor for ridge regression (either t(X.region.aw) or X.svd.mon.all.v[,eofs])


#effect of alpha and predictor type on prediction
#------------------------------------------------------------------------------------------

#alpha
#correlation.alpha=sapply(X=seq(0,1,0.1),FUN=function(alpha) cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,X.svd.mon.all.v[,1:100],lags=lags,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
#plot(seq(0,1,0.1),correlation.alpha,type="l",ylim=c(0,1),main=paste0("effect of alpha on correlation, ",month.abb[month]))

#effect of predictors (PSL vs PSL EOFs) on prediction
#correlation.psl=cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=lags,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s))

# ------------------------------------------------------------------------------------------ 
#correlation between monthly mean SST and elast.net predicted SST for specific lags, only 2 months
#------------------------------------------------------------------------------------------
months=c(1,7)
correlation.lag=lapply(X=months,FUN=function(month){
    sapply(X=lags,FUN=function(lag) cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=lag,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
})

#plot
{figname=paste0("01_elasticnet_lag_cor_",region,".jpg")
  jpeg(figname)
  plot(lags,correlation.lag[[1]],main="ridge regression lag correlation", 
       ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="correlation",type="l")
  lines(lags,correlation.lag[[2]],col=2)
  legend("topright",month.abb[months],pch="-",col=1:2,bty="y",cex = 1,pt.cex=2)
  dev.off()}

# ------------------------------------------------------------------------------------------ 
#correlation between monthly mean SST and elast.net predicted SST for specific lags, all months
#------------------------------------------------------------------------------------------
correlation.lag=list()
for (i in 1:12){
  print(i)
  Y.obs.cur.month=extract.lag.period.idx(Y.mean.aw,i,0)
  correlation.lag[[month.abb[i]]]=sapply(X=lags,FUN=function(lag) cor(Y.obs.cur.month[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=lag,train.years = train.years,pred.years=pred.years,month=i,alpha=alpha,s=s)))
}

#plot monthly
{figname=paste0("01_elasticnet_lag_cor_",paste(region,sep="_"),".pdf")
  pdf(figname, width = 6*4, height = 5.5*3)
  par(mfrow=c(4, 3))
  for (i in 1:12){
    plot(lags,correlation.lag[[i]],main=month.abb[i], 
         ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="linear correlation",type="l")
  }
  dev.off()
}

#plot seasonally
season.abb=c("DJF","MAM","JJA","SON")
plot.months=c(12,1,2,3,4,5,6,7,8,9,10,11)

for (i in 1:4){
  season=season.abb[i]
  {figname=paste0("01_elasticnet_lag_cor_",region,"_",season,".jpg")
    {jpeg(figname)
    cur.months=plot.months[(1:3)+(i-1)*3]
    plot(lags,correlation.lag[[cur.months[1]]],ylim=c(0,1),type="l",col=1,main=season,ylab="Correlation",xlab="lag in months")
    for (j in 2:3) lines(lags,correlation.lag[[cur.months[j]]],col=j)
    legend("topright",month.abb[cur.months],pch="-",col=1:3,bty="y",cex = 1,pt.cex=2)
    dev.off()}
  }
}

# ------------------------------------------------------------------------------------------   
#correlation between monthly mean sst and elast.net predicted SST for cumulative lags, only 2 months
#------------------------------------------------------------------------------------------------

months=c(1,7)
correlation.cum.lag=lapply(X=months,FUN=function(month){
  sapply(X=lags,FUN=function(lag) cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=0:lag,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
})

#plot
{figname=paste0("01_elasticnet_cumulative_lag_cor_",region,".jpg")
  jpeg(figname)
  plot(lags,correlation.cum.lag[[1]],main="ridge regression lag correlation", 
       ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="correlation",type="l")
  lines(lags,correlation.cum.lag[[2]],col=2)
  legend("topright",month.abb[months],pch="-",col=1:2,bty="y",cex = 1,pt.cex=2)
  dev.off()}

# ------------------------------------------------------------------------------------------ 
#correlation between monthly mean sst and elast.net predicted SST for cumulative lags, all months
#------------------------------------------------------------------------------------------------
correlation.cum.lag=list()
for (i in 1:12){
  print(i)
  Y.obs.cur.month=extract.lag.period.idx(Y.mean.aw,i,0)
  correlation.cum.lag[[month.abb[i]]]<-sapply(X=lags,FUN=function(lag) cor(Y.obs.cur.month[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=0:lag,train.years = train.years,pred.years=pred.years,month=i,alpha=alpha,s=s)))
}

#plot monthly
{figname=paste0("01_elasticnet_cumulative_lag_cor_",region,".pdf")
  pdf(figname, width = 6*4, height = 5.5*3)
  par(mfrow=c(4, 3))
  for (i in 1:12){
    plot(lags,correlation.cum.lag[[i]],main=month.abb[i], 
         ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="linear correlation",type="l")
  }
  dev.off()
}

#plot seasonally
season.abb=c("DJF","MAM","JJA","SON")
plot.months=c(12,1,2,3,4,5,6,7,8,9,10,11)

for (i in 1:4){
  season=season.abb[i]
  figname=paste0("01_elasticnet_cumulative_lag_cor_",region,"_",season,".jpg")
  {jpeg(figname)
    cur.months=plot.months[(1:3)+(i-1)*3]
    plot(lags,correlation.cum.lag[[cur.months[1]]],ylim=c(0,1),type="l",col=1,main=season,ylab="Correlation",xlab="time lag[months]")
    for (j in 2:3) lines(lags,correlation.cum.lag[[cur.months[j]]],col=j)
    legend("topright",month.abb[cur.months],pch="-",col=1:3,bty="y",cex = 1,pt.cex=2)
    dev.off()}
}



#-----------------------------------------------------------------------------------------------
#correlation of elasticnet regression depending on region
# ------------------------------------------------------------------------------------------ 
if (region=="Atlantic") subregions=c("Sargasso","Gulfstream","Spain")
if (region=="Pacific") subregions=c("East","West")

month=7
lags=0:12
train.years=1:1000
pred.years=1001:1998
s="lambda.min"
alpha=0.1

#define subregion, crop data to regional extent
correlation.cum.lag.subregions=lapply(X=subregions,FUN=function(subregion){
  print(subregion)
  Y.extent=define.region.extent(region,subregion)
  X.extent=Y.extent+40
  X.cur.region.RB = crop(X.data.RB, X.extent)
  X.cur.region.aw.RB=area.weight.RB(crop(X.data.RB,X.extent))
  Y.cur.mean.aw=mean.area.weighted.RB(crop(Y.data.RB,Y.extent))
  
  cur.predictor=t(values(X.cur.region.aw.RB)) #define predictor for ridge regression (either t(X.region.aw) or X.svd.mon.all.v[,eofs])
  correlation.cum.lag=sapply(X=lags,FUN=function(lag) cor(extract.lag.period.idx(Y.cur.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.cur.mean.aw,cur.predictor,lags=0:lag,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
})

if (region=="Pacific") correlation.cum.lag.subregions.pacific=correlation.cum.lag.subregions
if (region=="Atlantic") correlation.cum.lag.subregions.atlantic=correlation.cum.lag.subregions

setwd(figpath)
#plot
{figname=paste0("01_elasticnet_cumulative_lag_cor_subregions_",region,".jpg")
  jpeg(figname)
  plot(lags,correlation.cum.lag.subregions[[1]],main="ridge regression lag correlation, Jul", 
       ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="correlation",type="l")
  for (i in 2:length(subregions)) lines(lags,correlation.cum.lag.subregions[[i]],col=i)
  legend("topright",subregions,pch="-",col=1:length(subregions),bty="y",cex = 1,pt.cex=2)
  dev.off()}

setwd("/net/h2o/climphys1/rdaenzer/figures/elasticnet_regression/")
{jpeg("01_elasticnet_cumulative_lag_cor_allsubregions.jpg")
plot(lags,correlation.cum.lag.subregions.atlantic[[1]],main="ridge regression lag correlation, Jul", 
     ylim = c(0, 1), xlim = c(0, length(lags)),xlab="time lag [months]", ylab="correlation",type="l")
for (i in 2:3) lines(lags,correlation.cum.lag.subregions.atlantic[[i]],col=i)
for (i in 1:2) lines(lags,correlation.cum.lag.subregions.pacific[[i]],col=i+3)
legend("bottomright",c("Sargasso","Gulfstream","Spain","East Pacific","West Pacific"),pch="-",col=1:5,bty="y",cex = 1,pt.cex=2)
dev.off()}
