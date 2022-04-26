#----------------------------------------------------------------------------------------
# Dynamical Adjustment Script for Model Ensemble
#----------------------------------------------------------------------------------------

#structure of script
#0. read in data
#1. process data --> detrend
#2. plot some preliminary analysis plots
#3. train dataset on a region
#4. predict Y for the region
#5. compare prediction with observations
#6. subtract prediction from observations to obtain a ts with higher signal to noise ratio

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
# other functions to read:
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/convert.to.eurocentric.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/frenchcolormap.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/project_raster.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/gridcorts.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/filled.contour2.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/filled.contour3.R")


source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/02_elasticnet_regression.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/03_dynamical_adjustment.R")


# ------------------------------------------------------------------------------------------
# 0.c) Read in Data
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")
for (region in regions){
  subregion=NULL #e.g. North, Gulfstream, Sargasso, Spain
  
  print(paste0(region,", ",subregion))
  
  
  outpath="/net/h2o/climphys1/rdaenzer/output"
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region)
  setwd(path)
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  
  #read ncdf files into RB objects
  #control data
  control.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control")
  setwd(control.path)
  X.control.data.RB = brick(paste0(X, "_",region,".nc")) 
  Y.control.data.RB = brick(paste0(Y,"_",region,".nc"))
  X.control.anom.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
  Y.control.anom.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
  
  #scenario data
  X.ensemble.data=list()
  Y.ensemble.data=list()
  X.ensemble.anom.data=list()
  Y.ensemble.anom.data=list()
  
  ensemble.members=seq(580,980,20)
  for (i in 1:length(ensemble.members)) {
    member=ensemble.members[i]
    scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",member)
    setwd(scenario.path)
    X.ensemble.data[[i]]=brick(paste0(X, "_cesm",member,"_monthly_",region,".nc"))
    Y.ensemble.data[[i]]=brick(paste0(Y, "_cesm",member,"_monthly_remap_",region,".nc"))
    X.ensemble.anom.data[[i]]=brick(paste0(X, "_cesm",member,"_monthly_anom_",region,".nc"))
    Y.ensemble.anom.data[[i]]=brick(paste0(Y, "_cesm",member,"_monthly_remap_anom_",region,".nc"))
  }
  names(X.ensemble.data)=ensemble.members
  names(Y.ensemble.data)=ensemble.members
  names(X.ensemble.anom.data)=ensemble.members
  names(Y.ensemble.anom.data)=ensemble.members
  
  #ENSMEAN
  member="ENSMEAN"
  scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",member)
  setwd(scenario.path)
  X.ENSMEAN.data.RB=brick(paste0(X, "_cesm",member,"_monthly_",region,".nc"))
  Y.ENSMEAN.data.RB=brick(paste0(Y, "_cesm",member,"_monthly_remap_",region,".nc"))
  
  
  # ------------------------------------------------------------------------------------------ 
  # 01. Preprocessing: Define regions, crop data, area weighting and mean
  # ------------------------------------------------------------------------------------------
  #define subregion, crop data to regional extent
  Y.extent=define.region.extent(region,subregion)
  X.extent=Y.extent+40
  if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)
  
  #crop and area weigh data
  X.control.region.RB = crop(X.control.data.RB, X.extent) #crop and area weigh control data
  Y.control.region.RB = crop(Y.control.data.RB, Y.extent)
  X.control.anom.region.RB = crop(X.control.anom.data.RB, X.extent) #crop and area weigh control data
  Y.control.anom.region.RB = crop(Y.control.anom.data.RB, Y.extent)
  Y.control.mean.aw=mean.area.weighted.RB(Y.control.region.RB)
  Y.control.anom.mean.aw=mean.area.weighted.RB(Y.control.anom.region.RB)
  
  X.ensemble.region=list()
  Y.ensemble.region=list()
  X.ensemble.anom.region=list()
  Y.ensemble.anom.region=list()
  Y.ensemble.region.mean=c()
  Y.ensemble.anom.region.mean=c()
  
  for (i in 1:length(ensemble.members)) { #crop  and area weigh scenario data
    X.ensemble.region[[i]] = crop(X.ensemble.data[[i]], X.extent)
    Y.ensemble.region[[i]] = crop(Y.ensemble.data[[i]], Y.extent)
    X.ensemble.anom.region[[i]] = crop(X.ensemble.anom.data[[i]], X.extent)
    Y.ensemble.anom.region[[i]] = crop(Y.ensemble.anom.data[[i]], Y.extent)
    Y.ensemble.region.mean=rbind(Y.ensemble.region.mean,mean.area.weighted.RB(Y.ensemble.region[[i]]))
    Y.ensemble.anom.region.mean=rbind(Y.ensemble.anom.region.mean,mean.area.weighted.RB(Y.ensemble.anom.region[[i]]))
    
  }
  
  Y.ENSMEAN.region.RB=crop(Y.ENSMEAN.data.RB,Y.extent)
  Y.ENSMEAN.region.mean=mean.area.weighted.RB(Y.ENSMEAN.region.RB)
  
  #get dates
  date.control = as.Date(substring(text = names(X.control.region.RB), first = 2, last = 11), "%Y.%m.%d")
  date.scenario = as.Date(substring(text = names(X.ensemble.region[[1]]), first = 2, last = 11), "%Y.%m.%d")
  
  # ------------------------------------------------------------------------------------------ 
  # 02. Dynamical Adjustment of mean SST, whole ensemble
  # ------------------------------------------------------------------------------------------
  #train model using control data
  months=c(1,7)
  lags=1:12
  train.years=1001:2020
  pred.years=1851:2040
  
  for (month in months){
    Y.hat=list()
    Y.residual=c()
    for (i in 1:length(ensemble.members)){
      print(ensemble.members[i])
      Y.hat[[i]]=dynamical.adjustment.gridcell.monthly(Y.train=Y.control.anom.mean.aw,X.pred=X.ensemble.anom.region[[i]],X.train=X.control.anom.region.RB,train.years=train.years,pred.years=pred.years,
                                                   alpha=0.1,lags=lags,months=month,ret.cv.model=F)
      Y.residual=rbind(Y.residual,(Y.ensemble.region.mean[i,pred.idx] - Y.hat[[i]]))  
    }
    
    
    #save/read stored data instead of calculating
    setwd(outpath)
    #saveRDS(Y.hat,paste0("dynamical_adjustment_compare_ensemble_",region,"_",month.abb[month],".rds"))
    #saveRDS(Y.residual,paste0("dynamical_adjustment_compare_ensemble_residual",region,"_",month.abb[month],".rds"))
    Y.hat=readRDS(paste0("dynamical_adjustment_compare_ensemble_Yhat_",region,"_",month.abb[month],".rds"))
    Y.residual=readRDS(paste0("dynamical_adjustment_compare_ensemble_residual_",region,"_",month.abb[month],".rds"))
    
    # ------------------------------------------------------------------------------------------ 
    # 03.Calculate and visualize range of ensemble and dyn.adj. ensemble
    # ------------------------------------------------------------------------------------------
    setwd(figpath)
    
    #calculate mean and sd for shading --> maybe better?
    ensemble.mean=sapply(1:dim(Y.ensemble.region.mean)[2], FUN=function(idx) mean(Y.ensemble.region.mean[,idx]))
    ensemble.sd=sapply(1:dim(Y.ensemble.region.mean)[2], FUN=function(idx) sd(Y.ensemble.region.mean[,idx]))
    residual.mean=sapply(1:dim(Y.residual)[2], FUN=function(idx) mean(Y.residual[,idx]))
    residual.sd=sapply(1:dim(Y.residual)[2], FUN=function(idx) sd(Y.residual[,idx]))
    
    
    #plot with sd and mean
    {figname=paste0("dynamical_adjustment_ensemble_comparison_",region,"_",month.abb[month],".png")
    png(figname,width=500,height=350)
    par(mar=c(2,2,2,2),oma=c(0,0,0,0))
    plot(pred.years,Y.ENSMEAN.region.mean[pred.idx],type="l",ylab="mean SST",xlab="year",
         ylim=c(floor(min(Y.ENSMEAN.region.mean[pred.idx])),ceiling(max(Y.ENSMEAN.region.mean[pred.idx]))),main=paste(region,month.abb[month]),xaxt="n")
    polygon(c(pred.years,rev(pred.years)),c((ensemble.mean-ensemble.sd)[pred.idx],rev((ensemble.mean+ensemble.sd)[pred.idx])),col=rgb(0.3,0,1,0.1))
    polygon(c(pred.years,rev(pred.years)),c((residual.mean-residual.sd),rev((residual.mean+residual.sd))),col=rgb(1,0,0.3,0.3))
    lines(pred.years,Y.ENSMEAN.region.mean[pred.idx],type="l",lwd=2)
    lines(pred.years,(ensemble.mean-ensemble.sd)[pred.idx],type="l",col=rgb(0.3,0,1,1))
    lines(pred.years,(ensemble.mean+ensemble.sd)[pred.idx],type="l",col=rgb(0.3,0,1,1))
    lines(pred.years,residual.mean-residual.sd,type="l",col=rgb(1,0,0.3,1))
    lines(pred.years,residual.mean+residual.sd,type="l",col=rgb(1,0,0.3,1))
    axis(side=1,at=c(1850,1900,1950,2000,2040), labels=c(1850,1900,1950,2000,2040))
    legend("topleft",c("st. dev. original time series","st. dev. residual time series","ensemble mean"),pch="-",col=c(rgb(0.3,0,1,1),rgb(1,0,0.3,1),1),bty="y",cex = 1,pt.cex=5)
    dev.off()}
  }
}
# ------------------------------------------------------------------------------------------ 
# 04. calculate and visualize trends of individual ensemble members and residuals
# ------------------------------------------------------------------------------------------
#calculate trends
trend.years=1971:2020
trend.season="djf"
trend.month=7
trend.residual=c()
trend.original=c()
for (i in 1:length(ensemble.members)){
  #seasonally
  #cur.orig.ts=period.mean.ts(Y.ensemble.region.mean[i,],trend.years,months=trend.season)
  #cur.res.ts=period.mean.ts(Y.residual[i,],trend.years,months=trend.season)
  
  #monthlz
  cur.orig.ts=extract.period.ts(Y.ensemble.region.mean[i,],trend.years,trend.month)
  cur.res.ts=extract.period.ts(Y.residual[i,],trend.years,trend.month)
  years=as.numeric(format(as.Date(names(cur.orig.ts)), "%Y"))
  cur.orig.trend=lm(cur.orig.ts~years)$coefficients[2] * length(years)
  cur.residual.trend=lm(cur.res.ts~years)$coefficients[2] * length(years)
  trend.original=rbind(trend.original,cur.orig.trend)
  trend.residual=rbind(trend.residual,cur.residual.trend)

}

#calculate season trend
ens.seasonal.mean=period.mean.ts(Y.ENSMEAN.region.mean,trend.years,months=trend.season)
trend.ensemble=lm(ens.seasonal.mean~years)$coefficients[2] * length(years)

hist.orig=hist(trend.original)
hist.orig$counts=hist.orig$counts / sum(hist.orig$counts)

hist.res=hist(trend.residual)  
hist.res$counts=hist.res$counts / sum(hist.res$counts)

plot(hist.orig,xlim=c(0.5,1.5),ylim=c(0,0.4))
plot(hist.res,add=T,border=4)
