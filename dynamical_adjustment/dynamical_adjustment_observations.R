#field dyn.adj. testscript

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
region="Atlantic"
subregion=NULL

path=file.path("/net/h2o/climphys1/rdaenzer/data/observations")
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
#control data
X.obs.data.RB = brick(paste0("PSL_observations_mean_anom_remap_",region,".nc")) 
Y.obs.data.RB = brick(paste0("SST_observations_median_anom_remap_",region,".nc"))
Y.obs.data.RB=subset(Y.obs.data.RB,13:1980) #shorten to same length as X.data
# ------------------------------------------------------------------------------------------ 
# 01. Preprocessing: Define regions, crop data, area weighting and mean
# ------------------------------------------------------------------------------------------
#define subregion, crop data to regional extent
Y.extent=define.region.extent(region,subregion)
X.extent=Y.extent+40

if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)

#crop data
X.obs.region.RB = crop(X.obs.data.RB, X.extent)
Y.obs.region.RB = crop(Y.obs.data.RB, Y.extent)
Y.obs.region.mean=mean.area.weighted.RB(Y.obs.region.RB)
Y.obs.region.detr.mean=detrend.ts(Y.obs.region.mean)


# ------------------------------------------------------------------------------------------ 
# 02. Dynamical adjustment
# ------------------------------------------------------------------------------------------
#train model using control data
months=1:12
lags=1:12
train.years=1901:1970
pred.years=1971:2014

date.obs= as.Date(substring(text = names(X.obs.data.RB), first = 2, last = 11), "%Y.%m.%d")

train.idx=which(as.numeric(format(date.obs, "%Y")) %in% train.years & as.numeric(format(date.obs, "%m")) %in% months)
pred.idx=which(as.numeric(format(date.obs, "%Y")) %in% pred.years  & as.numeric(format(date.obs, "%m")) %in% months)

#create X.pred and X.train manually for using SST as predictor
#X.pred=cbind(t(values(X.scenario.aw.RB)), t(values(Y.scenario.aw.RB)))
#X.train=cbind(t(values(X.control.aw.RB)), t(values(Y.control.aw.RB)))


cv.model=dynamical.adjustment.gridcell.monthly(Y.train=Y.obs.region.detr.mean,X.pred=X.obs.region.RB,X.train=X.obs.region.RB,train.years=train.years,pred.years=pred.years,
                                               alpha=0.1,lags=lags,months=months,ret.cv.model=T)

Y.hat=cv.model[[1]]
Y.model=cv.model[[2]]


#evaluate model

# ------------------------------------------------------------------------------------------ 
# 03. Plotting
# ------------------------------------------------------------------------------------------
plot.months=c(1,7)
for (plot.month in plot.months) {
  plot.obs.idx=which(as.numeric(format(date.obs, "%Y")) %in% pred.years  & as.numeric(format(date.obs, "%m")) %in% plot.month)
  plot.pred.idx=plot.obs.idx-min(plot.obs.idx)+plot.month
  
  Y.dyn.adj.detr.mean=Y.obs.region.detr.mean[plot.obs.idx]-Y.hat[plot.pred.idx]
  Y.dyn.adj.mean=Y.obs.region.mean[plot.obs.idx]-Y.hat[plot.pred.idx]
  
  #plot
  setwd(figpath)
  {jpeg(paste0("dynamical_adjustment_detr_",region,"_",month.abb[plot.month],"_observations",".jpeg"),width=500,height=300)
    plot(date.scenario[plot.obs.idx],Y.obs.region.detr.mean[plot.obs.idx],type="l",xlab="year a.d.",ylab="mean SST")
    lines(date.scenario[plot.obs.idx],Y.hat[plot.pred.idx],col=3)
    lines(date.scenario[plot.obs.idx],Y.dyn.adj.detr.mean,col=2)
    legend("bottomright",c("observation","prediction","residual timeseries"),pch="-",col=c(1,3,2),bty="y",cex = 1,pt.cex=2)
    legend("bottomleft",paste("cor:", round(cor(Y.hat[plot.pred.idx],Y.obs.region.detr.mean[plot.obs.idx]),3)),cex = 1)
    dev.off()}
  
  {jpeg(paste0("dynamical_adjustment_observations_",region,"_",month.abb[plot.month],"_observations",".jpeg"),width=700,height=500)
    plot(date.scenario[plot.obs.idx],Y.obs.region.mean[plot.obs.idx],col=1,type='l',xlab="year",ylab="mean SST [K]")
    lines(date.scenario[plot.obs.idx],Y.dyn.adj.mean,col=2)
    legend("bottomright",c("scenario","residual timeseries"),pch="-",col=c(1,2),bty="y",cex = 1,pt.cex=2)
    dev.off()}
  
}

