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
  
  subregions="whole"
  #if (region=="Atlantic") subregions=c("North","Sargasso","Gulfstream","Spain")
  #if (region=="Pacific") subregions=c("East","West")
    
  
  outpath="/net/h2o/climphys1/rdaenzer/output"
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region)
  setwd(path)
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  #create arrays for dataframe-inputs
  df.region=c()
  df.member=c()
  df.month=c()
  df.cor.pred.scenario.anom=c()
  df.rmse.residuals.ensmean=c()
  df.rmse.scenario.ensmean=c()
  df.signal.to.noise.scenario=c()
  df.signal.to.noise.dyn.adj=c()
  
  #read ncdf files into RB objects
  #control data
  control.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control")
  setwd(control.path)
  X.control.anom.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
  Y.control.anom.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
  
  
  #read ensemble mean
  ensmean.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/ENSMEAN/")
  setwd(ensmean.path)
  X.ENSMEAN.data.RB=brick(paste0(X, "_cesmENSMEAN_monthly_",region,".nc"))
  Y.ENSMEAN.data.RB=brick(paste0(Y, "_cesmENSMEAN_monthly_remap_",region,".nc"))
  
  #read scenario data
  for (ensemble.member in seq(580,980,20)){
    scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",ensemble.member)
    setwd(scenario.path)
    X.scenario.anom.data.RB=brick(paste0(X, "_cesm",ensemble.member,"_monthly_anom_",region,".nc"))
    Y.scenario.anom.data.RB=brick(paste0(Y, "_cesm",ensemble.member,"_monthly_remap_anom_",region,".nc"))
    Y.scenario.data.RB=brick(paste0(Y, "_cesm",ensemble.member,"_monthly_remap_",region,".nc"))
    
    # ------------------------------------------------------------------------------------------ 
    # 01. Preprocessing: Define regions, crop data, area weighting and mean
    # ------------------------------------------------------------------------------------------
    #define subregion, crop data to regional extent
    for (subregion in subregions){
      if(subregion=="whole") subregion=NULL
      
      print(paste0(region,", ",subregion,", member ",ensemble.member))
    
      Y.extent=define.region.extent(region,subregion)
      X.extent=Y.extent+40
      
      all.region=region
      if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)
      
      #crop data
      X.control.anom.region.RB = crop(X.control.anom.data.RB, X.extent)
      Y.control.anom.region.RB = crop(Y.control.anom.data.RB, Y.extent)
      X.scenario.anom.region.RB = crop(X.scenario.anom.data.RB, X.extent)
      Y.scenario.anom.region.RB = crop(Y.scenario.anom.data.RB, Y.extent)
      Y.scenario.region.RB = crop(Y.scenario.data.RB, Y.extent)
      X.ENSMEAN.region.RB = crop(X.ENSMEAN.data.RB, X.extent)
      Y.ENSMEAN.region.RB = crop(Y.ENSMEAN.data.RB, Y.extent)
      
      #area weighting
      X.control.anom.aw.RB=area.weight.RB(X.control.anom.region.RB) #area weighting of predictor field
      X.scenario.anom.aw.RB=area.weight.RB(X.scenario.anom.region.RB)
      X.ENSMEAN.aw.RB=area.weight.RB(X.ENSMEAN.region.RB)
      
      #area weighted mean
      Y.control.anom.mean.aw=mean.area.weighted.RB(Y.control.anom.region.RB)
      Y.scenario.anom.mean.aw=mean.area.weighted.RB(Y.scenario.anom.region.RB)
      Y.scenario.mean.aw=mean.area.weighted.RB(Y.scenario.region.RB)
      Y.ENSMEAN.mean.aw=mean.area.weighted.RB(Y.ENSMEAN.region.RB)
      
      date.control = as.Date(substring(text = names(X.control.anom.region.RB), first = 2, last = 11), "%Y.%m.%d")
      date.scenario = as.Date(substring(text = names(X.scenario.anom.region.RB), first = 2, last = 11), "%Y.%m.%d")
      
      #detrend Scenario SST
      Y.detr.anom.scenario.mean.aw = detrend.ts(Y.scenario.anom.mean.aw)
      Y.detr.scenario.mean.aw = detrend.ts(Y.scenario.mean.aw)
      Y.detr.ENSMEAN.mean.aw=detrend.ts(Y.ENSMEAN.mean.aw)
      
      #spatial raster
      spatial.raster=raster(X.scenario.region.RB,layer=1)
      values(spatial.raster)=NA
      
      
  # ------------------------------------------------------------------------------------------ 
  # 02. Dynamical Adjustment of mean SST
  # ------------------------------------------------------------------------------------------
      #train model using control data
      months=1:12
      lags=1:12;s.lags=1:6;w.lags=1:12
      train.years=1001:2000
      pred.years=1851:2099
      
      #create X.pred and X.train manually for using SST as predictor
      #X.pred=cbind(t(values(X.scenario.aw.RB)), t(values(Y.scenario.aw.RB)))
      #X.train=cbind(t(values(X.control.aw.RB)), t(values(Y.control.aw.RB)))
      
      cv.model=dynamical.adjustment.gridcell.monthly(Y.train=Y.control.anom.mean.aw,X.pred=X.scenario.anom.aw.RB,X.train=X.control.anom.aw.RB,train.years=train.years,pred.years=pred.years,
                                          alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,ret.cv.model=T)
      
      
      setwd(outpath)
      saveRDS(cv.model,file=paste0("dynamical_adjustment_model_",region,"_member",ensemble.member,".rds"))
      #cv.model=readRDS(paste0("dynamical_adjustment_model_",region,"_member",ensemble.member,".rds"))
      Y.hat=cv.model[[1]]
      Y.model=cv.model[[2]]
      

# ------------------------------------------------------------------------------------------ 
# 03. Calculate Dynamically adjusted Time Series and quantitative values
# ------------------------------------------------------------------------------------------
      
      Y.dyn.adj.anom.mean=Y.hat
      Y.dyn.adj.mean=Y.hat

      for (month in 1:12) {
        idx=which(as.numeric(format(date.scenario, "%Y")) %in% pred.years  & as.numeric(format(date.scenario, "%m")) %in% month)
        Y.dyn.adj.anom.mean[idx-12]=Y.detr.anom.scenario.mean.aw[idx]-Y.hat[idx-12]
        Y.dyn.adj.mean[idx-12]=Y.scenario.mean.aw[idx]-Y.hat[idx-12]
        
        #calculate signal to noise
        stn.years=1971:2020

        #create dataframe data to store as CSV
        df.region=append(df.region, region)
        df.member=append(df.member, ensemble.member)
        df.month=append(df.month, month)
        df.cor.pred.scenario.anom=append(df.cor.pred.scenario.anom, round(cor(Y.hat[idx-12],Y.detr.anom.scenario.mean.aw[idx]),3))
        df.rmse.residuals.ensmean=append(df.rmse.residuals.ensmean, round(rmse(Y.ENSMEAN.mean.aw[idx],Y.dyn.adj.mean[idx-12]),3))
        df.rmse.scenario.ensmean=append(df.rmse.scenario.ensmean, round(rmse(Y.ENSMEAN.mean.aw[idx],Y.scenario.mean.aw[idx]),3))
        df.signal.to.noise.scenario=append(df.signal.to.noise.scenario, round(signal.to.noise(Y.scenario.mean.aw[idx],stn.years,month),3))
        df.signal.to.noise.dyn.adj=append(df.signal.to.noise.dyn.adj, round(signal.to.noise(Y.dyn.adj.mean[idx-12],stn.years,month),3))
      }
      
  
# ------------------------------------------------------------------------------------------ 
# 04. Plotting
# ------------------------------------------------------------------------------------------
      plot.months=c(1,7)
      plot.years=1861:2040
      
      for (plot.month in plot.months) {
        plot.idx=which(as.numeric(format(date.scenario, "%Y")) %in% plot.years  & as.numeric(format(date.scenario, "%m")) %in% plot.month)
        
        #plot
        setwd(figpath)
        ylim=c(floor(min(Y.detr.anom.scenario.mean.aw[plot.idx])),ceiling(max(Y.detr.anom.scenario.mean.aw[plot.idx])))
        {png(paste0("dynamical_adjustment_detr_",region,"_",month.abb[plot.month],"_member",ensemble.member,".png"),width=500,height=350)
        plot(pred.years,Y.detr.anom.scenario.mean.aw[plot.idx],ylim=ylim,type="l",xlab="year a.d.",ylab="mean SST anomalies [°C]",lwd=1.5,xaxt="n")
        axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
        lines(pred.years,Y.hat[plot.idx-12],col=3)
        lines(pred.years,Y.dyn.adj.anom.mean[plot.idx-12],col=2)
        legend("bottomright",c("original time series","prediction","residual time series"),pch="-",col=c(1,3,2),bty="y",cex = 1,pt.cex=2)
        dev.off()}
        
        ylim=c(floor(min(Y.scenario.mean.aw[plot.idx])),ceiling(max(Y.scenario.mean.aw[plot.idx])))
        {png(paste0("dynamical_adjustment_",region,"_",month.abb[plot.month],"_member",ensemble.member,".png"),width=500,height=350)
        plot(pred.years,Y.scenario.mean.aw[plot.idx],col=4,type='l',ylim=ylim,xlab="year",ylab="mean SST [K]",xaxt="n")
        axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
        lines(pred.years,Y.dyn.adj.mean[plot.idx-12],col=2)
        lines(pred.years,Y.ENSMEAN.mean.aw[plot.idx],col=1,lwd=1.4)
        legend("topleft",c("original time series","residual time series","ensemble mean"),pch="-",col=c(1,2,4),bty="y",cex = 1,pt.cex=2)
        dev.off()}
        
        }
    
    
      region=all.region
    }
  }
  
  csv.df=data.frame(df.region,df.member,df.month,df.cor.pred.scenario.anom, df.rmse.scenario.ensmean,
                    df.rmse.residuals.ensmean,df.signal.to.noise.scenario,df.signal.to.noise.dyn.adj)
  setwd(figpath)
  write.csv(csv.df,"dynamical_adjustment_values.csv")
}
#--------------------------------------------------------------------------------------------------
# 04. calculate and plot beta coefficients on map and as time series
#--------------------------------------------------------------------------------------------------
lag=0
idx=(lag*dim(values(X.scenario.aw.RB))[1]+1):((lag+1)*dim(values(X.scenario.aw.RB))[1])
beta1=Y.model[[plot.month]]$glmnet.fit$beta[idx,which(Y.model[[plot.month]]$lambda %in% Y.model[[plot.month]]$lambda.min)]
#plot with filled.contour
color.palette=colorRampPalette(c('blue','green','white',
                                 'yellow','red'),interpolate='spline')
levels=seq(-2e-4,2e-4,length.out=61)

#set lat and lon
lon=unique(coordinates(X.scenario.aw.RB)[,1])
lat=rev(unique(coordinates(X.scenario.aw.RB)[,2]))

mapmat=matrix(beta1,nrow=dim(X.scenario.aw.RB)[2])
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
{png(paste0("elasticnet_coefficients_",region,"_",month.abb[plot.month],"_lag",lag,".png"),width=500,height=350)
  
filled.contour(lon, lat, mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main=paste("elasticnet regression coefficients for ",month.abb[plot.month],", lag ",lag)),
               plot.axes={axis(1); axis(2);maps::map('world', add=TRUE);grid()},
               key.title=title(main="cor"))
dev.off()}



#calculate beta norm
cv.mean=matrix(NA,nrow=12,ncol=13)
for (month in months){
  for (i in 0:12){
    idx=(i*dim(values(X.scenario.aw.RB))[1]+1):((i+1)*dim(values(X.scenario.aw.RB))[1])
    beta1=Y.model[[month]]$glmnet.fit$beta[idx,which(Y.model[[month]]$lambda %in% Y.model[[month]]$lambda.min)]
    cv.mean[month,i+1]=norm(beta1,type="2")
  }
}

#plot for one month per season
months=c(1,4,7,10)  
figname=paste0("elasticnet_beta_norm_",region,".png")
{png(figname,width=500,height=350)
  plot(cv.mean[months[1],],ylim=c(0,0.001),type="l",col=1,main="beta2 norm",xlab="time lag [months]",ylab="beta2 norm")
  for (j in 2:4) lines(cv.mean[months[j],],col=j)
  legend("topright",month.abb[months],pch="-",col=1:4,bty="y",cex = 1,pt.cex=2)
  
dev.off()}
