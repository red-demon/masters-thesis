
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


# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
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
regions=c("Atlantic","Pacific")
  
  
for (region in regions){
  
  if (region=="Atlantic") subregions=c("Sargasso","Gulfstream","Spain")
  if (region=="Pacific") subregions=c("East","West","Central","NorthNP")
  
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
  for (subregion in subregions){
    if(subregion=="whole") subregion=NULL
    
    print(paste(region,subregion))
    
    #define subregion, crop data to regional extent
    Y.extent=define.region.extent(region,subregion)
    X.extent=Y.extent+40 #if aggregated, the predictor field should be the whole X-field
    
    all.region=region #to save main region for later
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
    setwd(figpath)
    
    #elasticnet regression for Y.mean using first 100 EOFs as predictors
    months=c(1,7)
    lags=0:24
    train.years=1:1000
    pred.years=1001:1998
    s="lambda.min"
    alpha=0.1
    
    #extract and define Y.obs and predictor
    predictor=t(values(X.region.aw.RB)) #define predictor for ridge regression (either t(X.region.aw) or X.svd.mon.all.v[,eofs])
    
    
    # ------------------------------------------------------------------------------------------ 
    #correlation between monthly mean SST and elast.net predicted SST for specific lags, only 2 months
    #------------------------------------------------------------------------------------------
    correlation.lag=lapply(X=months,FUN=function(month){
        sapply(X=lags,FUN=function(lag) cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=lag,train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
    })
    
    # ------------------------------------------------------------------------------------------   
    #correlation between monthly mean sst and elast.net predicted SST for cumulative lags, only 2 months
    #------------------------------------------------------------------------------------------------
    correlation.cum.lag=lapply(X=months,FUN=function(month){
      sapply(X=lags,FUN=function(lag) cor(extract.lag.period.idx(Y.mean.aw,month,0)[pred.years],elasticnet.gridcell.monthly(Y.mean.aw,predictor,lags=c(0:lag),train.years = train.years,pred.years=pred.years,month=month,alpha=alpha,s=s)))
    })
    
    #save/read correlations
    setwd(outpath)
    saveRDS(correlation.lag,file=paste0("glmnet_regression_correlation_",region,".rds"))
    saveRDS(correlation.cum.lag,file=paste0("glmnet_regression_correlation_cumulative_",region,".rds"))
    #correlation.lag=readRDS(file=paste0("glmnet_regression_correlation_",region,".rds"))
    #correlation.cum.lag=readRDS(file=paste0("glmnet_regression_correlation_cumulative_",region,".rds"))
    
    #------------------------------------------------------------------------------------------   
    # 02. plotting correlations
    #------------------------------------------------------------------------------------------
    setwd(figpath)
    #plot correlation between monthly mean sst and elast.net predicted SST for specific lags, only 2 months
    figname=paste0("01_elasticnet_lag_cor_",region,".png")
    {png(figname,width=500,height=350)
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(lags,correlation.lag[[1]],col=4,
           ylim = c(0, 1), xlim = c(0, length(lags)),axes=F,xlab="time lag [months]", ylab="correlation",type="l")
      lines(lags,correlation.lag[[2]],col=2)
      legend("topright",month.abb[months],lty=1,col=c(4,2),bty="y",cex = 1,pt.cex=2)
      axis(1,at=seq(0,max(lags),3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,max(lags),3),lty="dotted",col="darkgrey")
      dev.off()}
    
  
    #plot correlation between monthly mean sst and elast.net predicted SST for cumulative lags, only 2 months
    months=c(1,7)
    figname=paste0("01_elasticnet_cumulative_lag_cor_",region,".png")
    {png(figname,width=500,height=350)
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(lags,correlation.cum.lag[[1]], col=4,
           ylim = c(0, 1), xlim = c(0, length(lags)),axes=F,xlab="time lag [months]", ylab="correlation",type="l")
      lines(lags,correlation.cum.lag[[2]],col=2)
      legend("bottomright",month.abb[months],lty=1,col=c(4,2),bty="y",cex = 1,pt.cex=2)
      axis(1,at=seq(0,max(lags),3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,max(lags),3),lty="dotted",col="darkgrey")
      dev.off()
    }
    
    region=all.region
    
  }
}

    
    
    #plot lag maximum
    months=1:12
    max_lag_elasticnet=sapply(X=1:12,FUN=function(month)  which(correlation.cum.lag[[months[month]]] %in% max(correlation.cum.lag[[months[month]]])))
    figname=paste0("01_elasticnet_cumulative_max_lag_",region,".png")
    {png(figname,width=500,height=350)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(1:12,max_lag_elasticnet,ylim=c(0,24),axes=F,type="l",col=4,ylab="maximum lag",xlab="month")
      axis(1,at=seq(0,nr.lags,3))
      axis(2,at=seq(0,24,6))
      abline(h=seq(0,24,6),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      dev.off()}
    
    #plot cumulative lag correlation for 4 months
    months=c(1,4,7,10)
    figname=paste0("01_elasticnet_cumulative_lag_cor_JAJO_",region,".png")
    {png(figname,width=500,height=300)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(lags,correlation.cum.lag[[months[1]]],axes=F,ylim=c(0,1),type="l",col=1,ylab="correlation",xlab="time lag [months]")
      for (i in 2:length(months)) lines(x = lags, y = correlation.cum.lag[[months[i]]],col=i)
      legend("bottomright",month.abb[months],lty=1,col=1:length(months),bty="y",cex = 1,pt.cex=2)
      axis(1,at=seq(0,nr.lags,3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      dev.off()
    }
  

    
  
  
