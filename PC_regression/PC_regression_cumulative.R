#load packages and sources
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)

source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/additional_functions/filled_contour3.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/colormaps.R")

# ------------------------------------------------------------------------------------------ 
# 00. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")

for (region in regions){
  subregions="whole"
  #if (region=="Atlantic") subregions=c("North","Sargasso","Gulfstream","Spain")
  #if (region=="Pacific") subregions=c("East","West","Central","NorthNP")
  
  #set directories
  path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/PC_regression/",region)
  outpath="/net/h2o/climphys1/rdaenzer/output/"
  setwd(path)
  
  # ------------------------------------------------------------------------------------------ 
  # 00. Read processed data and SVDs
  # ------------------------------------------------------------------------------------------
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  
  #read ncdf files into RB objects
  X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
  Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
  
  #calculate total SVDs and store in separate objects
  X.region.svd=values(area.weight.RB(X.region.RB,svd=TRUE))
  X.svd=svd(X.region.svd)
  X.svd.u=X.svd$u
  X.svd.v=X.svd$v
  X.svd.d=X.svd$d
  X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 
  
  #calculate SVDs for each month separately
  {X.svd.mon=list()
    X.svd.mon.u=list()
    X.svd.mon.v=list()
    X.svd.mon.d=list()
    for (mon in 1:12){
      month=month.abb[mon]
      X.svd.mon.temp=svd(X.region.svd[,seq(mon,length(X.region.svd[1,]),12)])
      X.svd.mon.temp$d=X.svd.mon.temp$d^2/sum(X.svd.mon.temp$d^2)
      X.svd.mon[[month]] <- X.svd.mon.temp
      X.svd.mon.v[[month]]<-X.svd.mon.temp$v
    }
  }
  
  #store all monthly PC ts in one array
  X.svd.mon.all.v=matrix(NA,nrow=dim(X.svd.v)[1],ncol=dim(X.svd.mon.v[[1]])[2])
  for (i in 1:12){
    idx=seq(i,length(X.svd$v[,1]),12)
    X.svd.mon.all.v[idx,]=X.svd.mon.v[[i]]
  }
  for (subregion in subregions){
    if(subregion=="whole") subregion=NULL
    
    print(paste(region,subregion))
    
    Y.extent=define.region.extent(region,subregion)
    Y.region.RB = crop(Y.data.RB, Y.extent)
    
    #area weighting and centering
    
    Y.region.aw.RB=area.weight.RB(Y.region.RB)
    Y.mean.aw=mean.area.weighted.RB(Y.region.RB)
    

    

    
    all.region=region #to save region for later 
    if(!is.null(subregion)) region=subregion #for plotting names
    
    #-------------------------------------------------------------------------------------
    # 08. Cumulative Lagged SVD regression
    #-------------------------------------------------------------------------------------
    print("08")
    nr.lags=24
    lags=0:nr.lags
    train.years=1:1000
    pred.years=1001:1996
    all.eofs=c(5,20)
    
    for (eofs in all.eofs){
      correlation.cum.lag=sapply(X=1:12,FUN=function(month.idx) svd.regression.lag.cumulative(Y.ts=Y.mean.aw,X.svd=X.svd.mon.all.v,eofs=1:eofs,train.years=train.years,pred.years=pred.years,month.idx=month.idx,lags=lags))
      
      #for plotting all EOFs and both months in one plot
      if(eofs==5) cum.lag.5eofs=correlation.cum.lag
      if(eofs==20) cum.lag.20eofs=correlation.cum.lag
      
    }
    
    setwd(figpath)
    #plot cumulative lag correlation
    months=c(1,7)  
    figname=paste0("08_cumulative_lag_corr_",region,"_Jan_Jul_wholeSVD.png")
    {png(figname,width=500,height=350)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,3,0),bty="n")
      plot(lags,cum.lag.5eofs[,months[1]],axes=F,ylim=c(0,1),type="l",col=4,main=region,ylab="correlation",xlab="time lag [months]")
      lines(x = lags, y = cum.lag.5eofs[,months[2]],col=2)
      
      #20 EOFs
      lines(x = lags, y = cum.lag.20eofs[,months[1]],lty=2,col=4)
      lines(x = lags, y = cum.lag.20eofs[,months[2]],lty=2,col=2)
      
      legend("bottomright",c("Jan, 5 EOFs","Jan, 20 EOFs","Jul, 5 EOFs","Jul, 20 EOFs"),lty=c(1,2,1,2),col=c(4,4,2,2),bty="y",cex = 1,pt.cex=2)
      axis(1,at=seq(0,nr.lags,3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      dev.off()}
    
    region=all.region
    
  }
}

    #more plots
    #plot lag maximum
    max_lag_20eofs=sapply(X=1:12,FUN=function(month)  which(cum.lag.20eofs[,month] %in% max(cum.lag.20eofs[,month])))
    figname=paste0("08_cumulative_lag_corr_",region,"_max_lag.png")
    {png(figname,width=500,height=350)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,3,0),bty="n")
      plot(1:12,max_lag_20eofs,ylim=c(0,24),axes=F,type="l",col=4,main=region,ylab="maximum lag",xlab="month")
      axis(1,at=seq(0,nr.lags,3))
      axis(2,at=seq(0,24,6))
      abline(h=seq(0,24,6),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      dev.off()}
    
    #plot cumulative lag correlation for 4 months
    months=c(1,4,7,10)
    figname=paste0("08_cumulative_lag_corr_20eofs_",region,"_JAJO.png")
    {png(figname,width=500,height=300)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,3,0),bty="n")
      plot(lags,cum.lag.20eofs[,months[1]],axes=F,ylim=c(0,1),type="l",col=1,main=region,ylab="correlation",xlab="time lag [months]")
      for (i in 2:length(months)) lines(x = lags, y = cum.lag.20eofs[,months[i]],col=i)
      legend("bottomright",month.abb[months],lty=1,col=1:length(months),bty="y",cex = 1,pt.cex=2)
      axis(1,at=seq(0,nr.lags,3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      dev.off()}
    

