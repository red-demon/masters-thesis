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
  
  if (region=="Atlantic") subregions=c("Sargasso","Gulfstream","Spain")
  if (region=="Pacific") subregions=c("East","West","Central","NorthNP")

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
 
  
  for (subregion in subregions){
    if(subregion=="whole") subregion=NULL
    
    print(paste(region,subregion))
    
    Y.extent=define.region.extent(region,subregion)
    X.extent=Y.extent+40
    
    X.region.RB = crop(X.data.RB, X.extent)
    Y.region.RB = crop(Y.data.RB, Y.extent)
  
    #area weighting and centering
    
    X.region.aw.RB=area.weight.RB(X.region.RB)
    Y.region.aw.RB=area.weight.RB(Y.region.RB)
    
    Y.mean.aw=mean.area.weighted.RB(Y.region.RB)
    
    X.region.svd=values(area.weight.RB(X.region.RB,svd=TRUE))
    
    #calculate total SVDs and store in separate objects
    #X.svd=svd(X.region.svd)
    #X.svd.u=X.svd$u
    #X.svd.v=X.svd$v
    #X.svd.d=X.svd$d
    #X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 
    
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
        X.svd.mon.u[[month]]<-X.svd.mon.temp$u
        X.svd.mon.d[[month]]<-X.svd.mon.temp$d
        
       
      }
    }
    
    #store all monthly PC ts in one array
    X.svd.mon.all.v=matrix(NA,nrow=dim(X.svd.v)[1],ncol=dim(X.svd.mon.v[[1]])[2])
    for (i in 1:12){
      idx=seq(i,length(X.svd$v[,1]),12)
      X.svd.mon.all.v[idx,]=X.svd.mon.v[[i]]
    }

    
    
    #-----------------------------------------------------------------------------------
    # 05. calculate and plot ccf and acf
    #----------------------------------------------------------------------------------
    print("05")
    setwd(figpath)
    #worldmaps
    if(region=="Pacific") {map="world2";fig.height=350;fig.width=450}
    if(region=="Atlantic") {map="world";fig.height=350;fig.width=360}
    
    all.region=region #to save region for later 
    if(!is.null(subregion)) region=subregion #for plotting names
    
    
    #calculate lag correlation from SVDs calculated monthly
    eofs=100
    nr.lags=48
    lag.cor.mon.eof = lapply(X=1:eofs,FUN=function(eof) sapply(X = 1:12, FUN = function(mon.ix) ccf.monthly(Y.ts = Y.mean.aw,X.ts=X.svd.mon.all.v[,eof], mon.ix = mon.ix, nr.lags = nr.lags)))
    lag.cor.mon.sst =sapply(X = 1:12, FUN = function(mon.ix) ccf.monthly(Y.ts = Y.mean.aw,X.ts=Y.mean.aw, mon.ix = mon.ix, nr.lags = nr.lags))
    sign.level=1.96/sqrt((length(Y.mean.aw)-nr.lags)/12)
    
    #plot months in one plot
    months=c(1,4,7,10)  
    {figname=paste0("05_acf_",region,".png")
        {png(figname,width=500,height=350)
          par(oma=c(0,0,0,0),mar=c(5,4,0,0))
          plot(x = 0:nr.lags, y = abs(lag.cor.mon.sst[,months[1]]),axes=F,ylim=c(0,1),type="l",lwd=1.2,col=1,xlab = "time lag [months]", ylab = "correlation")
          for (j in 2:4) lines(x = 0:nr.lags, y = abs(lag.cor.mon.sst[,months[j]]),col=j,lwd=1.2)
          axis(1,at=seq(0,48,12))
          axis(2)
          abline(h=c(seq(0,1,0.2)),v=seq(0,48,12),lty="dotted",col="darkgrey")
          legend("topright",month.abb[months],lty=1,col=1:4,bty="y",cex = 1,pt.cex=2)
          dev.off()}
    }
    
    
    
    #plot EOFs in one plot
    months=c(1,7)
    eofs=1:5
    for (month in months){
      figname=paste0("05_ccf_acf_",month.abb[month],"_",region,".png")
        {png(figname,width=500,height=350)
          par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
          plot(x = 0:nr.lags, y = abs(lag.cor.mon.sst[,month]),ylim=c(0,1),axes=F,type="l",col=1,xlab = "time lag [months]", ylab = "correlation")
          for (j in eofs) lines(x = 0:nr.lags, y = abs(lag.cor.mon.eof[[j]][,month]),col=j+1)
          axis(1,at=seq(0,48,12))
          axis(2)
          abline(h=c(seq(0,1,0.2)),v=seq(0,48,12),lty="dotted",col="darkgrey")
          legend("topright",c("SST acf",paste0("PC ",eofs)),lty=1,col=1:6,bty="y",cex = 1,pt.cex=2)
          
          dev.off()}
      }
      
      #plot EOFs separately
      months=c(1,7)
      eofs=1:5
      for (month in months){
        for (j in eofs){
          figname=paste0("05_ccf_acf_EOF",j,"_",month.abb[month],"_",region,".png")
          png(figname,width=500,height=300)
          par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
          plot(x = 0:nr.lags, y = abs(lag.cor.mon.sst[,month]),ylim=c(0,1),axes=F,type="l",col=1,xlab = "time lag [months]", ylab = "correlation")
          lines(x = 0:nr.lags, y = abs(lag.cor.mon.eof[[j]][,month]),col=j+1)
          axis(1,at=seq(0,48,12))
          axis(2)
          abline(h=c(seq(0,1,0.2)),v=seq(0,48,12),lty="dotted",col="darkgrey")
          legend("topright",c("SST acf",paste0("PC ",j)),lty=1,col=c(1,j+1),bty="y",cex = 1,pt.cex=2)
          
          dev.off()}
      
    }
    
    
    #-----------------------------------------------------------------------------------
    # 06. How does correlation change with increasing EOF as predictors? 
    #------------------------------------------------------------------------------------
    print("06")
    
    #choose period for the prediction
    months=c(1,7)
    for( month in months){
      idx=seq(month,length(Y.mean.aw),12)
      train.idx=idx[1:1000]#indices for training step
      pred.idx=idx[1001:2000] #indices for prediction
      
      #correlate Y.mean with all PC time series
      correlation.eofs=sapply(X=1:length(X.svd.mon.all.v[1,]),FUN=function(eof) cor(Y.mean.aw[idx],X.svd.mon.all.v[idx,eof]))
      
      #correlate modelled Y.mean with PC-Regression predictedd Y.mean using increasing #EOFs
      #how many eofs to test?
      eofs=300
      correlation.svd=rep(NA,eofs)
      predictor=c()
      correlation.month=list()
      #regress temp on SVD PCs
      for (i in seq(1,eofs)){
        predictor=cbind(predictor,X.svd.mon.all.v[train.idx,i]) #add new eof to predictor-array
        model.eof=lm(Y.mean.aw[train.idx]~predictor) #regression based on predictor-array
        
        #test predictedd data on subset of data
        if(i==1) Y.pred= X.svd.mon.all.v[pred.idx,1:i] * model.eof$coefficients[2:(i+1)] # multiply subset of eofs with coefficients to obtain predictions
        if(i>1) Y.pred= X.svd.mon.all.v[pred.idx,1:i] %*% model.eof$coefficients[2:(i+1)] # multiply subset of eofs with coefficients to obtain predictions
        Y.pred= model.eof$coefficients[1]+Y.pred #add intercept to predictions
        correlation.svd[i]=cor(Y.pred,Y.mean.aw[pred.idx]) # calculate correlation and store in correlation-array
      
      }
      if (month==1) jan=correlation.svd
      if (month==7) jul=correlation.svd

      #plot correlation in dependence of EOFs
      setwd(figpath)
      figname=paste0("06_cor_numbereofs_",month.abb[month],"_",region,".png")
      {png(figname,width=500,height=350)
        par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
        plot(seq(1,eofs),abs(correlation.svd[1:eofs]), type='l',axes=F,ylim=c(0,1),
             ylab="correlation",xlab="# Principal Components")
        #lines(cumsum(correlation.eofs[1:100]^2),type="l",main="cumulative squared correlation of mean SST and PSL EOFs",ylab="pearson correlation,squared",xlab="EOF")
        lines(abs(correlation.eofs[1:eofs]),type="l",ylab="squared correlation",xlab="Principal Components",col=2)
        axis(1,at=seq(0,300,50))
        axis(2)
        abline(h=c(seq(0,1,0.2)),v=seq(0,300,50),lty="dotted",col="darkgrey")
        legend("topright",c("predicted vs. modelled SST","modelled SST vs. PCs"),lty=1,col=1:2,bty="y",cex = 1,pt.cex=2)
        dev.off()}
    }
    
    
    #plot correlation in dependence of EOFs for two months
    figname=paste0("06_cor_numbereofs_Jan_Jul_",region,".png")
    {png(figname,width=500,height=350)
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(seq(1,eofs),jan[1:eofs], type='l',axes=F,ylim=c(0,1), main="",
           ylab="correlation",xlab="Principal Components",col=4)
      #lines(cumsum(correlation.eofs[1:100]^2),type="l",main="cumulative squared correlation of mean SST and PSL PCs",ylab="pearson correlation,squared",xlab="EOF")
      lines(jul[1:eofs],type="l",col=2)
      axis(1,at=seq(0,300,50))
      axis(2)
      abline(h=c(seq(0,1,0.2)),v=seq(0,300,50),lty="dotted",col="darkgrey")
      legend("topright",c("Jan","Jul"),lty=1,col=c(4,2),bty="y",cex = 1,pt.cex=2)
      dev.off()}
    
    #-----------------------------------------------------------------------------------------
    # 07. Lagged SVD Regression
    #------------------------------------------------------------------------------------------
    print("07")
    
    period.names=month.abb
    nr.lags=24
    lags=seq(0,nr.lags)
    all.eofs=c(5,20)
    
    for (eofs in all.eofs){
      #create empty matrix for storing lagged correlations with dimensions lags x periods
      correlation.lagged.svd=matrix(NA,nrow=length(lags),ncol=length(period.names),dimnames=list(lags,period.names))
      
      for (i in 1:length(period.names)){
        
        start.idx=i+nr.lags
        Y.idx=seq(start.idx,(length(Y.mean.aw)-12),12)
        Y.mean.mon=Y.mean.aw[Y.idx]
        
        for (lag in lags){
          #iterate trough lag from month 1 to 12
          
          X.idx=Y.idx-lag
          
          #choose period for the prediction and training
          train.years=seq(1,1000) #indices for training step
          pred.years=seq(1001,length(Y.mean.mon)) #indices for prediction
          
          #svd regression
          Y.pred=svd.regression(Y.mean.mon,X.svd.mon.all.v[X.idx,],predictors=eofs,train.years=train.years,pred.years=pred.years)
          
          #calculate correlation betw. Y.pred and Y.mean for prediction Period
          correlation.lagged.svd[which(lags %in% lag),i]=cor(Y.mean.mon[pred.years],Y.pred)
        }
      }
      
      if (eofs==5) lag.cor.5eofs=correlation.lagged.svd
      if (eofs==20) lag.cor.20eofs=correlation.lagged.svd
      
    }
      
    #plot January and July
    months=c(1,7)
    for (month in months){
      {figname=paste0("07_pc_regr_lag_corr_",paste(Y,region,month.abb[month],sep="_"),".png")
        png(figname,width=500,height=350)
        par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
        plot(lags,lag.cor.mon.sst[,month][lags+1],ylim=c(0,1),type="l",axes=F,col=1,
             ylab="correlation",xlab="time lag [months]")    
        lines(lags,lag.cor.5eofs[,month],type="l",col=2)    
        lines(lags,lag.cor.20eofs[,month],type="l",col=4)    
        
        axis(1,at=seq(0,nr.lags,3))
        axis(2)
        abline(h=seq(0,1,0.2),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
        legend("topright",c("SST  autocorrelation","lagged regression, PC 1-5","lagged regression, PC 1-20"),lty=1,col=c(1,2,4),bty="y",cex = 1,pt.cex=2)
        
        dev.off()
      }
      
      
    }

    
  
    
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
    
    
    months=c(1,7)  
    figname=paste0("08_cumulative_lag_corr_",region,"_Jan_Jul.png")
    {png(figname,width=500,height=350)
      #5 EOFs
      par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
      plot(lags,cum.lag.5eofs[,months[1]],axes=F,ylim=c(0,1),type="l",col=4,ylab="correlation",xlab="time lag [months]")
      lines(x = lags, y = cum.lag.5eofs[,months[2]],col=2)
      
      #20 EOFs
      lines(x = lags, y = cum.lag.20eofs[,months[1]],lty=2,col=4)
      lines(x = lags, y = cum.lag.20eofs[,months[2]],lty=2,col=2)
      
      axis(1,at=seq(0,nr.lags,3))
      axis(2)
      abline(h=seq(0,1,0.2),v=seq(0,nr.lags,3),lty="dotted",col="darkgrey")
      legend("bottomright",c("Jan, PC 1-5","Jan, PC 1-20","Jul, PC 1-5","Jul, PC 1-20"),lty=c(1,2,1,2),col=c(4,4,2,2),bty="y",cex = 1,pt.cex=2)
      dev.off()}
  
    
    region=all.region
    
  }
}
    
