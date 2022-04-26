#----------------------------------------------------------------------------------------
# Dynamical Adjustment Script for Plotting
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
# 16.12.2020

# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------


library(fields)
library(ncdf4)
library(rworldmap)
library(hydroGOF)
library(glmnet)
library(SpatialEpi)
library(raster)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------

# other functions to read:
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/02_elasticnet_regression.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/03_dynamical_adjustment.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/gridcorts.R")


#

# ------------------------------------------------------------------------------------------
# 1.a) define regions, ensemble members and paths
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")
ensemble.members=seq(580,980,20)



X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST

for (region in regions){
    
  if (region=="Atlantic") subregions=c("North","Sargasso","Gulfstream","Spain")
  if (region=="Pacific") subregions=c("East","West","Central","NorthNP")

  outpath="/net/h2o/climphys1/rdaenzer/output"
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region)
  
  #define subregion, crop data to regional extent
  for (subregion in subregions){
    if(subregion=="whole") subregion=NULL
    
    print(paste0(region,", ",subregion))


    #define subregion, crop data to regional extent
    Y.extent=define.region.extent(region,subregion)
    X.extent=Y.extent+40
    
    # ------------------------------------------------------------------------------------------
    # 1.b) read and process control data
    # ------------------------------------------------------------------------------------------
    print("01. read and process data")
    
    #read control data
    control.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control")
    setwd(control.path)
    X.control.anom.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
    Y.control.anom.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))
    
    #crop and area weigh data
    X.control.anom.region.RB = area.weight.RB(crop(X.control.anom.data.RB, X.extent)) #crop and area weigh control data
    Y.control.anom.region.RB = crop(Y.control.anom.data.RB, Y.extent)
    Y.control.anom.region.mean=mean.area.weighted.RB(Y.control.anom.region.RB)
    
    # ------------------------------------------------------------------------------------------
    # 1.c) read and process scenario data
    # ------------------------------------------------------------------------------------------
    #ENSMEAN
    member="ENSMEAN"
    scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",member)
    setwd(scenario.path)
    Y.ENSMEAN.anom.data.RB=brick(paste0(Y, "_cesm",member,"_monthly_remap_anom_",region,".nc"))
    Y.ENSMEAN.anom.region.RB=crop(Y.ENSMEAN.anom.data.RB,Y.extent)
    Y.ENSMEAN.anom.region.mean=mean.area.weighted.RB(Y.ENSMEAN.anom.region.RB)
    
    #set up lists for scenariodata
    X.ensemble.anom.data=list(); X.ensemble.anom.region=list()
    Y.ensemble.anom.data=list(); Y.ensemble.anom.region=list(); 
    Y.ensemble.anom.region.mean=c()
  
    
    #read scenario data
    for (i in (1:length(ensemble.members))) {
      #read scenario data
      member=ensemble.members[i]
      print(member)
      scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",member)
      setwd(scenario.path)
      X.ensemble.anom.data=brick(paste0(X, "_cesm",member,"_monthly_anom_",region,".nc"))
      Y.ensemble.anom.data=brick(paste0(Y, "_cesm",member,"_monthly_remap_anom_",region,".nc"))
      
      #crop and area weigh scenario data
      X.ensemble.anom.region[[i]] = area.weight.RB(crop(X.ensemble.anom.data, X.extent))
      Y.ensemble.anom.region[[i]] = crop(Y.ensemble.anom.data, Y.extent)
      Y.ensemble.anom.region.mean=rbind(Y.ensemble.anom.region.mean,mean.area.weighted.RB(Y.ensemble.anom.region[[i]]))
    }
    
    # ------------------------------------------------------------------------------------------ 
    # 02. Calculate Dynamically adjusted Time Series and quantitative values
    # ------------------------------------------------------------------------------------------
    all.region=region
    if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)
    
    setwd(outpath)
    Y.hat.control=readRDS(file=paste0("dynamical_adjustment_Yhat_control_",region,".rds"))
    #Y.hat.ensemble=readRDS(file=paste0("dynamical_adjustment_Yhat_ensemble_",region,".rds"))
    #Y.hat.ensemble.detr=readRDS(file=paste0("dynamical_adjustment_Yhat_detr_ensemble",region,".rds"))
      
      
    date.control = as.Date(substring(text = names(Y.control.anom.region.RB), first = 2, last = 11), "%Y.%m.%d")
    date.scenario = as.Date(substring(text = names(Y.ENSMEAN.anom.region.RB), first = 2, last = 11), "%Y.%m.%d")
    
    pred.months=1:12
    pred.years=1851:2099
    
    Y.residual.control=t(sapply(X=1:length(Y.hat.control),FUN=function(x) unlist(Y.hat.control[[x]])))
    #Y.residual.ensemble=t(sapply(X=1:length(Y.hat.ensemble),FUN=function(x) unlist(Y.hat.control[[x]])))
    #Y.residual.ensemble.detr=t(sapply(X=1:length(Y.hat.ensemble.detr),FUN=function(x) unlist(Y.hat.control[[x]])))
    
    for (i in (1:length(ensemble.members))) {
      idx=which(as.numeric(format(date.scenario, "%Y")) %in% pred.years)
      Y.residual.control[i,idx-12]=Y.ensemble.anom.region.mean[i,idx]-Y.hat.control[[i]][idx-12]
      #Y.residual.ensemble[i,idx-12]=Y.ensemble.anom.region.mean[i,idx]-Y.hat.ensemble[[i]][idx-12]
      #Y.residual.ensemble.detr[i,idx-12]=Y.ensemble.anom.region.mean[i,idx]-Y.hat.ensemble.detr[[i]][idx-12]
      
    }
          
    
    # ------------------------------------------------------------------------------------------ 
    # 03. Plotting
    # ------------------------------------------------------------------------------------------
    #plot 3 different model types: control, ensemble, ensemble.detr
    
    #diff_sd=c()
    
    for (i in seq(1,1,1)){
      if (i==1) {Y.hat=Y.hat.control; Y.residual=Y.residual.control;model.type="control_model"}
      if (i==2) {Y.hat=Y.hat.ensemble;Y.residual=Y.residual.ensemble;model.type="ensemble_model"}
      if (i==3) {Y.hat=Y.hat.ensemble.detr;Y.residual=Y.residual.ensemble.detr;model.type="ensemble_detr_model"}
      
      print(model.type)
      
      fig.width=450;fig.height=400
      plot.months=c(1,7)
      plot.years=1861:2040
      
      if (model.type=="control_model") {
        for (i in 1:length(ensemble.members)) {
          for (plot.month in plot.months) {
            plot.idx=which(as.numeric(format(date.scenario, "%Y")) %in% plot.years  & as.numeric(format(date.scenario, "%m")) %in% plot.month)
            
            #plot prediction and residual and original ts
            setwd(figpath)
            ylim=c(-1,2);yticks=seq(-1,2,0.5); if(!is.null(subregion)) ylim=c(-1.25,2.25); if(!is.null(subregion)) yticks=seq(-1,2,0.5)
            {png(paste0("dynamical_adjustment_prediction_",model.type,"_",region,"_",month.abb[plot.month],"_member",ensemble.members[i],".png"),width=fig.width,height=fig.height)
              par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
              plot(plot.years,Y.ensemble.anom.region.mean[i,plot.idx],ylim=ylim,type="l",axes=F,
                   xlab="year a.d.",ylab="mean SST anomaly [°C]",lwd=1.5,xaxt="n",col=1)
              lines(plot.years,Y.hat[[i]][plot.idx-12],col=3)
              lines(plot.years,Y.residual[i,plot.idx-12],col=2)
              axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
              axis(2,at=yticks)
              abline(h=yticks,v=seq(1860,2040,20),lty="dotted",col="darkgrey")
              legend("topright",paste("member",i),bty="n",cex=1)
              legend("topleft",c("original time series","prediction","residual time series"),lty=1,col=c(1,3,2),bty="y",cex = 1,pt.cex=2)
              dev.off()}
            
            #plot ensmean, residual and original ts
            {png(paste0("dynamical_adjustment_",model.type,"_",region,"_",month.abb[plot.month],"_member",ensemble.members[i],".png"),width=fig.width,height=fig.height)
              par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
              plot(plot.years,Y.ensemble.anom.region.mean[i,plot.idx],col=1,type='l',axes=F,
                   ylim=ylim,xlab="year a.d.",ylab="mean SST anomaly [°C]",xaxt="n")
              lines(plot.years,Y.residual[i,plot.idx-12],col=2)
              lines(plot.years,Y.ENSMEAN.anom.region.mean[plot.idx],col=4,lwd=2)
              axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
              axis(2,at=yticks)
              abline(h=yticks,v=seq(1860,2040,20),lty="dotted",col="darkgrey")
              legend("topright",paste("member",i),bty="n",cex=1)
              legend("topleft",c("original time series","residual time series","ensemble mean"),lty=1,col=c(1,2,4),bty="y",cex = 1,pt.cex=2)
              dev.off()}
            
           }
          }
        }
      
      #plot 3 ensemble members in 1 plot
      for (plot.month in plot.months){
        plot.idx=which(as.numeric(format(date.scenario, "%Y")) %in% plot.years  & as.numeric(format(date.scenario, "%m")) %in% plot.month)
        
        #plot prediction and residual and original ts
        setwd(figpath)
        ylim=c(-0.75,2);yticks=seq(-0.5,2,0.5); if(!is.null(subregion)) ylim=c(-1.25,2.5); if(!is.null(subregion)) yticks=seq(-1,2.5,0.5)
        {png(paste0("dynamical_adjustment_prediction_",model.type,"_",region,"_",month.abb[plot.month],"_3members.png"),width=fig.width,height=fig.height*2)
          par(mfrow=c(3,1),oma=c(5,0,0,0),mar=c(0,5,0,0),bty="n",cex=1.15)
          for (i in 1:3){
            xlab=""; if(i==3) xlab="year a.d."; ylab=""; if(i==2) ylab="mean SST anomaly [°C]"
            plot(plot.years,Y.ensemble.anom.region.mean[i,plot.idx],ylim=ylim,type="l",axes=F,
                 xlab=xlab,ylab=ylab,lwd=1.5,col=1)
            lines(plot.years,Y.hat[[i]][plot.idx-12],col=3)
            lines(plot.years,Y.residual[i,plot.idx-12],col=2)
            if (i==3) axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
            axis(2,at=seq(-0.5,2,0.5))
            abline(h=yticks,v=seq(1860,2040,20),lty="dotted",col="darkgrey")
            legend("topright",paste("member",i),bty="n",cex=1.5)
            if (i==1) legend("topleft",c("original time series","prediction","residual time series"),lty=1,col=c(1,3,2),bty="y",cex = 1,pt.cex=2)
          }
        dev.off()}
        
        #plot ensmean, residual and original ts
        {png(paste0("dynamical_adjustment_",model.type,"_",region,"_",month.abb[plot.month],"_3members.png"),width=fig.width,height=2*fig.height)
          par(mfrow=c(3,1),oma=c(5,0,0,0),mar=c(0,5,0,0),bty="n",cex=1.15)
          for (i in 1:3){
            xlab=""; if(i==3) xlab="year a.d."; ylab=""; if(i==2) ylab="mean SST anomaly [°C]"
            plot(plot.years,Y.ensemble.anom.region.mean[i,plot.idx],col=1,type='l',axes=F,
               ylim=ylim,xlab=xlab,ylab=ylab,lwd=1.5)
            lines(plot.years,Y.residual[i,plot.idx-12],col=2)
            lines(plot.years,Y.ENSMEAN.anom.region.mean[plot.idx],col=4,lwd=2)
            if (i==3) axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
            axis(2,at=seq(-0.5,2,0.5))
            abline(h=yticks,v=seq(1860,2040,20),lty="dotted",col="darkgrey")
            legend("topright",paste("member",i),bty="n")
            if (i==1) legend("topleft",c("original time series","residual time series","ensemble mean"),lty=1,col=c(1,2,4),bty="y",cex = 1,pt.cex=2)
            }
        dev.off()}
      }
      
      # ------------------------------------------------------------------------------------------ 
      # 04.Calculate and visualize range of ensemble and dyn.adj. ensemble
      # ------------------------------------------------------------------------------------------
      setwd(figpath)
      ###TODO: CHANGE MATRIX UNLIST TO NEW INDICES FROM MATRICES
      #calculate mean and sd for shading --> maybe better?
      ensemble.mean=sapply(1:dim(Y.ensemble.anom.region.mean)[2], FUN=function(idx) mean(Y.ensemble.anom.region.mean[,idx]))
      ensemble.sd=sapply(1:dim(Y.ensemble.anom.region.mean)[2], FUN=function(idx) sd(Y.ensemble.anom.region.mean[,idx]))
      residual.control.mean=sapply(1:dim(Y.residual)[2], FUN=function(idx) mean(Y.residual[,idx]))
      residual.control.sd=sapply(1:dim(Y.residual)[2], FUN=function(idx) sd(Y.residual[,idx]))
      mean.residual.sd=mean(residual.control.sd)
      mean.ensemble.sd=mean(ensemble.sd)
      
      #plot with sd and mean
      plot.months=c(1,7)
      plot.years=c(1861:2040)
      
      for (month in plot.months) {
        ylim=c(-0.6,1.7);yticks=seq(-0.5,1.5,0.5); if(!is.null(subregion)) ylim=c(-1,2); if(!is.null(subregion)) yticks=seq(-1,2,0.5)
        plot.idx=which(as.numeric(format(date.scenario, "%Y")) %in% plot.years  & as.numeric(format(date.scenario, "%m")) %in% month)
        #diff_sd=append(diff_sd,(mean.ensemble.sd[plot.idx] - mean.residual.sd[plot.idx])/mean.ensemble.sd[plot.idx])

        {figname=paste0("dynamical_adjustment_ensemble_comparison_",model.type,"_",region,"_",month.abb[month],".png")
          png(figname,width=fig.width,height=fig.height)
          par(oma=c(0,0,0,0),mar=c(5,4,0,0),bty="n")
          plot(plot.years,Y.ENSMEAN.anom.region.mean[plot.idx],type="l",axes=F,ylab="mean SST anomaly [°C]",xlab="year",
               ylim=ylim,xaxt="n")
          polygon(c(plot.years,rev(plot.years)),c((ensemble.mean-ensemble.sd)[plot.idx],rev((ensemble.mean+ensemble.sd)[plot.idx])),col=rgb(0.1,0.1,.1,0.2))
          polygon(c(plot.years,rev(plot.years)),c((residual.control.mean-residual.control.sd)[plot.idx-12],rev((residual.control.mean+residual.control.sd)[plot.idx-12])),col=rgb(1,0,0.3,0.3))
          lines(plot.years,Y.ENSMEAN.anom.region.mean[plot.idx],type="l",lwd=2,col=4)
          lines(plot.years,(ensemble.mean-ensemble.sd)[plot.idx],type="l",col=1)
          lines(plot.years,(ensemble.mean+ensemble.sd)[plot.idx],type="l",col=1)
          lines(plot.years,(residual.control.mean-residual.control.sd)[plot.idx-12],type="l",col=rgb(1,0,0.3,1))
          lines(plot.years,(residual.control.mean+residual.control.sd)[plot.idx-12],type="l",col=rgb(1,0,0.3,1))
          axis(side=1,at=seq(1860,2040,20),labels=seq(1860,2040,20))
          axis(2,at=yticks)
          abline(h=yticks,v=seq(1860,2040,20),lty="dotted",col="darkgrey")
          legend("topleft",c("SD original time series","SD residual time series","ensemble mean"),lty=1,col=c(1,2,4),bty="y",cex = 1,pt.cex=5)
          dev.off()}
        }
      
      
      
      
      
      # ------------------------------------------------------------------------------------------ 
      # 05. calculate quantitative values
      # ------------------------------------------------------------------------------------------
      #calculate signal to noise
      stn.years=1971:2020
      stn.month=1:12
      
      df.months=c()
      df.member=c()
      
      df.cor.pred.scenario.control=c()

      
      df.rmse.scenario.ensmean=c()
      df.rmse.residuals.ensmean.control=c()

      
      df.signal.to.noise.scenario=c()
      df.signal.to.noise.dyn.adj.control=c()

      
      for (member in 1:length(ensemble.members)){
        for (month in c(1,7)){
          #create dataframe data to store as CSV
          df.months=append(df.months,month)
          df.member=append(df.member,member)
          idx=which(as.numeric(format(date.scenario, "%Y")) %in% plot.years  & as.numeric(format(date.scenario, "%m")) %in% month)
          
          df.cor.pred.scenario.control=append(df.cor.pred.scenario.control, round(cor(Y.hat[[member]][idx-12],detrend.ts(Y.ensemble.anom.region.mean[member,idx])),3))
          
          df.rmse.scenario.ensmean=append(df.rmse.scenario.ensmean, round(rmse(Y.ENSMEAN.anom.region.mean[idx],Y.ensemble.anom.region.mean[member,idx]),3))
          df.rmse.residuals.ensmean.control=append(df.rmse.residuals.ensmean.control, round(rmse(Y.ENSMEAN.anom.region.mean[idx],Y.residual[member,idx-12]),3))
       
          df.signal.to.noise.scenario=append(df.signal.to.noise.scenario, round(signal.to.noise(Y.ensemble.anom.region.mean[member,idx],stn.years,month),3))
          df.signal.to.noise.dyn.adj.control=append(df.signal.to.noise.dyn.adj.control, round(signal.to.noise(Y.residual[member,idx-12],stn.years,month),3))
        
        }
      }
      
      csv.df=data.frame(df.member,df.months,
                        df.cor.pred.scenario.control, 
                        df.rmse.scenario.ensmean,df.rmse.residuals.ensmean.control,
                        df.signal.to.noise.scenario,df.signal.to.noise.dyn.adj.control)
      setwd(figpath)
      write.csv(csv.df,paste0("dynamical_adjustment_values_",model.type,"_",region,".csv"))
      
    }
    #names(diff_sd)=c("control.model.jan","control.model.jul","ensemble.model.jan","ensemble.model.jul","ensemble.detr.model.jan","ensemble.detr.model.jul")
    #write.csv(diff_sd,paste0("change_in_sd_",region,".csv"))
    
    region=all.region
  }
}
