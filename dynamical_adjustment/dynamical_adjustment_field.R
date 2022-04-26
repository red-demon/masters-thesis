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

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
# other functions to read:
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/gridcorts.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/02_elasticnet_regression.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/03_dynamical_adjustment.R")

source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/colormaps.R")

# ------------------------------------------------------------------------------------------
# 1.) Read in Data and preprocess
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")
subregion=NULL
for (region in regions){
  print(region)
  outpath="/net/h2o/climphys1/rdaenzer/output"
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region)
  
  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  #define subregion, crop data to regional extent
  Y.extent=define.region.extent(region,subregion)
  X.extent=Y.extent+40
  
  #read ncdf files into RB objects
  #control data
  control.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control")
  setwd(control.path)
  X.control.data.RB = raster::brick(paste0(X, "_",region,"_anom.nc")) 
  Y.control.data.RB = raster::brick(paste0(Y,"_",region,"_anom.nc"))
  
  #crop data
  X.control.region.RB = raster::crop(X.control.data.RB, X.extent)
  Y.control.region.RB = raster::crop(Y.control.data.RB, Y.extent)
  
  #scenario data
  ensemble.members=seq(580,980,20)
  X.scenario.region.RB=list()
  Y.scenario.region.RB=list()
  for(i in 1:length(ensemble.members)){
    scenario.path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/scenario/",ensemble.members[i])
    setwd(scenario.path)
    X.scenario.data.RB=raster::brick(paste0(X, "_cesm",ensemble.members[i],"_monthly_anom_",region,".nc"))
    Y.scenario.data.RB=raster::brick(paste0(Y, "_cesm",ensemble.members[i],"_monthly_remap_anom_",region,".nc"))
    X.scenario.region.RB[[i]] = (raster::crop(X.scenario.data.RB, X.extent))
    Y.scenario.region.RB[[i]] = (raster::crop(Y.scenario.data.RB, Y.extent))
    if (Y=="TREFHT") Y.scenario.data.RB=raster::brick(paste0(Y, "_cesm",ensemble.members[i],"_monthly_anom_",region,".nc"))
  }
  
  #spatial raster
  spatial.raster=raster(Y.scenario.region.RB[[1]],layer=1)
  values(spatial.raster)=NA
  
  # ------------------------------------------------------------------------------------------ 
  # 02. Dynamical Adjustment of whole grid
  # ------------------------------------------------------------------------------------------
  print("02 DA of grid")
  months=c(1,7)
  lags=1:12;s.lags=1:6;w.lags=1:12
  train.years=1001:2000
  pred.years=1861:2040
  
  
  for (month in months){
    
    #train  and save models (aggregated and non-aggregated)
    #cv.model=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,X.RB=X.scenario.region.RB[[1]],X.train.RB=X.control.region.RB,
    #                                       train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
     #                                     months=month,alpha=0.1,s="lambda.min",x.domain=20,y.domain=20,ret.cv.model=T)
    #cv.model.aggregated=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,X.RB=aggregate(X.scenario.region.RB[[1]],2),X.train.RB=aggregate(X.control.region.RB,2),
    #                                                train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
    #                                               months=month,alpha=0.1,s="lambda.min",x.domain=40,y.domain=40,ret.cv.model=T)
    
    setwd(outpath)
    #saveRDS(cv.model,file=paste0("field_dynamical_adjustment_cvmodel",region,"_",month.abb[month],".rds"))
    #saveRDS(cv.model.aggregated,file=paste0("field_dynamical_adjustment_cvmodel",region,"_",month.abb[month],"_aggregated.rds"))
    
    cv.model=readRDS(file=paste0("field_dynamical_adjustment_cvmodel",region,"_",month.abb[month],".rds"))
    cv.model.aggregated=readRDS(file=paste0("field_dynamical_adjustment_cvmodel",region,"_",month.abb[month],"_aggregated.rds"))
    
    #run dynamical adjustment prediction for ensembles
    for (i in 1:length(ensemble.members)){
      
      ensemble.member=ensemble.members[i]
      print(ensemble.member)
      X.RB=X.scenario.region.RB[[i]]

      Y.hat.RB=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,model.cv.glmnet=cv.model,X.RB=X.RB,
                                             train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
                                             months=month,alpha=0.1,s="lambda.min",x.domain=20,y.domain=20)
      Y.hat.RB.aggregated=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,model.cv.glmnet=cv.model.aggregated,X.RB=aggregate(X.RB,2),
                                                        train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
                                                        months=month,alpha=0.1,s="lambda.min",x.domain=40,y.domain=40)
      
      #extract model data and compare to predictions
      Y.model.RB=extract.period.RB(Y.scenario.region.RB[[i]],years=pred.years, months=month)
      Y.model.RB.detr=detrend.RB(Y.model.RB)
      test.cor = gridcorts(rasterstack = brick(list(Y.hat.RB, Y.model.RB.detr)), method = "pearson", type = "corel")
      test.cor.aggregated = gridcorts(rasterstack = brick(list(Y.hat.RB.aggregated, Y.model.RB.detr)), method = "pearson", type = "corel")
      
      setwd(outpath)
      saveRDS(test.cor,file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,".rds"))
      saveRDS(test.cor.aggregated,file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,"_aggregated.rds"))
      
    }
    
    #run prediction for model, aggregated and non-aggregated
    #ensemble.member="control"
    #print(ensemble.member)
    #X.RB=X.control.region.RB
    #run dynamical adjustment
    #Y.hat.RB=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,model.cv.glmnet=cv.model,X.RB=X.RB,
    #                                       train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
    #                                       months=month,alpha=0.1,s="lambda.min",x.domain=20,y.domain=20)
    #Y.hat.RB.aggregated=dynamical.adjustment.field.RB(Y.RB=Y.control.region.RB,model.cv.glmnet=cv.model.aggregated,X.RB=aggregate(X.RB,2),
    #                                                  train.years=train.years,pred.years=pred.years,lags=NULL,w.lags=w.lags,s.lags=s.lags,
    #                                                  months=month,alpha=0.1,s="lambda.min",x.domain=40,y.domain=40)
    
    #extract model data and compare to predictions
    #Y.model.RB=extract.period.RB(Y.control.region.RB,years=pred.years, months=month)
    #test.cor = gridcorts(rasterstack = brick(list(Y.hat.RB, Y.model.RB)), method = "pearson", type = "corel")
    #test.cor.aggregated = gridcorts(rasterstack = brick(list(Y.hat.RB.aggregated, Y.model.RB)), method = "pearson", type = "corel")
    
    #setwd(outpath)
    #saveRDS(test.cor,file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,".rds"))
    #saveRDS(test.cor.aggregated,file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,"_aggregated.rds"))
  }
}
# ------------------------------------------------------------------------------------------ 
# 03. plotting of field
# ------------------------------------------------------------------------------------------

for(region in regions){
  for(month in c(1,7)){
    df.test.cor=c()
    df.test.cor.aggregated=c()
    
    for (ensemble.member in c(ensemble.members)){
      setwd(outpath)
      print(ensemble.member)
      test.cor=readRDS(file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,".rds"))
      if (ensemble.member=="control") test.cor.aggregated= readRDS(file=paste0("test_cor_",region,"_",month.abb[month],"_member",ensemble.member,"_aggregated.rds"))
      
      if(region=="Pacific") {map="world2";fig.height=350;fig.width=510}
      if(region=="Atlantic") {map="world";fig.height=350;fig.width=450}
      
      #plot and save
      setwd(file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region))
      figname=paste0(Y,"_field_dynamical_adjustment_",region,"_",month.abb[month],"_member",ensemble.member,".png")
      {png(figname,height=fig.height,width=fig.width)
        par(mar=c(2,2,1,6),oma=c(0,0,0,0))
        image(test.cor,col=color.YlOrBr,zlim=c(0,1))
        maps::map(map,add=T,fill=T)
        image.plot(zlim=c(0,1),col=color.YlOrBr,legend.only=T)
        dev.off()}
      
      if(ensemble.member=="control"){
        #plot and save
        setwd(file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region))
        figname=paste0(Y,"_field_dynamical_adjustment_",region,"_",month.abb[month],"_member",ensemble.member,"_aggregated.png")
        {png(figname,height=fig.height,width=fig.width)
          par(mar=c(2,2,1,6),oma=c(0,0,0,0))
          image(test.cor.aggregated,col=color.YlOrBr,zlim=c(0,1))
          maps::map(map,add=T,fill=T)
          image.plot(zlim=c(0,1),col=color.YlOrBr,legend.only=T)
          dev.off()}
      }
      
      df.test.cor=cbind(df.test.cor,values(test.cor))
      df.test.cor.aggregated=cbind(df.test.cor.aggregated,values(test.cor.aggregated))
    }
    
    colnames(df.test.cor)=c("control",ensemble.members); colnames(df.test.cor.aggregated)=c("control",ensemble.members)
    test.cor.mean=rowMeans(df.test.cor[,2:length(df.test.cor[1,])])
    test.cor.mean.aggregated=rowMeans(df.test.cor.aggregated[,2:length(df.test.cor.aggregated[,1])])
    
    
    #plot and save mean test cor
    spatial.raster=test.cor
    setwd(file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region))
    figname=paste0(Y,"_field_dynamical_adjustment_",region,"_",month.abb[month],"mean_correlation.png")
    values(spatial.raster)=test.cor.mean
    {png(figname,height=fig.height,width=fig.width)
      par(mar=c(2,2,1,6),oma=c(0,0,0,0))
      image(spatial.raster,col=color.YlOrBr,zlim=c(0,1))
      maps::map(map,add=T,fill=T)
      image.plot(zlim=c(0,1),col=color.YlOrBr,legend.only=T)
      dev.off()}
    
    #plot and save mean test cor aggregated
    spatial.raster=test.cor.aggregated
    setwd(file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region))
    figname=paste0(Y,"_field_dynamical_adjustment_",region,"_",month.abb[month],"mean_correlation_aggregated.png")
    values(spatial.raster)=test.cor.mean.aggregated
    {png(figname,height=fig.height,width=fig.width)
      par(mar=c(2,2,1,6),oma=c(0,0,0,0))
      image(spatial.raster,col=color.YlOrBr,zlim=c(0,1))
      maps::map(map,add=T,fill=T)
      image.plot(zlim=c(0,1),col=color.YlOrBr,legend.only=T)
      dev.off()}
    
  }
}
  

