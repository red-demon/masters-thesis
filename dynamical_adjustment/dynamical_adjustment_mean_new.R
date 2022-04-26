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
# 16.12.2020


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


#

# ------------------------------------------------------------------------------------------
# 1.a) define regions, ensemble members and paths
# ------------------------------------------------------------------------------------------
regions=c("Atlantic","Pacific")

for (region in regions){
  
  if (region=="Atlantic") subregions=c("Sargasso","Gulfstream","Spain")
  if (region=="Pacific") subregions=c("East","West","Central","NorthNP")
  ensemble.members=seq(580,980,20)
  
  outpath="/net/h2o/climphys1/rdaenzer/output"
  figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/dynamical_adjustment",region)

  X="PSL" #predictor variable
  Y="SST" #target variable: TREFHT or SST
  
  #define subregion, crop data to regional extent
  for (subregion in subregions){
    if(subregion=="whole") subregion=NULL
    if(!(subregion %in% subregions)) subregion=NULL
    
    print(paste0(region,", ",subregion))
    
    #define subregion, crop data to regional extent
    Y.extent=define.region.extent(region,subregion)
    X.extent=Y.extent+40
    
    # ------------------------------------------------------------------------------------------
    # 1.b) read and process control data
    # ------------------------------------------------------------------------------------------
    print("01. read and process data")
    
    #read control data
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
    #set up lists for data
    X.ensemble.anom.data=list(); X.ensemble.anom.region=list()
    Y.ensemble.anom.data=list(); Y.ensemble.anom.region=list(); 
    Y.ensemble.anom.region.mean=list()
    X.detr.ensemble.anom.region=list()
    Y.detr.ensemble.anom.region.mean=list()
  
    #read scenario data
    for (i in 1:length(ensemble.members)) {
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
      Y.ensemble.anom.region.mean[[i]]=mean.area.weighted.RB(Y.ensemble.anom.region[[i]])
      
      #detrend scenario data
      #X.detr.ensemble.anom.region[[i]]=detrend.RB(X.ensemble.anom.region[[i]])
      #Y.detr.ensemble.anom.region.mean[[i]] = detrend.ts(Y.ensemble.anom.region.mean[[i]])
    }
    
    # ------------------------------------------------------------------------------------------ 
    # 02. Dynamical Adjustment of mean SST
    # ------------------------------------------------------------------------------------------
    print("02. Dynamical Adjustment")
    
    all.region=region
    if(!is.null(subregion)) region=subregion #subregion is now region (important for plot names)
    
    #train model using control data
    #----------------------------------------------------------
    months=1:12
    lags=1:12;s.lags=1:6;w.lags=1:12
    train.years=1001:2000
    pred.years=1851:2099
    
    Y.hat=list()
    cv.model=dynamical.adjustment.gridcell.monthly(Y.train=Y.control.anom.region.mean,X.train=X.control.anom.region.RB,train.years=train.years,pred.years=pred.years,
                                                   alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,ret.cv.model=T)
   
    for (i in 1:length(ensemble.members)){
      print(ensemble.members[i])
      Y.hat[[i]]=dynamical.adjustment.gridcell.monthly(model.cv.glmnet=cv.model,X.pred=X.ensemble.anom.region[[i]],train.years=train.years,pred.years=pred.years,
                                                       alpha=0.1,lags=NULL,w.lags=w.lags,s.lags=s.lags,months=months)
    }
    
    setwd(outpath)
    saveRDS(cv.model,file=paste0("dynamical_adjustment_model_control_",region,".rds"))
    saveRDS(Y.hat,file=paste0("dynamical_adjustment_Yhat_control_",region,".rds"))
    
    
    #train model using detrended or no-detrended scenario data
    #----------------------------------------------------------
    #months=1:12
    #lags=1:12;s.lags=1:6;w.lags=1:12
    #train.years=1851:2000
    #pred.years=1851:2099
    
    
    #cv.model=list();cv.model.detr=list()
    #Y.hat=list();Y.hat.detr=list()
    
    #for (i in 1:length(ensemble.members)){
    #  print(ensemble.members[i])
    #  j=(i%%length(ensemble.members))+1 #current member + 1, member 21 --> 1
      
    #  cv.model[[j]]=dynamical.adjustment.gridcell.monthly(Y.train=Y.ensemble.anom.region.mean[[j]],X.pred=X.ensemble.anom.region[[i]],X.train=X.ensemble.anom.region[[j]],train.years=train.years,pred.years=pred.years,
    #                                           alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,ret.cv.model=T,detrend.X=F,detrend.Y=T)
    #  Y.hat[[i]]=dynamical.adjustment.gridcell.monthly(model.cv.glmnet=cv.model[[j]],X.pred=X.ensemble.anom.region[[i]],train.years=train.years,pred.years=pred.years,
    #                                                      alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,detrend.X=F,detrend.Y=T)
  
     # cv.model.detr[[j]]=dynamical.adjustment.gridcell.monthly(Y.train=Y.ensemble.anom.region.mean[[j]],X.pred=X.ensemble.anom.region[[i]],X.train=X.ensemble.anom.region[[j]],train.years=train.years,pred.years=pred.years,
       #                                     alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,ret.cv.model=T,detrend.X=T,detrend.Y=T)
      #Y.hat.detr[[i]]=dynamical.adjustment.gridcell.monthly(model.cv.glmnet=cv.model.detr[[j]],X.pred=X.ensemble.anom.region[[i]],train.years=train.years,pred.years=pred.years,
        #                                                       alpha=0.1,lags=NULL,s.lags=s.lags,w.lags=w.lags,months=months,detrend.X=T,detrend.Y=T)
    #}
  
    #setwd(outpath)
    #saveRDS(cv.model,file=paste0("dynamical_adjustment_model_ensemble_",region,".rds"))
    #saveRDS(cv.model.detr,file=paste0("dynamical_adjustment_model_detr_ensemble",region,".rds"))
    #saveRDS(Y.hat,file=paste0("dynamical_adjustment_Yhat_ensemble_",region,".rds"))
    #saveRDS(Y.hat.detr,file=paste0("dynamical_adjustment_Yhat_detr_ensemble",region,".rds"))
  
    #set region back to all.region for next subregion round
    region=all.region
  }
}
  
  
