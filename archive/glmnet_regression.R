
# ----------------------------------------------------------------------------------------------
# Elastic Net Regression script:
# Dynamically adjust regional SST or TREFHT mean temperatures using PSL as predictor with an
# elastic net regression.
# Data are read in using the ncdf4 Package
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 1.7.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------

library(ncdf4)
library(raster)
library(chron)


# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------

setwd("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")

#Read in data
PSL = nc_open("PSL_Europe_2000y_DJF-seasonal_anom.nc")
PSL_EUROPE=ncvar_get(PSL)
TREFHT = nc_open("TREFHT_Europe_2000y_DJF-seasonal_anom.nc")
TREFHT_EUROPE= ncvar_get(TREFHT)
lon=PSL$var$PSL$dim[[1]]$vals
lat=PSL$var$PSL$dim[[2]]$vals
Time<- PSL$dim$time$vals
tmdy=as.Date(Time,origin=as.Date("0001-01-01"))

#define extent
crop=TRUE
lat_min=30
lon_min=-30
lat_max=70
lon_max=30

#crop data
if (crop) 
{
  lon_ind=which(lon>lon_min&lon<lon_max)
  lat_ind=which(lat>lat_min&lat<lat_max)
  
  #create grid from lat_regional & lon_regional
  lon_region=lon[lon_ind]
  lat_region=lat[lat_ind]
  #  precst=matrix(0,nrow=length(lon_region),ncol=length(lat_region))
  #  grid=cbind(lat_region, lon_region, precst)
  
  #crop PSL and TREFHT to desired region
  PSL_region=PSL_EUROPE[lon_ind,lat_ind,]
  TREFHT_region=TREFHT_EUROPE[lon_ind,lat_ind,]
}else{
  #create grid from lat and lon
  precst=matrix(0,nrow=length(lon),ncol=length(lat))
  grid=cbind(lat, lon, precst)
  
  PSL_region=PSL_EUROPE
  TREFHT_region=TREFHT_EUROPE
}

#calculate means of temperature and PSL
TREFHT_mean=colMeans(TREFHT_matrix)
PSL_mean=colMeans(PSL_matrix)
# ------------------------------------------------------------------------------------------ 
# 2. predict temperature using glmnet 
# ------------------------------------------------------------------------------------------

#convert to 2d matrices

# dimensions of the arrays for converting them to matrices
length=length(PSL_region[1,,1])*length(PSL_region[,1,1])
width=length(PSL_region[1,1,])

#initialize 2d matrices
PSL_matrix=matrix(NA,nrow=length,ncol=width)
TREFHT_matrix=matrix(NA,nrow=length,ncol=width)

#convert data to 2d matrices with rows=gridpoints and columns=observations
for (i in seq(1,width)) {
  PSL_matrix[,i]=as.vector(PSL_region[,,i]) # converting the 3d array into a 2d matrix. the algorithm goes trough the matrix column-wise. 
  TREFHT_matrix[,i]=as.vector(TREFHT_region[,,i]) # as columns are latitude, this means that we create a vector latitude-wise
} 

#prediction using glmnet
train.idx=seq(1,1000)
pred.idx=seq(1001,2000)
model.glmnet=glmnet(t(PSL_matrix[,train.idx]),TREFHT_mean[train.idx])
TREFHT_hat=predict(model.glmnet,t(PSL_matrix[,pred.idx]))
plot(TREFHT_mean[pred.idx],TREFHT_hat[,100])
cor(TREFHT_mean[pred.idx],TREFHT_hat[,100])

#prediction using cv.glmnet
correlation=rep(NA,10)
train.idx=seq(1,1000)
pred.idx=seq(1001,2000)
alpha=seq(0,1,0.1)

for (i in seq(1,length(alpha))){
  print("round:")
  print(i)
  model.cv.glmnet=cv.glmnet(t(PSL_matrix[,train.idx]),TREFHT_mean[train.idx], alpha=alpha[i])
  TREFHT_hat=predict(model.cv.glmnet,t(PSL_matrix[,pred.idx]))
  correlation[i]=cor(TREFHT_mean[pred.idx],TREFHT_hat)
}
plot(TREFHT_mean[pred.idx],TREFHT_hat)
plot(alpha,correlation, type='l')
plot(model.cv.glmnet)
