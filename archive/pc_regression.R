# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment PC Regression simple script:
# Dynamically adjust regional SST or TREFHT mean temperatures using PSL as predictor.
# Data are read in using the NCDF4 Package
# ----------------------------------------------------------------------------------------------

#------------------
# 0. Read packages
#------------------

library(ncdf4)

#------------------
# 1. Read data
#------------------

region="Atlantic"
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable
season="DJF" # "_DJF", "_MAM", "_JJA", "_SON" or ""

#Read in data
data= nc_open(paste0(X, "_",region, "_anom.nc"))
X.data=ncvar_get(data)
data = nc_open(paste0(Y, "_",region, "_anom.nc"))
Y.data= ncvar_get(data)
lon=data$var$SST$dim[[1]]$vals
lat=data$var$SST$dim[[2]]$vals

#define extent
lat.min=20
lon.min=-70
lat.max=70
lon.max=0
lon.ind=which(lon>lon.min&lon<lon.max)
lat.ind=which(lat>lat.min&lat<lat.max)
lon.region=lon[lon.ind]
lat.region=lat[lat.ind]
X.region=X.data[lon.ind,lat.ind,]
Y.region=Y.data[lon.ind,lat.ind,]



#Convert data arrays to 2d Matrices
# dimensions of the arrays for converting them to matrices
length=length(X.region[1,,1])*length(X.region[,1,1])
width=length(X.region[1,1,])

#initialize 2d matrices
X.matrix=matrix(NA,nrow=length,ncol=width)
Y.matrix=matrix(NA,nrow=length,ncol=width)

#convert data to 2d matrices with rows=gridpoints and columns=observations
for (i in seq(1,width)) {
  X.matrix[,i]=as.vector(X.region[,,i]) # converting the 3d array into a 2d matrix. the algorithm goes trough the matrix column-wise. 
  Y.matrix[,i]=as.vector(Y.region[,,i]) # as columns are latitude, this means that we create a vector latitude-wise
} 

#calculate region mean
X.mean=colMeans(X.matrix,na.rm=TRUE)
Y.mean=colMeans(Y.matrix,na.rm=TRUE)

#plot means vs. means and regression
plot(X.mean,Y.mean,main="Mean SST vs. Mean PLS",
     xlab="PSL anomalies, mean [Pa]", ylab="SST anomalies, mean [K]")
model.mean=lm(Y.mean~X.mean) #regress means
abline(model.mean$coefficients,col="red")


#----------------
# 2. SVD
#----------------
#compute SVDs of predictor data (X)


#adjust area size
X.area.weighted=matrix(NA,nrow=length(X.matrix[,1]),ncol=length(X.matrix[1,]))
phi=lat.region/180*pi
for (i in seq(1,length(lat.region))) {
  idx=((i-1)*length(lon.region)+1):(i*length(lon.region))
  X.area.weighted[idx,]=X.matrix[idx,]*sqrt(cos(phi[i]))
}

#calculate SVD
X.svd=svd(X.area.weighted) #finding SVD of the matrix
X.svd.u=X.svd$u #spatial EOFs of the initial matrix
X.svd.v=X.svd$v # PCs or time series of the EOFs
X.svd.d=X.svd$d # variance of each eof

# standardize variance
X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 


#--------------
# 3. EOF plotting
#--------------
#plot EOFs and PCs of the predictor variable


#plot target variable

plot(1999:2999,Y.mean[1000:2000],type="l", main="Mean SST Anomalies",
     xlab="year", ylab="Mean SST anomalies [K]")

#plot variances of first x EOFs
plot(X.svd.d[1:10],type='l') #plot first 10 eof variances to see which ones are relevant

#Plot First 3 EOFs and PCs
rgb.palette=colorRampPalette(c('blue','green','white',
                               'yellow','red'),interpolate='spline')
int=seq(-0.15,0.15,length.out=61)

#plot EOF1:
mapmat=matrix(X.svd.u[,1],nrow=length(lon.region))
mapmat=mapmat[, seq(length(mapmat[1,]),1)]
filled.contour(lon.region, lat.region, -mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="1st EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,1], type='l',main="1. PC of PSL Data",
     xlab="index", ylab="1. PC")
acf(X.svd.v[,1]) #autocorrelation

#plot EOF2:
mapmat=matrix(X.svd.u[,2],nrow=length(lon.region))
mapmat=mapmat[, seq(length(mapmat[1,]),1)]
filled.contour(lon.region, lat.region, -mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="2nd EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,2], type='l',main="2. PC of PSL Data",
     xlab="index", ylab="2. PC")
acf(X.svd.v[,1]) #autocorrelation

#plot EOF3:
mapmat=matrix(X.svd.u[,3],nrow=length(lon.region))
mapmat=mapmat[, seq(length(mapmat[1,]),1)]
filled.contour(lon.region, lat.region, -mapmat, color.palette=rgb.palette, levels=int,
               plot.title=title(main="3rd EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,3], type='l',main="3. PC of PSL Data",
     xlab="index", ylab="3. PC")
acf(X.svd.v[,1]) #autocorrelation


# ------------------------------------------------------------------------------------------ 
# 4. Regression 
# ------------------------------------------------------------------------------------------
#choose period for the prediction
train.idx=seq(1,12000) #indices for training step
pred.idx=seq(12001,24000) #indices for prediction

#regress temp on SVD PCs
predictor=cbind(X.svd.v[train.idx,1],X.svd.v[train.idx,2],X.svd.v[train.idx,3],X.svd.v[train.idx,4],X.svd.v[train.idx,5])
model.eof=lm(Y.mean[train.idx]~predictor)
plot(model.eof$fitted.values,Y.mean[train.idx],main="Target anomalies vs. Fitted Values from PC Regression",
     xlab="Fitted values from EOF regression [K]", ylab="Mean Y anomalies [K]")
model.fit=lm(Y.mean[train.idx]~model.eof$fitted.values)
abline(model.fit$coefficients,col="red")
cor(model.eof$fitted.values,Y.mean[train.idx])

#test predicted data on subset of data
Y.pred= X.svd.v[pred.idx,1:5] %*% model.eof$coefficients[2:6]
Y.pred= model.eof$coefficients[1]+Y.pred

plot(Y.pred,Y.mean[pred.idx],main="Observed vs. predicted SST anomalies from regression with a subset of the PCs",
     xlab="Fitted values from EOF regression [K]", ylab="Observed Mean SST anomalies [K]")
model.fit=lm(Y.mean[pred.idx]~Y.pred)
abline(model.fit$coefficients,col="red")
correlation=cor(Y.pred,Y.mean[pred.idx])
text(0.1, -0.7, labels = paste("cor=",correlation),col="red")



# ------------------------------------------------------------------------------------------ 
# 5. How does correlation change with increasing EOF as predictors? 
# ------------------------------------------------------------------------------------------
#choose period for the prediction
train.idx=seq(1,12000) #indices for training step
pred.idx=seq(12001,24000) #indices for prediction

#how many eofs to test?
rm(predictor)
rm(model.eof)
eofs=500
correlation.numbereofs=rep(NA,eofs)
predictor=X.svd.v[train.idx,1]


#first iteration regression
model.eof=lm(Y.mean[train.idx]~predictor)
#test predicted data on subset of data
Y.pred = X.svd.v[pred.idx,1] * model.eof$coefficients[2] #predict target variable from regression coefficients for new period
Y.pred = model.eof$coefficients[1]+Y.pred #add intercept
correlation.numbereofs[1]=cor(Y.pred,Y.mean[pred.idx]) #store correlation between observations and predicted values in array

#regress temp on SVD PCs
for (i in seq(2,eofs)){
  print("round:") 
  print(i)
  predictor=cbind(predictor,X.svd.v[train.idx,i]) #add new eof to predictor-array
  model.eof=lm(Y.mean[train.idx]~predictor) #regression based on predictor-array
  
  #test predicted data on subset of data
  Y.pred= X.svd.v[pred.idx,1:i] %*% model.eof$coefficients[2:(i+1)] # multiply subset of eofs with coefficients to obtain predictions
  Y.pred= model.eof$coefficients[1]+Y.pred #add intercept to predictions
  correlation.numbereofs[i]=cor(Y.pred,Y.mean[pred.idx]) # calculate correlation and store in correlation-array
}

#plot correlation in dependence of EOFs
plot(seq(1,eofs),correlation.numbereofs[1:eofs], type='l', main="Correlation between predicted and observed SST for different #EOFS",
     xlab="# EOFs", ylab="Correlation")
correlation.numbereofs[eofs]

