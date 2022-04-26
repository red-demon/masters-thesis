#commands for plotting
# ------------------------------------------------------------------------------------------ 
# 3 plotting 
# ------------------------------------------------------------------------------------------
#plot means vs. means and regression
plot(X.mean,Y.mean,main=paste0("Mean ",Y," vs. Mean ",X),
     xlab=paste0(X," mean anomalies"), ylab=paste0(Y," mean anomalies"))
model.mean=lm(Y.mean~X.mean) #regress means
abline(model.mean$coefficients,col="red")

plotting.range=seq(1:50)

#plot target variable
plot(Y.mean,type="l", main=paste0("Mean ",Y," anomalies, ",season),
     xlab=paste0("year a.d."), ylab=paste0("Mean ", Y," anomalies"))


#plot variances
plot(X.svd.d[1:10],type='l') #plot first 10 eof variances to see which ones are relevant


#Plot First 3 EOFs and PCs
color.palette=colorRampPalette(c('blue','green','white',
                                 'yellow','red'),interpolate='spline')
levels=seq(-0.15,0.15,length.out=61)

#plot EOF1:
mapmat=matrix(X.svd.u[,1],nrow=X.region.dim[2]) #create matrix from 1st PC Time Series on 2d-Grid
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(X.lon, X.lat, -mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="1st EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,1], type='l',main="1. Principal Component",
     xlab="index", ylab="1. PC")
acf(X.svd.v[,1]) #autocorrelation

#plot EOF2:
mapmat=matrix(X.svd.u[,2],nrow=X.region.dim[2]) #create matrix from 1st PC Time Series on 2d-Grid
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(X.lon, X.lat, -mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="2nd EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,2], type='l',main="2. Principal Component",
     xlab="index", ylab="2. PC")
acf(X.svd.v[,1]) #autocorrelation

#plot EOF3:
mapmat=matrix(X.svd.u[,3],nrow=X.region.dim[2]) #create matrix from 1st PC Time Series on 2d-Grid
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(X.lon, X.lat, -mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="3rd EOF"),
               plot.axes={axis(1); axis(2);map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))
#plot PC1:
plot(X.svd.v[1:1000,3], type='l',main="3. Principal Component",
     xlab="index", ylab="3. PC")
acf(X.svd.v[,1]) #autocorrelation



