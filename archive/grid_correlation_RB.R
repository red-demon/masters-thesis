# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment PC Regression simple script:
# Dynamically adjust regional SST or TREFHT mean temperatures using PSL as predictor.
# Data are read in using the Raster Package
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini
# 29.7.2020


# ------------------------------------------------------------------------------------------
# 0.a) Read Relevant packages
# ------------------------------------------------------------------------------------------
library(raster)
library(ncdf4)
library(fields)
library(rworldmap)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/00_data_processing_RB.R")
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/tools/01_data_analysis_RB.R")


# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
region=("Atlantic")

print(region)
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
figpath="/net/h2o/climphys1/rdaenzer/figures/"
outpath="/net/h2o/climphys1/rdaenzer/output/"
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))

# ------------------------------------------------------------------------------------------ 
# 2. Data processing
# ------------------------------------------------------------------------------------------

#crop data to smaller region
if (region=="Atlantic"){
  Y.extent = extent(c(-90, 20,20, 70))
  X.extent = Y.extent + 40
} else if (region=="Pacific") {
  Y.extent = extent(c(140,240, 20, 70))
  X.extent = Y.extent + 40
}

#crop data
X.region.RB = crop(X.data.RB, X.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)

#area weighted RB objects
X.region.area.weighted.RB=area.weighted.RB(X.region.RB)
Y.region.area.weighted.RB=area.weighted.RB(Y.region.RB)
Y.mean=cellStats(Y.region.area.weighted.RB,mean)

#---------------------------------------------------------------------------------------------
# 3. calculate gridwise correlation and plot
#---------------------------------------------------------------------------------------------

#calculate gridpoint correlation between X and Y.mean for all months
X.correlation.RB=cor.mean.gridcells.RB(Y.mean,X.region.area.weighted.RB)

#calculate gridpoint correlation between Y and Y.mean for all months
Y.correlation.RB=cor.mean.gridcells.RB(Y.mean,Y.region.area.weighted.RB)

#calculate gridpoint lag correlation betw. X and Y.mean by month
X.lag.cor.mon.RB = sapply(X=1:12, FUN=function(mon.ix) cor.mean.gridcells.monthly.RB(Y.ts = Y.mean,X.RB=X.region.area.weighted.RB, mon.ix = mon.ix, nr.lags = 12))
names(X.lag.cor.mon.RB)=month.abb    

#calculate gridpoint lag correlation betw.Y and Y.mean by month
Y.lag.cor.mon.RB = sapply(X=1:12, FUN=function(mon.ix) cor.mean.gridcells.monthly.RB(Y.ts = Y.mean,X.RB=Y.region.area.weighted.RB, mon.ix = mon.ix, nr.lags = 12))
names(Y.lag.cor.mon.RB)=month.abb                        



#---------------------------------------------------------------------------------------------
# 4.plotting
#---------------------------------------------------------------------------------------------
setwd(figpath)

#set map
if (region=="Pacific") map="world2"
if (region=="Atlantic") map="world"

#plot contemporaneous correlation map
#--------------------------------------------------------------------------------------------
#convert RB to dataframe
df=data.frame(coordinates(correlation.RB),values(correlation.RB))
colnames(df)=c("lon","lat","correlation")

#plot with ggplot2
xlim=X.extent[1:2]
ylim=X.extent[3:4]
ggplot(data=df, aes(x=lon, y=lat, fill=correlation)) +
  geom_tile()+ scale_fill_distiller(palette = "RdYlBu") + 
  coord_map(projection="mercator",xlim = xlim,ylim=ylim) + borders(map, xlim = xlim, ylim =ylim)+
  ggtitle("Gridwise correlation of PSL and mean SST")


#Plot gridwise correlations for SST and PSL for whole years with filled.contour
{figname=paste0("03_cor_gridwise_",region,".pdf")
pdf(figname)
  
color.palette=colorRampPalette(c('blue','green','white',
                                 'yellow','red'),interpolate='spline')
levels=seq(-1,1,length.out=61)

#plot X.correlation
lon=unique(coordinates(X.correlation.RB)[,1])
lat=unique(coordinates(X.correlation.RB)[,2])[dim(X.correlation.RB)[1]:1]
mapmat=matrix(values(X.correlation.RB),nrow=dim(X.correlation.RB)[2])
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(lon, lat, mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="Gridwise correlation of Atlantic PSL and mean SST"),
               plot.axes={axis(1); axis(2);maps::map(map, add=TRUE);grid()},
               key.title=title(main="Scale"))

#plot Y.correlation
lon=unique(coordinates(Y.correlation.RB)[,1])
lat=unique(coordinates(Y.correlation.RB)[,2])[dim(Y.correlation.RB)[1]:1]
mapmat=matrix(values(Y.correlation.RB),nrow=dim(Y.correlation.RB)[2])
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(lon, lat, mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="Gridwise correlation of Atlantic SST and mean SST"),
               plot.axes={axis(1); axis(2);maps::map(map, add=TRUE);grid()},
               key.title=title(main="Scale"))
dev.off()
}


#plot lag correlation maps for SST
#--------------------------------------------------------------------------------------------
#plotting

for (i in 1:12){
  #calculate monthly lag correlation
  print(month.abb[i])
  figname=paste0("03_lag_cor_SST_gridwise_",region,"_",month.abb[i],".pdf")
  pdf(figname,width = 6*4, height = 5.5*4)
  par(mfrow=c(4, 4))
  pal=colorRampPalette(c("blue","white", "red"))
  
  for (j in 1:13) {
    lag=13-j
    print(lag)
    plot(Y.lag.cor.mon.RB[[i]][[j]],col=pal(60),zlim=c(-0.5,0.5),xlab = "lon", ylab = "lon", main=paste0("Correlation of ",month.abb[i]," mean SST and SST at lag ",lag))
    lines(coasts1)
    lines(coasts2)
    }
  dev.off()
}




#plot lag correlation maps for for PSL
#--------------------------------------------------------------------------------------------
for (i in 1:12){
  #calculate monthly lag correlation
  print(month.abb[i])
  figname=paste0("03_lag_cor_PSL_gridwise_",region,"_",month.abb[i],".pdf")
  pdf(figname,width = 6*4, height = 5.5*4)
  par(mfrow=c(4, 4))
  pal=colorRampPalette(c("blue","white", "red"))
  
  for (j in 1:13) {
    lag=13-j
    print(lag)
    pal=colorRampPalette(c("blue","white", "red"))
    plot(X.lag.cor.mon.RB[[i]][[j]],col=pal(60),zlim=c(-0.5,0.5), xlab = "lon", ylab = "lon", main=paste0("Correlation of ",month.abb[i]," mean SST and PSL at lag ",lag))
    lines(coasts1)
    lines(coasts2)
  }
  dev.off()
}

#color.palette=colorRampPalette(c('blue','green','white',
#                                 'yellow','red'),interpolate='spline')
#levels=seq(-1,1,length.out=61)

#set lat and lon
#lon=unique(coordinates(X.region.area.weighted.RB)[,1])
#lat=unique(coordinates(X.region.area.weighted.RB)[,2])[dim(X.region.area.weighted.RB)[1]:1]

#for (i in 1:12){
  #calculate monthly lag correlation
#  print(month.abb[i])
#  figname=paste0("03_lag_cor_PSL_gridwise_",region,"_",month.abb[i],".pdf")
#  pdf(figname)
#  for (j in 1:13){
#    lag=13-j
#    print(lag)
#    mapmat=matrix(values(X.lag.cor.mon.RB[[i]][[j]]),nrow=dim(X.lag.cor.mon.RB[[i]][[j]])[2])
#    mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
#    filled.contour(lon, lat, mapmat, color.palette=color.palette, levels=levels,
#                   plot.title=title(main=paste0("Correlation of ",month.abb[i]," mean SST and PSL at lag ",lag)),
#                   plot.axes={axis(1); axis(2);maps::map(map, add=TRUE);grid()},
#                   key.title=title(main="cor"))
#    
#  }
#  dev.off()
#}

