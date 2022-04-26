## ----------------------------------------------------------------------------------------------
# Dynamical Adjustment EOF plotting simple script:
# Plot EOFS of PSL
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
source("/net/h2o/climphys1/rdaenzer/code/dynamical_adjustment/PC_regression/svd_regression.R")

# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
regions=c("Pacific","Atlantic")
period.type="seasonal" #seasonal or monthly

region="Atlantic"

print(region)
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/",region)
figpath="/net/h2o/climphys1/rdaenzer/figures"
outpath="/net/h2o/climphys1/rdaenzer/output"
setwd(path)

X="PSL" #predictor variable
Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))

# ------------------------------------------------------------------------------------------ 
# 2. Define regions, crop data
# ------------------------------------------------------------------------------------------

#crop data to smaller region
if (region=="Atlantic"){
  Y.extent = extent(c(-70, 0, 0, 70))
  X.extent = Y.extent + 40
} else if (region=="Pacific") {
  Y.extent = extent(c(140,240, 0, 70))
  X.extent = Y.extent + 40
}

X.region.RB = crop(X.data.RB, X.extent)
Y.region.RB = crop(Y.data.RB, Y.extent)
X.region.coordinates=coordinates(X.region.RB)
X.lat=unique(X.region.coordinates[,2])[dim(X.region.RB)[1]:1]
X.lon=unique(X.region.coordinates[,1])

Y.region= values(Y.region.RB)

#-------------------------------------------------------------------------------------
# 4. SVD
#-----------------------------------------------------------------------------
#calculate SVD
#adjust area size
area.vec=values(raster::area(X.region.RB))
X.region.area.weighted = values(X.region.RB) * area.vec^2

#calculate SVD
X.svd=svd(X.region.area.weighted) #finding SVD of the matrix
X.svd.u=X.svd$u #spatial EOFs of the initial matrix
X.svd.v=X.svd$v # PCs or time series of the EOFs
X.svd.d=X.svd$d # variance of each eof

# standardize variance
X.svd.d=X.svd.d^2/sum(X.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4 
  
#---------------------------------------------------------------------------------
# 4. Plotting
#---------------------------------------------------------------------------------

#Plot First 3 EOFs and PCs
color.palette=colorRampPalette(c('blue','green','white',
                                 'yellow','red'),interpolate='spline')
levels=seq(-0.15,0.15,length.out=61)

#plot EOF1:
mapmat=matrix(X.svd.u[,1],nrow=dim(X.region.RB)[1]) #create matrix from 1st PC Time Series on 2d-Grid
mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
filled.contour(X.lon, X.lat, -mapmat, color.palette=color.palette, levels=levels,
               plot.title=title(main="1st EOF"),
               plot.axes={axis(1); axis(2);maps::map('world', add=TRUE);grid()},
               key.title=title(main="Scale"))


