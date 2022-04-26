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

region="Atlantic"

#set directories
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
figpath=file.path("/net/h2o/climphys1/rdaenzer/figures/PC_regression/",region)
outpath="/net/h2o/climphys1/rdaenzer/output/"
setwd(path)

# ------------------------------------------------------------------------------------------ 
# 00. Read processed data and SVDs
# ------------------------------------------------------------------------------------------

Y="SST" #target variable: TREFHT or SST


#read ncdf files into RB objects
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))

Y.extent=define.region.extent(region,subregion)
X.extent=Y.extent+40
Y.region.RB = crop(Y.data.RB, X.extent)


#area weighting and centering

#calculate total SVDs and store in separate objects
Y.region.svd=values(area.weight.RB(Y.region.RB,svd=TRUE))

NONA.idx = which(apply(X = Y.region.svd, MARGIN=1, FUN=function(x) !all(is.na(x))))             # Get NA values in Y raster
Y.RB.values.nona = Y.region.svd[NONA.idx,]      

Y.svd=svd(Y.RB.values.nona)
Y.svd.u=Y.svd$u
Y.svd.v=Y.svd$v
Y.svd.d=Y.svd$d
Y.svd.d=Y.svd.d^2/sum(Y.svd.d^2) #standardizing variance as written in Shen Tutorial, Chapter 4

#---------------------------------------------------------------------------------
# 03. Plotting leading EOFs and corresponding PCs
#---------------------------------------------------------------------------------
print("03")
setwd(figpath)
#worldmaps
if(region=="Pacific") {map="world2";fig.height=350;fig.width=510}
if(region=="Atlantic") {map="world";fig.height=350;fig.width=410}

#plot total EOFs 1-5:
Y.spatial.raster=raster(Y.region.RB,1)
for (i in 1:20){
figname=paste0("03_EOF_",paste(i,X,region,sep="_"),".png")
#png(figname,height=fig.height,width=fig.width)
par(mar=c(2,2,2,2),oma=c(0,0,0,0))
values(Y.spatial.raster)[NONA.idx]=Y.svd$u[,5]
plot(Y.spatial.raster,col=color.warm.cold,zlim=c(-0.2,0.2),main=paste0(i,". EOF"));maps::map(map,add=T,interior=F)

#dev.off()
}
