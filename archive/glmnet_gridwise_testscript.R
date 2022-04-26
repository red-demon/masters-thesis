
# ----------------------------------------------------------------------------------------------
# Dynamical Adjustment simple script:
# ----------------------------------------------------------------------------------------------

# Rafael Bonafini, adapted from S. Sippel
# 20.8.2020


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

library(foreach)
library(doParallel)
library(bigmemory)
registerDoParallel(cores = 4)

# ------------------------------------------------------------------------------------------
# 0.b) Read Relevant Functions
# ------------------------------------------------------------------------------------------
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_00_preprocessing_4DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_00_preprocessing_SPACETIME_4DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_01_GRIDCELL_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_01_SPACE_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/dynamical_adjustment_elasticnet/_02_ANALYSISFUN_DYNADJ_RB.R")

# other functions to read:
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/convert.to.eurocentric.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/frenchcolormap.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/project_raster.R")
source("/net/h2o/climphys1/rdaenzer/code/code_sippels/tools/gridcorts.R")

#set directories
path=file.path("/net/h2o/climphys1/rdaenzer/data/cesm122/monthly_control/")
figpath="/net/h2o/climphys1/rdaenzer/figures/"
outpath="/net/h2o/climphys1/rdaenzer/output/"
setwd(path)

# ------------------------------------------------------------------------------------------ 
# 1. Read model data from CESM1.2.2 
# ------------------------------------------------------------------------------------------
setwd(path)
region="Pacific"

#read ncdf files into RB objects
X.data.RB = brick(paste0(X, "_",region,"_anom.nc")) 
Y.data.RB = brick(paste0(Y,"_",region,"_anom.nc"))

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
Y.region.RB= crop(Y.data.RB,Y.extent)

#area weighted RB objects
X.region.area.weighted.RB=area.weighted.RB(X.region.RB)
Y.region.area.weighted.RB=area.weighted.RB(Y.region.RB)


# ------------------------------------------------------------------------------------------ 
# 2. TRAIN MODEL FOR AN ENTIRE TEMPERATURE FIELD:
# ------------------------------------------------------------------------------------------ 

# run dynamical adjustment on the whole field:
Y.region.area.weighted.hat.RB = train.dyn.adj.elastnet.annual.RB(Y.RB = Y.region.area.weighted.RB, X.RB = X.region.area.weighted.RB, train.years = 1000:2000, train.months = 1, add.mon = 1, alpha = 0.2, nfolds = 10, s = "lambda.1se", 
                                                  lags.to.aggregate = list(1), n.pc = 10, nr.cores = 1, x.domain = 20, y.domain = 20)

Y.region.area.weighted.pred.RB = extract.period.RB(RB = Y.region.area.weighted.hat.RB, years = 2001:3000, months = 1)
Y.region.area.weighted.obs.RB = extract.period.RB(RB = Y.region.area.weighted.RB, years = 2001:3000, months = 1)


test.cor = gridcorts(rasterstack = brick(list(Y.region.area.weighted.pred.RB, Y.region.area.weighted.obs.RB)), method = "pearson", type = "corel")
setwd(outpath)
saveRDS(test.cor,paste0("04_glmnet_gridwise_corr_",region,".RDS"))

## Plot correlation between prediction and true values:
setwd(figpath)
pdf(paste0("04_glmnet_gridwise_corr_",region,".pdf"))
plot(test.cor)
lines(coastsCoarse)
dev.off()

#plot with filled.contour
color.palette=colorRampPalette(c('blue','green','white',
                                 'yellow','red'),interpolate='spline')
levels=seq(-1,1,length.out=61)

#set lat and lon
lon=unique(coordinates(test.cor)[,1])
lat=unique(coordinates(test.cor)[,2])[dim(test.cor)[1]:1]

#calculate monthly lag correlation
  print(month.abb[i])
  figname=paste0("04_glmnet_gridwise_corr_",region,"_",month.abb[i],".pdf")
  pdf(figname)
  mapmat=matrix(values(test.cor),nrow=dim(test.cor)[2])
  mapmat=mapmat[, seq(length(mapmat[1,]),1)] 
  filled.contour(lon, lat, mapmat, color.palette=color.palette, levels=levels,
                 plot.title=title(main="Gridwise correlation of pred. SST (elasticnet) and obs. SST"),
                 plot.axes={axis(1); axis(2);maps::map('world2', add=TRUE);grid()},
                 key.title=title(main="cor"))
  
  dev.off()

