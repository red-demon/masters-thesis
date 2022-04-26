# ---------------------------------------------------------
# Function Repository for Dynamical Adjustment Elasticnet Regression 
# Rafael Bonafini 
# 04.08.2020
#-------------------------------------------------------------

require(zoo)
require(xts)
require(raster)


# remove NA from RB
# to do: create an XTS as returned object
rm.NA.from.RB <-function(RB) {
  
  RB.values=values(RB)
  
  NONA.idx = which(apply(X = RB.values, MARGIN=1, FUN=function(x) !all(is.na(x))))          
  RB.values.nona = RB.values[NONA.idx,]                                                  
  RB.coord.nona = coordinates(RB)[NONA.idx,]
  ret=list(RB.values.nona,RB.coord.nona)
  names(ret)=list("values","coordinates")
  return(ret)
}


#define regions for data processing
define.region.extent <-function(region,subregion=NULL){
  #region
  atlantic.extent=extent(c(-70,0,20,60))
  #atlantic subregions
  north.extent=extent(c(-60,-10,50,60))
  gulfstream.extent=extent(c(-70,-40,35,45))
  sargasso.extent=extent(c(-70,-50,27.5,32.5))
  spain.extent=extent(c(-30,-10,20,35))
  
  #region
  pacific.extent=extent(c(140,240,20,60))
  #pacific subregions
  west.extent=extent(c(140,180,20,40))
  east.extent=extent(c(210,240,20,40))
  central.extent=extent(c(180,210,20,40))
  northnp.extent=extent(c(160,230,45,60))
  
  #
  #crop data to smaller region
  if (region=="Atlantic"){
    if(is.null(subregion)) Y.extent = atlantic.extent
    else if (subregion=="North") Y.extent=north.extent
    else if(subregion=="Gulfstream") Y.extent=gulfstream.extent
    else if(subregion=="Sargasso") Y.extent=sargasso.extent
    else if(subregion=="Spain") Y.extent=spain.extent
  } else if (region=="Pacific"){
    if(is.null(subregion)) Y.extent = pacific.extent
    else if(subregion=="West") Y.extent = west.extent
    else if(subregion=="East") Y.extent = east.extent
    else if(subregion=="Central") Y.extent = central.extent
    else if(subregion=="NorthNP") Y.extent = northnp.extent
  }
  return(Y.extent)
}

#function to extract time lagged monthly or seasonal predictor data from monthly data
#used with extract.period.RB function or period.mean.RB functio
#---------------------------------------------------------------------------------
extract.lag.period.RB <- function(RB,years=NULL,month=NULL,lag=NULL) {
  #RB: raster::brick object
  #month: reference month for lag
  #lag: time lag in number of months (e.g. 1)

    #get years
  date.seq = as.Date(substring(text = names(RB), first = 2, last = 11), "%Y.%m.%d")
  if(is.null(years)) years = unique(as.numeric(format(date.seq, "%Y")))
  if(is.null(month)) month = unique(as.numeric(format(date.seq, "%m")))
  
  
  if (!is.numeric(month)) stop("'month' requires numeric argument if period='monthly'.")
    
    #shift years if necessary
    res.month=month-lag
    
    if (res.month<1) {
      shift.year=-floor((res.month-1)/12) # years shifted if res.month is 0 or negative (= change of year)
    } else {
      shift.year=0
    }
    
    # define years and period to be extracted
    years=years[1:(length(years)-shift.year)]
    months=((res.month+11)%%12)+1 #values go only from 1 to 12

    time.idx = which(as.numeric(format(date.seq, "%Y")) %in% years & as.numeric(format(date.seq, "%m")) %in% months)
    
    lag.RB=subset(RB, time.idx)
    
  return(lag.RB)
}

#function to extract time lagged months from monthly data from xts-object
#---------------------------------------------------------------------------------
extract.lag.period.ts <- function(ts,month,lag=0) {
  #ts: ts object
  #month: reference month for lag
  #lag: time lag in number of months (e.g. 1)
  
  #get years
  date.seq = as.Date(substring(text = names(ts), first = 2, last = 11), "%Y.%m.%d")
  years = unique(as.numeric(format(date.seq, "%Y")))
  
    #shift years if necessary
  res.month=month-lag
  
  if (res.month<1) {
    shift.year=-floor((res.month-1)/12) # years shifted if res.month is 0 or negative (= change of year)
  } else {
    shift.year=0
  }
  
  # define years and period to be extracted
  years=years[1:(length(years)-shift.year)]
  months=((res.month+11)%%12)+1 #values go only from 1 to 12
  
  time.idx = which(as.numeric(format(date.seq, "%Y")) %in% years & as.numeric(format(date.seq, "%m")) %in% months)
  
  #extract period for X
  lag.ts=ts[time.idx]
  
  return(lag.ts)
}

#fix by changing start.idx to something else
#extract month from array using indices, not dates
extract.lag.period.idx <- function(matrix,month.idx,lag=0) {

  start.idx=month.idx+24
  
  if (is.null(dim(matrix))) { #check if matrix is just 1-dimensional array or 2-dimensional matrix
    Y.idx=seq(start.idx, length(matrix), 12)
    X.idx=Y.idx-lag
    matrix.lag=matrix[X.idx]
  } else {
    Y.idx=seq(start.idx, length(matrix[1,]), 12)
    X.idx=Y.idx-lag
    matrix.lag=matrix[,X.idx]
  }
  
  return(matrix.lag)
}
  
#extract period from RB
extract.period.RB <- function(RB,years=NULL,months=NULL) {
  #RB: raster::brick object
  #month: reference month for lag

  #get years
  date.seq = as.Date(substring(text = names(RB), first = 2, last = 11), "%Y.%m.%d")
  if(is.null(years)) years = unique(as.numeric(format(date.seq, "%Y")))
  if(is.null(month)) months = unique(as.numeric(format(date.seq, "%m")))
  
  time.idx = which(as.numeric(format(date.seq, "%Y")) %in% years & as.numeric(format(date.seq, "%m")) %in% months)
  
  new.RB=subset(RB, time.idx)
  
  return(new.RB)
}

#extract period from RB
extract.period.ts <- function(ts,years=NULL,months=NULL) {
  #ts:time series
  #years: years to extract
  #month: reference month for lag
  
  if(!is.matrix(ts)&!is.vector(ts)) ts=values(ts)
  
  #get years
  date.seq = as.Date(substring(text = names(ts), first = 2, last = 11), "%Y.%m.%d")
  if(is.null(years)) years = unique(as.numeric(format(date.seq, "%Y")))
  if(is.null(month)) months = unique(as.numeric(format(date.seq, "%m")))
  
  time.idx = which(as.numeric(format(date.seq, "%Y")) %in% years & as.numeric(format(date.seq, "%m")) %in% months)

  new.ts=ts[time.idx]
  names(new.ts)=date.seq[time.idx]
  
  return(new.ts)
}
  
# Function to calculate area-weighted mean on Rasterbrick:
# --------------------------------------------------------
fldmean.RB <- function(RB, w = "area", mask = NULL, maskvalue = NA, ret.xts = T) {
  
  # check for potential mask:
  if (!is.null(mask)) RB = mask(RB, mask = mask, maskvalue = maskvalue)
  
  RB.values = values(RB)
  area.vec = values(raster::area(RB))
  
  if (w == "area") {
    RB.mean.ts = apply(X = RB.values, MARGIN = 2, FUN=function(x) {
      weighted.mean(x, w = area.vec, na.rm = T)
    })  
  } else if (w == "none") {
    RB.mean.ts = cellStats(RB, stat="mean", na.rm = T)
  }
  
  if (ret.xts == T) {
    ret.ts = xts(x = RB.mean.ts, order.by = as.Date(substring(names(RB), first = 2, last = 11), "%Y.%m.%d"))
  } else if (ret.xts == F) {
    ret.ts = RB.mean.ts
  }
  
  return(ret.ts)
}


# Get period mean from Rasterbrick:
# --------------------------------------------------------
# RB = EOBS_TG_anom
# years = 1960:2017
# months = c(12, 1, 2)
# RB.values = values(RB)
period.mean.RB <- function(RB, years, months = "annual") {
  
  date.seq = as.Date(substring(text = names(RB), first = 2, last = 11), "%Y.%m.%d")
  years.seq = as.numeric(format(date.seq, "%Y"))
  months.seq = as.numeric(format(date.seq, "%m"))
  years.idx = which(unique(years.seq) %in% years) #only the selected years are returned
  
  if (months == "annual") {
    ret.RB = brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-07-01", sep=""))
    
  } else if (months == "djf"){
    months=c(12,1,2)
    ret.RB= brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year & months.seq %in% c(1,2) | years.seq == cur.year-1 & months.seq ==12)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-01-15", sep=""))
    
  } else if (months == "mam") {
    months=c(3,4,5)
    ret.RB= brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year & months.seq %in% months)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-04-15", sep=""))
    
  } else if (months == "jja") {
    months=c(6,7,8)
    ret.RB= brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year & months.seq %in% months)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-07-15", sep=""))
    
  } else if (months == "son") {
    months=c(9,10,11)
    ret.RB= brick(sapply(X = unique(years.seq), FUN=function(cur.year) mean(subset(RB, which(years.seq == cur.year & months.seq %in% months)))))
    names(ret.RB) = as.Date(paste(unique(years.seq), "-10-15", sep=""))
    
  } 
  
  return(subset(ret.RB,years.idx))
  
  
  # extract.ts.from.RB(TGhat_LDT_pred20CR.ERAI_ridge_pres.sfc)
}

#period (e.g. seasonal mean)
period.mean.ts <- function(ts, years, months = "annual") {
  
  date.seq = as.Date(substring(text = names(ts), first = 2, last = 11), "%Y.%m.%d")
  years.seq = as.numeric(format(date.seq, "%Y"))
  months.seq = as.numeric(format(date.seq, "%m"))
  years.idx = which(unique(years.seq) %in% years) #only the selected years are returned
  
  if (months == "annual") {
    ret.ts = sapply(X = unique(years.seq), FUN=function(cur.year) mean(ts[which(years.seq == cur.year)]))
    names(ret.ts) = as.Date(paste(unique(years.seq), "-07-01", sep=""))
    
  } else if (months == "djf"){
    months=c(12,1,2)
    ret.ts= sapply(X = unique(years.seq), FUN=function(cur.year) mean(ts[which(years.seq == cur.year & months.seq %in% c(1,2) | years.seq == cur.year-1 & months.seq ==12)]))
    names(ret.ts) = as.Date(paste(unique(years.seq), "-01-15", sep=""))
    
  } else if (months == "mam") {
    months=c(3,4,5)
    ret.ts= sapply(X = unique(years.seq), FUN=function(cur.year) mean(ts[which(years.seq == cur.year & months.seq %in% months)]))
    names(ret.ts) = as.Date(paste(unique(years.seq), "-04-15", sep=""))
    
  } else if (months == "jja") {
    months=c(6,7,8)
    ret.ts= sapply(X = unique(years.seq), FUN=function(cur.year) mean(ts[which(years.seq == cur.year & months.seq %in% months)]))
    names(ret.ts) = as.Date(paste(unique(years.seq), "-07-15", sep=""))
    
  } else if (months == "son") {
    months=c(9,10,11)
    ret.ts= sapply(X = unique(years.seq), FUN=function(cur.year) mean(ts[which(years.seq == cur.year & months.seq %in% months)]))
    names(ret.ts) = as.Date(paste(unique(years.seq), "-10-15", sep=""))
    
  } 
  
  return(ret.ts[years.idx])
}

#area weighted RB
area.weight.RB <- function(RB,svd=FALSE){
  RB.cos=cos(coordinates(RB)[,2]*2*pi/360)
  
  #where data arrays are always NA, areaweighting vector also has to be NA
  na.idx=which(is.na(values(RB)[,1]))
  for (i in na.idx) if(all(is.na(values(RB)[i,]))) RB.cos[i]=NA
  
  aw.RB=RB
  
  #area weighting
  if(svd==F){
  values(aw.RB)=values(RB) * RB.cos
  }
  
  #area weighting and centering if area weighting for eofs
  if(svd==T){
  RB.mean = rowMeans(values(RB)) 
  values(aw.RB)=(values(RB) - RB.mean) * sqrt(RB.cos) #area weighting and centering for EOFs
  }
  return(aw.RB)
}

#area weighted mean of RB, returns ts
mean.area.weighted.RB <- function(RB){
  RB.cos=cos(coordinates(RB)[,2]*2*pi/360)
  
  #where data arrays are always NA, areaweighting vector also has to be NA
  na.idx=which(is.na(values(RB)[,1]))
  for (i in na.idx) if(all(is.na(values(RB)[i,]))) RB.cos[i]=NA
  
  aw=values(RB) * RB.cos
  mean.aw=colSums(aw,na.rm=TRUE)/sum(RB.cos,na.rm=TRUE)
  return(mean.aw)
}

#detrend a single time series
detrend.ts <- function(ts,span=0.75,degree=2,return.trend=F){
  if (all(is.na(ts))) return(ts)
  ts.loess=loess(ts ~ c(1:length(ts)), span = span, degree = degree,na.action=na.exclude)
  ts.loess.predict=predict(ts.loess)
  detrended.ts=ts-ts.loess.predict
  
  if (return.trend==F) return(detrended.ts)
  if (return.trend==T) return(list(detrended.ts,ts.loess.predict))
}

#detrend each gridcell in a RB and return RB
detrend.RB <- function(RB,span=0.75,degree=2){
  RB.values=values(RB)
  matrix.detrended=apply(X=RB.values,MARGIN=1,FUN=function(x) detrend.ts(x,span=span,degree=degree,return.trend=F))
  values(RB)=t(matrix.detrended)
  return(RB)
}
