#----------------------------------------------------------------------------------------
#function for calculating correlation between mean SST and PSL or SSTs at each gridpoint
#with possibility for lagged data
#-----------------------------------------------------------------------------------------
# Rafael Bonafini
# 17.8.2020
#-----------------------------------------------------------------------------------------


ccf.monthly <- function(Y.ts,X.ts, mon.ix=1, nr.lags = 48) {
  #calculates monthly correlations between two univariate time series
  start.ix = mon.ix + nr.lags
  Y.mon.ts = Y.ts[seq(start.ix, length(Y.ts), 12)]
  lag.cor = sapply(X = start.ix:mon.ix, FUN=function(ix) cor(Y.mon.ts, X.ts[seq(ix, length(X.ts), 12)][1:length(Y.mon.ts)]))
  # plot(0:48, lag.cor)
  return(lag.cor)
}

ccf.monthly.gridcell <- function(Y.RB,X.RB,month=1,nr.lags=48) {
  # (1) Subset relevant years for Y in training sample and KICK out NA grid cells:
  # -----------------------------------------------
  Y.RB.values = values(Y.RB)               # Get Y values (only!) for training time 
  NONA.idx = which(apply(X = Y.RB.values, MARGIN=1, FUN=function(x) !all(is.na(x))))             # Get NA values in Y raster
  Y.RB.values.nona = Y.RB.values[NONA.idx,]                                                     # New Y dataset with NA's removed
  Y.RB.coord.nona = coordinates(Y.RB)[NONA.idx,]                                                # Coordinates of Y dataset.
  
  if(isS4(X.RB)){
    X.RB.values = values(X.RB)
    X.RB.values.nona = X.RB.values[NONA.idx,]
    X.RB.coord.nona = coordinates(X.RB)[NONA.idx,]                                                # Coordinates of Y dataset.
  
    rm(Y.RB.values)
    rm(X.RB.values)
    
    # 3.2 Define small raster object to use in parallel loops:
    X.coords = coordinates(X.RB)
    X.template = raster(subset(X.RB, 1)) #raster with same extent as X.RB to store coordinates
    values(X.template) <- 1:(dim(X.coords)[1])
    
    grid.correlation=sapply(X=1:(dim(Y.RB.values.nona)[1]),FUN=function(i){ 
      print(i)
      # (b) Subset/Define new raster based on present data point:
  
      cur.correlation= ccf.monthly(Y.RB.values.nona[i,],X.RB.values.nona[i,],mon.ix=month,nr.lags=nr.lags)
      return(cur.correlation)
      
    })
  }
  
  #if correlation is with a vector, not a rasterbrick (e.g. with 1. EOF of SLP field)
  if(is.vector(X.RB)){
    grid.correlation=sapply(X=1:(dim(Y.RB.values.nona)[1]),FUN=function(i){ 
      print(i)
      cur.correlation= ccf.monthly(Y.RB.values.nona[i,],X.RB,mon.ix=month,nr.lags=nr.lags)
      return(cur.correlation)
      
    })
  }
  
  ## (3) Derive NEW Yhat raster and return:
  ## ---------------------------------------------------------
  grid.correlation.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, dim(grid.correlation)[1]))
  values(grid.correlation.RB)[NONA.idx,] = t(grid.correlation)    # Put predictions into Yhat raster raster (!!)
  names(grid.correlation.RB) = paste0("lag",0:nr.lags)
  return(grid.correlation.RB)
}
  


#calculate correlation between Y.mean and each gridcell of X.RB
#returns raster
cor.mean.gridcells<- function(Y.ts,X.ts){
  #Y.ts: mean time series of Y, 1 dimensional
  # X.ts: timeseries with N gridpoints and K observations
  ncell=length(X.ts[,1])
  cor.gridcells=sapply(X=1:ncell, FUN=function(gridcell) cor(Y.ts,X.ts[gridcell,]))
  return(cor.gridcells)
}

#calculate monthly correlation between Y.mean and each gridcell of X.RB
cor.lag.mean.gridcells.monthly <- function(Y.ts,X.ts, mon.ix=1, nr.lags = 12){
  ncell=length(X.ts[,1])
  start.ix = mon.ix + nr.lags
  Y.mon.ts = Y.ts[seq(start.ix, length(Y.ts), 12)]
  cor.gridcells = sapply(X=start.ix:mon.ix,FUN=function(ix) sapply(X=1:ncell,FUN=function(gridcell) cor(Y.mon.ts, X.ts[gridcell,seq(ix, length(X.ts[1,]), 12)][1:length(Y.mon.ts)])))

  return(t(cor.gridcells))
}

#calculate monthly correlation function
monthly.cor.fun <- function(Y.ts,X.ts, mon.ix=1, nr.lags = 12) {
  start.ix = mon.ix + nr.lags
  Y.mon.ts = Y.ts[seq(start.ix, length(Y.ts), 12)]
  lag.cor = sapply(X = start.ix:mon.ix, FUN=function(ix) cor(Y.mon.ts, X.ts[seq(ix, length(X.ts), 12)][1:length(Y.mon.ts)]))
  # plot(0:48, lag.cor)
  return(lag.cor)
}

#calculate signal to noise ratio of a time series
signal.to.noise <- function(ts,years,months){
  date = as.Date(substring(text = names(ts), first = 2, last = 11), "%Y.%m.%d")
  if(is.null(years)) years = as.numeric(format(date, "%Y"))
  if(is.null(months)) months = as.numeric(format(date, "%m"))
  
  idx=which(as.numeric(format(date, "%Y")) %in% years  & as.numeric(format(date, "%m")) %in% months)
  stn.detrend=detrend.ts(ts[idx],return.trend=T)
  stn=diff(range(stn.detrend[[2]]))/sqrt(var(stn.detrend[[1]]))
  return(stn)
}

# ------------------------------------------------------------------------------------------
# SVD regression: calculation of SVDs of a 2d-matrix (datapoints as rows, timesteps as cols)
# as predictor and a mean time series as predictant
# uses SVD function
#-------------------------------------------------------------------------------------------

svd.regression <- function(Y.ts,X.svd.v,predictors=1,train.years,pred.years){
  #predictors: numbers of PCs used for the linear regression
  #train.years: timesteps used to train the linear model
  #pred.years: timesteps where Y.hat is predicted
  
  #train.years and pred.years have to be passed
  if (missing(train.years) | missing(pred.years)){
    stop("Arguments for train.years or pred.years are missing.")
  }
  
 
  # create predictor array
  predictor=X.svd.v[train.years,1]
  if (predictors>1){
    for (i in seq(2,predictors)) {
      predictor=cbind(predictor,X.svd.v[train.years,i])
    }
  }
  
  #regress temp on SVD PCs
  model.eof=lm(Y.ts[train.years]~predictor)
  
  #test predicted data on subset of data
  Y.hat= X.svd.v[pred.years,1:predictors] %*% model.eof$coefficients[2:(predictors+1)]
  Y.hat= model.eof$coefficients[1]+Y.hat
  return(Y.hat)
}

# Lagged SVD regression: Pass X SVD and Y mean time series and use contemporary and lagged PCs
# for a prediction of Y
#predictors: numbers of PCs used for the linear regression
#train.years: timesteps used to train the linear model
#pred.years: timesteps where Y.hat is predicted
#month.idx: index of predicted month in Y
#nr.lags: nr of lags, has to be multiple of 12

svd.regression.lag.cumulative <- function(Y.ts,X.svd,eofs=1:5,train.years=1:1000,pred.years=1001:2000,month.idx=1,lags=0:12){
  #months: array of lags or months to include in the predictor
  #type: 'lag' or 'months', indicates whether the predictors consist of lags or fixed months
  
  
  #train.years and pred.years have to be passed
  if (missing(train.years) | missing(pred.years)){
    stop("Arguments for train.years or pred.years are missing.")
  }

  # create predictor array
  start.idx=month.idx+max(lags)
  Y.idx=seq(start.idx, length(Y.ts), 12)
  Y.mon=Y.ts[Y.idx]
  predictor=c()
  correlation.cumulative=array(NA,length(lags))
  for (i in 1:length(lags)) {
    lag=lags[i]
    X.idx=Y.idx-lag
    X.predictor=X.svd[X.idx,eofs]
    predictor=cbind(predictor,X.predictor)
    model.eof=lm(Y.mon[train.years]~predictor[train.years,])
    Y.hat= predictor %*% model.eof$coefficients[-1]
    Y.hat= model.eof$coefficients[1]+Y.hat
    
    correlation.cumulative[i]=cor(Y.hat[pred.years],Y.mon[pred.years])
  }
  

  
  return(correlation.cumulative)
}

