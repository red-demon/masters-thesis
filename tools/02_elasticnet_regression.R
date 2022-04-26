# ---------------------------------------------------------
# Function Repository for Elasticnet regression 
# Rafael Bonafini 
# 2.10.2020
#-------------------------------------------------------------


#elasticnet model for monthly data on 1 Gridcell or region mean, using lagged data as predictors
elasticnet.gridcell.monthly <- function(Y.ts,X.ts,train.years=1:1000,pred.years=1001:1998,lags=0:12,month,alpha=0.1,s="lambda.1se",ret.cv.model=F){
#Y.ts=univariate timeseries of a predictant variable
#X.ts=multivariate timeseries of a predictor variable (format: PxN, so rows=observations, cols=variables)
  Y.mon=extract.lag.period.idx(Y.ts,month=month,lag=0)
  predictor=c()
  for (lag in lags) {
    predictor=rbind(predictor, extract.lag.period.idx(t(X.ts),month=month,lag=lag))
  }
  

  
  
  predictor=t(predictor)
  
  tic(paste(month,max(lags)))
  model.cv.glmnet=cv.glmnet(predictor[train.years,],Y.mon[train.years], alpha=alpha)
  Y.hat=predict(model.cv.glmnet,predictor[pred.years,],s=s)
  toc()
  
  if (ret.cv.model==F) return(Y.hat)
  if (ret.cv.model==T) return(model.cv.glmnet)
}  


