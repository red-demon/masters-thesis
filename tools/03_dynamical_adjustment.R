# ---------------------------------------------------------
# Function Repository for Dynamical Adjustment 
# Rafael Bonafini 
# 2.10.2020
#-------------------------------------------------------------



#elasticnet model for monthly data on 1 Gridcell or region mean, using lagged data as predictors
dynamical.adjustment.gridcell.monthly <- function(Y.train=NULL,X.pred=NULL,X.train=NULL,model.cv.glmnet=NULL,detrend.X=F,detrend.Y=F,train.years=NULL,pred.years=NULL,lags=NULL,w.lags=NULL,s.lags=NULL,months=1,alpha=0.1,s="lambda.min",ret.cv.model=F){
  #Y.ts=univariate timeseries of a predictant variable
  #X.pred=multivariate timeseries of a predictor variable (format: PxN, so rows=observations, cols=variables)
  
  # a) no model is provided, return Y.hat
  #-----------------------------------------------
  if(is.null(model.cv.glmnet) & ret.cv.model==F){
    # CHECK IF Separate dataset of TRAIN Years are given (can EITHER GIVE train.years OR train on entire Y time):
    if (is.null(train.years)) train.years = unique(as.numeric(format(as.Date(substring(text = names(X.train), first = 2, last = 11), "%Y"))))
    if (is.null(pred.years)) pred.years = unique(as.numeric(format(as.Date(substring(text = names(X.pred), first = 2, last = 11), "%Y"))))
    
    # CHECK IF X.train and X.pred are matrices or RB
    if (!is.matrix(X.pred)) X.pred=t(values(X.pred))
  
    #use X.pred as training dataset if X.train is NULL
    if (is.null(X.train)) X.train=X.pred  
    
    # Check if X.train is Matrix and convert to matrix if not
    if (!is.null(X.train) & !is.matrix(X.train)) X.train=t(values(X.train))
  
    #create list to store the cv.glmnet models for each month
    model.cv.glmnet=list()
    
    #array to store the monthly predictions Y.hat, should have the length of X.train truncatet to pred.years
    date.seq=as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d")
    Y.date.seq=date.seq[which(as.numeric(format(date.seq, "%Y")) %in% pred.years & as.numeric(format(date.seq, "%m")) %in% months)]
    Y.date.names=names(X.pred[,1][which(as.numeric(format(date.seq, "%Y")) %in% pred.years & as.numeric(format(date.seq, "%m")) %in% months)])
    Y.hat=rep(NA,length(Y.date.names))
    
    for (cur.month in months){  
      print(cur.month)
      
      #different lags for summer and winter 
      if(!is.null(w.lags) & cur.month %in% c(12,1,2,3,4,5)) lags=w.lags
      if(!is.null(s.lags) & cur.month %in% c(6,7,8,9,10,11)) lags=s.lags
      
      # define Y indices for training
      Y.train.idx = which(as.numeric(format(as.Date(substring(text = names(Y.train), first = 2, last = 11),"%Y.%m.%d"), "%m")) %in% cur.month & 
                            as.numeric(format(as.Date(substring(text = names(Y.train), first = 2, last = 11),"%Y.%m.%d"), "%Y")) %in% train.years)    # WHICH INDICES ARE IN PRESENT MONTHS AND IN train.years?
      
      # add (lag) predictors and define X and Y for the model
      #detrend X?
      if(detrend.X==T){ #detrend 
      X=apply(X.train[Y.train.idx,],MARGIN=2,FUN=detrend.ts) #add concurrent field as predictor
        if(!is.null(lags)){
          for (lag in lags) {
            X=cbind(X, apply(X.train[Y.train.idx-lag,],MARGIN=2,FUN=detrend.ts))
          }
        } 
      } else { #no detrend
            X=X.train[Y.train.idx,] #add concurrent field as predictor
            if(!is.null(lags)){
              for (lag in lags) {
                X=cbind(X, X.train[Y.train.idx-lag,])
              }
            }
      }
       #detrend Y?
      if(detrend.Y==F) Y=Y.train[Y.train.idx]
      if(detrend.Y==T) Y=detrend.ts(Y.train[Y.train.idx])
      
      
      #train model with predictor
      model.cv.glmnet[[cur.month]]=cv.glmnet(X,Y, alpha=alpha)
      
      #predict Yhat using X.pred only if not the model is returned
      X.pred.idx = which(as.numeric(format(as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d"),"%m")) %in% cur.month & 
                           as.numeric(format(as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d"), "%Y")) %in% pred.years)    # WHICH INDICES ARE IN PRESENT MONTHS AND IN train.years?
      
      if(detrend.X==T){
        predictor=apply(X.pred[X.pred.idx,],MARGIN=2,FUN=detrend.ts) #add concurrent field as predictor
        if(!is.null(lags)){
          for (lag in lags) {
            predictor=cbind(predictor, apply(X.pred[X.pred.idx-lag,],MARGIN=2,FUN=detrend.ts))
          }
        } 
      } else { 
        predictor=X.pred[X.pred.idx,]
          if(!is.null(lags)){
            for (lag in lags) {
              predictor=cbind(predictor, X.pred[X.pred.idx-lag,])
            }
          }
        }
      Y.hat[seq(which(months %in% cur.month),length(Y.hat),length(months))]=predict(model.cv.glmnet[[cur.month]],predictor,s=s)
    
    }
    
    #create names for Y.hat
    names(Y.hat)=Y.date.names
    
    return(Y.hat)
  }
  
  # b) no model provided, return cv.model but no Y.hat
  #-----------------------------------------------
  if(is.null(model.cv.glmnet) & ret.cv.model==T){
    # CHECK IF Separate dataset of TRAIN Years are given (can EITHER GIVE train.years OR train on entire Y time):
    if (is.null(train.years)) train.years = unique(as.numeric(format(as.Date(substring(text = names(X.train), first = 2, last = 11), "%Y"))))

    # Check if X.train is Matrix and convert to matrix if not
    if (!is.null(X.train) & !is.matrix(X.train)) X.train=t(values(X.train))
    
    #create list to store the cv.glmnet models for each month
    model.cv.glmnet=list()
    
    for (cur.month in months){  
      print(cur.month)
      
      #different lags for summer and winter 
      if(!is.null(w.lags) & cur.month %in% c(12,1,2,3,4,5)) lags=w.lags
      if(!is.null(s.lags) & cur.month %in% c(6,7,8,9,10,11)) lags=s.lags
      
      # define Y indices for training
      Y.train.idx = which(as.numeric(format(as.Date(substring(text = names(Y.train), first = 2, last = 11),"%Y.%m.%d"), "%m")) %in% cur.month & 
                            as.numeric(format(as.Date(substring(text = names(Y.train), first = 2, last = 11),"%Y.%m.%d"), "%Y")) %in% train.years)    # WHICH INDICES ARE IN PRESENT MONTHS AND IN train.years?
      
      # add (lag) predictors and define X and Y for the model
      #detrend?
      if(detrend.X==T){ #detrend 
        X=apply(X.train[Y.train.idx,],MARGIN=2,FUN=detrend.ts) #add concurrent field as predictor
        if(!is.null(lags)){
          for (lag in lags) {
            X=cbind(X, apply(X.train[Y.train.idx-lag,],MARGIN=2,FUN=detrend.ts))
          }
        } 
        } else { #no detrend
          X=X.train[Y.train.idx,] #add concurrent field as predictor
          if(!is.null(lags)){
            for (lag in lags) {
              X=cbind(X, X.train[Y.train.idx-lag,])
            }
          }
        }
      #detrend Y?
      if(detrend.Y==F) Y=Y.train[Y.train.idx]
      if(detrend.Y==T) Y=detrend.ts(Y.train[Y.train.idx])
      
        #train model with predictor
        model.cv.glmnet[[cur.month]]=cv.glmnet(X,Y, alpha=alpha)
      
    }
    
    return(model.cv.glmnet)
  }
  
  # c) if cv model is given, it can be used to make the prediction
  #-----------------------------------------------
  if(!is.null(model.cv.glmnet)){
    
    # check if pred.years are given and create otherwise
    if (is.null(pred.years)) pred.years = unique(as.numeric(format(as.Date(substring(text = names(X.pred), first = 2, last = 11), "%Y"))))
    
    # CHECK IF X.train and X.pred are matrices or RB
    if (!is.matrix(X.pred)) X.pred=t(values(X.pred))
    
    #array to store the monthly predictions Y.hat, should have the length of X.train truncatet to pred.years
    date.seq=as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d")
    Y.date.seq=date.seq[which(as.numeric(format(date.seq, "%Y")) %in% pred.years & as.numeric(format(date.seq, "%m")) %in% months)]
    Y.date.names=names(X.pred[,1][which(as.numeric(format(date.seq, "%Y")) %in% pred.years & as.numeric(format(date.seq, "%m")) %in% months)])
    Y.hat=rep(NA,length(Y.date.names))
    
    for (cur.month in months){  
      print(cur.month)
      
      #different lags for summer and winter 
      if(!is.null(w.lags) & cur.month %in% c(12,1,2,3,4,5)) lags=w.lags
      if(!is.null(s.lags) & cur.month %in% c(6,7,8,9,10,11)) lags=s.lags
      
      #predict Yhat using X.pred
      X.pred.idx = which(as.numeric(format(as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d"),"%m")) %in% cur.month & 
                           as.numeric(format(as.Date(substring(text = names(X.pred[,1]), first = 2, last = 11), "%Y.%m.%d"), "%Y")) %in% pred.years)    # WHICH INDICES ARE IN PRESENT MONTHS AND IN train.years?
      
      #detrend?
      if(detrend.X==T){ #detrend
        predictor=apply(X.pred[X.pred.idx,],MARGIN=2,FUN=detrend.ts) #add concurrent field as predictor
        if(!is.null(lags)){
          for (lag in lags) {
            predictor=cbind(predictor, apply(X.pred[X.pred.idx-lag,],MARGIN=2,FUN=detrend.ts))
          }
        } 
      } else { #no detrend
          predictor=X.pred[X.pred.idx,]
          if(!is.null(lags)){
            for (lag in lags) {
              predictor=cbind(predictor, X.pred[X.pred.idx-lag,])
          }
        }
      }
          
      Y.hat[seq(which(months %in% cur.month),length(Y.hat),length(months))]=predict(model.cv.glmnet[[cur.month]],predictor,s=s)
    }
    
    #create ts from Y.hat or just assign date names to Y.hat vector
    #Y.hat=xts(x=Y.hat,order.by=Y.date.seq)
    names(Y.hat)=Y.date.names
    return(Y.hat)
  }
  
} 


#dynamical adjustment of whole field

dynamical.adjustment.field.RB <- function(Y.RB,X.RB,X.train.RB=NULL,model.cv.glmnet=NULL,detrend.X=F,detrend.Y=F,train.years=NULL,pred.years=NULL,lags=NULL,s.lags=NULL,w.lags=NULL,months=1,alpha=0.1,s="lambda.min",x.domain=20,y.domain=20,ret.cv.model=F) {
  
  # (1) Subset relevant years for Y in training sample and KICK out NA grid cells:
  # -----------------------------------------------
  if(is.null(train.years)) train.years = unique(as.numeric(format(as.Date(substring(names(Y.RB), first = 2, last = 11), format = "%Y.%m.%d"), "%Y")))
  Y.RB.values = t(values(Y.RB))                 # Get Y values (only!) for training time 
  Y.RB.date = as.Date(substring(names(Y.RB), first = 2, last = 11), format = "%Y.%m.%d")        # Get dates for training entie y time series
  Y.RB.date = Y.RB.date[which(as.numeric(format(Y.RB.date, "%Y")) %in% train.years)]            # Get only dates fpr training time series
  
  NONA.idx = which(apply(X = Y.RB.values, MARGIN=2, FUN=function(x) !all(is.na(x))))             # Get NA values in Y raster
  Y.RB.values.nona = Y.RB.values[,NONA.idx]                                                     # New Y dataset with NA's removed
  Y.RB.coord.nona = coordinates(Y.RB)[NONA.idx,]                                                # Coordinates of Y dataset.
  
  # (2) Sanity CHECKS and (possible) early return: 
  # ---------------------------------------------------------------------
  if (dim(Y.RB.values.nona)[2] == 0) {    # are all values NA in Y.RB.values.nona?
    Yhat.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, nlayers(X.RB)))
    values(Yhat.RB) = NA; names(Yhat.RB) = names(X.RB);  # ALL VALUES SHOULD BE NA !
    return(Yhat.RB)
  }
  # remove unnecessary variables to avoid exceeding memory:
  rm(Y.RB.values)
  

  
  # (3) RUN DYNAMICAL ADJUSTMENT:
  ## ---------------------------------------------------------------------------------------

    # Define small raster object to use in parallel loops:
  X.coords = coordinates(X.RB)
  X.template = raster(subset(X.RB, 1)) #raster with same extent as X.RB to store coordinates
  values(X.template) <- 1:(dim(X.coords)[1])
  
  # Start the clock!
  # ptm <- proc.time()
  Xpred=t(values(X.RB))
  Xtrain=Xpred
  Ytrain=Y.RB.values.nona
  if(!is.null(X.train.RB)) Xtrain=t(values(X.train.RB))
  
  #3a) return normal Yhat.Rb
  #----------------------------------------------
  if(ret.cv.model==F & is.null(model.cv.glmnet)){
  Yhat=sapply(X=1:(dim(Y.RB.values.nona)[2]),FUN=function(i){ 
    print(i)
    # (b) Subset/Define new raster based on present data point:
    cur.extent = extent(c(Y.RB.coord.nona[i,1] - x.domain, Y.RB.coord.nona[i,1] + x.domain, Y.RB.coord.nona[i,2] - y.domain, Y.RB.coord.nona[i,2] + y.domain))
    cur.X.idx = values(crop(X.template, cur.extent))
    cur.X.RB = Xpred[,cur.X.idx]
    cur.Xtrain = NULL
    if (!is.null(X.train.RB)) cur.Xtrain = Xtrain[,cur.X.idx]
    
    cur.Yhat= dynamical.adjustment.gridcell.monthly(Y.train=Ytrain[,i],X.pred=cur.X.RB,X.train=cur.Xtrain,train.years=train.years,pred.years=pred.years,
                                                 lags=lags,s.lags=s.lags,w.lags=w.lags,months=months,alpha=alpha,s=s,detrend.X=detrend.X,detrend.Y=detrend.Y)
    return(cur.Yhat)

  })
  
  ## (4) Derive NEW Yhat raster and return:
  ## ---------------------------------------------------------
  Yhat.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, dim(Yhat)[1]))
  values(Yhat.RB)[NONA.idx,] = t(Yhat)    # Put predictions into Yhat raster raster (!!)
  names(Yhat.RB) = names(Yhat[,1])
  
  return(Yhat.RB)
  }
  
  #3b)return cv.model list for each gridcell
  #----------------------------------------
  if(ret.cv.model==T){
    model.cv.glmnet=lapply(X=1:(dim(Y.RB.values.nona)[2]),FUN=function(i){ 
      print(i)
      # (b) Subset/Define new raster based on present data point:
      cur.extent = extent(c(Y.RB.coord.nona[i,1] - x.domain, Y.RB.coord.nona[i,1] + x.domain, Y.RB.coord.nona[i,2] - y.domain, Y.RB.coord.nona[i,2] + y.domain))
      cur.X.idx = values(crop(X.template, cur.extent))
      cur.X.RB = Xpred[,cur.X.idx]
      cur.Xtrain = NULL
      if (!is.null(X.train.RB)) cur.Xtrain = Xtrain[,cur.X.idx]
      
      cur.cv.model= dynamical.adjustment.gridcell.monthly(Y.train=Ytrain[,i],X.pred=cur.X.RB,X.train=cur.Xtrain,train.years=train.years,pred.years=pred.years,
                                                      lags=lags,s.lags=s.lags,w.lags=w.lags,months=months,alpha=alpha,s=s,ret.cv.model=T,detrend.Y=detrend.Y,detrend.X=detrend.X)
      return(cur.cv.model)
      
    })
    
    
    return(model.cv.glmnet)
  }
  
  #3c) use cv.model list for each gridcell
  #--------------------------------------
  if(!is.null(model.cv.glmnet)){
    Yhat=sapply(X=1:(dim(Y.RB.values.nona)[2]),FUN=function(i){ 
      print(i)
      # (b) Subset/Define new raster based on present data point:
      cur.extent = extent(c(Y.RB.coord.nona[i,1] - x.domain, Y.RB.coord.nona[i,1] + x.domain, Y.RB.coord.nona[i,2] - y.domain, Y.RB.coord.nona[i,2] + y.domain))
      cur.X.idx = values(crop(X.template, cur.extent))
      cur.X.RB = Xpred[,cur.X.idx]
      cur.Xtrain = NULL
      if (!is.null(X.train.RB)) cur.Xtrain = Xtrain[,cur.X.idx]
      
      cur.Yhat= dynamical.adjustment.gridcell.monthly(X.pred=cur.X.RB,model.cv.glmnet=model.cv.glmnet[[i]],train.years=train.years,pred.years=pred.years,
                                                      lags=lags,s.lags=s.lags,w.lags=w.lags,months=months,alpha=alpha,s=s,ret.cv.model=F,detrend.Y=detrend.Y,detrend.X=detrend.X)
      return(cur.Yhat)
      
    })
    
    
    ## (4) Derive NEW Yhat raster and return:
    ## ---------------------------------------------------------
    Yhat.RB = subset(brick(raster(subset(Y.RB, 1))), rep(1, dim(Yhat)[1]))
    values(Yhat.RB)[NONA.idx,] = t(Yhat)    # Put predictions into Yhat raster raster (!!)
    names(Yhat.RB) = names(Yhat[,1])
    
    return(Yhat.RB)
  }
}



