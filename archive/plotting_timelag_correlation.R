#plotting of timelag correlations
region="Atlantic"
Y="SST"

type="seasonal" #seasonal or monthly
outpath="/net/h2o/climphys1/rdaenzer/output/"
filename=paste("pc_regr_lagcor",Y,region,type,"anom.RDS",sep="_")
file=file.path(outpath,filename)

figpath="/net/h2o/climphys1/rdaenzer/figures/"
setwd(figpath)


timelag.correlation=readRDS(file)

if (type=="monthly") period.names=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
if (type=="seasonal") period.names=c("djf","mam","jja","son")
lag=1:length(period.names)


  
{figname=paste0("02_pc_regr_lag_corr_",paste(Y,region,type,period.names[i],sep="_"),".pdf")
pdf(figname, width = 6*4, height = 5.5*3)
par(mfrow=c(4, 3))
for (i in length(period.names)){
  plot(lag,timelag.correlation[,i],main=paste(type, "lag correlations of pred. vs. obs. mean SST",period.names[i], "anomalies",sep=" "), 
       xlab="time lag", ylab="linear correlation",type="l")}
dev.off()
}

