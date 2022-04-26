setwd(paste0("/net/h2o/climphys1/rdaenzer/figures/basic_analysis"))

#for maps
{png("test.png",height=2*300,width=2*350)
par(mfrow=c(2,2),oma=c(4,3,1,5),mar=c(0,0,2,1),cex=1.1)
for (i in 1:4){
  raster.plot=subset(RB,i)
  image(raster.plot,xlab="",ylab="",axes=F,col=color.high.low)
  title(main=paste("plot",i),line=0.2)
  grid(col="black")
  axis(1,labels=(if(i>2) labels=TRUE else labels=FALSE))
  axis(2,labels=(if(i%%2==1) labels=TRUE else labels=FALSE))
  lines(coastsCoarse)
}
title(main="test plot",outer=T,line=0)
title(xlab="lat",ylab="lon",outer=T,line=2)

par(mfrow=c(1,1),oma=c(4,0,2,0),mar=c(0,0,0,0))
image.plot(z=pretty(c(-600,600),20),legend.only=T,add=T,col=color.high.low,horizontal=F,legend.shrink=0.7)
dev.off()}

#for maps
{png("test2.png",height=300,width=2*300)
  par(mfrow=c(1,2),oma=c(4,3,1,5),mar=c(0,0,2,1))
  for (i in 1:2){
    raster.plot=subset(RB,i)
    image(raster.plot,xlab="",ylab="",axes=F,col=color.high.low)
    title(main=paste("plot",i),line=0.2)
    grid(col="black")
    axis(1,labels=T)
    axis(2,labels=(if(i%%2==1) labels=TRUE else labels=FALSE))
    lines(coastsCoarse)
  }
  title(main="test plot",outer=T,line=0)
  title(xlab="lat",ylab="lon",outer=T,line=2)
  
  par(mfrow=c(1,1),oma=c(2,0,2,0),mar=c(0,0,0,0))
  image.plot(z=pretty(c(-600,600),20),legend.only=T,add=T,col=color.high.low,horizontal=F,legend.shrink=0.7)
  dev.off()}

#for line plots
x=1:100
y=2*rnorm(100,0,0.01)*x + rnorm(100,0,6)
y2=3*rnorm(100,0,0.01)*x + rnorm(100,0,6)
y3=4*rnorm(100,0,0.01)*x + rnorm(100,0,6)
y4=5*rnorm(100,0,0.01)*x + rnorm(100,0,6)

{png("test_line.png")
  par(mfrow= c(2, 2),oma=c(5,4,0,0),mar=c(0,0,0.5,1),bty="n") #bty="n" to remove frame and have nicer axes  
  plot(x,y,type="l",axes=F, panel.first=grid());text(10,10,label="panel1")
  plot(x,y2,type="l",xaxt="n",mar=c(0,2,0,0), panel.first=grid());title(sub="panel1")
  plot(x,y3,type="l", panel.first=grid());legend("topleft",legend="panel1")
  plot(x,y4,type="l", panel.first=grid())
  dev.off()}

df=data.frame(x,y,y2,y3,y4)
list=list(y,y2,y3,y4)
#create dataframe from list
df2=data.frame(sapply(list,FUN=function(x) (unlist(x))))
