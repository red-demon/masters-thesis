#plot frame with lower and upper left and right corners preset
lx=150; ux=250; ly=20; uy=60
lines(c(lx,ux),c(ly,ly),col="red") + lines(c(ux,ux),c(ly,uy),col="red") + lines(c(ux,lx),c(uy,uy),col="red")+lines(c(lx,lx),c(uy,ly),col="red")
coastsCoarse

#convert RB to dataframe
df=data.frame(coordinates(X.region.RB),X.svd.u[,1])
colnames(df)=c("lon","lat","svd")

#ggplot2 plotting
#plot dataframe by using "world2" data
world_map=map_data("world2",xlim = c(120, 260),ylim=c(0,70),interior=FALSE)
ggplot() +
  geom_polygon(data=world_map,aes(x=long,y=lat,group=group)) +
  geom_raster(data=df, aes(x=lon, y=lat, fill=svd), alpha=0.8)+
  scale_fill_distiller(palette = "RdYlBu")

#plot dataframe by using coord_map
xlim=X.extent[1:2]
ylim=X.extent[3:4]
ggplot(data=df, aes(x=lon, y=lat, fill=svd)) +
  geom_tile()+ scale_fill_distiller(palette = "RdYlBu") + 
  coord_map(projection="mercator",xlim = xlim,ylim=ylim) + borders("world2", xlim = xlim, ylim =ylim)


ggplot(data=df,aes(x=lon,y=lat,fill=svd))+geom_tile


#raster plots
plot(X.region.RB,1,useRaster=T,cl=color.palette(100))
maps::map("world2",add=T)
