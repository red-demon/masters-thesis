library(RColorBrewer)

#diverging
color.RdBu=rev(colorRampPalette(brewer.pal(n=9, name="RdBu"))(100))
color.RdYlBu=rev(colorRampPalette(brewer.pal(n=9, name="RdYlBu"))(100))
color.RdYlGn=rev(colorRampPalette(brewer.pal(n=9, name="RdYlGn"))(21))
color.spectral=rev(colorRampPalette(brewer.pal(n=9, name="Spectral"))(100))
color.spectral.white=rev(colorRampPalette(c( "#F46D43", "#FDAE61", "#FEE08B", "#FFFFFF", "#E6F598", "#ABDDA4","#66C2A5"))(100)) #Spectral, but with white in the middle


#highlow
color.YlBu=colorRampPalette(brewer.pal(n=5, name="YlGnBu"))(100)
color.BrBG=colorRampPalette(brewer.pal(n=5, name="BrBG"))(100)
color.Bu=colorRampPalette(brewer.pal(n=5, name="Blues"))(100)
color.Gn=colorRampPalette(brewer.pal(n=5, name="Greens"))(100)

color.YlGn=colorRampPalette(brewer.pal(n=9, name="YlGn"))(100)
color.YlGnBu=colorRampPalette(c("#E6F598", "#ABDDA4","#66C2A5", "#3288BD","#5E4FA2"))(100) #Spectral, but only half of it
color.YlOrBr=colorRampPalette(brewer.pal(n=7, name="YlOrBr"))(100)
color.viridis=viridis::viridis(100)

#main colors
color.high.low=color.viridis
color.high.low.negative=color.RdYlGn
color.warm.cold=color.RdYlBu

