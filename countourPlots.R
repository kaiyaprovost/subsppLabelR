library(raster)
myras = "/Users/kprovost/Documents/Research/subsppLabelR/Phainopepla nitens/Phainopepla nitens nitens_raw_raster.tif"
ras = raster(myras)

quantiles = as.numeric(quantile(values(ras),seq(0.5,1,0.2)))
quanmap = ras
values(quanmap) = 0.3
values(quanmap)[values(ras)>quantiles[1]] = 0.50
values(quanmap)[values(ras)>quantiles[2]] = 0.70
values(quanmap)[values(ras)>quantiles[3]] = 0.90

plot(quanmap)

myras2 = "/Users/kprovost/Documents/Research/subsppLabelR/Phainopepla nitens/Phainopepla nitens lepida_raw_raster.tif"
ras2 = raster(myras2)

quantiles2 = as.numeric(quantile(values(ras2),seq(0.5,1,0.2)))
quanmap2 = ras
values(quanmap2) = 0.3
values(quanmap2)[values(ras2)>quantiles2[1]] = 0.50
values(quanmap2)[values(ras2)>quantiles2[2]] = 0.70
values(quanmap2)[values(ras2)>quantiles2[3]] = 0.90

pal1 <- colorRampPalette(c(rgb(0,1,1,1),rgb(0,0,1,1)),alpha=T)
pal2 <- colorRampPalette(c(rgb(1,1,0,1),rgb(1,0,0,1)),alpha=T)

contour1 = rasterToContour(quanmap)
contour2 = rasterToContour(quanmap2)
quanmap[values(quanmap==0.3)] = NA
quanmap2[values(quanmap2==0.3)] = NA

par(mfrow=c(1,2))
plot(quanmap,col=pal1(3))
#plot(contour1,add=T,lty=2,lwd=2)
plot(quanmap2,col=pal2(3),add=F)
#plot(contour2,add=T)

legend("bottomleft",
       legend=c(1,2,3,4,5,6),
       fill=c(pal1(3),pal2(3)),
       lty=c(2,2,2,1,1,1),
       lwd=c(2,2,2,1,1,1))
