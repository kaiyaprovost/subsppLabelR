library(subsppLabelR)
library(e1071)

## import data I already have

df2 <- read.delim("~/Work/GitHub/subsppLabelR/Phainopepla nitens/LOCALITIES/Phainopepla nitens_0.9_subspplabelR_RAW.txt")
df = df2[,c("longitude","latitude","subspecies")]
#df = unique(df) -- NO
df = df[df$subspecies!="unknown",]
df$subspecies = as.factor(df$subspecies)
## quick error removal
df = df[df$latitude<=45,]
df = df[df$longitude<=-90,]

library(raster)

#r=raster(nrows=(90*2),ncols=(180*2),xmn=-180,xmx=180,ymn=-90,ymx=90)
sampling_bias1 = raster("C:/Users/kaiya/Documents/Work/density_sampling_maps_land/land/0.50/Aves_number_occ_per_cell.tif")
sampling_bias2 = raster("C:/Users/kaiya/Documents/Work/density_sampling_maps_land/land/0.50/Aves_number_species_per_cell.tif")
sb1ps = sampling_bias1 / sampling_bias2 ## assumes all species are equally probable to be detected in a square

plot(sampling_bias1,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(sampling_bias2,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(sb1ps,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

## https://datadryad.org/stash/dataset/doi:10.5061/dryad.zw3r2287z

r = sb1ps
values(r) = NA

## do this process for each subspecies?
cells=cellFromXY(r,df[,c("longitude","latitude")])
celltab = table(cells)
r2 = r
r2[as.numeric(names(celltab))] = celltab
plot(r2,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
r0 = r
r0[as.numeric(names(celltab))] = 1
plot(r0,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

bias0 = sb1ps*r0
bias0_sum = sum(values(bias0),na.rm=T)
bias0_scl = bias0/bias0_sum*100

plot(bias0,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(bias0_scl,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

## can we multiply the scaled by the number of points in our dataset?
#expected_pts = bias0_scl*sum(values(r2),na.rm=T)
#plot(expected_pts,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
#r3 = r2/expected_pts
#plot(r3,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
## or maybe its a more or less than thing?
#r4 = (r2-expected_pts)/expected_pts
#plot(r4,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

cells_ssp1=cellFromXY(r2,df[df$subspecies=="lepida",c("longitude","latitude")])
cells_ssp2=cellFromXY(r2,df[df$subspecies=="nitens",c("longitude","latitude")])
celltab_ssp1 = table(cells_ssp1)
celltab_ssp2 = table(cells_ssp2)
r20 = r2
values(r20)[!(is.na(values(r20)))]=0
r2_ssp1 = r2
r2_ssp2 = r2
r2_ssp1[as.numeric(names(celltab_ssp1))] = celltab_ssp1
r2_ssp2[as.numeric(names(celltab_ssp2))] = celltab_ssp2
plot(r2_ssp1,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(r2_ssp2,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

r2_ssp1_rel = r2_ssp1/sum(values(r2_ssp1),na.rm=T)*100
r2_ssp2_rel = r2_ssp2/sum(values(r2_ssp2),na.rm=T)*100
plot(r2_ssp1_rel,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(r2_ssp2_rel,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")


## i think these relative values are what need to be scaled by sample size
r2_ssp1_rel_bias = (r2_ssp1_rel*bias0_scl)
r2_ssp2_rel_bias = (r2_ssp2_rel*bias0_scl)
plot(values(r2_ssp1_rel),values(bias0_scl))
abline(0,1,col="red")

par(mfrow=c(1,2))
plot(r2_ssp1_rel_bias,ylim=c(10,50),xlim=c(-130,-90),colNA="grey",col=rwbramp(100),zlim=c(-0.1,0.1))
plot(r2_ssp2_rel_bias,ylim=c(10,50),xlim=c(-130,-90),colNA="grey",col=rwbramp(100),zlim=c(-0.1,0.1))


detectSpatialOutliers()
mgdad = function(localities,epsilon=0.0001) {
  ## Multivariate Gaussian Distribution anomaly detection
  m=nrow(localities) ## the number of rows
  n = ncol(localities) ## the number of columns, aka dimensions
  mu = colMeans(localities,na.rm=T) ## takes the mean of each column
  Sigma = cov(localities) ## get the covariance matrix
  centered <- caret::preProcess(localities, method = "center") ## generates a scaling with which to center the values
  space_mu <- as.matrix(predict(centered, localities)) ## converts the data to the centered values
  A = (2 * pi) ^ (-n / 2) * det(Sigma) ^ (-0.5)
  B = exp(-0.5 * rowSums((space_mu %*% MASS::ginv(Sigma)) * space_mu))
  p_x = A * B ## this is the probability of each value coming from the same normal distribution
  anomalies = which(p_x < epsilon) ## anything where the probability is less than epsilon is considered too improbable
  return(anomalies)
}
df = df[complete.cases(df),]
df_u = unique(df)
outliers=mgdad(localities=df[,c("longitude","latitude")])
outliers_u=mgdad(localities=(df_u[,c("longitude","latitude")]))
plot(df$longitude,df$latitude)
points(df$longitude[outliers],df$latitude[outliers],col="red")
points(df_u$longitude[outliers_u],df_u$latitude[outliers_u],col="blue")


rwbramp = colorRampPalette(c("blue","white","red"))

par(mfrow=c(1,2))
plot(r2_ssp1_rel,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(r2_ssp2_rel,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
rx=(r2_ssp1_rel - r2_ssp2_rel)/(r2_ssp1_rel+r2_ssp2_rel)
plot(rx,ylim=c(10,50),xlim=c(-130,-90),colNA="grey",zlim=c(-1,1),
     col=rwbramp(100))

expected_pts_ssp1 = bias0_scl*sum(values(r2_ssp1),na.rm=T)
expected_pts_ssp2 = bias0_scl*sum(values(r2_ssp2),na.rm=T)
plot(expected_pts_ssp1,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(expected_pts_ssp2,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

r4_ssp1 = (r2_ssp1-expected_pts_ssp1)/expected_pts_ssp1
r4_ssp2 = (r2_ssp2-expected_pts_ssp2)/expected_pts_ssp2
plot(r4_ssp1,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")
plot(r4_ssp2,ylim=c(10,50),xlim=c(-130,-90),colNA="grey")

plot(r2_ssp1)

## need to rescale the data based on the sampling bias map?
## clip the sampling bias map to my data
## divide the clipped sampling bias map by the total sum of pts on clipped map


## might need to thin the data
#df$longitude = round(df$longitude)
#df$latitude = round(df$latitude)
#df = unique(df)
plot(df$longitude,df$latitude,col=as.numeric(as.factor(df$subspecies)),
     pch=as.numeric(as.factor(df$subspecies)),
     cex=as.numeric(as.factor(df$subspecies)))

## after error removal not better so also thin
svmfit = svm(subspecies~longitude+latitude,data=df,
             kernel="radial",scale=F,cost=100)
## reached max iterations? taking a while?
## linear took a long time and didn't work
## radial did some weird stuff
## polynomial was also weird
## sigmoid didn't work

svmfit
plot(x=svmfit,data=df)
