library(spocc)
library(knor)
library(sp)
library(EMCluster)
library(viridis)
library(MASS)
library(raster)

Env = raster::stack(list.files(
  path='/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  pattern="\\.bil$",
  full.names=T))
ext = raster::extent(c(-125, -60, 10, 50))
Env = raster::crop(Env, ext)
plot(Env[[1]], col=viridis::viridis(99))

occ1 = occ("Phainopepla nitens",limit=500,has_coords=T,
          from=c("gbif","bison","inat","ebird","ecoengine","vertnet"))
occ2 = occ("Phainopepla nitens nitens",limit=500,has_coords=T,
          from=c("gbif","bison","inat","ebird","ecoengine","vertnet"))
occ3 = occ("Phainopepla nitens lepida",limit=500,has_coords=T,
          from=c("gbif","bison","inat","ebird","ecoengine","vertnet"))

## calculate center of the subspecies
## then get the 90% of points closest to the center

## get centroid of all points
## plot
## take 90% points closest to center
## plot
## draw polygon?

## change occ to points for each
subspp = "UNKNOWN"
occ1d = data.frame(occ2df(occ1))
loc1all = (na.omit(occ1d[,2:3]))
loc1 = unique(loc1all)
loc1lab = cbind(loc1,rep(subspp,length(loc1[,1])))
colnames(loc1lab) = c(colnames(loc1),"subspecies")
plot(loc1)
cent_x1 = mean(loc1$longitude)
cent_y1 = mean(loc1$latitude)
centroid1 = cbind(cent_x1,cent_y1)
points(centroid1,col="red")

## for nitens
subspp = "nitens"
occ2d = data.frame(occ2df(occ2))
loc2all = (na.omit(occ2d[,2:3]))
loc2 = unique(loc2all)
loc2lab = cbind(loc2,rep(subspp,length(loc2[,1])))
colnames(loc2lab) = c(colnames(loc2),"subspecies")
plot(loc2)
cent_x = mean(loc2$longitude)
cent_y = mean(loc2$latitude)
centroid2 = cbind(cent_x,cent_y)
points(centroid2,col="red")

## for lepida
subspp = "lepida"
occ3d = data.frame(occ2df(occ3))
loc3all = na.omit(occ3d[,2:3])
loc3 = unique(loc3all)
loc3lab = cbind(loc3,rep(subspp,length(loc3[,1])))
colnames(loc3lab) = c(colnames(loc3),"subspecies")
plot(loc3)

## full matrix
locmerge = unique(rbind(loc1lab,loc2lab,loc3lab))

#min convex poly
mcp <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}

## with all
plot(loc1,pch=22,col="black",bg="black")
points(loc2,col="red",pch="+")
points(loc3,col="cyan",pch="x")
plot(mcp(loc3),add=T,border="cyan")
plot(mcp(loc2),add=T,border="red")

## with occ2
MCP.thinlocs2 = mcp(loc2)
plot(loc2)
cent_x2 = mean(loc2$longitude)
cent_y2 = mean(loc2$latitude)
centroid2 = cbind(cent_x2,cent_y2)
points(centroid2,col="red")
plot(MCP.thinlocs2, add=T)
## rasterize it?
r2 <- raster(MCP.thinlocs2)
res(r2) = res(r2)/1
r2 = rasterize(MCP.thinlocs2,r2)
plot(r2)
plot(MCP.thinlocs2,add=T)
nc2 <- rasterize(coordinates(loc2), r2, fun='count', background=0)
plot(nc2)
plot(MCP.thinlocs2, add=TRUE)

## with occ 3
MCP.thinlocs3 = mcp(loc3)
plot(loc3)
cent_x3 = mean(loc3$longitude)
cent_y3 = mean(loc3$latitude)
centroid3 = cbind(cent_x3,cent_y3)
points(centroid3,col="red")
plot(MCP.thinlocs3, add=T)
## rasterize it?
r3 <- raster(MCP.thinlocs3)
res(r3) = res(r3)/1
r3 = rasterize(MCP.thinlocs3,r3)
plot(r3,col=viridis::viridis(99))
plot(MCP.thinlocs3,add=T)
nc3 <- rasterize(coordinates(loc3), r3, fun='count', background=0)
plot(nc3,col=viridis::viridis(99))
plot(MCP.thinlocs3, add=TRUE)

## points that are in the other polygon
pt3 = loc3
coordinates(pt3) = ~longitude+latitude
proj4string(pt3) = CRS("+proj=longlat")
proj4string(MCP.thinlocs3) = CRS("+proj=longlat")
pt2 = loc2
coordinates(pt2) = ~longitude+latitude
proj4string(pt2) = CRS("+proj=longlat")
proj4string(MCP.thinlocs2) = CRS("+proj=longlat")
pt1 = loc1
coordinates(pt1) = ~longitude+latitude
proj4string(pt1) = CRS("+proj=longlat")
#proj4string(MCP.thinlocs1) = CRS("+proj=longlat")

## which of subspp 2 is in subspp 3 polygon
#(over(pt2,MCP.thinlocs3))
inOtherPolygon2 = loc2[which(!is.na(over(pt2,MCP.thinlocs3))),]
inCorrectPolygon2 = loc2[which(is.na(over(pt2,MCP.thinlocs3))),]

## which of subspp 3 is in subspp 2 polygon
#(over(pt3,MCP.thinlocs2))
inOtherPolygon3 = loc3[which(!is.na(over(pt3,MCP.thinlocs2))),]
inCorrectPolygon3 = loc3[which(is.na(over(pt3,MCP.thinlocs2))),]

plot(locmerge[,1:2],col="grey",pch=16)
points(inCorrectPolygon3,col="darkcyan",pch=4)
points(inCorrectPolygon2,col="cyan",pch=3)
points(inOtherPolygon3,col="darkred",pch=4)
points(inOtherPolygon2,col="red",pch=3)

range = c(-125, -60, 10, 50)
ext = raster::extent(range)
w1 = matrix(1,3,3)
x1 = kde2d(loc1$longitude, loc1$latitude, #n=length(loc1$longitude)/20,
           lims=range)
r1 = raster(x1)
quan1 = quantile(r1[r1],0.95)
r1q = r1
r1q[r1q <= quan1] = NA; plot(r1q)
#r1[r1 <= 0.000] <- NA
w3 = matrix(1,3,3)
x3 = kde2d(loc3$longitude, loc3$latitude, #n=length(loc3$longitude)/20,
           lims=range)
r3 = raster(x3)
quan3 = quantile(r3[r3],0.95)
r3q = r3
r3q[r3q <= quan3] = NA; plot(r3q)
#r3[r3 <= 0.000] <- NA
w2 = matrix(1,3,3)
x2 = kde2d(loc2$longitude, loc2$latitude, #n=length(loc2$longitude)/20,
           lims=range)
r2 = raster(x2)
quan2 = quantile(r2[r2],0.95)
r2q = r2
r2q[r2q <= quan2] = NA; plot(r2q)

plot(Env[[1]], col="black")
plot(r1q,add=T,col=viridis::viridis(99))
plot(r3q,add=T,col=viridis::viridis(99))
plot(r2q,add=T,col=viridis::viridis(99))

#r2[r2 <= 0.000] <- NA
#r = aggregate(r,fact=2)

## get polygons
qpoly1 = r1q
qpoly1[!is.na(qpoly1)] = 1
qpoly1 <- disaggregate(rasterToPolygons(qpoly1, fun=NULL, na.rm=T,dissolve=T))
plot(qpoly1)

qpoly2 = r2q
qpoly2[!is.na(qpoly2)] = 1
qpoly2 <- disaggregate(rasterToPolygons(qpoly2, fun=NULL, na.rm=T,dissolve=T))
plot(qpoly2)

qpoly3 = r3q
qpoly3[!is.na(qpoly3)] = 1
qpoly3 <- disaggregate(rasterToPolygons(qpoly3, fun=NULL, na.rm=T,dissolve=T))
plot(qpoly3)

plot(Env[[1]], col="black")
plot(qpoly1,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(qpoly2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(qpoly3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)

## which of subspp 2 is in subspp 3 polygon
proj4string(qpoly2) = CRS("+proj=longlat")
proj4string(qpoly3) = CRS("+proj=longlat")
qpoly2 = spTransform(qpoly2,CRS("+proj=longlat"))
qpoly3 = spTransform(qpoly3,CRS("+proj=longlat"))

#(over(pt2,qpoly3))
inOtherPolygon2 = loc2[which(!is.na(over(pt2,qpoly3))),]
inCorrectPolygon2 = loc2[which(is.na(over(pt2,qpoly3))),]

## which of subspp 3 is in subspp 2 polygon
#(over(pt3,qpoly2))
inOtherPolygon3 = loc3[which(!is.na(over(pt3,qpoly2))),]
inCorrectPolygon3 = loc3[which(is.na(over(pt3,qpoly2))),]

plot(Env[[1]], col="black")
points(locmerge[,1:2],col="darkgrey",pch=16)
plot(qpoly1,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(qpoly2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(qpoly3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)
points(inCorrectPolygon3,col="darkcyan",pch=4)
points(inCorrectPolygon2,col="darkred",pch=3)
points(inOtherPolygon3,col="cyan",pch=4)
points(inOtherPolygon2,col="red",pch=3)


plot(gIntersection(qpoly2,qpoly3,byid=T))
badList_i = c()
badList_j = c()
for (i in range(1,length(qpoly2))){
  for(j in range(1,length(qpoly3))) {
    intersect = gIntersection(qpoly2[i,],qpoly3[j,],byid=T)
    #print(length(intersect))
    if (length(intersect) != 0){
      plot(Env[[1]], col="black")
      plot(qpoly2[i,],border="red",add=T)
      plot(qpoly3[j,],border="cyan",add=T)
      areaInt = gArea(intersect)
      area2 = gArea(qpoly2[i,])
      area3 = gArea(qpoly3[j,])

      if (areaInt >= area2) {
        #print("remove area 2")
        badList_i = c(badList_i,i)
      }
      if (areaInt >= area3) {
        #print("remove area 3")
        badList_j = c(badList_j,j)
      }

    }


    #invisible(readline(prompt="Press [enter] to continue"))
  }
}

fullspp = qpoly1
subspp2 = qpoly2
subspp3 = qpoly3

if(length(badList_i) > 0) {
  subspp2 = subspp2[-badList_i,]
}
if(length(badList_j) > 0) {
  subspp3 = subspp3[-badList_j,]
}

plot(Env[[1]], col="black")
plot(fullspp,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(subspp2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(subspp3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)

## redo points to see which are where
## which of subspp 2 is in subspp 3 polygon
proj4string(subspp2) = CRS("+proj=longlat")
proj4string(subspp3) = CRS("+proj=longlat")
subspp2 = spTransform(subspp2,CRS("+proj=longlat"))
subspp3 = spTransform(subspp3,CRS("+proj=longlat"))

#(over(pt2,subspp3))
inOtherPolygon2 = loc2[which(!is.na(over(pt2,subspp3))),]
inCorrectPolygon2 = loc2[which(!is.na(over(pt2,subspp2))),]
x = loc2[which(is.na(over(pt2,subspp2))),]
ptx = x
coordinates(ptx) = ~longitude+latitude
proj4string(ptx) = CRS("+proj=longlat")
inNoPolygon2 = x[which(is.na(over(ptx,subspp3))),]

## which of subspp 3 is in subspp 2 polygon
#(over(pt3,subspp2))
inOtherPolygon3 = loc3[which(!is.na(over(pt3,subspp2))),]
inCorrectPolygon3 = loc3[which(!is.na(over(pt3,subspp3))),]
x = loc3[which(is.na(over(pt3,subspp3))),]
ptx = x
coordinates(ptx) = ~longitude+latitude
proj4string(ptx) = CRS("+proj=longlat")
inNoPolygon3 = x[which(is.na(over(ptx,subspp2))),]

plot(Env[[1]], col="black")
points(locmerge[,1:2],col="darkgrey",pch=16)
#plot(fullspp,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(subspp2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(subspp3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)
points(inCorrectPolygon3,col="darkcyan",pch=4)
points(inCorrectPolygon2,col="darkred",pch=3)
points(inOtherPolygon3,col="cyan",pch=4)
points(inOtherPolygon2,col="red",pch=3)
points(inNoPolygon2,col="pink",pch=3)
points(inNoPolygon3,col="blue",pch=4)

## remove points in the wrong polygon
cleanloc2 = loc2[which(!(rownames(loc2) %in% rownames(inOtherPolygon2))),]
cleanloc3 = loc3[which(!(rownames(loc3) %in% rownames(inOtherPolygon3))),]

plot(Env[[1]], col="black")
points(locmerge[,1:2],col="darkgrey",pch=16)
#plot(fullspp,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(subspp2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(subspp3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)
points(cleanloc3,col="darkcyan",pch=4)
points(cleanloc2,col="darkred",pch=3)

cleanloc2lab = cleanloc2
cleanloc2lab$subspecies = loc2lab$subspecies[1]
cleanloc2lab$how_subspp_determined = "GBIF"
cleanloc3lab = cleanloc3
cleanloc3lab$subspecies = loc3lab$subspecies[1]
cleanloc3lab$how_subspp_determined = "GBIF"

## get points that are unlabeled and put in polygons
pt1 = loc1
coordinates(pt1) = ~longitude+latitude
proj4string(pt1) = CRS("+proj=longlat")
unlabeledInPolygon3 = loc1[which(!is.na(over(pt1,subspp3))),]
unlabeledInPolygon2 = loc1[which(!is.na(over(pt1,subspp2))),]

plot(Env[[1]], col="black")
points(locmerge[,1:2],col="darkgrey",pch=16)
#plot(fullspp,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(subspp2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(subspp3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)
points(unlabeledInPolygon3,col="darkcyan",pch=4)
points(unlabeledInPolygon2,col="darkred",pch=3)

addSubspp2 = loc1[which((rownames(loc1) %in% rownames(unlabeledInPolygon2))),]
addSubspp2$subspecies = cleanloc2lab$subspecies[1]
addSubspp2$how_subspp_determined = "assumed"
addSubspp3 = loc1[which((rownames(loc1) %in% rownames(unlabeledInPolygon3))),]
addSubspp3$subspecies = cleanloc3lab$subspecies[1]
addSubspp3$how_subspp_determined = "assumed"

plot(Env[[1]], col="black")
points(locmerge[,1:2],col="darkgrey",pch=16)
#plot(fullspp,add=T,col=rgb(1,1,1,0),border="grey",lwd=7)
plot(subspp2,add=T,col=rgb(1,0,0,0),border="red",lwd=4)
plot(subspp3,add=T,col=rgb(0,1,1,0),border="cyan",lwd=1)
points(addSubspp3,col="darkcyan",pch=4)
points(addSubspp2,col="darkred",pch=3)

cleanlocMerge = unique(rbind(cleanloc2lab,cleanloc3lab,addSubspp2,addSubspp3))
plot(Env[[1]], col="grey")
points(cleanlocMerge,col=cleanlocMerge$subspecies,
       pch=cleanlocMerge$how_subspp_determined)


# f <- function(X) max(X, na.rm=TRUE)
# localmax <- focal(r, w, fun = f, pad=TRUE, padValue=NA)
# r2 <- r==localmax
# maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
# image(x,asp=1)
# points(maxXY)

## calculate based on eucledian distance only
distance = pointDistance(p1=c(cent_x,cent_y), p2=loc2,
              lonlat=T, allpairs=FALSE)
hist(distance)
quantile(distance,seq(0,1,0.05))
distloc = cbind(loc2,cbind(distance))
names(distloc)=c(names(loc2),"distance")
distorder = distloc[order(distloc$distance),]

loc2trim = distorder[distorder$distance<=quantile(distorder$distance,0.9),]
plot(loc2trim$longitude,loc2trim$latitude)
cent_x = mean(loc2trim$longitude)
cent_y = mean(loc2trim$latitude)
centroid1 = cbind(cent_x,cent_y)
points(centroid2,col="red")
points(centroid1,col="blue")

## can you do this, make a polygon around the points, pull anything inside the polygon?
## add a buffer just in case?

## kmeans
means2 = Kmeans(as.matrix(loc2),centers=2,tolerance=1e-10,
                iter.max=10,dist.type="eucl")
clus2 = means2$cluster
plot(loc2,col=clus2)

## testing EM cluster -- unsupervised
{
 library(EMCluster, quiet = TRUE)
 set.seed(1234)
 x <- da1$da
 TC <- da1$class
 n <- nrow(x)
 p <- ncol(x)
 k <- 10
 ret.em <- init.EM(x, nclass = k, method = "em.EM")
 ret.Rnd <- init.EM(x, nclass = k, method = "Rnd.EM", EMC = .EMC.Rnd)
 ret.Rndp <- init.EM(x, nclass = k, method = "Rnd.EM", EMC = .EMC.Rndp)
 ret.svd <- emgroup(x, nclass = k)
 par(mfrow = c(2, 2))
 plotem(ret.em, x, main = "em")
 plotem(ret.Rnd, x, main = "Rnd")
 plotem(ret.Rndp, x, main = "Rnd+")
 plotem(ret.svd, x, main = "svd")
 ret.all <-
   cbind(
       c(ret.em$llhdval, ret.Rnd$llhdval, ret.Rndp$llhdval,
               ret.svd$llhdval),
       c(RRand(ret.em$class, TC)$adjRand,
               RRand(ret.Rnd$class, TC)$adjRand,
               RRand(ret.Rndp$class, TC)$adjRand,
               RRand(ret.svd$class, TC)$adjRand)
     )
 rownames(ret.all) <- c("em", "Rnd", "Rnd+", "svd")
 colnames(ret.all) <- c("logL", "adjR")
 ret.all
}

 emClustPicker = function(maxk=2,loc,iter=100,label){
  #maxk=3
  kval = seq(1,maxk,1)
   #k = 2
   #loc=loc2
   #loc=loc2trim[,1:2]
   #iter=10000
   a = lapply(kval,function(x,d=loc,numRandomIter=iter){
     #print(x)
     k = as.character(x)
     ret.em <- init.EM(d, nclass = x, method = "em.EM",min.n.iter=numRandomIter)
     ret.Rnd <- init.EM(d, nclass = x, method = "Rnd.EM", EMC = .EMC.Rnd,min.n.iter=numRandomIter)
     ret.Rndp <- init.EM(d, nclass = x, method = "Rnd.EM", EMC = .EMC.Rndp,min.n.iter=numRandomIter)
     ret.init <- emcluster(d, emobj, assign.class = TRUE)
     ret.svd <- emgroup(d, nclass = x)
     png(paste("/Users/kprovost/Documents/Classes/Spatial Bioinformatics/",label,"SubspeciesDifferences_k",
               as.character(k),".png",sep=""))
     par(mfrow=c(2,2))
     plotem(ret.em, d, main = "em")
     plotem(ret.Rnd, d, main = "Rnd")
     plotem(ret.Rndp, d, main = "Rnd+")
     plotem(ret.svd, d, main = "svd")
     dev.off()

     aics = rbind(em.aic(d,ret.em),
              em.aic(d,ret.Rnd),
              em.aic(d,ret.Rndp),
              em.aic(d,ret.svd))
     rownames(aics) = c("em","Rnd","Rnd+","svd")
     colnames(aics) = k
     ## check if all of these are the same
     if(x!=1){
       groupdif=FALSE
       chk = cbind(d,ret.em$class,ret.Rnd$class,
                   ret.Rndp$class,ret.svd$class)

       for(i in 1:x){
         ch = chk[chk$`ret.em$class`==i,]
         l = sum(sapply(3:6, function(x){length(unique(ch[,x]))}))
         if(l!=4){
           groupdif=TRUE
           #print(paste("group ",i," is NOT same across models for k=",x,sep=""))
         }
       }
       if(groupdif==FALSE){
         #print("yay")
         bestAIC = 1
         suffix = "allsameAIC"
       }
       else {
         bestAIC = which(aics==min(aics))
         suffix = "lowestAIC"
       }
     }
     else {
       bestAIC=1
       suffix = "k1"
     }
     bestmodel = paste(rownames(aics)[bestAIC])
     finalaic = aics[bestAIC]
     returnaic = cbind(as.numeric(finalaic),bestmodel,suffix)
     colnames(returnaic) = c("best AIC","best model","howpicked")
     return(returnaic)
     })
  a = do.call(rbind,a)
  rownames(a) = paste("k",as.character(1:maxk),sep="")
  #a
  return(a)
}

 loc2clust = emClustPicker(maxk=2,loc=loc2,iter=100,label="Phainopepla_nitens_nitens")
 loc3clust = emClustPicker(maxk=2,loc=loc3,iter=100,label="Phainopepla_nitens_lepida")



 kval = seq(1,5,1)
 test = sapply(kval,function(x,y=loc2){
   emobj = simple.init(y,nclass=x)
   return(em.aic(y,emobj))
 })
deltatest = test-min(test)
ktest = which(test==min(test))

x = loc2
k = ktest

 ret.em <- init.EM(x, nclass = 2, method = "em.EM")
 ret.Rnd <- init.EM(x, nclass = k, method = "Rnd.EM", EMC = .EMC.Rnd)
 emobj <- simple.init(x, nclass = k)
 ret.init <- emcluster(x, emobj, assign.class = TRUE)

 par(mfrow = c(2, 2))
 plotem(ret.em, x)
 plotem(ret.Rnd, x)
 plotem(ret.init, x)

 em.aic(x,emobj)

 ## EM -- supervised
 library(EMCluster, quiet = TRUE)
 set.seed(1234)
 x <- da1$da
 TC <- da1$class
 lab <- da1$class
 n <- nrow(x)
 p <- ncol(x)
 k <- 10
 nk <- floor(k * 2 / 3)
 tmp <- names(sort(table(lab), decreasing = TRUE))[1:nk]
 lab[!(lab %in% tmp)] <- 0
 for(i in 1:nk){
   tmp.id <- lab == tmp[i]
   id <- sample(which(tmp.id), sum(tmp.id) * 0.5, replace = FALSE)
   tmp.id[id] <- FALSE
   lab[tmp.id] <- 0
 }
 index.lab <- sort(unique(lab))
 lab <- sapply(lab, function(i) which(index.lab == i) - 1)
 ret.em <- init.EM(x, nclass = k, lab = lab, method = "em.EM")
 ret.Rnd <- init.EM(x, nclass = k, lab = lab, method = "Rnd.EM",
                    EMC = .EMC.Rnd)
 ret.Rndp <- init.EM(x, nclass = k, lab = lab, method = "Rnd.EM",
                     EMC = .EMC.Rndp)
 par(mfrow = c(2, 2))
 plotem(ret.em, x, main = "em")
 plotem(ret.Rnd, x, main = "Rnd")
 plotem(ret.Rndp, x, main = "Rnd+")
 ret.all <-
   cbind(
     c(ret.em$llhdval, ret.Rnd$llhdval, ret.Rndp$llhdval),
     c(RRand(ret.em$class, TC, lab = lab)$adjRand,
       RRand(ret.Rnd$class, TC, lab = lab)$adjRand,
       RRand(ret.Rndp$class, TC, lab = lab)$adjRand)
   )
 rownames(ret.all) <- c("em", "Rnd", "Rnd+")
 colnames(ret.all) <- c("logL", "adjR")
 ret.all
