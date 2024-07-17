library(subsppLabelR)

nearestNeighborCleanup = function(good_points,bad_points,k=21){
  pr <- class::knn(train=good_points[,c("longitude","latitude")],
                   test=bad_points[,c("longitude","latitude")],
                   cl=good_points[,c("assigned")],
                   k=k,
                   use.all=T)
  bad_points$KNN = pr
  return(bad_points)
}

good_points = read.table("~/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspplabelR_loc_good.txt",sep="\t",header=T)
bad_points = read.table("~/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspplabelR_loc_suspect.txt",sep="\t",header=T)

plot(good_points[,c("longitude","latitude")],
     col=as.numeric(as.factor(good_points[,c("assigned")])))
points(bad_points[,c("longitude","latitude")],
     col="goldenrod",cex=0.5,pch=0)

newpoints = nearestNeighborCleanup(good_points,bad_points,21)

points(newpoints[,c("longitude","latitude")],
     col=as.numeric(as.factor(newpoints[,c("KNN")])),pch=0)
