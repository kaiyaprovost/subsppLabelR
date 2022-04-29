library(subsppLabelR)

good_points = read.table("~/Phainopepla_nitens_subspplabelR_loc_good.txt",sep="\t",header=T)
bad_points = read.table("~/Phainopepla_nitens_subspplabelR_loc_suspect.txt",sep="\t",header=T)

plot(good_points[,c("longitude","latitude")],
     col=as.numeric(as.factor(good_points[,c("assigned")])))

library(class)
pr <- knn(train=good_points[,c("longitude","latitude")],
          test=bad_points[,c("longitude","latitude")],
          cl=good_points[,c("assigned")],
          k=20,
          use.all=T)
