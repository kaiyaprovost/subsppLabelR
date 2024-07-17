library("plyr")
library("dplyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("caret")
library("h2o")
library("doFuture")
library(AppliedPredictiveModeling)
library(RColorBrewer)
library(ENMTools)
library(raster)
library(sp)
library(ENMeval)
library(rJava)
library(ecospat)

datalist = list.files(path="~/bioacoustics/GBIF/",full.names=T,recursive=T,
                      pattern="_occurrences_subspplabelR.txt$")

for(datafile in datalist){

#datafile="~/Work/Amphispiza-quinquestriata_occurrences_subspplabelR.txt"

## need to add a check in here -- remove unknown points with exact same lat/long as a labeled point

## get environment 
## right now just gonna test with R
if(file.exists("~/bioacoustics/GBIF/worldclim.tif")){
  r1 = stack("~/bioacoustics/GBIF/worldclim.tif")
} else {
  r1 <- getData("worldclim",var="bio",res=10)
  writeRaster(r1,"~/bioacoustics/GBIF/worldclim.tif",format="GTiff")
}

if(file.exists("~/bioacoustics/GBIF/worldclim_1degree.tif")){
  r2 = stack("~/bioacoustics/GBIF/worldclim_1degree.tif")
} else {
  ## convert to 1 degree
  r2 = aggregate(r1,fact=6)
  writeRaster(r2,"~/bioacoustics/GBIF/worldclim_1degree.tif",format="GTiff")
  bg = r2[[1]]
}
## thin to 1 per grid cell 

## import the data
data <- read.csv(datafile,sep=" ")
data = data[,c("longitude","latitude")]
## assign to pixel values 

## thin to one value per pixel 
## first tested with 1 degree, now testing with 10 arc min
cells=cellFromXY(r1,data)
data$cells = cells
localities = as.matrix(data[!duplicated(cells),1:2])

## generate background points

bg1 = (ENMTools::background.points.buffer(localities, radius = 200000,n = 100*nrow(localities), mask = r1[[1]]))
plot(r1[[1]])
points(bg1)
points(localities,col="red",cex=0.3)
bgcells=cellFromXY(r1,bg1)
extracted=extract(r1,bgcells)
backgroundpts = cbind(bg1,extracted)

e.mx = ENMevaluate(occs=localities,bg = backgroundpts[,1:2], envs = r1, 
                   tune.args=list(fc=c("L", "LQ", "H"),rm=seq(0.5,4,0.5)), 
                   parallel=F, numCores=1,algorithm = "maxent.jar",
                   partitions = "block")
res <- eval.results(e.mx)
#opt.aicc <- res %>% filter(delta.AICc == 0)
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
#mod.seq <- eval.models(e.mx)[[opt.seq$tune.args]]
#plot(mod.seq, type = "cloglog")
pred.seq <- eval.predictions(e.mx)[[opt.seq$tune.args]] ## THIS IS THE PREDICTED ENM
plot(pred.seq) 
writeRaster(pred.seq,paste(datafile,"_predicted_raster_enmeval.asc",sep=""),format="ascii")
write.table(res,paste(datafile,"_predicted_results_enmeval.txt",sep=""),sep="\t",row.names = F)

#mod.null <- ENMnulls(e.mx, mod.settings = list(fc = "LQ", rm = 5), no.iter = 100)

## https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html#select
}
