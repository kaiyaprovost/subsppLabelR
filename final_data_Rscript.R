## maybe calculate the density map based on points per grid cell?
## add a "method" param?

## TODO: figure out alaska

detach("package:subsppLabelR", unload = TRUE)
devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
library(subsppLabelR)

EBIRD_KEY = "f49839r87f7g"

## TODO: add support for when too few points are given

## phainopepla
if(!(file.exists("~/Phainopepla_nitens_subspplabelR_RAW.txt"))){
nitens_listFromSubspeciesOcc = subspeciesOccQuery(spp="Phainopepla nitens",
  subsppList=c("lepida","nitens"),pointLimit=10000,
  c("gbif","inat","bison","vertnet"))
nitens_labeledLoc = labelSubspecies(subsppOccList=nitens_listFromSubspeciesOcc)
head(nitens_labeledLoc)
write.table(nitens_labeledLoc,"~/Phainopepla_nitens_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  nitens_labeledLoc = read.table("~/Phainopepla_nitens_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!(file.exists("~/Phainopepla_nitens_subspplabelR_loc_good.txt"))){
nitens = subsppLabelR::databaseToAssignedSubspecies(spp="Phainopepla nitens",
                                                    subsppList = c("lepida","nitens"),
                                                    pointLimit=10000,dbToQuery=c("gbif","inat","bison","vertnet"),
                                                    quantile=0.95,
                                                    #xmin=-125,xmax=-60,ymin=10,ymax=55,
                                                    plotIt=T,
                                                    #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                    datafile = nitens_labeledLoc,
                                                    outputDir="~/")
nitens$loc_suspect
nitens$loc_good
nitens$pol
write.table(nitens$loc_suspect,"~/Phainopepla_nitens_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(nitens$loc_good,"~/Phainopepla_nitens_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}

## sinuatus
if(!(file.exists("~/Cardinalis_sinuatus_subspplabelR_RAW.txt"))){
  sinuatus_listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
                                                    subsppList = c("sinuatus","fulvescens","peninsulae"),
                                                    pointLimit=10000,
                                                    c("gbif","inat","bison","vertnet"))
  sinuatus_labeledLoc = labelSubspecies(subsppOccList=sinuatus_listFromSubspeciesOcc)
  head(sinuatus_labeledLoc)
  write.table(sinuatus_labeledLoc,"~/Cardinalis_sinuatus_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  sinuatus_labeledLoc = read.table("~/Cardinalis_sinuatus_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("~/Cardinalis_sinuatus_subspecies_subspplabelR_loc_good.txt")){
sinuatus = subsppLabelR::databaseToAssignedSubspecies(spp="Cardinalis sinuatus",
                                                      subsppList = c("sinuatus","fulvescens","peninsulae"),
                                                      pointLimit=10000,dbToQuery=c("gbif","inat","bison","vertnet"),
                                                      quantile=0.95, ## works with 0 but does not work well 
                                                      #xmin=-130,xmax=-60,ymin=10,ymax=60,
                                                      plotIt=T,
                                                      datafile = sinuatus_labeledLoc,
                                                      #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                      outputDir="~/")
write.table(sinuatus$loc_suspect,"~/Cardinalis_sinuatus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(sinuatus$loc_good,"~/Cardinalis_sinuatus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}

## melodia
if(!(file.exists("~/Melospiza_melodia_subspplabelR_RAW.txt"))){
  melodia_listFromSubspeciesOcc = subspeciesOccQuery(spp="Melospiza melodia",
                                                     subsppList = c("adusta","amaka","atlantica","beata","caurina",
                                                                    "clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella",
                                                                    "goldmani","gouldii","graminea","heermanni","inexspectata","insignis",
                                                                    "juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia",
                                                                    "merrilli","mexicana","micronyx","montana","morphna","pectoralis",
                                                                    "pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka",
                                                                    "santaecrucis","villai","yuriria","zacapu"),
                                                      pointLimit=10000,
                                                      c("gbif","inat","bison","vertnet"))
  melodia_labeledLoc = labelSubspecies(subsppOccList=melodia_listFromSubspeciesOcc)
  head(melodia_labeledLoc)
  write.table(melodia_labeledLoc,"~/Melospiza_melodia_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  melodia_labeledLoc = read.table("~/Melospiza_melodia_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("~/Melospiza_melodia_subspecies_subspplabelR_loc_good.txt")){
melodia = subsppLabelR::databaseToAssignedSubspecies(spp="Melospiza melodia",
                                                     subsppList = c("adusta","amaka","atlantica","beata","caurina",
                                                                    "clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella",
                                                                    "goldmani","gouldii","graminea","heermanni","inexspectata","insignis",
                                                                    "juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia",
                                                                    "merrilli","mexicana","micronyx","montana","morphna","pectoralis",
                                                                    "pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka",
                                                                    "santaecrucis","villai","yuriria","zacapu"),
                                                     pointLimit=10000,dbToQuery=c("gbif","inat","bison","vertnet"),
                                                     quantile=0.95,
                                                     #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                     plotIt=T,
                                                     datafile=melodia_labeledLoc,
                                                     #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                     outputDir="~/")
write.table(melodia$loc_suspect,"~/Melospiza_melodia_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(melodia$loc_good,"~/Melospiza_melodia_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}
mg = melodia$loc_good
mean_lon = aggregate(mg$longitude~mg$assigned,FUN=function(x){mean(x,na.rm=T)})
sd_lon = aggregate(mg$longitude~mg$assigned,FUN=function(x){sd(x,na.rm=T)})
mean_lat = aggregate(mg$latitude~mg$assigned,FUN=function(x){mean(x,na.rm=T)})
sd_lat = aggregate(mg$latitude~mg$assigned,FUN=function(x){sd(x,na.rm=T)})
means=merge(mean_lon,mean_lat); colnames(means) = c("name","meanlon","meanlat")
sds=merge(sd_lon,sd_lat); colnames(sds) = c("name","sdlon","sdlat")
means_sds = merge(means,sds)
plot(means_sds[,2:3],pch=means_sds[,1],col=as.numeric(as.factor(means_sds[,1])))
segments(x0=means_sds[,2]-means_sds[,4], means_sds[,3], x1 = means_sds[,2]+means_sds[,4], y1 = means_sds[,3],
         col=as.numeric(as.factor(means_sds[,1])))
segments(x0=means_sds[,2], means_sds[,3]-means_sds[,5], x1 = means_sds[,2], y1 = means_sds[,3]+means_sds[,5],
         col=as.numeric(as.factor(means_sds[,1])))

## californianus 
if(!(file.exists("~/Geococcyx_californianus_subspplabelR_RAW.txt"))){
  californianus_listFromSubspeciesOcc = subspeciesOccQuery(spp="Geococcyx californianus",
                                                     pointLimit=10000,
                                                     dbToQuery=c("gbif","inat","bison","vertnet"))
  californianus_labeledLoc = labelSubspecies(subsppOccList=californianus_listFromSubspeciesOcc)
  head(californianus_labeledLoc)
  write.table(californianus_labeledLoc,"~/Geococcyx_californianus_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  californianus_labeledLoc = read.table("~/Geococcyx_californianus_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("~/Geococcyx_californianus_subspecies_subspplabelR_loc_good.txt")){
californianus = subsppLabelR::databaseToAssignedSubspecies(spp="Geococcyx californianus",
                                                           pointLimit=10000,dbToQuery=c("gbif","inat","bison","vertnet"),
                                                           quantile=0.95,
                                                           #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                           datafile=californianus_labeledLoc,
                                                           plotIt=T,
                                                           #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                           outputDir="~/")
write.table(californianus$loc_suspect,"~/Geococcyx_californianus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(californianus$loc_good,"~/Geococcyx_californianus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}
