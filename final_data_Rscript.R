## TODO: update so that the density map gets plotted onto the environmental raster of choice?
## maybe calculate the density map based on points per grid cell?
## add a "method" param?

detach("package:subsppLabelR", unload = TRUE)
devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
library(subsppLabelR)

EBIRD_KEY = "f49839r87f7g"

## phainopepla
## TODO: add ebird support, currently not working
## TODO: add support for when too few points are given

if(!(file.exists("~/Phainopela_nitens_subspplabelR_RAW.txt"))){
nitens_listFromSubspeciesOcc = subspeciesOccQuery(spp="Phainopepla nitens",
  subsppList=c("lepida","nitens"),pointLimit=10000,
  c("gbif","inat","bison","vertnet"))
nitens_labeledLoc = labelSubspecies(subsppOccList=nitens_listFromSubspeciesOcc)
head(nitens_labeledLoc)
write.table(nitens_labeledLoc,"~/Phainopela_nitens_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  nitens_labeledLoc = read.table("~/Phainopela_nitens_subspplabelR_RAW.txt",sep="\t",header=T)
}


if(!(file.exists("~/Phainopela_nitens_subspplabelR_loc_good.txt"))){
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
write.table(nitens$loc_suspect,"~/Phainopela_nitens_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(nitens$loc_good,"~/Phainopela_nitens_subspplabelR_loc_good.txt",sep="\t",row.names = F)
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

## californianus 

## this does not work at the "matching subspecies" stage
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




##
