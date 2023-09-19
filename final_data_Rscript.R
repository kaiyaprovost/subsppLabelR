## maybe calculate the density map based on points per grid cell?
## add a "method" param?

## TODO: figure out alaska

detach("package:subsppLabelR", unload = TRUE)
#devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
install.packages("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",repos=NULL,type="source",force=T)
library(subsppLabelR)

## parameters
pointLimit=10000
quant=0.95
dbToQuery=c("gbif","inat","bison","vertnet")
#EBIRD_KEY = "f49839r87f7g"

## TODO: add support for when too few points are given

## DONE
{

## phainopepla
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_RAW.txt"))){
nitens_listFromSubspeciesOcc = subspeciesOccQuery(spp="Phainopepla nitens",
  subsppList=c("lepida","nitens"),pointLimit=pointLimit,
  c("gbif","inat","bison","vertnet"))
nitens_labeledLoc = labelSubspecies(subsppOccList=nitens_listFromSubspeciesOcc)
head(nitens_labeledLoc)
write.table(nitens_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  nitens_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_loc_good.txt"))){
  nitens_labeledLoc$longitude = as.numeric(nitens_labeledLoc$longitude)
  nitens_labeledLoc$latitude = as.numeric(nitens_labeledLoc$latitude)
  nitens_labeledLoc$subspecies = as.factor(nitens_labeledLoc$subspecies)

  nitens = subsppLabelR::databaseToAssignedSubspecies(spp="Phainopepla nitens",
                                                    subsppList = c("lepida","nitens"),
                                                    pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                    quantile=quant,
                                                    #xmin=-125,xmax=-60,ymin=10,ymax=55,
                                                    plotIt=T,
                                                    #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                    datafile = nitens_labeledLoc,
                                                    outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
nitens$loc_suspect
nitens$loc_good
nitens$pol
write.table(nitens$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(nitens$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Phainopepla_nitens_0.95/Phainopepla_nitens_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}

## sinuatus
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspplabelR_RAW.txt"))){
  sinuatus_listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
                                                    subsppList = c("sinuatus","fulvescens","peninsulae"),
                                                    pointLimit=pointLimit,
                                                    c("gbif","inat","bison","vertnet"))
  sinuatus_labeledLoc = labelSubspecies(subsppOccList=sinuatus_listFromSubspeciesOcc)
  head(sinuatus_labeledLoc)
  write.table(sinuatus_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  sinuatus_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspecies_subspplabelR_loc_good.txt")){
sinuatus = subsppLabelR::databaseToAssignedSubspecies(spp="Cardinalis sinuatus",
                                                      subsppList = c("sinuatus","fulvescens","peninsulae"),
                                                      pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                      quantile=quant, ## works with 0 but does not work well
                                                      #xmin=-130,xmax=-60,ymin=10,ymax=60,
                                                      plotIt=T,
                                                      datafile = sinuatus_labeledLoc,
                                                      #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                      outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
write.table(sinuatus$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(sinuatus$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Cardinalis_sinuatus_0.95/Cardinalis_sinuatus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}

## leucophrys
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspplabelR_RAW.txt"))){
  leucophrys_listFromSubspeciesOcc = subspeciesOccQuery(spp="Zonotrichia leucophrys",
                                                      subsppList = c("leucophrys","gambelii","nuttalli",
                                                                     "pugetensis","oriantha"),
                                                      pointLimit=pointLimit,
                                                      c("gbif","inat","bison","vertnet"))
  leucophrys_labeledLoc = labelSubspecies(subsppOccList=leucophrys_listFromSubspeciesOcc)
  head(leucophrys_labeledLoc)
  write.table(leucophrys_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  leucophrys_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspecies_subspplabelR_loc_good.txt")){
  leucophrys = subsppLabelR::databaseToAssignedSubspecies(spp="Zonotrichia leucophrys",
                                                          subsppList = c("leucophrys","gambelii","nuttalli",
                                                                         "pugetensis","oriantha"),                                                        pointLimit=10000,dbToQuery=dbToQuery,
                                                        quantile=quant, ## works with 0 but does not work well
                                                        #xmin=-130,xmax=-60,ymin=10,ymax=60,
                                                        plotIt=T,
                                                        datafile = leucophrys_labeledLoc,
                                                        #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                        outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(leucophrys$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(leucophrys$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Zonotrichia_leucophrys_0.95/Zonotrichia_leucophrys_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## melodia
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspplabelR_RAW.txt"))){
  melodia_listFromSubspeciesOcc = subspeciesOccQuery(spp="Melospiza melodia",
                                                     subsppList = c("adusta","amaka","atlantica","beata","caurina",
                                                                    "clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella",
                                                                    "goldmani","gouldii","graminea","heermanni","inexspectata","insignis",
                                                                    "juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia",
                                                                    "merrilli","mexicana","micronyx","montana","morphna","pectoralis",
                                                                    "pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka",
                                                                    "santaecrucis","villai","yuriria"),
                                                      pointLimit=pointLimit,
                                                      c("gbif","inat","bison","vertnet"))
  melodia_labeledLoc = labelSubspecies(subsppOccList=melodia_listFromSubspeciesOcc)
  head(melodia_labeledLoc)
  write.table(melodia_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  melodia_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspecies_subspplabelR_loc_good.txt")){
melodia = subsppLabelR::databaseToAssignedSubspecies(spp="Melospiza melodia",
                                                     subsppList = c("adusta","amaka","atlantica","beata","caurina",
                                                                    "clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella",
                                                                    "goldmani","gouldii","graminea","heermanni","inexspectata","insignis",
                                                                    "juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia",
                                                                    "merrilli","mexicana","micronyx","montana","morphna","pectoralis",
                                                                    "pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka",
                                                                    "santaecrucis","villai","yuriria","zacapu"),
                                                     pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                     quantile=quant,
                                                     #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                     plotIt=T,
                                                     datafile=melodia_labeledLoc,
                                                     #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                     outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
write.table(melodia$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(melodia$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melospiza_melodia_0.95/Melospiza_melodia_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}
# mg = melodia$loc_good
# mean_lon = aggregate(mg$longitude~mg$assigned,FUN=function(x){mean(x,na.rm=T)})
# sd_lon = aggregate(mg$longitude~mg$assigned,FUN=function(x){sd(x,na.rm=T)})
# mean_lat = aggregate(mg$latitude~mg$assigned,FUN=function(x){mean(x,na.rm=T)})
# sd_lat = aggregate(mg$latitude~mg$assigned,FUN=function(x){sd(x,na.rm=T)})
# means=merge(mean_lon,mean_lat); colnames(means) = c("name","meanlon","meanlat")
# sds=merge(sd_lon,sd_lat); colnames(sds) = c("name","sdlon","sdlat")
# means_sds = merge(means,sds)
# plot(means_sds[,2:3],pch=means_sds[,1],col=as.numeric(as.factor(means_sds[,1])))
# segments(x0=means_sds[,2]-means_sds[,4], means_sds[,3], x1 = means_sds[,2]+means_sds[,4], y1 = means_sds[,3],
#          col=as.numeric(as.factor(means_sds[,1])))
# segments(x0=means_sds[,2], means_sds[,3]-means_sds[,5], x1 = means_sds[,2], y1 = means_sds[,3]+means_sds[,5],
#          col=as.numeric(as.factor(means_sds[,1])))

## californianus
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspplabelR_RAW.txt"))){
  californianus_listFromSubspeciesOcc = subspeciesOccQuery(spp="Geococcyx californianus",
                                                     pointLimit=pointLimit,
                                                     dbToQuery=c("gbif","inat","bison","vertnet"))
  californianus_labeledLoc = labelSubspecies(subsppOccList=californianus_listFromSubspeciesOcc)
  head(californianus_labeledLoc)
  write.table(californianus_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  californianus_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspecies_subspplabelR_loc_good.txt")){
californianus = subsppLabelR::databaseToAssignedSubspecies(spp="Geococcyx californianus",
                                                           pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                           quantile=quant,
                                                           #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                           datafile=californianus_labeledLoc,
                                                           plotIt=T,
                                                           restrictNominate = F,
                                                           #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                           outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
write.table(californianus$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(californianus$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_californianus_0.95/Geococcyx_californianus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## flaviceps
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspplabelR_RAW.txt"))){
  flaviceps_listFromSubspeciesOcc = subspeciesOccQuery(spp="Auriparus flaviceps",
                                                       subsppList = c("acaciarum","flaviceps","hidalgensis","lamprocephalus","ornatus","sinaloae"),
                                                           pointLimit=pointLimit,
                                                           dbToQuery=c("gbif","inat","bison","vertnet"))
  flaviceps_labeledLoc = labelSubspecies(subsppOccList=flaviceps_listFromSubspeciesOcc)
  head(flaviceps_labeledLoc)
  write.table(flaviceps_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  flaviceps_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspecies_subspplabelR_loc_good.txt")){
  flaviceps = subsppLabelR::databaseToAssignedSubspecies(spp="Auriparus flaviceps",
                                                         subsppList = c("acaciarum","flaviceps","hidalgensis","lamprocephalus","ornatus","sinaloae"),
                                                             pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                             quantile=quant,
                                                             #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                             datafile=flaviceps_labeledLoc,
                                                             plotIt=T,
                                                             restrictNominate = F,
                                                             #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                             outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(flaviceps$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(flaviceps$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Auriparus_flaviceps_0.95/Auriparus_flaviceps_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}

## bilineata
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspplabelR_RAW.txt"))){
  bilineata_listFromSubspeciesOcc = subspeciesOccQuery(spp="Amphispiza bilineata",
                                                       subsppList = c("bangsi","bilineata","belvederei","cana","deserticola","carmenae",
                                                                      "grisea","opuntia","pacifica","tortugae"),
                                                       pointLimit=pointLimit,
                                                       dbToQuery=c("gbif","inat","bison","vertnet"))
  bilineata_labeledLoc = labelSubspecies(subsppOccList=bilineata_listFromSubspeciesOcc)
  head(bilineata_labeledLoc)
  write.table(bilineata_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  bilineata_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspecies_subspplabelR_loc_good.txt")){
  bilineata = subsppLabelR::databaseToAssignedSubspecies(spp="Amphispiza bilineata",
                                                         subsppList = c("bangsi","bilineata","belvederei","cana","deserticola","carmenae",
                                                                        "grisea","opuntia","pacifica","tortugae"),
                                                         pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                         quantile=quant,
                                                         #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                         datafile=bilineata_labeledLoc,
                                                         plotIt=T,
                                                         restrictNominate = F,
                                                         #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                         outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(bilineata$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(bilineata$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Amphispiza_bilineata_0.95/Amphispiza_bilineata_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## brunneicapillus
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspplabelR_RAW.txt"))){
  brunneicapillus_listFromSubspeciesOcc = subspeciesOccQuery(spp="Campylorhynchus brunneicapillus",
                                                       subsppList = c("affinis","brunneicapillus","anthonyi","bryanti","guttatus","sandiegensis",
                                                                      "seri"),
                                                       pointLimit=pointLimit,
                                                       dbToQuery=c("gbif","inat","bison","vertnet"))
  brunneicapillus_labeledLoc = labelSubspecies(subsppOccList=brunneicapillus_listFromSubspeciesOcc)
  head(brunneicapillus_labeledLoc)
  write.table(brunneicapillus_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  brunneicapillus_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspecies_subspplabelR_loc_good.txt")){
  brunneicapillus = subsppLabelR::databaseToAssignedSubspecies(spp="Campylorhynchus brunneicapillus",
                                                         subsppList = c("affinis","brunneicapillus","anthonyi","bryanti","guttatus","sandiegensis",
                                                                        "seri"),
                                                         pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                         quantile=quant,
                                                         #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                         datafile=brunneicapillus_labeledLoc,
                                                         plotIt=T,
                                                         restrictNominate = F,
                                                         #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                         outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(brunneicapillus$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(brunneicapillus$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Campylorhynchus_brunneicapillus_0.95/Campylorhynchus_brunneicapillus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## bellii
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspplabelR_RAW.txt"))){
  bellii_listFromSubspeciesOcc = subspeciesOccQuery(spp="Vireo bellii",
                                                       subsppList = c("medius","bellii","arizonae","pusillus"),
                                                       pointLimit=pointLimit,
                                                       dbToQuery=c("gbif","inat","bison","vertnet"))
  bellii_labeledLoc = labelSubspecies(subsppOccList=bellii_listFromSubspeciesOcc)
  head(bellii_labeledLoc)
  write.table(bellii_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  bellii_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspecies_subspplabelR_loc_good.txt")){
  bellii = subsppLabelR::databaseToAssignedSubspecies(spp="Vireo bellii",
                                                         subsppList = c("medius","bellii","arizonae","pusillus"),
                                                         pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                         quantile=quant,
                                                         #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                         datafile=bellii_labeledLoc,
                                                         plotIt=T,
                                                         restrictNominate = F,
                                                         #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                         outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(bellii$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(bellii$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Vireo_bellii_0.95/Vireo_bellii_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}



## crissale
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspplabelR_RAW.txt"))){
  crissale_listFromSubspeciesOcc = subspeciesOccQuery(spp="Toxostoma crissale",
                                                      subsppList = c("coloradense","crissale","dumosum","trinitatis"),
                                                      pointLimit=pointLimit,
                                                      dbToQuery=c("gbif","inat","bison","vertnet"))
  crissale_labeledLoc = labelSubspecies(subsppOccList=crissale_listFromSubspeciesOcc)
  head(crissale_labeledLoc)
  write.table(crissale_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  crissale_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspecies_subspplabelR_loc_good.txt")){
  crissale = subsppLabelR::databaseToAssignedSubspecies(spp="Toxostoma crissale",
                                                        subsppList =  c("coloradense","crissale","dumosum","trinitatis"),
                                                        pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                        quantile=quant,
                                                        #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                        datafile=crissale_labeledLoc,
                                                        plotIt=T,
                                                        restrictNominate = F,
                                                        #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                        outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(crissale$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(crissale$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_crissale_0.95/Toxostoma_crissale_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## curvirostre
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspplabelR_RAW.txt"))){
  curvirostre_listFromSubspeciesOcc = subspeciesOccQuery(spp="Toxostoma curvirostre",
                                                      subsppList = c("celsum","curvirostre","insularum","maculatum",
                                                                     "oberholseri","occidentale","palmeri"),
                                                      pointLimit=pointLimit,
                                                      dbToQuery=c("gbif","inat","bison","vertnet"))
  curvirostre_labeledLoc = labelSubspecies(subsppOccList=curvirostre_listFromSubspeciesOcc)
  head(curvirostre_labeledLoc)
  write.table(curvirostre_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  curvirostre_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspecies_subspplabelR_loc_good.txt")){
  curvirostre = subsppLabelR::databaseToAssignedSubspecies(spp="Toxostoma curvirostre",
                                                        subsppList =   c("celsum","curvirostre","insularum","maculatum",
                                                                         "oberholseri","occidentale","palmeri"),
                                                        pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                        quantile=quant,
                                                        #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                        datafile=curvirostre_labeledLoc,
                                                        plotIt=T,
                                                        restrictNominate = F,
                                                        #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                        outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(curvirostre$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(curvirostre$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Toxostoma_curvirostre_0.95/Toxostoma_curvirostre_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## fusca
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspplabelR_RAW.txt"))){
  fusca_listFromSubspeciesOcc = subspeciesOccQuery(spp="Melozone fusca",
                                                         subsppList = c("campoi","fusca","intermedia","jamesi","mesata","mesoleuca","perpallida","potosina","texana","toroi"),
                                                         pointLimit=pointLimit,
                                                         dbToQuery=c("gbif","inat","bison","vertnet"))
  fusca_labeledLoc = labelSubspecies(subsppOccList=fusca_listFromSubspeciesOcc)
  head(fusca_labeledLoc)
  write.table(fusca_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  fusca_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspecies_subspplabelR_loc_good.txt")){
  fusca = subsppLabelR::databaseToAssignedSubspecies(spp="Melozone fusca",
                                                           subsppList =   c("campoi","fusca","intermedia","jamesi","mesata","mesoleuca","perpallida","potosina","texana","toroi"),
                                                           pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                           quantile=quant,
                                                           #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                           datafile=fusca_labeledLoc,
                                                           plotIt=T,
                                                           restrictNominate = F,
                                                           #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                           outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(fusca$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(fusca$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Melozone_fusca_0.95/Melozone_fusca_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}


## melanura
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspplabelR_RAW.txt"))){
  melanura_listFromSubspeciesOcc = subspeciesOccQuery(spp="Polioptila melanura",
                                                    subsppList = c("curtata","melanura","lucida"),
                                                    pointLimit=pointLimit,
                                                    dbToQuery=c("gbif","inat","bison","vertnet"))
  melanura_labeledLoc = labelSubspecies(subsppOccList=melanura_listFromSubspeciesOcc)
  head(melanura_labeledLoc)
  write.table(melanura_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  melanura_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspecies_subspplabelR_loc_good.txt")){
  melanura = subsppLabelR::databaseToAssignedSubspecies(spp="Polioptila melanura",
                                                      subsppList = c("curtata","melanura","lucida"),
                                                      pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                      quantile=quant,
                                                      #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                      datafile=melanura_labeledLoc,
                                                      plotIt=T,
                                                      restrictNominate = F,
                                                      #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                      outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(melanura$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(melanura$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Polioptila_melanura_0.95/Polioptila_melanura_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}
}

## try at species level
## Geococcyx
dir.create("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/")
if(!(file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspplabelR_RAW.txt"))){
  Geococcyx_listFromSubspeciesOcc = subspeciesOccQuery(spp="Geococcyx",
                                                           subsppList = c("californianus","velox"),
                                                           pointLimit=pointLimit,
                                                           dbToQuery=c("gbif","inat","bison","vertnet"))
  Geococcyx_labeledLoc = labelSubspecies(subsppOccList=Geococcyx_listFromSubspeciesOcc)
  head(Geococcyx_labeledLoc)
  write.table(Geococcyx_labeledLoc,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspplabelR_RAW.txt",sep="\t",row.names = F)
} else {
  Geococcyx_labeledLoc = read.table("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspplabelR_RAW.txt",sep="\t",header=T)
}
if(!file.exists("C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspecies_subspplabelR_loc_good.txt")){
  californianus = subsppLabelR::databaseToAssignedSubspecies(spp="Geococcyx",
                                                             subsppList = c("californianus","velox"),
                                                             pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                             quantile=quant,
                                                             #xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                             datafile=Geococcyx_labeledLoc,
                                                             plotIt=T,
                                                             restrictNominate = F,
                                                             #bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                             outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
  write.table(californianus$loc_suspect,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
  write.table(californianus$loc_good,"C:/Users/kaiya/Documents/GitHub/subsppLabelR/Geococcyx_0.95/Geococcyx_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)
}



## now get the files

occfiles=list.files(path="C:/Users/kaiya/Documents/GitHub/subsppLabelR/",
           pattern="_subspplabelR_loc_good.txt",recursive=T,full.names = T)
wcdata = raster::getData(name="worldclim",download=F,path="~/",var="bio",res=10)

occ=read.table(occfiles[9],header=T)
occ_clean = subsppLabelR::localitiesToNicheMath(Env=wcdata,loc=occ,species=basename(occfiles[9]))


