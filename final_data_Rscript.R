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
                                                      quantile=0.95,
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
                                                           subsppList = c(""),
                                                           pointLimit=10000,dbToQuery=c("gbif","inat","bison","vertnet"),
                                                           quantile=0.95,
                                                           xmin=-180,xmax=-60,ymin=0,ymax=90,
                                                           plotIt=T,
                                                           bgLayer=raster::raster(ext=extent(c(xmin,xmax,ymin,ymax)),nrow=100,ncol=100,vals=0),
                                                           outputDir="~/")
write.table(californianus$loc_suspect,"~/Geococcyx_californianus_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(californianus$loc_good,"~/Geococcyx_californianus_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)




##


spp = "Melospiza melodia"
subsppList = c("adusta","amaka","atlantica","caurina",
               "clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella",
               "goldmani","gouldii","graminea","heermanni","inexspectata","insignis",
               "juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia",
               "merrilli","mexicana","micronyx","montana","morphna",
               "pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka",
               "villai","yuriria","zacapu")

pointLimit=2000
dbToQuery=c("gbif")

print("Downloading species occurrences")
listFromMelospizaOcc = subsppLabelR::subspeciesOccQuery(
  spp = spp,
  subsppList = subsppList,
  pointLimit = pointLimit,
  dbToQuery = dbToQuery
)
## label the data by subspeces
print("Labeling data by subspecies")
melospizaLoc = subsppLabelR::labelSubspecies(subsppOccList = listFromMelospizaOcc)
melospizaLoc = melospizaLoc[!(is.na(melospizaLoc$longitude)),]
melospizaLoc = melospizaLoc[!(is.na(melospizaLoc$latitude)),]

head(melospizaLoc)

max_long = max(melospizaLoc$longitude,na.rm=T)
min_long = min(melospizaLoc$longitude,na.rm=T)
max_lat = max(melospizaLoc$latitude,na.rm=T)
min_lat = min(melospizaLoc$latitude,na.rm=T)

cols = c("red","darkred","purple","mediumpurple","goldenrod","darkgoldenrod","green","darkgreen",
         "blue","darkblue","cyan","darkcyan","black")
palette(cols)
#palette(sample(colors()))

ext = raster::extent(c(min_long,max_long,min_lat,max_lat))
bgLayer = raster::raster(ext=ext,nrow=100,ncol=100,vals=0)

subsppNames = unique(melospizaLoc$subspecies)

raster::plot(bgLayer,
             col = "grey",
             colNA = "darkgrey",
             main = spp)
points(melospizaLoc$longitude,
       melospizaLoc$latitude,
       col = as.numeric(as.factor(melospizaLoc$subspecies)),
       pch = as.numeric(as.factor(melospizaLoc$subspecies)) %% 25)
legend(
  "top",
  legend = as.factor(unique(melospizaLoc$subspecies)),
  bty = "n",
  col = as.numeric(as.factor(unique(melospizaLoc$subspecies))),
  pch = as.numeric(as.factor(unique(melospizaLoc$subspecies))) %% 25)
print("Starting anomaly detection for whole species")

purged_list = c()
kept_list = c()
list_of_anomalies = c()
epsilon = 1e-6

for (i in 0:length(c(subsppNames))) {
  if (i == 0) {
    #print("full")
    #print(nrow(melospizaLoc))
    detectedLocs = subsppLabelR::detectSpatialOutliers(localities = melospizaLoc, epsilon = epsilon)
  }
  else {
    name = subsppNames[[i]]
    if (name != "unknown") {
      #print(name)
      subset = melospizaLoc[melospizaLoc$subspecies == name,]
      #print(nrow(subset))
      detectedLocs = subsppLabelR::detectSpatialOutliers(localities = subset, epsilon = epsilon)
    }
  }
  purged = detectedLocs[[1]]
  anomalies = detectedLocs[[2]]
  kept = detectedLocs[[3]] ## DO NOT KEEP THIS
  #print(length(anomalies))
  
  purged_list = rbind(purged_list, purged)
  list_of_anomalies = c(list_of_anomalies, anomalies)
  kept_list = rbind(kept_list, kept)
}
rows_purged = sort(unique(as.integer(rownames(purged_list))))
head(rows_purged)

print(paste(
  "Removing",
  length(rows_purged),
  "of",
  length(melospizaLoc[, 1]),
  "detected anomalies"
))
removed = melospizaLoc[(rows_purged),]
melospizaLoc = melospizaLoc[-(rows_purged),]

raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = paste("Anomalies")
)
points(melospizaLoc$longitude,
       melospizaLoc$latitude,
       col = "lightgrey",
       pch = 0)
points(removed$longitude,
       removed$latitude,
       col = as.factor(removed$subspecies))
legend(
  "top",
  legend = as.factor(unique(removed$subspecies)),
  pch = 1,
  bty = "n",
  col = as.factor(unique(removed$subspecies))
)

## removing single individual subspecies
print("Removing single-individual subspecies")

for (sub in unique(melospizaLoc$subspecies)) {
  #print(sub)
  rows = (nrow(melospizaLoc[melospizaLoc$subspecies == sub,]))
  if (rows <= 1) {
    print(sub)
    melospizaLoc = melospizaLoc[melospizaLoc$subspecies != sub,]
  }
}
subsppNames = unique(melospizaLoc$subspecies)

max_long = max(melospizaLoc$longitude,na.rm=T)
min_long = min(melospizaLoc$longitude,na.rm=T)
max_lat = max(melospizaLoc$latitude,na.rm=T)
min_lat = min(melospizaLoc$latitude,na.rm=T)
#palette(c("red","blue","black"))
palette(cols)

ext = raster::extent(c(min_long,max_long,min_lat,max_lat))
bgLayer = raster::raster(ext=ext,nrow=100,ncol=100,vals=0)
raster::plot(bgLayer,
             col = "grey",
             colNA = "darkgrey",
             main = spp)
points(melospizaLoc$longitude,
       melospizaLoc$latitude,
       col = as.numeric(as.factor(melospizaLoc$subspecies)),
       pch = as.numeric(as.factor(melospizaLoc$subspecies)))
legend(
  "topright",
  legend = as.factor(unique(melospizaLoc$subspecies)),
  bty = "n",
  col = as.numeric(as.factor(unique(melospizaLoc$subspecies))),
  pch = as.numeric(as.factor(unique(melospizaLoc$subspecies))) %% 25)


xmax = max_long 
xmin = min_long 
ymax = max_lat
ymin = min_lat
quantile=0.95

print("Building species kernel density maps")

## UPDATE THIS CODE
densityRasters = lapply(subsppNames, function(subspp) {
  print(subspp)
  #subspp="amaka"
  locs = melospizaLoc[melospizaLoc$subspecies == subspp,]
  #print(head(locs))
  dens = NULL
  ## TODO: NEED TO CHANGE THIS LINE IN THE PACKAGE
  try({
    dens = subsppLabelR::subspeciesDensityMap(
      localities = locs,
      quantile = quantile,
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    )
  })
  if (is.null(dens) || 
      (length((raster::unique(dens,na.last=NA)))) == 0 ||
      max(values(dens),na.rm=T)==min(values(dens),na.rm=T)) {
    dens = NA
  } else {
    dens = crop(dens,ext)
  }
  
  names(dens) = subspp
  return(dens)
})
#print("done density")
names(densityRasters) = subsppNames

## remove failed ones
densityRasters = densityRasters[!(is.na(densityRasters))]
subsppNames = names(densityRasters)

par(mfrow=c(1,3))
for (i in 1:length(densityRasters)) {
  print(i)
  name = names(densityRasters)[[i]]
  raster::plot(
    bgLayer,
    col = "grey",
    colNA = "darkgrey",
    legend=F,
    main = paste("Density, subspp:", name)
  )
  raster::plot(densityRasters[[i]],
               add = T,
               col = viridis::viridis(99))
}


## convert to polygons
print("Converting density maps to polygons")
## can't handle if there's no data in previous step

## this throws errors if you end up with a density plot that falls outside of normal boundaries
densityPolygons = lapply(densityRasters, function(dens) {
  print(names(dens))
  dens = crop(dens,ext,snap="in")
  densPol = subsppLabelR::densityMapToPolygons(densityMap = dens)
  return(densPol)
})

raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = paste("Polygons, untrimmed")
)
for (i in 1:length(densityPolygons)) {
  raster::plot(
    densityPolygons[[i]],
    add = T,
    border = i,
    lwd = ((3 * i) / 10)
  )
}
legend(
  "top",
  legend = names(densityPolygons),
  bty = "n",
  fill = rgb(0, 0, 0, 0),
  border = 1:length(densityPolygons)
)


## check overlaps between polygons

print("Checking Overlaps of Polygons and Removing Overlaps")

polygonTrimmer_temp = function(polygonList, namesList) {
  newPolygonList = polygonList
  for (slotA in 1:length(namesList)) {
    for (slotB in 1:length(namesList)) {
      if (namesList[[slotA]] != "unknown" &&
          namesList[[slotB]] != "unknown" && slotA != slotB) {
        print(paste(slotA,slotB,sep=" "))
        print(paste(namesList[[slotA]],"with",namesList[[slotB]],sep=" "))
        polA = newPolygonList[[slotA]]
        polB = newPolygonList[[slotB]]
        
        print(class(polA))
        print(class(polB))
        
        # plot(bg,col="grey",colNA="darkgrey")
        # plot(polA,add=T,border="cyan",lwd=7)
        # plot(polB,add=T,border="red",lwd=4)
        # #invisible(readline(prompt="Press [enter] to continue"))
        # if(!is.null(raster::intersect(polA,polB))){
        #   plot(raster::intersect(polA,polB),add=T,lwd=1,border="black")
        # }
        
        
        ## this throws the warnings
        polygonsToRemove = (flagPolygonOverlap2(polA, polB)) ######### CHANGED
        print("CALLING NEW FUNCTION")
        subsppA_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
        subsppB_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
        overlapToRemove_subsppA = polygonsToRemove$subsppA_intToRemove
        overlapToRemove_subsppB = polygonsToRemove$subsppB_intToRemove
        
        subsppA_densityPolygon_trim = trimPolygonsByOverlap(polygon = polA,
                                                            idList = subsppA_polygonsToRemove,
                                                            intList = overlapToRemove_subsppA)
        subsppB_densityPolygon_trim = trimPolygonsByOverlap(polygon = polB,
                                                            idList = subsppB_polygonsToRemove,
                                                            intList = overlapToRemove_subsppB)
        
        ## this is acting up and I do not know why
        ## the problem with this is that it is removing the entirety of the subspecies???
        ## something is going wrong here 
        subsppA = namesList[[slotA]]
        subsppB = namesList[[slotB]]
        if(!(is.null(subsppA_densityPolygon_trim))){newPolygonList[[slotA]] = subsppA_densityPolygon_trim}
        if(!(is.null(subsppB_densityPolygon_trim))){newPolygonList[[slotB]] = subsppB_densityPolygon_trim}
        names(newPolygonList) = names(polygonList)
        print(names(newPolygonList))
        
        #plot(bg,col="grey",colNA="darkgrey")
        #plot(subsppA_densityPolygon_trim,add=T,border="cyan",lwd=7)
        #plot(subsppB_densityPolygon_trim,add=T,border="red",lwd=4)
        #plot(intersect(polA,polB),add=T,lwd=1,border="black")
        
      }
    }
  }
  return(newPolygonList)
}

densityPolygons_trim1 = polygonTrimmer_temp(polygonList = densityPolygons, namesList =
                                                       subsppNames)



for (i in 1:length(densityPolygons_trim1)) {
  name = names(densityPolygons_trim1)[[i]]
  raster::plot(
    bgLayer,
    col = "grey",
    colNA = "darkgrey",
    main = paste("Polygon, subspp:", name)
  )
  raster::plot(densityPolygons_trim1[[i]],
               add = T,
               col = viridis::viridis(99))
}

raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = paste("Polygon, subspp:", name)
)
for (i in 1:length(densityPolygons_trim1)) {
  name = names(densityPolygons_trim1)[[i]]
  raster::plot(
    densityPolygons_trim1[[i]],
    add = T,
    border = i,
    lwd = ((3 * i) / 3)
  )
}
legend(
  "top",
  legend = names(densityPolygons_trim1),
  bty = "n",
  fill = rgb(0, 0, 0, 0),
  border = 1:length(densityPolygons_trim1)
)


densityPolygons_trim = densityPolygons_trim1

print("Locating points relative to polygons")

polyLocations = melospizaLoc

for (slotA in 1:length(subsppNames)) {
  for (slotB in 1:length(subsppNames)) {
    if (subsppNames[[slotA]] != "unknown" &&
        subsppNames[[slotB]] != "unknown" && slotA != slotB) {
      print(paste(slotA,slotB,sep=" "))
      polyLocations = subsppLabelR::locatePolygonPoints(
        test_points = polyLocations,
        polygonA = densityPolygons_trim[[slotA]],
        polygonB = densityPolygons_trim[[slotB]],
        nameA = subsppNames[[slotA]],
        nameB = subsppNames[[slotB]],
        setcoord = T
      )
      
    }
    
  }
  
}

print ("Cleaning up duplicate columns")
colsToDelete = c()

#print(polyLocations)

for (colNumA in 5:length(colnames(polyLocations))) {
  for (colNumB in 6:length(colnames(polyLocations))) {
    if (colNumA < colNumB) {
      print(paste("compare",colNumA,colNumB,sep=" "))
      if (identical(polyLocations[[colNumA]], polyLocations[[colNumB]])) {
        #print("identical, deleting")
        colsToDelete = c(colsToDelete, colNumB)
      }
    }
  }
}
if (!(is.null(colsToDelete))) {
  polyLocations = polyLocations[,-colsToDelete]
}

head(polyLocations)


print("Matching subspecies")
checked = subsppLabelR::subspeciesMatchChecker(locfile = polyLocations, subsppNames =
                                                 subsppNames)
loc_suspect = checked$suspect
loc_good = checked$good
pol = densityPolygons_trim

write.table(loc_suspect,"~/Melospiza_melodia_subspecies_subspplabelR_loc_suspect.txt",sep="\t",row.names = F)
write.table(loc_good,"~/Melospiza_melodia_subspecies_subspplabelR_loc_good.txt",sep="\t",row.names = F)


par(mfrow=c(1,1))

raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = "Unknown or Suspect Poimts"
)
points(loc_suspect$longitude,loc_suspect$latitude)
raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = "Good Points by Subspecies"
)
points(loc_good$longitude,loc_good$latitude,col=as.numeric(as.factor(loc_good$assigned)))

raster::plot(
  bgLayer,
  col = "grey",
  colNA = "darkgrey",
  main = "Final Polygons by Subspecies"
)
lapply(1:length(pol),FUN=function(i){
  plot(pol[[i]],add=T,border=as.numeric(as.factor(names(pol)))[i])
})
