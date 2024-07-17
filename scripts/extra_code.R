# if(method=="polygon") {
#print("endplot1")
## convert to polygons
print("Converting density maps to polygons")
## can't handle if there's no data in previous step
densityPolygons = lapply(densityRasters, function(dens) {
  #print(names(dens))
  densPol = NULL
  try({densPol = densityMapToPolygons(densityMap = dens)})
  return(densPol)
})

## optionally restrict the nominate
if(restrictNominate==T){
  print("Restricting the nominate")
  polygons_notnom = densityPolygons[!(names(densityPolygons) %in% c(nominateSubspecies,"unknown"))]
  #polygons_notnom = densityPolygons[!(names(densityPolygons) %in% c(nominateSubspecies))]
  fullpoly=raster::bind(polygons_notnom)
  if(length(fullpoly)==1){fullpoly=fullpoly[[1]]}
  densityPolygons[[nominateSubspecies]]=rgeos::gDifference(densityPolygons[[nominateSubspecies]], fullpoly)
}

#print(densityPolygons)
if (plotIt == T) {
  for (i in 1:length(densityPolygons)) {
    name = names(densityPolygons)[[i]]
    png(paste("RawDensityPolygon_", spp, " ", name,quant, ".png", sep = ""))
    #pdf(paste("RawDensityPolygon_", spp, " ", name,quant, ".pdf", sep = ""))
    raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
    #sp::spplot(densityPolygons[[i]],add = T,col = "red")
    plot(densityPolygons[[i]],add=T,col="red")
    dev.off()
  }
  png(paste("RawDensityPolygon_", spp, quant," ALL.png", sep = ""))
  #pdf(paste("RawDensityPolygon_", spp, quant," ALL.pdf", sep = ""))
  raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = "Polygon, ALL")
  cols = c("black","red", "blue", "green",
           "cyan", "magenta", "pink", "white",
           "purple", "orange", "yellow", "sienna",
           "thistle", "palegreen", "powderblue", "aquamarine",
           "violet", "mediumslateblue", "lightsalmon", "lightblue")
  for (i in 1:length(densityPolygons)) {
    name = names(densityPolygons)[[i]]
    #sp::spplot(densityPolygons[[i]],add = T,border = cols[i],lwd = ((3 * i) / 3))
    plot(densityPolygons[[i]],add = T,border = cols[i],lwd = ((3 * i) / 3))
  }
  legend("top",legend = names(densityPolygons),bty = "n",fill = rgb(0, 0, 0, 0),border = cols)
  dev.off()
}
## check overlaps between polygons
print("Checking Overlaps of Polygons and Removing Overlaps")
## remove polygons that are completely within other polygon
## TODO: what about things that are in neither polgon?
## TODO: what about things that are in both?
## there is a bug -- if one subspp range is entirely subsumed within another polygon,
## will delete that subspecies. no bueno
## TODO: nominate subspecies special case

densityPolygons_trim1 = polygonTrimmer2(polygonList = densityPolygons, namesList = subsppNames)
if (plotIt == T) {
  for (i in 1:length(densityPolygons_trim1)) {
    name = names(densityPolygons_trim1)[[i]]
    png(paste("TrimDensityPolygon_", spp, " ", name,quant, ".png", sep = ""))
    raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
    raster::plot(densityPolygons_trim1[[i]],add = T,col = viridis::viridis(99))
    dev.off()
  }
  png(paste("TrimDensityPolygon_", spp, quant," ALL.png", sep = ""))
  raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
  cols = c( "black", "red", "blue", "green", "cyan", "magenta",
            "pink", "white", "purple", "orange", "yellow", "sienna",
            "thistle", "palegreen", "powderblue", "aquamarine", "violet", "mediumslateblue",
            "lightsalmon", "lightblue")
  for (i in 1:length(densityPolygons_trim1)) {
    print(i)
    name = names(densityPolygons_trim1)[[i]]
    raster::plot(densityPolygons_trim1[[i]],add = T,border = cols[i],lwd = ((3 * i) / 3))
  }
  legend("top",legend = names(densityPolygons_trim1),bty = "n",fill = rgb(0, 0, 0, 0),border = cols)
  dev.off()
}
#densityPolygons_trim = closestPolygonFunction(listOfPolygons=densityPolygons_trim1)
##TODO: get the closestPolygonFunction working
densityPolygons_trim = densityPolygons_trim1
## this isn't working!!!!!!
##TODO: figure out how to remove small polygons that are closer to other subspp than their own
## remove points that are in wrong polygon
## if labeled and in wrong polygon, unlabel
## label unlabeled points based on polygons
print("Locating points relative to polygons")
## TODO: this is hanging
polyLocations = labeledLoc
## this is taking a long time
## we are gonna try something else
## iterate through each polygon
#print(densityPolygons_trim)
#print(subsppNames)
for(polygonSlot in 1:length(subsppNames)){
  if (subsppNames[[polygonSlot]] != "unknown"){
    print(polygonSlot)
    print(subsppNames[[polygonSlot]])
    try({
      print(densityPolygons_trim[[polygonSlot]])
      polygonA=densityPolygons_trim[[polygonSlot]]
      polygonA = as(polygonA, "SpatialPolygonsDataFrame")
      nameA = subsppNames[[polygonSlot]]
      print(polygonA) ## this needs to be a spatial polygon dataframe
      polyLocations = locatePolygonPoints2(test_points = polyLocations,
                                           polygonA = polygonA,
                                           nameA = nameA,
                                           setcoord = T)
    })
  }
}
#
# for (slotA in 1:length(subsppNames)) {
#   for (slotB in 1:length(subsppNames)) {
#     if (subsppNames[[slotA]] != "unknown" && subsppNames[[slotB]] != "unknown" && slotA != slotB) {
#       print(paste(slotA,slotB,sep=" "))
#       polyLocations = locatePolygonPoints(test_points = polyLocations,
#                                           polygonA = densityPolygons_trim[[slotA]],polygonB = densityPolygons_trim[[slotB]],
#                                           nameA = subsppNames[[slotA]],nameB = subsppNames[[slotB]],
#                                           setcoord = T)
#     }
#   }
# }
print ("Cleaning up duplicate columns")
## this does not work with only one species
colsToDelete = c()
#print(polyLocations)
print(length(colnames(polyLocations)))
if(length(colnames(polyLocations)) > 6){
  for (colNumA in 5:length(colnames(polyLocations))) {
    for (colNumB in 6:length(colnames(polyLocations))) {
      if (colNumA < colNumB) {
        #print(paste("compare",colNumA,colNumB,sep=" "))
        if (identical(polyLocations[[colNumA]], polyLocations[[colNumB]])) {
          #print("identical, deleting")
          colsToDelete = c(colsToDelete, colNumB)
        }
      }
    }
  }
  if (!(is.null(colsToDelete))) {
    #print("is null cols")
    #print(colsToDelete)
    #print(names(polyLocations))
    #print(head(polyLocations))
    polyLocations = polyLocations[,-colsToDelete]
    #print("success")
  }
  # }
  # }

  # if(method=="density") {
  #
  #   density_stack = stack(densityRasters[names(densityRasters)!="unknown"])
  #   density_max = max(density_stack,na.rm=T)
  #   density_sum = sum(density_stack,na.rm=T)
  #   max_density_list = lapply(1:length(densityRasters),FUN=function(i){
  #     ras = densityRasters[[i]]
  #     #plot(ras,colNA="grey")
  #     ras2 = density_max-ras
  #     values(ras2)[values(ras2)!=0]=NA
  #     values(ras2)[values(ras2)==0]=i
  #     #plot(ras2,colNA="grey")
  #     return(ras2)
  #   })
  #   names(max_density_list) = names(densityRasters)
  #   max_density_stack = stack(max_density_list[names(max_density_list)!="unknown"])
  #   #plot(max_density_stack)
  #   max_density = sum(max_density_stack,na.rm=T)
  #   plot(max_density)
  #
  # }



  ## TODO: see if you can label unlabeled points based on polygons or nearest neighbor
  ## or nearest neighbor
  ## output the new data
  ## KNN algorithm, k nearest neighbors
  ## plot (optional)
  ## workflow
  ## go through data and look for instances where:
  ## something is assigned to multiple subspecies
  ## or subspecies assignment a priori does not match final
  ## TODO: consider putting this in the other script file
  ## not working right now
  print("Matching subspecies")
  checked = subspeciesMatchChecker(locfile = polyLocations, subsppNames =
                                     subsppNames)
  #print("c1")
  checked_suspect = checked$suspect
  #print("c2")
  checked_good = checked$good
  #print("done")
  ## return nice clean data
  print("Warning: no valid definition for subspecies given!") ## this is a joke
  return(list(loc_suspect = checked_suspect,loc_good = checked_good,pol = densityPolygons_trim))
}
# Env = raster::stack(list.files(
#   path='/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
#   pattern="\\.bil$",
#   full.names=T))
# ext = raster::extent(c(-125,-60,10,50))
# Env = raster::crop(Env, ext)
# bg = Env[[1]]
# phainopeplaNitens = databaseToAssignedSubspecies(spp="Phainopepla nitens",
#     subsppList=c("nitens","lepida"),
#     pointLimit=500,
#     dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#    quant=0.95,
#    xmin=-125,
#    xmax=-60,
#    ymin=10,
#    ymax=50,
#    plotIt=T,
#    bgLayer=bg,
#    outputDir="/Users/kprovost/Documents/Classes/Spatial Bioinformatics/project/")
#
# cardinalisSinuatus = databaseToAssignedSubspecies(spp="Cardinalis sinuatus",
# subsppList = c("sinuatus","fulvescens","peninsulae"),
#    pointLimit=500,
#    dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#    quant=0.95,
#    xmin=-125,
#    xmax=-60,
#    ymin=10,
#    ymax=50,
#    plotIt=T,
#    bgLayer=bg,
#    outputDir="/Users/kprovost/Documents/Classes/Spatial Bioinformatics/project/")
