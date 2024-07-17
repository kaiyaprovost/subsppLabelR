#' @import raster
#' @import MASS
#' @import spocc
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
NULL
#' Species Occurrences to Subspecies Occurrences
#'
#' This function is a wrapper for the package and takes a species with a list of subspecies, and then
#' assigns points to subspecies. Points that have problems with assignment are listed as "suspect"
#' and points without problems are listed as "good." Also returns polygons used to make subspecies
#' assignments.
#'
#' @param spp Genus and species to query, as string
#' @param subsppList Strings of subspecies to query
#' @param pointLimit Maximum point limit to return for each database -- see spocc::occ
#' @param dbToQuery List of databases to search through -- see spocc::occ
#' @param quant quant for density, below which points are removed. E.g., if set to 0.95, removes 95% least dense squares.
#' @param xmin Minimum longitude extent to clip rasters to
#' @param xmax Maximum longitude extent to clip rasters to
#' @param ymin Minimum latitute extent to clip rasters to
#' @param ymax Maximum latitude extent to clip rasters to
#' @param plotIt Whether to generate plots
#' @param bgLayer A background layer for generating plots
#' @param outputDir What directory to output to
#' @param datafile if already ran and saved output from spocc:occ, put file here -- default NULL
#' @param epsilon Parameter for anomaly detection
#'
#' @export
#' @examples
#'
#' Env = raster::stack(list.files(path='~/wc2-5/',pattern="\\.bil$",full.names=T))
#' ext = raster::extent(c(-125,-60,10,50))
#' Env = raster::crop(Env, ext)
#' bg = Env[[1]]
#' phainopeplaNitens = databaseToAssignedSubspecies(spp="Phainopepla nitens",
#'    subsppList=c("nitens","lepida"),
#'    pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#'    quant=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#'    outputDir="~/project/")
#' suspect_occurrences = phainopeplaNitens$loc_suspect,
#' good_occurrences = phainopeplaNitens$loc_good,
#' subspecies_polygons = phainopeplaNitens$pol
databaseToAssignedSubspecies = function(spp,
                                        subsppList,
                                        pointLimit,
                                        dbToQuery,
                                        method="polygon", ## methods: polygon, density
                                        quant = 0.95,
                                        xmin = NULL,
                                        xmax = NULL,
                                        ymin = NULL,
                                        ymax = NULL,
                                        plotIt = F,
                                        bgLayer = NULL,
                                        outputDir,
                                        datafile = NULL,
                                        epsilon = 1e-6,
                                        restrictNominate=F,
                                        cleanup_nominate=T,
                                        num_digits_latlong=2,
                                        cells_per_bgLayer=50,
                                        ...) {
  ## TODO: allow to begin from any step?
  setwd(outputDir)
  ##TODO: give option to supplement these data with data from other sources
  #library(dplyr)
  #library(spocc)
  #library(MASS)
  #library(raster)
  #library(sp)
  ## get the species and subspecies
  #print(is.null(datafile))
  if (is.null(datafile)) {
    print("Downloading species occurrences")
    ## TODO: add progress bar
    listFromSubspeciesOcc = subspeciesOccQuery(spp = spp,subsppList = subsppList,pointLimit = pointLimit,dbToQuery = dbToQuery)
    ## label the data by subspeces
    print("Labeling data by subspecies")
    ## why is this failing here
    labeledLoc = labelSubspecies(subsppOccList = listFromSubspeciesOcc)
    ## TODO: add a check where you make sure that the species is correct
    ## EXPORT THE OCCURRENCE DATA!
    ## NOTE: you don't need to do this. it is identical to merging goodsubspp and suspectsubspp and only taking first four cols
    ## (name	longitude	latitude	subspecies)
    #print("Exporting Occurrence Data")
    #write.table(labeledLoc,paste(paste("OccurrenceDatabase",spp,paste(subspecies,collapse=" "),sep="_"),".occ",sep=""),
    #  quote=FALSE,sep="\t",row.names=FALSE)
  } else if (!(is.null(datafile))) {
    if (class(datafile) == "character") {
      ## check whether given an object or a string
      print("Uploading datafile")
      labeledLoc = read.csv(datafile, sep = "\t", stringsAsFactors = F)
      print("Extracting datafile relevant cols")
      labeledLoc = labeledLoc[, c("name","longitude","latitude","subspecies")]
    } else {
      labeledLoc = datafile[, c("name", "longitude", "latitude", "subspecies")]
    }
  }

  nominateSubspecies = strsplit(spp," ")[[1]][2]

  print("Cleaning bad lat/longs")
  labeledLoc = labeledLoc[!(is.na(labeledLoc$longitude)),] ## fine
  labeledLoc = labeledLoc[!(is.na(labeledLoc$latitude)),]
  labeledLoc = labeledLoc[labeledLoc$latitude<=90,]
  labeledLoc = labeledLoc[labeledLoc$latitude>=-90,]
  labeledLoc = labeledLoc[labeledLoc$longitude<=180,]
  labeledLoc = labeledLoc[labeledLoc$longitude>=-180,]
  labeledLoc = labeledLoc[!(is.na(labeledLoc$latitude)),]
  labeledLoc = labeledLoc[!(is.na(labeledLoc$longitude)),]
  print(paste("Rounding lat/longs to",num_digits_latlong,"decimal places",sep=" "))
  labeledLoc$latitude = round(labeledLoc$latitude,num_digits_latlong)
  labeledLoc$longitude = round(labeledLoc$longitude,num_digits_latlong)
  labeledLoc = unique(labeledLoc)

  if(cleanup_nominate==T){
    print("RELABELING NOMINATE AFTER CLEANUP")
    good_nominate_rows=which(grepl(paste(nominateSubspecies,nominateSubspecies,sep=" "),labeledLoc$name))
    labeled_nominate_rows = which(labeledLoc$subspecies==nominateSubspecies)
    nominate_rows_to_keep = intersect(good_nominate_rows,labeled_nominate_rows)
    to_relabel=labeled_nominate_rows[!(labeled_nominate_rows %in% nominate_rows_to_keep)]
    labeledLoc$subspecies[to_relabel] = "unknown"
    labeledLoc = unique(labeledLoc)
  }

  print("Check xy maxmin")
  if(is.null(xmin)) { xmin = as.numeric(min(as.numeric(labeledLoc$longitude),na.rm=T)) } ## fine?
  if(is.null(xmax)) { xmax = as.numeric(max(as.numeric(labeledLoc$longitude),na.rm=T)) }
  if(is.null(ymin)) { ymin = as.numeric(min(as.numeric(labeledLoc$latitude),na.rm=T)) }
  if(is.null(ymax)) { ymax = as.numeric(max(as.numeric(labeledLoc$latitude),na.rm=T)) }
  print("Cleaning bgLayer")
  if(is.null(bgLayer)){
    print(paste(xmin,xmax,ymin,ymax))
    ext = raster::extent(c(as.numeric(xmin),as.numeric(xmax),
                           as.numeric(ymin),as.numeric(ymax)))
    print(ext)
    bgLayer = raster::raster(ext=ext,nrow=cells_per_bgLayer,ncol=cells_per_bgLayer,vals=0)
    print(bgLayer)
  }
  subsppNames = unique(labeledLoc$subspecies)
  if (plotIt == T) {
    png(paste("Labeled occurences", spp,quant, ".png"))
    #print("Plotting")
    ## TODO: make this work again it doesn't
    raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = spp)
    ## need to make sure not factors and plotting numeric
    points(labeledLoc$longitude,labeledLoc$latitude,col = as.factor(labeledLoc$subspecies))
    legend("top",
           legend = as.factor(unique(labeledLoc$subspecies)),
           pch = 1,bty = "n",col = as.factor(unique(labeledLoc$subspecies)))
    dev.off()
  }
  print("Starting anomaly detection for whole species")
  list_of_anomalies = c()
  detectedLocs = detectSpatialOutliers(localities = labeledLoc, epsilon = epsilon)
  anomalies = as.numeric(detectedLocs)
  list_of_anomalies = c(list_of_anomalies, anomalies)
  rows_purged = sort(unique(as.integer(list_of_anomalies)))
  if (length(rows_purged) > 0) {
    print(paste("Removing",length(rows_purged),"detected anomalies of",length(labeledLoc[, 1]),"rows"))
    removed = labeledLoc[list_of_anomalies,]
    labeledLoc = labeledLoc[-(list_of_anomalies),]
    if (plotIt == T) {
      png(paste("AnomaliesRemoved_", spp,quant,".png", sep = ""))
      raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Anomalies"))
      points(labeledLoc$longitude,labeledLoc$latitude,
             col = "lightgrey",pch = 0)
      points(removed$longitude,removed$latitude,
             col = as.factor(removed$subspecies))
      legend("top",
             legend = as.factor(unique(removed$subspecies)),
             pch = 1,bty = "n",
             col = as.factor(unique(removed$subspecies)))
      dev.off()
    }
  } else {
    print("No anomalies found at whole species level")
  }
  ## removing single individual subspecies
  print("Removing single-individual subspecies")
  for (sub in unique(labeledLoc$subspecies)) {
    #print(sub)
    rows = (nrow(labeledLoc[labeledLoc$subspecies == sub,]))
    if (rows <= 1) {
      print(sub)
      labeledLoc = labeledLoc[labeledLoc$subspecies != sub,]
    }
  }


  print("Starting anomaly detection for each subspecies")

  good_subset = NULL
  bad_subset = NULL
  for(i_name in subsppNames){
    subset = labeledLoc[labeledLoc$subspecies==i_name,]
    anomalies=detectSpatialOutliers(subset)
    if(length(anomalies)>=1){
      good_subset = rbind(good_subset,subset[-anomalies,])
      bad_subset = rbind(bad_subset,subset[anomalies,])
    } else {
      good_subset = rbind(good_subset,subset)
    }
  }
  labeledLoc = good_subset
  removed2 = bad_subset

  if (nrow(removed2) > 0) {
    print(paste("Removing",nrow(removed2),"detected subspecies anomalies"))
    if (plotIt == T) {
      png(paste("SubspeciesAnomaliesRemoved_", spp,quant,".png", sep = ""))
      raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Anomalies"))
      points(labeledLoc$longitude,labeledLoc$latitude,
             col = "lightgrey",pch = 0)
      points(removed2$longitude,removed2$latitude,
             col = as.numeric(as.factor(removed2$subspecies)))
      legend("top",
             legend = as.factor(unique(removed2$subspecies)),
             pch = 1,bty = "n",
             col = as.factor(unique(removed2$subspecies)))
      dev.off()
    }
  } else {
    print("No anomalies found at subspecies level")
  }

  subsppNames = unique(labeledLoc$subspecies)
  ## clean up the bgLayer again in case it needs to be smaller
  print("Check xy maxmin 2nd time")
  xmin2 = as.numeric(min(as.numeric(labeledLoc$longitude),na.rm=T))
  xmax2 = as.numeric(max(as.numeric(labeledLoc$longitude),na.rm=T))
  ymin2 = as.numeric(min(as.numeric(labeledLoc$latitude),na.rm=T))
  ymax2 = as.numeric(max(as.numeric(labeledLoc$latitude),na.rm=T))
  ## make sure the new boundaries aren't outside your given boundaries
  xmin2 = max(xmin,xmin2)
  xmax2 = min(xmax,xmax2)
  ymin2 = max(ymin,ymin2)
  ymax2 = min(ymax,ymax2)
  print("Removing points outside of bounds")
  labeledLoc = labeledLoc[labeledLoc$latitude<=ymax2,]
  labeledLoc = labeledLoc[labeledLoc$latitude>=ymin2,]
  labeledLoc = labeledLoc[labeledLoc$longitude<=xmax2,]
  labeledLoc = labeledLoc[labeledLoc$longitude>=xmin2,]


  print("Cleaning bgLayer 2nd time")
  #print(paste(xmin2,xmax2,ymin2,ymax2))
  ext2 = raster::extent(c(as.numeric(xmin2),as.numeric(xmax2),
                          as.numeric(ymin2),as.numeric(ymax2)))
  #print(ext)
  if(nrow(bgLayer)==cells_per_bgLayer & ncol(bgLayer)==cells_per_bgLayer) {
    bgLayer = raster::raster(ext=ext2,nrow=cells_per_bgLayer,ncol=cells_per_bgLayer,vals=0)
  } else {
    bgLayer = raster::crop(bgLayer,ext=ext2)
  }
  #print(bgLayer)
  ## to reduce error take only subspecies within main density
  ## clean up the polygons so that if grouping way out in middle of nowhere, get rid of it
  ## remove points that fall within the other subspecies' polygon
  ## and account for data being poor
  ## build the density of the points
  ## remove all but 95% (or quant) most dense cells
  ## TODO: some subspecies come back with the entire range as their range due to the way the distribution is
  ## need to generate a raster file to crop the density maps to
  #print("x1")
  total_range = bgLayer
  #print("x2")
  raster::values(total_range)[!is.na(raster::values(bgLayer))] = NA
  #print("x3")
  cells <- raster::cellFromXY(total_range, as.matrix(labeledLoc[,c("longitude","latitude")]));
  #print("x4")
  celltable = table(cells)
  #print("x5")
  total_range[as.numeric(names(celltable))] = celltable
  #print("x6")
  if (plotIt == T) {
    png(paste("FullDistribution", spp,quant,".png", sep = ""))
    raster::plot(total_range,colNA = "darkgrey",main = paste("Distribution"))
    dev.off()
  }

  # if(method %in% c("raw","biasraw")) {
  #   labeledLoc_cells = labeledLoc
  #   cells2 <- raster::cellFromXY(bgLayer, as.matrix(labeledLoc_cells[,c("longitude","latitude")]));
  #   labeledLoc_cells$cells = cells2
  #   #labeledLoc_cells
  #   subspp_cell_table = as.matrix(table(labeledLoc_cells$cells,labeledLoc_cells$subspecies))
  #
  #   x=lapply(1:ncol(subspp_cell_table),FUN=function(i){
  #     mylayer = bgLayer
  #     mylayer[as.numeric(rownames(subspp_cell_table))] = subspp_cell_table[,i]
  #     values(mylayer)[values(mylayer)==0] = NA
  #     return(mylayer)
  #   })
  #   names(x) = colnames(subspp_cell_table)
  #   x_st = stack(x)
  #
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
  #
  #
  # }
  #
  # if(method=="biasraw") {
  #   ## divide the raw
  # }

  # if(method %in% c("polygon","density")) {
  print("Building species kernel density maps")
  ## TODO: add something to increase quants dynamically here? so if 0.95 does not return all the valid subspecies, reduce to 0.90 etc?
  densityRasters = lapply(subsppNames, function(subspp) {
    print(subspp)
    #subspp="relicta"
    locs = labeledLoc[labeledLoc$subspecies == subspp,]
    #print(head(locs))
    dens = subspeciesDensityMap(localities = locs,quant = quant,xmin = xmin,xmax = xmax,
                                ymin = ymin,ymax = ymax, total_range=total_range,subspp=subspp,
                                spp=spp,outputDir=outputDir,raw_raster=T)
    if (is.null(dens)) { dens = NA }
    if ((length((raster::unique(dens,na.last=NA)))) <= 0) { dens = NA }
    names(dens) = subspp
    return(dens)
  })
  #print("done density")
  names(densityRasters) = subsppNames
  ## remove failed ones
  densityRasters = densityRasters[!(is.na(densityRasters))]
  subsppNames = names(densityRasters)
  #print("start plot1")
  if (plotIt == T) {
    for (i in 1:length(densityRasters)) {
      #print("plotforloop")
      name = names(densityRasters)[[i]]
      png(paste("DensityRaster_", spp, " ", name,quant, ".png", sep = ""))
      raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Density, subspp:", name))
      raster::plot(densityRasters[[i]],add = T,col = viridis::viridis(99))
      dev.off()
    }
  }


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
    densityPolygons[[nominateSubspecies]]=sf::st_Difference(densityPolygons[[nominateSubspecies]], fullpoly)
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

  densityPolygons_trim1 = polygonTrimmer(polygonList = densityPolygons, namesList = subsppNames)
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
        polyLocations = locatePolygonPoints(test_points = polyLocations,
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


}

