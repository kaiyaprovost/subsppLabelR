#' @import raster
#' @import MASS
#' @import spocc
#' @import rgeos
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
  #library(rgeos)
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
}

