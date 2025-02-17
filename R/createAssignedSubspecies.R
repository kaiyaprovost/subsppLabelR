#' @import raster
#' @import MASS
#' @import spocc
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
#' @import AppliedPredictiveModeling
#' @import doFuture
#' @import h2o
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
#' @param method Whether to calculate overlaps with rasters or polygons
#' @param quant quant for density, below which points are removed. E.g., if set to 0.95, removes 95 percent least dense squares.
#' @param xmin Minimum longitude extent to clip rasters to
#' @param xmax Maximum longitude extent to clip rasters to
#' @param ymin Minimum latitute extent to clip rasters to
#' @param ymax Maximum latitude extent to clip rasters to
#' @param plotIt Whether to generate plots
#' @param bgLayer A background layer for generating plots
#' @param outputDir What directory to output to
#' @param datafile if already ran and saved output from spocc:occ, put file here -- default NULL
#' @param epsilon Parameter for anomaly detection
#' @param spp_epsilon Parameter for anomly detection specifically for species
#' @param subspp_epsilon Parameter for anomaly detection specifically for subspecies
#' @param restrictNominate Whether or not to restrict the nominate and remove extra points
#' @param cleanup_nominate Whether or not to clean up the nominate
#' @param num_digits_latlong The number of digits to round the coordinates to
#' @param cells_per_bgLayer How many cells to include per background layer
#' @param downloadOnly Whether to run as download occurrences only
#'
#' @export
#' @examples
#'
#' Env = raster::stack(list.files(path = '~/wc2-5/',pattern = "\\.bil$",full.names = T))
#' ext = raster::extent(c(-125, -60, 10, 50))
#' Env = raster::crop(Env, ext)
#' bg = Env[[1]]
#' phainopeplaNitens = createAssignedSubspecies(spp = "Phainopepla nitens",#' subsppList=c("nitens", "lepida"),#' pointLimit=500,dbToQuery=c("gbif", "bison", "inat", "ebird", "ecoengine", "vertnet"),#' quant=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,#' outputDir="~/project/")
#' suspect_occurrences = phainopeplaNitens$loc_suspect,#' good_occurrences = phainopeplaNitens$loc_good,#' subspecies_polygons = phainopeplaNitens$pol
createAssignedSubspecies = function(spp, subsppList, pointLimit, dbToQuery, method = "raster", quant = 0.95, xmin = -180, xmax = 180, ymin = -90, ymax = 90, plotIt = F, bgLayer = NULL, outputDir, datafile = NULL, epsilon = 1e-6, spp_epsilon = epsilon, subspp_epsilon = epsilon, restrictNominate = F, cleanup_nominate = T, num_digits_latlong = 2, cells_per_bgLayer = 50, downloadOnly = FALSE, ...) {
  setwd(outputDir)
  print(is.null(datafile))
  if (is.null(datafile)){
    print("Downloading species occurrences")
    listFromSubspeciesOcc = subspeciesOccQuery(spp = spp, subsppList = subsppList, pointLimit = pointLimit, dbToQuery = dbToQuery)
    print("Labeling data by subspecies")
    labeledLoc = labelSubspecies(subsppOccList = listFromSubspeciesOcc)
    print("Exporting Occurrence Data")
  } else if (!(is.null(datafile))){
    if (class(datafile) == "character"){
      print("Uploading datafile")
      labeledLoc = read.csv(datafile, sep = "\t", stringsAsFactors = F)
      print("Extracting datafile relevant cols")
      labeledLoc = labeledLoc[, c("name", "longitude", "latitude", "subspecies")]
    } else {labeledLoc = datafile[, c("name", "longitude", "latitude", "subspecies")]}
  }
  nominateSubspecies = strsplit(spp, " ")[[1]][2]
  print("Cleaning bad lat/longs")
  labeledLoc = labeledLoc[!(is.na(labeledLoc$longitude)), ] 
  labeledLoc = labeledLoc[!(is.na(labeledLoc$latitude)), ]
  labeledLoc = labeledLoc[labeledLoc$latitude <= 90, ]
  labeledLoc = labeledLoc[labeledLoc$latitude >= -90, ]
  labeledLoc = labeledLoc[labeledLoc$longitude <= 180, ]
  labeledLoc = labeledLoc[labeledLoc$longitude >= -180, ]
  labeledLoc = labeledLoc[!(is.na(labeledLoc$latitude)), ]
  labeledLoc = labeledLoc[!(is.na(labeledLoc$longitude)), ]
  print(paste("Rounding lat/longs to", num_digits_latlong,"decimal places", sep = " "))
  labeledLoc$latitude = round(labeledLoc$latitude, num_digits_latlong)
  labeledLoc$longitude = round(labeledLoc$longitude, num_digits_latlong)
  labeledLoc = unique(labeledLoc)
  print("Removing points outside of bounds")
  labeledLoc = labeledLoc[labeledLoc$latitude <= ymax, ]
  labeledLoc = labeledLoc[labeledLoc$latitude >= ymin, ]
  labeledLoc = labeledLoc[labeledLoc$longitude <= xmax, ]
  labeledLoc = labeledLoc[labeledLoc$longitude >= xmin, ]
  if (cleanup_nominate == T){
    print("RELABELING NOMINATE AFTER CLEANUP")
    good_nominate_rows = which(grepl(paste(nominateSubspecies, nominateSubspecies, sep = " "), labeledLoc$name))
    labeled_nominate_rows = which(labeledLoc$subspecies == nominateSubspecies)
    nominate_rows_to_keep = intersect(good_nominate_rows, labeled_nominate_rows)
    to_relabel = labeled_nominate_rows[!(labeled_nominate_rows %in% nominate_rows_to_keep)]
    labeledLoc$subspecies[to_relabel] = "unknown"
    labeledLoc = unique(labeledLoc)
  }
  if(downloadOnly==TRUE) {return(list(labeledLoc = labeledLoc, loc_suspect = NULL, loc_good = NULL))}
  
  max_long = max(labeledLoc$longitude, na.rm = T)
  min_long = min(labeledLoc$longitude, na.rm = T)
  max_lat = max(labeledLoc$latitude, na.rm = T)
  min_lat = min(labeledLoc$latitude, na.rm = T)
  if (is.null(bgLayer)) {
    ext = raster::extent(c(min_long, max_long, min_lat, max_lat))
    print(ext)
    bgLayer = raster::raster(ext = ext, nrow = cells_per_bgLayer, ncol = cells_per_bgLayer, vals = 0)
    print(bgLayer)
  }
  subsppNames = unique(labeledLoc$subspecies)
  if (plotIt == T){
    png(paste("Labeled occurences", spp, quant, ".png"))
    raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = spp)
    points(labeledLoc$longitude, labeledLoc$latitude, col = as.factor(labeledLoc$subspecies))
    legend("top", legend = as.factor(unique(labeledLoc$subspecies)), pch = 1, bty = "n", col = as.factor(unique(labeledLoc$subspecies)))
    dev.off()
  }
  print("Starting anomaly detection for whole species")
  list_of_anomalies = subsppLabelR::detectSpatialOutliers(localities = labeledLoc, epsilon = spp_epsilon)
  list_of_anomalies_names = names(list_of_anomalies)
  rows_purged = sort(unique(list_of_anomalies_names))
  print("Starting anomaly detection for each subspecies")
  list_of_anomalies_sub = c()
  for (i in 1:length(c(subsppNames))){
    name = subsppNames[[i]]
    if (name != "unknown") {
      print(name)
      subset = labeledLoc[labeledLoc$subspecies == name, ]
      print(nrow(subset))
      anomalies = subsppLabelR::detectSpatialOutliers(localities = subset, epsilon = subspp_epsilon)
      print(length(anomalies))
      list_of_anomalies_sub = c(list_of_anomalies_sub, anomalies)
    }
  }
  list_of_anomalies_sub_names = names(list_of_anomalies_sub)
  rows_purged_sub = sort(unique(as.integer(list_of_anomalies_sub_names)))
  rows_purged = sort(unique(c(rows_purged, rows_purged_sub)))
  print(paste("Removing", length(rows_purged), "of", length(labeledLoc[, 1]), "detected anomalies"))
  removed = labeledLoc[(rownames(labeledLoc)%in% rows_purged), ]
  labeledLoc = labeledLoc[!(rownames(labeledLoc)%in% rows_purged), ]
  removed$subspecies = as.factor(removed$subspecies)
  levels(removed$subspecies) = levels(as.factor(labeledLoc$subspecies))
  if (plotIt == T){
    png(paste("SubspeciesAnomaliesRemoved_", spp, spp_epsilon, "_", subspp_epsilon, "_", quant, ".png", sep = ""))
    raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = paste("Anomalies"))
    points(labeledLoc$longitude, labeledLoc$latitude, col = "lightgrey", pch = 0)
    points(removed$longitude, removed$latitude, col = as.factor(removed$subspecies))
    legend("top", legend = as.factor(unique(removed$subspecies)), pch = 1, bty = "n", col = as.factor(unique(removed$subspecies)))
    dev.off()
  }
  print("Removing single-individual subspecies")
  for (sub in unique(labeledLoc$subspecies)){
    print(sub)
    rows = (nrow(labeledLoc[labeledLoc$subspecies == sub, ]))
    if (rows <= 1){
      print(sub)
      labeledLoc = labeledLoc[labeledLoc$subspecies != sub, ]
    }
  }
  subsppNames = unique(labeledLoc$subspecies)
  max_long = max(labeledLoc$longitude, na.rm = T)
  min_long = min(labeledLoc$longitude, na.rm = T)
  max_lat = max(labeledLoc$latitude, na.rm = T)
  min_lat = min(labeledLoc$latitude, na.rm = T)
  print("Cleaning bgLayer 2nd time")
  ext2 = raster::extent(c(min_long, max_long, min_lat, max_lat))
  bgLayer = raster::raster(ext = ext2, nrow = cells_per_bgLayer, ncol = cells_per_bgLayer, vals = 0)
  print("Building species kernel density maps")
  xmax = max_long
  xmin = min_long
  ymax = max_lat
  ymin = min_lat
  densityRasters = lapply(subsppNames, function(subspp){
    print(subspp)
    locs = labeledLoc[labeledLoc$subspecies == subspp, ]
    dens = subspeciesDensityMap(localities = locs, quant = quant, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, total_range = bgLayer, subspp = subspp, spp = spp, outputDir = outputDir)
    if (is.null(dens)){dens = NA}
    if ((length((raster::unique(dens, na.last = NA)))) <= 0) {dens = NA}
    names(dens) = subspp
    return(dens)
  })
  names(densityRasters) = subsppNames
  densityRasters = densityRasters[!(is.na(densityRasters))]
  subsppNames = names(densityRasters)
  if (plotIt == T){
    for (i in 1:length(densityRasters)){
      name = names(densityRasters)[[i]]
      png(paste("DensityRaster_", spp, " ", name, quant, ".png", sep = ""))
      raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = paste("Density, subspp:", name))
      raster::plot(densityRasters[[i]], add = T, col = viridis::viridis(99))
      dev.off()
    }
  }
  print("Outputting density rasters")
  densityStack = stack(densityRasters)
  densityStackFile = paste("DensityRaster_", spp, "_", quant, ".tif", sep = "")
  names(densityStack) = names(densityRasters)
  writeRaster(densityStack, densityStackFile, format = "GTiff", overwrite = T, suffix = "names")
  if (method == "raster"){
    print("Removing overlapping raster sections from density map")
    validRasters = densityRasters
    if (length(validRasters)>= 3) {
      for (rasterA_i in 2:length(validRasters)) {
        for (rasterB_i in 2:length(validRasters)) {
          if (rasterA_i > rasterB_i) {
            densA = validRasters[[rasterA_i]]
            densB = validRasters[[rasterB_i]]
            validRasterList = densityRasterRemoveIntersection(densA =densA, densB = densB, verbose = F)
            validRasters[[rasterA_i]] = validRasterList[[1]]
            validRasters[[rasterB_i]] = validRasterList[[2]]
          }
        }
      }
    }
    if (plotIt == T){
      for (i in 1:length(validRasters)){
        name = names(validRasters)[[i]]
        png(paste("ValidRaster_", spp, " ", name, quant, ".png", sep = ""))
        raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = paste("Density, subspp:", name))
        raster::plot(validRasters[[i]], add = T, col = viridis::viridis(99))
        dev.off()
      }
    }
    print("Locating points relative to rasters")
    polyLocations = labeledLoc
    for (slotA in 1:length(subsppNames)) {
      if (subsppNames[[slotA]] != "unknown") {
        polyLocations = subsppLabelR::locateRasterPoints(test_points = polyLocations, rasterA = validRasters[[slotA]], name = subsppNames[[slotA]])
      }
    }
  } 
  if (method == "polygon") {
    print("Converting density maps to polygons")
    densityPolygons = lapply(densityRasters, function(dens){
      print(names(dens))
      densPol = NULL 
      try({densPol = densityMapToPolygons(densityMap = dens)})
      return(densPol)
    })
    if (restrictNominate == T){
      print("Restricting the nominate")
      polygons_notnom = densityPolygons[!(names(densityPolygons) %in% c(nominateSubspecies, "unknown"))]
      fullpoly = raster::bind(polygons_notnom)
      if (length(fullpoly) == 1) {fullpoly = fullpoly[[1]]}
      densityPolygons[[nominateSubspecies]] = sf::st_Difference(densityPolygons[[nominateSubspecies]], fullpoly)
    }
    print(densityPolygons)
    if (plotIt == T) {
      for (i in 1:length(densityPolygons)){
        name = names(densityPolygons)[[i]]
        png(paste("RawDensityPolygon_", spp, " ", name, quant, ".png", sep = ""))
        raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
        plot(densityPolygons[[i]],add = T,col = "red")
        dev.off()
      }
      png(paste("RawDensityPolygon_", spp, quant, " ALL.png", sep = ""))
      raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = "Polygon, ALL")
      cols = c("black", "red", "blue", "green", "cyan", "magenta", "pink", "white", "purple", "orange", "yellow", "sienna", "thistle", "palegreen", "powderblue", "aquamarine", "violet", "mediumslateblue", "lightsalmon", "lightblue")
      for (i in 1:length(densityPolygons)) {
        name = names(densityPolygons)[[i]]
        plot(densityPolygons[[i]],add = T,border = cols[i],lwd =((3 * i)  / 3))
      }
      legend("top", legend = names(densityPolygons), bty = "n", fill = rgb(0, 0, 0, 0), border = cols)
      dev.off()
    }
    print("Checking Overlaps of Polygons and Removing Overlaps")
    densityPolygons_trim1 = polygonTrimmer(polygonList = densityPolygons, namesList = subsppNames)
    if (plotIt == T) {
      for (i in 1:length(densityPolygons_trim1)){
        name = names(densityPolygons_trim1)[[i]]
        png(paste("TrimDensityPolygon_", spp, " ", name, quant, ".png", sep = ""))
        raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = paste("Polygon, subspp:", name))
        raster::plot(densityPolygons_trim1[[i]], add = T, col = viridis::viridis(99))
        dev.off()
      }
      png(paste("TrimDensityPolygon_", spp, quant, " ALL.png", sep = ""))
      raster::plot(bgLayer, col = "grey", colNA = "darkgrey", main = paste("Polygon, subspp:", name))
      cols = c("black", "red", "blue", "green", "cyan", "magenta", "pink", "white", "purple", "orange", "yellow", "sienna", "thistle", "palegreen", "powderblue", "aquamarine", "violet", "mediumslateblue", "lightsalmon", "lightblue")
      for (i in 1:length(densityPolygons_trim1)) {
        print(i)
        name = names(densityPolygons_trim1)[[i]]
        raster::plot(densityPolygons_trim1[[i]], add = T, border = cols[i], lwd = ((3 * i) / 3))
      }
      legend("top", legend = names(densityPolygons_trim1),bty = "n", fill = rgb(0, 0, 0, 0), border = cols)
      dev.off()
    }
    densityPolygons_trim = densityPolygons_trim1
    print("Locating points relative to polygons")
    polyLocations = labeledLoc
    print(densityPolygons_trim)
    print(subsppNames)
    for (slotA in 1:length(subsppNames)) {
      if (subsppNames[[slotA]] != "unknown") {polyLocations = subsppLabelR::locatePolygonPoints(test_points = polyLocations, polygonA = densityPolygons_trim[[slotA]], name = subsppNames[[slotA]])}
    }
  }
  print("Matching subspecies")
  print(head(polyLocations))
  checked = subspeciesMatchChecker(locfile = polyLocations, subsppNames = subsppNames)
  checked_suspect = checked$suspect
  print(head(checked_suspect))
  checked_good = checked$good
  print(head(checked_good))
  print("Warning: no valid definition for subspecies given!")
  ## this is a joke
  return(list(labeledLoc = labeledLoc, loc_suspect = checked_suspect,loc_good = checked_good))
} 
