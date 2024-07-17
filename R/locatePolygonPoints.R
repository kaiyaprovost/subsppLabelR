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
#' Point to Polygon Locator
#'
#' This function takes a two polygons and a list of points
#' and determines whether the points fall within the polygons, then
#' adds a column for each polygon's name with a boolean of whether
#' the point overlaps.
#'
#' @param test_points The points to check where they fall
#' @param polygonA The first polygon
#' @param polygonB The second polygon
#' @param setcoord Whether to set a coordinate reference system to the points and polygons
#' @param crs A coordinate reference system to assign to the polygons and the test points if setcoord = TRUE
#' @param nameA The (subspecies) name associated with the first polygon
#' @param nameB The (subspecies) name associated with the second polygon
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polyLocations = labeledLoc
#' polyLocations = locatePolygonPoints(test_points=polyLocations,polygonA=densPol_sin,
#'    polygonB=densPol_ful,crs="+proj=longlat +ellps=WGS84",setcoord = TRUE,nameA="sinuatus",nameB="fulvescens")
locatePolygonPoints = function(test_points,
                               polygonA,
                               polygonB,
                               crs = "+proj=longlat +ellps=WGS84",
                               setcoord = TRUE,
                               nameA = "A",
                               nameB = "B") {
  ## this function determines whether points overlap with polygon A, polygon B, neither, or both
  ## and then adds columns to label them
  #library(dplyr)
  pts = test_points
  pts$longitude = as.numeric(pts$longitude)
  pts$latitude = as.numeric(pts$latitude)
  if (setcoord == T) {
    sp::coordinates(pts) = ~ longitude + latitude
    sp::proj4string(pts) = sp::CRS(crs)
    polygonA = sp::spTransform(polygonA, crs)
    polygonB = sp::spTransform(polygonB, crs)
  }
  inpolygonA = test_points[which(!is.na(sp::over(pts, polygonA))),]
  inpolygonB = test_points[which(!is.na(sp::over(pts, polygonB))),]
  notInpolygonA = test_points[which(is.na(sp::over(pts, polygonA))),]
  notInpolygonB = test_points[which(is.na(sp::over(pts, polygonB))),]
  inBothPolygons = dplyr::intersect(inpolygonA, inpolygonB)
  inNeitherPolygon = dplyr::intersect(notInpolygonA, notInpolygonB)
  onlypolygonA = dplyr::intersect(inpolygonA, notInpolygonB)
  onlypolygonB = dplyr::intersect(inpolygonB, notInpolygonA)
  #print("testpoint colnames")
  #print(colnames(test_points))
  #print("colnames: inboth, inneither, onlyA, onlyB")
  colnames(inBothPolygons) = colnames(test_points)
  colnames(inNeitherPolygon) = colnames(test_points)
  colnames(onlypolygonA) = colnames(test_points)
  colnames(onlypolygonB) = colnames(test_points)
  #print(colnames(inBothPolygons))
  #print(colnames(inNeitherPolygon))
  #print(colnames(onlypolygonA))
  #print(colnames(onlypolygonB))
  #print("names to add")
  #print(paste(nameA,nameB))
  inBothPolygons_1 = cbind(inBothPolygons,
                           testcol1 = rep(1, length(inBothPolygons[, 1])),
                           testcol2 = rep(1, length(inBothPolygons[, 1])))#,rep("both",length(inBothPolygons[,1]))
  #)
  #print(head(inBothPolygons_1))
  inNeitherPolygon_1 = cbind(
    inNeitherPolygon,
    testcol1 = rep(0, length(inNeitherPolygon[, 1])),
    testcol2 = rep(0, length(inNeitherPolygon[, 1]))
  )#,rep("neither",length(inNeitherPolygon[,1]))
  #)
  onlypolygonA_1 = cbind(onlypolygonA,
                         testcol1 = rep(1, length(onlypolygonA[, 1])),
                         testcol2 = rep(0, length(onlypolygonA[, 1])))#,rep(nameA,length(onlypolygonA[,1]))
  #)
  onlypolygonB_1 = cbind(onlypolygonB,
                         testcol1 = rep(0, length(onlypolygonB[, 1])),
                         testcol2 = rep(1, length(onlypolygonB[, 1])))#,rep(nameB,length(onlypolygonB[,1]))
  #)
  #print("AHHHHHHH")
  # length_colnames = length(colnames(inBothPolygons))
  # to_replace_A = length_colnames-1
  # to_replace_B = length_colnames
  #
  # inBothPolygons_1 <- setNames(inBothPolygons_1, c(colnames(inBothPolygons[,1:(length_colnames-2)]),nameA,nameB))
  colnames(inBothPolygons_1)[which(colnames(inBothPolygons_1) == "testcol1")] = paste(nameA) ## WHY IS THIS SETTING TO 2 AND 3 INSTEAD OF 2 AND 1 ETCCCCCC
  colnames(inBothPolygons_1)[which(colnames(inBothPolygons_1) == "testcol2")] = paste(nameB)
  colnames(inNeitherPolygon_1)[which(colnames(inNeitherPolygon_1) == "testcol1")] = paste(nameA)
  colnames(inNeitherPolygon_1)[which(colnames(inNeitherPolygon_1) == "testcol2")] = paste(nameB)
  colnames(onlypolygonA_1)[which(colnames(onlypolygonA_1) == "testcol1")] = paste(nameA)
  colnames(onlypolygonA_1)[which(colnames(onlypolygonA_1) == "testcol2")] = paste(nameB)
  colnames(onlypolygonB_1)[which(colnames(onlypolygonB_1) == "testcol1")] = paste(nameA)
  colnames(onlypolygonB_1)[which(colnames(onlypolygonB_1) == "testcol2")] = paste(nameB)
  #colnames(inBothPolygons_1)[to_replace_A] = nameA
  #colnames(inBothPolygons_1)[to_replace_B] = nameB
  # colnames(inBothPolygons_1) = c(colnames(inBothPolygons),nameA,nameB)#,"assigned_subspecies"
  # #)
  # colnames(inNeitherPolygon_1) = c(colnames(inNeitherPolygon),nameA,nameB)#,"assigned_subspecies"
  # #)
  # colnames(onlypolygonA_1) = c(colnames(onlypolygonA),nameA,nameB)#,"assigned_subspecies"
  # #)
  # colnames(onlypolygonB_1) = c(colnames(onlypolygonB),nameA,nameB)#,"assigned_subspecies"
  # #)
  #print("colnames post adding")
  #print(colnames(inBothPolygons_1))
  #print(colnames(inNeitherPolygon_1))
  #print(colnames(onlypolygonA_1))
  #print(colnames(onlypolygonB_1))
  toReturn = rbind(inBothPolygons_1,
                   inNeitherPolygon_1,
                   onlypolygonA_1,
                   onlypolygonB_1)
  return(toReturn)
}
