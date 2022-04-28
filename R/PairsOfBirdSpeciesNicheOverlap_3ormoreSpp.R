#' @import raster
#' @import MASS
#' @import spocc
#' @import rgeos
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import rebird
NULL

#' Pull Subspecies Occurrences
#'
#' This function uses spocc::occ to query GBIF and other
#' repositories for one species and N number of
#' subspecies, specified by the user. Returns a list of
#' occurrence record objects.
#' Note: currently only tested for N=2 and N=3 subspecies.
#'
#' @param spp Genus and species to query, as string
#' @param subsppList Strings of subspecies to query
#' @param pointLimit Maximum point limit to return for each database -- see spocc::occ
#' @param dbToQuery List of databases to search through -- see spocc::occ
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
subspeciesOccQuery = function(spp,
                              subsppList = NULL,
                              pointLimit = 500,
                              dbToQuery = c("gbif", "bison", "inat", "ecoengine", "vertnet"),
                              ...) {
  ## TODO: support for no subspecies given 
  ## this function uses spocc to query for one species and multiple subspecies
  #library(spocc)
  
  print(paste("Getting Species: ", spp))
  
  sppOcc = spocc::occ(query = spp,limit = pointLimit,has_coords = T,from = dbToQuery,...)
  
  if(is.null(subsppList)) {
    print("WARNING: No subspecies given -- not doing subspecies labeling")
    subSppListOcc = NULL
  } else {
    subSppListOcc = lapply(subsppList, function(x) {
      print(paste("     Getting Subspecies: ", x))
      return(spocc::occ(query = paste(spp, x, sep = " "),
                        limit = pointLimit,has_coords = T,from = dbToQuery))
    })
    names(subSppListOcc) = subsppList
    
    print(sppOcc)
    print(subSppListOcc)
  }
  
  toReturn = list(sppOcc, subSppListOcc)
  names(toReturn) = c("unknown", "labeled")
  return(toReturn)
  
}

#' Convert Occurence Records to Labeled Dataframe
#'
#' This function converts a labeled list of subspecies occurences
#' into a dataframe of occurences, with a column for subspecies
#' It is called within labelSubspecies()
#'
#' @param subsppOccList_object A list of subspecies occurences from subspeciesOccQuery() or similar
#' @param subsppOccList_name Names of the subspecies represented in subsppOccList_object
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' sppLocLab = occ2df_subspeciesLabels(subsppOccList_object=listFromSubspeciesOcc[[1]],
#'    subsppOccList_name=names(listFromSubspeciesOcc)[1])
occ2df_subspeciesLabels = function(subsppOccList_object,
                                   subsppOccList_name) {
  ## thus function turns an occ object into a dataframe with a column for subspecies
  ## TODO: make it optional to do the "unique" thing for future processing
  ## TODO: does not work if zero records
  
  sppDf = data.frame(spocc::occ2df(subsppOccList_object))
  
  if (nrow(sppDf) <= 0) {
    sppLocLab = sppDf
    print("THIS SUBSPECIES HAS ZERO RECORDS")
  } else {
    sppLoc = unique(na.omit(sppDf[, 1:3]))
    sppLocLab = sppDf
    sppLocLab$subspecies = subsppOccList_name
    
  }
  
  ## TODO: add a way to make this output the dataframe
  
  #print("return")
  return(sppLocLab)
}

#' Convert Occurence List to Labeled Dataframe
#'
#' This function converts a labeled list of subspecies occurences
#' into a dataframe of occurences, with a column for subspecies, using occ2df_subspeciesLabels().
#' The latter function works only on a single subspecies whereas this function
#' performs the task for all subspecies in the list.
#'
#' @param subsppOccList A list of subspecies occurences from subspeciesOccQuery() or similar
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' subsppNames = unique(labeledLoc$subspecies)
labelSubspecies = function(subsppOccList) {
  ## this function takes a list of three taxa and labels them with subspecies information
  ## TODO: turn this into a function to give multiple subspecies and return it
  ## TODO: doesn't work if one subspp has zero records
  sppOcc = subsppOccList[[1]]
  name_sppOcc = names(subsppOccList)[1]
  #print("Giving occ2df labels")
  sppLocLab = occ2df_subspeciesLabels(subsppOccList_object = sppOcc,
                                      subsppOccList_name = name_sppOcc)
  
  labeledOcc = subsppOccList[[2]]
  if(is.null(labeledOcc)) {
    print("No subspecies present -- returning for species only")
  } else {
    #print(paste("Length labeledOcc:",length(labeledOcc)))
    for (occ in 1:length(labeledOcc)) {
      print(names(labeledOcc)[[occ]])
      subsppLoc = occ2df_subspeciesLabels(subsppOccList_object = labeledOcc[[occ]],subsppOccList_name = names(labeledOcc)[[occ]])
      #print("check1")
      sppLocLab = rbind(sppLocLab, subsppLoc)
      #print("check2")
    }
  }
  return(sppLocLab)
}

#' Generate Subspecies Density Maps
#'
#' This function uses 2D kernel density estimation to make a density raster
#' of points for each subspecies. It then removes all but the most
#' dense cells as definied by quantile.
#'
#' @param localities Labeled localities as generated from labelSubspecies() or subspeciesOccQuery()
#' @param quantile Quantile for density, below which points are removed. E.g., if set to 0.95, removes 95% least dense squares.
#' @param xmin Minimum longitude extent to clip to
#' @param xmax Maximum longitude extent to clip to
#' @param ymin Minimum latitute extent to clip to
#' @param ymax Maximum latitude extent to clip to
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' dens = subspeciesDensityMap(localities=locs,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
subspeciesDensityMap = function(localities,
                                quantile = 0.95,
                                xmin = NULL,
                                xmax = NULL,
                                ymin = NULL,
                                ymax = NULL,
                                bgLayer) {
  ## this function uses kernel density to make a raster that will then be used to filter
  ## the data to remove all but the 5% (or 1-quantile) most dense cells
  ## TODO: allow for subspecies-specific quantiles 
  
  #library(MASS)
  #library(raster)
  if(is.null(xmin)) { xmin = min(localities$longitude,na.rm=T) }
  if(is.null(xmax)) { xmax = max(localities$longitude,na.rm=T) }
  if(is.null(ymin)) { ymin = min(localities$latitude,na.rm=T) }
  if(is.null(ymax)) { ymax = max(localities$latitude,na.rm=T) }
  
  range = c(xmin, xmax, ymin, ymax)
  ext = raster::extent(range)
  
  w1 = matrix(1, 3, 3)
  ## generate the two dimensional kernel density estimation
  
  if (nrow(localities) == 1) {
    print("NOT ENOUGH LOCALITIES")
    return(NULL)
  }
  
  density = NULL
  try({density = MASS::kde2d(as.numeric(localities$longitude),as.numeric(localities$latitude),
                             lims = range,n = 50)})
  
  if(is.null(density)){
    print("NULL")
    return(NULL)
  } else {
  ## convert to raster
  densRas = raster::raster(density)
  ## take the top percentile of the points, only the densest areas
  quan = quantile(densRas[densRas], quantile)
  densRas_trim = densRas
  densRas_trim[densRas_trim <= quan] = NA
  #plot(densRas_trim,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  
  ## crop to the existing layers
  densRas_trim = raster::crop(densRas_trim,bgLayer)
  
  return(densRas_trim)
  }
  
}

#' Convert Density Maps to Polygons
#'
#' This function converts density maps (or other raster files) to a polygon
#'
#' @param densityMap Raster of densities from subspeciesDensityMap() or similar
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
densityMapToPolygons = function(densityMap) {
  ## this function converts density maps to polygons
  ## will work on other kinds of polygons as well
  #library(sp)
  #library(raster)
  polygon = densityMap
  polygon[!is.na(polygon)] = 1
  polygon <- sp::disaggregate(raster::rasterToPolygons(polygon,fun = NULL,na.rm = T,dissolve = T))
  #plot(polygon)
  return(polygon)
}

#' Flag Bad Polygon Overlaps
#'
#' This function takes two polygons representing subspecies ranges for Subspecies-A and Subspecies-B
#' and flags for removal intersections between the polygons. Called within polygonTrimmer()
#'
#' Checks for intersections in two ways:
#'
#' 1) Any polygon features that are completely subsumed within the other
#' subspecies range are flagged for removal, e.g., if Subspecies-A feature is completely subsumed by
#' Subspecies-B, then the Subspecies-A feature is flagged for removal in its entirety.
#' 2) Any intersections between polygons that do not comprise entire polygon features are flagged
#' for removal such that the intersection is flagged for removal from the larger polygon, but
#' not the smaller polygon. If polygons are the same size, intersection is removed from both.
#' E.g., Subspecies-A has a larger range than Subspecies-B. An intersection between them is removed
#' only from Subspecies-A.
#'
#' @param subsppPoly1 Polygon for first subspecies to compare
#' @param subsppPoly2 Polygon for second subspecies to compare
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polygonsToRemove = (flagPolygonOverlap(subsppPoly1=densPol_sin,
#'    subsppPoly2=densPol_ful))
#' sin_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
#' ful_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
#' overlapToRemove_sin = polygonsToRemove$subsppA_intToRemove
#' overlapToRemove_ful = polygonsToRemove$subsppB_intToRemove
flagPolygonOverlap = function(subsppPoly1 = polA,
                              subsppPoly2 = polB) {
  ##function(subsppPoly1,subsppPoly2)
  ## this function checks for overlaps between polygons
  ## TODO: remove polygon if not touching another of same spp
  ## but also closer to polygon of other species
  #library(rgeos)
  #library(raster)
  
  ## there is a bug -- if one subspp range is entirely subsumed within another polygon,
  ## will delete that subspecies. no bueno
  
  badList_subsppA_features = c()
  badList_subsppB_features = c()
  overlapsToRemove_subsppA = c()
  overlapsToRemove_subsppB = c()
  
  for (feature_subsppA in (1:length(subsppPoly1))) {
    ## get the features within the subspecies polygon
    for (feature_subsppB in (1:length(subsppPoly2))) {
      ## get the features within the subspecies polygon
      ## check areas
      totArea1 = rgeos::gArea(subsppPoly1) ## the whole area of the subspecies
      totArea2 = rgeos::gArea(subsppPoly2) ## the whole area of the subspecies
      area1 = rgeos::gArea(subsppPoly1[feature_subsppA,]) ## the area of the single feature
      area2 = rgeos::gArea(subsppPoly2[feature_subsppB,]) ## the area of the single feature
      
      # subsppPoly1$totalArea = totArea1
      # subsppPoly2$totalArea = totArea2
      # subsppPoly1[feature_subsppA,]$featureArea = area1
      # subsppPoly2[feature_subsppB]$featureArea = area2
      
      if (rgeos::gIntersects(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])) {
        #print("INTERSECTS")
        testInt = rgeos::gIntersection(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
        #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
        #print(length(intersect))
        
        
        ## if they overlap
        
        #plot(subsppPoly1[feature_subsppA,],border="red",add=T)
        #plot(subsppPoly2[feature_subsppB,],border="cyan",add=T)
        testArea = rgeos::gArea(testInt)
        
        if (testArea > 0) {
          #print("OVERLAP AREA NOT ZERO")
          intersect = raster::intersect(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
          #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
          areaInt = rgeos::gArea(intersect)
          
          ## if they overlap
          
          area1percent = areaInt / totArea1 ## the area of the feature as a percent of the area of the subspecies
          area2percent = areaInt / totArea2
          
          if (areaInt >= area1 || areaInt >= area2) {
            #print("SUBSUMED")
            
            if (areaInt >= area1) {
              ## if the overlap entirely subsumes an area
              
              
              ## check if the area is the entire subspecies range
              if (area1percent != 1) {
                ## remove the area
                #print("SUBSPP1 SUBSUMED")
                badList_subsppA_features = c(badList_subsppA_features, feature_subsppA)
              }
              
              ## TODO: change to check for density?
              
            }
            if (areaInt >= area2) {
              ## if the overlap entirely subsumes an area
              
              ## check if the area is the entire subspecies range
              if (area2percent != 1) {
                ## remove the area
                ## TODO: change to check for density?
                #print("SUBSPP2 SUBSUMED")
                badList_subsppB_features = c(badList_subsppB_features, feature_subsppB)
              }
              
            }
          }
          else {
            #print("NOT SUBSUMED")
            if (area1percent <= area2percent) {
              #print("AREA1 IS LARGER")
              #print("remove area1")
              overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
            }
            else if (area1percent >= area2percent) {
              #print("AREA2 IS LARGER")
              #print("remove area1")
              overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
            }
          }
        }
      }
    }
  }
  toReturn = list(
    subsppApoly_toremove = badList_subsppA_features,
    subsppBpoly_toremove = badList_subsppB_features,
    subsppA_intToRemove = overlapsToRemove_subsppA,
    subsppB_intToRemove = overlapsToRemove_subsppB
  )
  return(toReturn)
}

#' Flag Bad Polygon Overlaps
#'
#' This function takes two polygons representing subspecies ranges for Subspecies-A and Subspecies-B
#' and flags for removal intersections between the polygons. Called within polygonTrimmer()
#'
#' Checks for intersections in two ways:
#'
#' 1) Any polygon features that are completely subsumed within the other
#' subspecies range are flagged for removal, e.g., if Subspecies-A feature is completely subsumed by
#' Subspecies-B, then the Subspecies-A feature is flagged for removal in its entirety.
#' 2) Any intersections between polygons that do not comprise entire polygon features are flagged
#' for removal such that the intersection is flagged for removal from the larger polygon, but
#' not the smaller polygon. If polygons are the same size, intersection is removed from both.
#' E.g., Subspecies-A has a larger range than Subspecies-B. An intersection between them is removed
#' only from Subspecies-A.
#'
#' @param subsppPoly1 Polygon for first subspecies to compare
#' @param subsppPoly2 Polygon for second subspecies to compare
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polygonsToRemove = (flagPolygonOverlap(subsppPoly1=densPol_sin,
#'    subsppPoly2=densPol_ful))
#' sin_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
#' ful_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
#' overlapToRemove_sin = polygonsToRemove$subsppA_intToRemove
#' overlapToRemove_ful = polygonsToRemove$subsppB_intToRemove
flagPolygonOverlap2 = function(subsppPoly1 = polA,
                               subsppPoly2 = polB,
                               crs = "+proj=longlat +ellps=WGS84") {
  ##function(subsppPoly1,subsppPoly2)
  ## this function checks for overlaps between polygons
  ## TODO: remove polygon if not touching another of same spp
  ## but also closer to polygon of other species
  #library(rgeos)
  #library(raster)
  
  ## there is a bug -- if one subspp range is entirely subsumed within another polygon,
  ## will delete that subspecies. no bueno
  
  badList_subsppA_features = c()
  badList_subsppB_features = c()
  overlapsToRemove_subsppA = c()
  overlapsToRemove_subsppB = c()
  
  ## NOT RIGHT BECAUSE NOT ADDING TO CORRECT SPOTS
  
  for (feature_subsppA in (1:length(subsppPoly1))) {
    ## get the features within the subspecies polygon
    for (feature_subsppB in (1:length(subsppPoly2))) {
      ## get the features within the subspecies polygon
      ## check areas
      
      ## Warning message:
      ## In RGEOSMiscFunc(spgeom, byid, "rgeos_area") :
      ##  Spatial object is not projected; GEOS expects planar coordinates
      
      
      
      totArea1 = rgeos::gArea(subsppPoly1) ## the whole area of the subspecies
      totArea2 = rgeos::gArea(subsppPoly2) ## the whole area of the subspecies
      area1 = rgeos::gArea(subsppPoly1[feature_subsppA,]) ## the area of the single feature
      area2 = rgeos::gArea(subsppPoly2[feature_subsppB,]) ## the area of the single feature
      
      # subsppPoly1$totalArea = totArea1
      # subsppPoly2$totalArea = totArea2
      # subsppPoly1[feature_subsppA,]$featureArea = area1
      # subsppPoly2[feature_subsppB]$featureArea = area2
      
      if (rgeos::gIntersects(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])) {
        #print("INTERSECTS")
        testInt = rgeos::gIntersection(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
        #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
        #print(length(intersect))
        
        
        ## if they overlap
        
        #plot(subsppPoly1[feature_subsppA,],border="red",add=T)
        #plot(subsppPoly2[feature_subsppB,],border="cyan",add=T)
        testArea = rgeos::gArea(testInt)
        
        if (testArea > 0) {
          #print("OVERLAP AREA NOT ZERO")
          intersect = raster::intersect(subsppPoly1[feature_subsppA,], subsppPoly2[feature_subsppB,])
          #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
          areaInt = rgeos::gArea(intersect)
          
          ## if they overlap
          
          area1percent = areaInt / totArea1 ## the area of the feature as a percent of the area of the subspecies
          area2percent = areaInt / totArea2
          
          ## NEW IF STATEMENTS
          if (area1 > area2) {
            ## if A big B small -- area1 vs area2
            if (testArea < area2) {
              ## if B not subsumed, remove from A -- testArea
              #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
              
            }
            else if (testArea >= area2) {
              ## else if B subsumed
              if (area2percent == 1) {
                ## if B = 100% of subspecies, remove from A -- area2percent
                #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
                overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
                
              }
              else if (area2percent != 1) {
                ## else if B not 100% of subspecies, remove from B
                badList_subsppB_features = c(badList_subsppB_features, feature_subsppB)
              }
            }
          }
          else if (area1 < area2) {
            ## else if A small B big
            if (testArea < area1) {
              ## if A not subsumed, remove from B -- area1percent
              #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
              
            }
            else if (testArea >= area1) {
              ## else if A subsumed
              if (area1percent == 1) {
                ## if A = 100% of subspecies, remove from B
                #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
                overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
                
              }
              else if (area1percent != 1) {
                ## else if A not 100% of subspecies, remove from A
                badList_subsppA_features = c(badList_subsppA_features, feature_subsppA)
              }
            }
          }
          else if (area1 == area2) {
            ## else if A = B
            ## check total areas and remove larger
            if (totArea1 > totArea2) {
              ## if total A is greater than total B, remove A
              #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
              
            }
            else if (totArea1 < totArea2) {
              ## if total B is greater than total A, remove B
              #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
              
            }
            else if (totArea1 == totArea2) {
              ## if they are the exact same size
              if (testArea < area1 &&
                  testArea < area2) {
                ## if neither subsumed, remove arbitrary
                flip = sample(1:2, 1)
                if (flip == 1) {
                  ## flip coin remove A
                  #badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
                  overlapsToRemove_subsppA = c(overlapsToRemove_subsppA, intersect)
                  
                }
                else if (flip == 2) {
                  ## flip coin remove B
                  #badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
                  overlapsToRemove_subsppB = c(overlapsToRemove_subsppB, intersect)
                  
                }
              }
              else if (testArea == area1 &&
                       testArea == area2) {
                ## if both are subsumed, keep both
                print("EXACT MATCH, KEEPING BOTH ")
              }
            }
          }
        }
      }
    }
    
    
    
    
    toReturn = list(
      subsppApoly_toremove = badList_subsppA_features,
      subsppBpoly_toremove = badList_subsppB_features,
      subsppA_intToRemove = overlapsToRemove_subsppA,
      subsppB_intToRemove = overlapsToRemove_subsppB
    )
    return(toReturn)
  }
}

## END







#' Check for Subspecies Contiguity
#'
#' This function takes multiple subspecies ranges as polygons and checks whether
#' 1) polygon features are adjacent to their own subspecies, and if not, 2)
#' whether polygon features are closer to their own subspecies than to another.
#' If a subspecies feature does not meet either criteria, it is eliminated.
#'
#' WARNING!!! Function is not working. Do not use.
#'
#' @param listOfPolygons List of subspecies range polygons.
#'
#' @export
#' @examples
#'
closestPolygonFunction = function(listOfPolygons) {
  ## TODO: is broken
  ##TODO: figure out how to do so that you check if polygons are between other polygons
  ## for now just remove if touching non-self boundary
  ## check each feature of each polygon
  newListPoly = listOfPolygons
  polyNames = names(newListPoly)
  
  for (subspp_i in 1:length(newListPoly)) {
    if (polyNames[subspp_i] != "unknown") {
      print(polyNames[subspp_i])
      wholePol = newListPoly[[subspp_i]]
      wholePol = sp::disaggregate(wholePol)
      
      if (length(wholePol) <= 1) {
        print("only one feature, skipping")
      }
      else {
        wholePol$area = rgeos::gArea(wholePol, byid = T)
        print(wholePol$area)
        
        ## sort features by size
        wholePol = wholePol[rev(order(wholePol$area)),]
        #print(wholePol)
        
        for (features_j in 1:length(wholePol)) {
          feat2check = wholePol[features_j,]
          if (rgeos::gArea(feat2check) < max(wholePol$area)) {
            ## if not the largest chunk, check if touching other chunks
            int = (rgeos::gIntersection(feat2check, wholePol[-features_j,]))
            if (is.null(int)) {
              ## diagonal is fine
              ## if not touching, check what polygon is shortest distance to
              polyListForCheck = newListPoly
              polyListForCheck[[subspp_i]] = wholePol[-features_j,]
              #plot(newListPoly[[subspp_i]])
              #plot(polyListForCheck[[subspp_i]])
              #polyListForCheck = polyListForCheck[names(polyListForCheck)!="unknown"]
              dists = sapply(1:length(polyListForCheck), function(i) {
                x = ((
                  rgeos::gDistance(feat2check, polyListForCheck[[i]])
                ))
                names(x) = names(polyListForCheck)[[i]]
                return(x)
              })
              dists[which(names(dists) == "unknown")] = 1e9
              print(paste(
                "Testing ",
                names(polyListForCheck)[subspp_i],
                " against ",
                paste(names(which(
                  dists == min(dists)
                )), collapse = ", "),
                sep = ""
              ))
              if (!(polyNames[subspp_i] %in% names(which(dists == min(
                dists
              ))))) {
                print("Closest is not same subspecies")
                if (min(dists) == 0) {
                  print("Touching other subspecies, removing")
                  removed = wholePol[-features_j,]
                  newListPoly[[subspp_i]] = removed
                }
              }
            }
          }
        }
      }
    }
  }
  return(newListPoly)
}

#' Remove Overlapping Polygons
#'
#' This function takes polygon features flagged for removal from flagPolygonOverlap()
#' and removes them from the polygon. Called within polygonTrimmer().
#'
#' @param polygon Polygon to remove flagged features from.
#' @param idList List of flagged features to completely remove.
#' @param intList List of intersections to remove.
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polygonsToRemove = (flagPolygonOverlap(subsppPoly1=densPol_sin,
#'    subsppPoly2=densPol_ful))
#' sin_polygonsToRemove = polygonsToRemove$subsppApoly_toremove
#' ful_polygonsToRemove = polygonsToRemove$subsppBpoly_toremove
#' overlapToRemove_sin = polygonsToRemove$subsppA_intToRemove
#' overlapToRemove_ful = polygonsToRemove$subsppB_intToRemove
#' sin_densityPolygon_trim = trimPolygonsByOverlap(polygon=densPol_sin,
#'    idList = sin_polygonsToRemove,intList=overlapToRemove_sin)
#' ful_densityPolygon_trim = trimPolygonsByOverlap(polygon=densPol_ful,
#'    idList = ful_polygonsToRemove,intList=overlapToRemove_ful)
trimPolygonsByOverlap = function(polygon,
                                 idList = NULL,
                                 intList = NULL) {
  ## function(polygon,idList=NULL,intList=NULL)
  ## it then removes polygons that are completely within another polygon as ID'd by flagPolygonOverlap()
  ## TODO: what about things that are in neither polygon?
  ## TODO: what about things that are in both?
  
  polytrim = polygon
  if (is.null(idList)) {
    if (is.null(intList)) {
      #print("no overlaps to remove")
      polytrim = rgeos::gUnion(polytrim, polytrim)
    }
    else {
      for (int in 1:length(intList)) {
        polytrim = rgeos::gDifference(polytrim, intList[[int]])
      }
      #print("removing differences from larger shapefile")
    }
  }
  else {
    polytrim = polygon[-idList,]
    if (is.null(intList)) {
      #print("removing subsumed polygons")
      polytrim = rgeos::gUnion(polytrim, polytrim)
    }
    else {
      for (int in 1:length(intList)) {
        polytrim = rgeos::gDifference(polytrim, intList[[int]])
      }
      #print("removing subsumed polygons and differences from larger shape")
    }
    
  }
  return(polytrim)
}

#' Polygon Trimmer
#'
#' This function takes a list of polygons as from densityMapToPolygons() or similar
#' and removes portions of the polygons as determined by flagPolygonOverlap() and
#' trimPolygonsByOverlap()
#'
#' @param polygonList A list of polygons to check
#' @param namesList The subspecies (or other) names associated with polygon list
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#'
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' densityPolygons = list(sinuatus=densPol_sin,fulvescens=densPol_ful)
#' densityPolygons_trim = polygonTrimmer(polygonList=densityPolygons,
#'    namesList=c("sinuatus","fulvescens"))
#'
polygonTrimmer = function(polygonList, namesList) {
  newPolygonList = polygonList
  for (slotA in 1:length(namesList)) {
    for (slotB in 1:length(namesList)) {
      if (namesList[[slotA]] != "unknown" &&
          namesList[[slotB]] != "unknown" && slotA != slotB) {
        #print(paste(slotA,slotB,sep=" "))
        #print(paste(namesList[[slotA]],"with",namesList[[slotB]],sep=" "))
        polA = newPolygonList[[slotA]]
        polB = newPolygonList[[slotB]]
        
        #print(class(polA))
        #print(class(polB))
        
        # plot(bg,col="grey",colNA="darkgrey")
        # plot(polA,add=T,border="cyan",lwd=7)
        # plot(polB,add=T,border="red",lwd=4)
        # #invisible(readline(prompt="Press [enter] to continue"))
        # if(!is.null(raster::intersect(polA,polB))){
        #   plot(raster::intersect(polA,polB),add=T,lwd=1,border="black")
        # }
        
        
        ## this throws the warnings
        polygonsToRemove = (flagPolygonOverlap2(polA, polB)) ######### CHANGED
        #print("CALLING NEW FUNCTION")
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
        
        ## TODO: fix this
        subsppA = namesList[[slotA]]
        subsppB = namesList[[slotB]]
        if(!(is.null(subsppA_densityPolygon_trim))){newPolygonList[[slotA]] = subsppA_densityPolygon_trim}
        if(!(is.null(subsppB_densityPolygon_trim))){newPolygonList[[slotB]] = subsppB_densityPolygon_trim}
        names(newPolygonList) = names(polygonList)
        #print(names(newPolygonList))
        
        #plot(bg,col="grey",colNA="darkgrey")
        #plot(subsppA_densityPolygon_trim,add=T,border="cyan",lwd=7)
        #plot(subsppB_densityPolygon_trim,add=T,border="red",lwd=4)
        #plot(intersect(polA,polB),add=T,lwd=1,border="black")
        
      }
    }
  }
  return(newPolygonList)
}

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
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,
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

#' Check Subspecies Matches and Labels
#'
#' This function takes the locations of points relative to polygons from locatePolygonPoints() or similar
#' and checks where points are assigned. Points that are in 1) no polygons, 2) multiple polygons, or
#' 3) have a subspecies label that does not match their assignment are flagged as "suspect".
#' Any points that are in single polygons and either 1) were originally unlabeled or 2) have a subspecies
#' label that matches their assignment are flagged as "good".
#'
#' @param locfile Points with subspecies assignments and boolean locations with respect to polygons.
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polyLocations = labeledLoc
#' polyLocations = locatePolygonPoints(test_points=polyLocations,polygonA=densPol_sin,polygonB=densPol_ful,crs="+proj=longlat +ellps=WGS84",setcoord = TRUE,nameA="sinuatus",nameB="fulvescens")
#' checked = subspeciesMatchChecker(locfile = polyLocations)
#' checked_suspect = checked$suspect
#' checked_good = checked$good
subspeciesMatchChecker = function(locfile, subsppNames) {
  #print("a")
  locWithSubspecies = locfile
  #print("b")
  #subsppNames = names(locWithSubspecies[5:length(names(locWithSubspecies))])
  #print("c")
  print(subsppNames)
  if("unknown" %in% subsppNames){
    numSub = length(subsppNames)-1
  } else {
    numSub = length(subsppNames)
  }
  
  if(numSub>1){
  
  #print("d")
  numPoints = length(rownames(locWithSubspecies))
  #print("e")
  lastSubsppCol = length(colnames(locWithSubspecies))
  #print("f")
  subsppAssignCol = locWithSubspecies[, (1+lastSubsppCol-numSub):lastSubsppCol]
  #print(head(locWithSubspecies))
  # for (colnum in 5:length(colnames(locWithSubspecies))){
  #   print(paste("colnum",colnum))
  #   print(head(locWithSubspecies[,colnum]))
  #   num_for_name = as.integer(colnames(locWithSubspecies)[colnum])
  #   print(class(num_for_name))
  #   print(paste("numforname",num_for_name))
  #   print(subsppNames)
  #   name_to_replace = subsppNames[num_for_name]
  #   print(paste("nametoreplace",name_to_replace))
  #   names(locWithSubspecies)[colnum] = name_to_replace
  #   print(names(locWithSubspecies))
  # }
  
  #print("g")
  subsppPriorCol = locWithSubspecies$subspecies
  #print("h")
  locWithSubspecies$numSubsppGroups <-
    rowSums(subsppAssignCol, na.rm = TRUE)
  #print("i")
  locWithSubspecies$assigned = NA
  #print("j")
  
  notassigned = locWithSubspecies[which(locWithSubspecies$numSubsppGroups ==
                                          0),]
  #print("k")
  if (nrow(notassigned) != 0) {
    notassigned$assigned = "none"
  }
  #print("l")
  
  multigroup = locWithSubspecies[which(locWithSubspecies$numSubsppGroups >
                                         1),]
  #print("m")
  if (nrow(multigroup) != 0) {
    multigroup$assigned = "multiple"
  }
  #print("n")
  
  singlegroup = locWithSubspecies[which(locWithSubspecies$numSubsppGroup ==
                                          1),]
  #print("o")
  
  ## check whether mismatch between apriori and not
  suspectpoints = rbind(multigroup, notassigned)
  goodpoints = data.frame()
  #print("p")
  
  for (column in 5:lastSubsppCol) {
    #print(column)
    name = colnames(singlegroup)[column]
    #print(name)
    assignedThat = singlegroup[singlegroup[, column] == 1,]
    
    ## breaking if none assigned
    
    if (nrow(assignedThat) == 0) {
      print("NONE ASSIGNED:")
      print(name)
    } else {
      ## subspp should be unk or the subspp
      assignedThat$assigned = name
      
    }
    okaynames = c("unknown", name)
    
    wrong1 = assignedThat[which(!(assignedThat$subspecies %in% okaynames)),]
    right1 = assignedThat[which(assignedThat$subspecies %in% okaynames),]
    
    
    
    notassignedThat = singlegroup[singlegroup[, column] == 0,]
    
    wrong2 = notassignedThat[which(notassignedThat$subspecies == name),]
    right2 = notassignedThat[which(notassignedThat$subspecies != name),]
    
    ## DO NOT ADD RIGHT/WRONG 2, NOT ENOUGH INFORMATION YET
    suspectpoints = rbind(suspectpoints, wrong1)
    goodpoints = rbind(goodpoints, right1)
  }
  #print("q")
  suspectpoints = unique(suspectpoints)
  goodpoints = unique(goodpoints)
  #print("r")
  } else {
    suspectpoints = data.frame()
    goodpoints = locWithSubspecies
  }
  
  return(list(suspect = suspectpoints,
              good = goodpoints))
}




#' Outlier Detection
#'
#' This does outlier detection on points.
#'
#' @param localities A list of localities to check for anomalies
#' @param epsilon An value with which to flag anomalies with probability less than epsilon
#'
#' @export
#' @examples
#'
#'
detectSpatialOutliers = function(localities = locs,
                                 epsilon = 0.0001) {
  ## use MASS to do linear algebra
  ## TODO: why is this removing central anomalies? 
  
  m = length(localities[, 1])
  
  if (m == 1) {
    anomalies = 0
    purged = c()
    kept = localities
    
    return(list(purged, anomalies, kept))
  }
  
  lat = as.numeric(localities$latitude)
  lon = as.numeric(localities$longitude)
  space = cbind(lat, lon)
  n = length(space[1,])
  mu = colMeans(space)
  Sigma = cov(space) ## this does the above but faster
  
  ## SOURCE: https://datascience-enthusiast.com/R/anomaly_detection_R.html
  ## And also the Corusera course Machine Learning by A. Ng
  ## this uses a somewhat different probability density function
  centered <- caret::preProcess(space, method = "center")
  space_mu <- as.matrix(predict(centered, space))
  
  #sigma2=diag(var(X2))
  #sigma2
  #sigma2=diag(sigma2)
  #sigma2 ## this is the same as Sigma on the diagonals, but different on the off-diags
  ## it doesn't matter which you use here you get the same answer
  
  A = (2 * pi) ^ (-n / 2) * det(Sigma) ^ (-0.5)
  B = exp(-0.5 * rowSums((space_mu %*% MASS::ginv(Sigma)) * space_mu))
  p_x = A * B
  
  #P_x = (1 /( ((2*pi)^(n/2) * (norm(Sigma) ^ (1/2))))) * exp (-1/2*(space-mu) %*% MASS::ginv(Sigma) %*% t(space-mu))
  
  p_x = cbind(p_x)
  colnames(p_x) = "anomaly"
  
  locs_2 = cbind(localities, p_x)
  
  anomalies = which(locs_2$anomaly < epsilon)
  purged = localities[anomalies,]
  kept = localities[-(anomalies),]
  
  return(list(purged, anomalies, kept))
  
}



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
#' @param quantile Quantile for density, below which points are removed. E.g., if set to 0.95, removes 95% least dense squares.
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
#'    quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#'    outputDir="~/project/")
#' suspect_occurrences = phainopeplaNitens$loc_suspect,
#' good_occurrences = phainopeplaNitens$loc_good,
#' subspecies_polygons = phainopeplaNitens$pol
databaseToAssignedSubspecies = function(spp,
                                        subsppList,
                                        pointLimit,
                                        dbToQuery,
                                        quantile = 0.95,
                                        xmin = NULL,
                                        xmax = NULL,
                                        ymin = NULL,
                                        ymax = NULL,
                                        plotIt = F,
                                        bgLayer = NULL,
                                        outputDir,
                                        datafile = NULL,
                                        epsilon = 1e-6,
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
  
  print("Cleaning bad lat/longs")
  labeledLoc = labeledLoc[!(is.na(labeledLoc$longitude)),] ## fine
  labeledLoc = labeledLoc[!(is.na(labeledLoc$latitude)),]
  labeledLoc = labeledLoc[labeledLoc$latitude<=90,]
  labeledLoc = labeledLoc[labeledLoc$latitude>=-90,]
  labeledLoc = labeledLoc[labeledLoc$longitude<=180,]
  labeledLoc = labeledLoc[labeledLoc$longitude>=-180,]  
  
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
    bgLayer = raster::raster(ext=ext,nrow=100,ncol=100,vals=0)
    print(bgLayer)
  }
  
  subsppNames = unique(labeledLoc$subspecies)
  
  if (plotIt == T) {
    png(paste("Labeled occurences", spp, ".png"))
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
  
  purged_list = c()
  kept_list = c()
  list_of_anomalies = c()
  
  for (i in 0:length(c(subsppNames))) {
    if (i == 0) {
      #print("full")
      #print(nrow(labeledLoc))
      detectedLocs = detectSpatialOutliers(localities = labeledLoc, epsilon = epsilon)
    }
    else {
      name = subsppNames[[i]]
      if (name != "unknown") {
        #print(name)
        subset = labeledLoc[labeledLoc$subspecies == name,]
        #print(nrow(subset))
        detectedLocs = detectSpatialOutliers(localities = subset, epsilon = epsilon)
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
  
  if (length(rows_purged) > 0) {
    print(paste("Removing",length(rows_purged),"of",length(labeledLoc[, 1]),"detected anomalies"))
    removed = labeledLoc[(rows_purged),]
    labeledLoc = labeledLoc[-(rows_purged),]
    
    if (plotIt == T) {
      png(paste("AnomaliesRemoved_", spp, ".png", sep = ""))
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
    print("No anomalies found")
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
  subsppNames = unique(labeledLoc$subspecies)
  
  ## clean up the bgLayer again in case it needs to be smaller 
  print("Check xy maxmin 2nd time")
  xmin2 = as.numeric(min(as.numeric(labeledLoc$longitude),na.rm=T))
  xmax2 = as.numeric(max(as.numeric(labeledLoc$longitude),na.rm=T))
  ymin2 = as.numeric(min(as.numeric(labeledLoc$latitude),na.rm=T))
  ymax2 = as.numeric(max(as.numeric(labeledLoc$latitude),na.rm=T))
  
  print("Cleaning bgLayer 2nd time")
  print(paste(xmin2,xmax2,ymin2,ymax2))
  ext2 = raster::extent(c(as.numeric(xmin2),as.numeric(xmax2),
                          as.numeric(ymin2),as.numeric(ymax2)))
  print(ext)
  if(nrow(bgLayer)==100 & ncol(bgLayer)==100) {
    bgLayer = raster::raster(ext=ext2,nrow=100,ncol=100,vals=0)
  } else {
    bgLayer = raster::crop(bgLayer,ext=ext2)
  }
  
  print(bgLayer)
  
  ## to reduce error take only subspecies within main density
  ## clean up the polygons so that if grouping way out in middle of nowhere, get rid of it
  ## remove points that fall within the other subspecies' polygon
  ## and account for data being poor
  ## build the density of the points
  ## remove all but 95% (or quantile) most dense cells
  
  print("Building species kernel density maps")
  
  ## TODO: add something to increase quantiles dynamically here? so if 0.95 does not return all the valid subspecies, reduce to 0.90 etc? 
  densityRasters = lapply(subsppNames, function(subspp) {
    print(subspp)
    #subspp="relicta"
    locs = labeledLoc[labeledLoc$subspecies == subspp,]
    #print(head(locs))
    dens = subspeciesDensityMap(localities = locs,quantile = quantile,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax, bgLayer=bgLayer)
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
      png(paste("DensityRaster_", spp, " ", name, ".png", sep = ""))
      plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Density, subspp:", name))
      plot(densityRasters[[i]],add = T,col = viridis::viridis(99))
      dev.off()
    }
  }
  #print("endplot1")
  
  ## convert to polygons
  print("Converting density maps to polygons")
  ## can't handle if there's no data in previous step
  
  densityPolygons = lapply(densityRasters, function(dens) {
    #print(names(dens))
    densPol = densityMapToPolygons(densityMap = dens)
    return(densPol)
  })
  
  #print(densityPolygons)
  
  if (plotIt == T) {
    for (i in 1:length(densityPolygons)) {
      name = names(densityPolygons)[[i]]
      png(paste("RawDensityPolygon_", spp, " ", name, ".png", sep = ""))
      raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
      sp::spplot(densityPolygons[[i]],add = T,col = viridis::viridis(99))
      dev.off()
    }
    
    png(paste("RawDensityPolygon_", spp, " ALL.png", sep = ""))
    raster::plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
    cols = c("black","red", "blue", "green",
             "cyan", "magenta", "pink", "white",
             "purple", "orange", "yellow", "sienna",
             "thistle", "palegreen", "powderblue", "aquamarine",
             "violet", "mediumslateblue", "lightsalmon", "lightblue")
    for (i in 1:length(densityPolygons)) {
      name = names(densityPolygons)[[i]]
      sp::spplot(densityPolygons[[i]],add = T,border = cols[i],lwd = ((3 * i) / 3))
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
  
  densityPolygons_trim1 = polygonTrimmer(polygonList = densityPolygons, namesList = subsppNames)
  
  if (plotIt == T) {
    for (i in 1:length(densityPolygons_trim1)) {
      name = names(densityPolygons_trim1)[[i]]
      png(paste("TrimDensityPolygon_", spp, " ", name, ".png", sep = ""))
      plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
      plot(densityPolygons_trim1[[i]],add = T,col = viridis::viridis(99))
      dev.off()
    }
    
    png(paste("TrimDensityPolygon_", spp, " ALL.png", sep = ""))
    plot(bgLayer,col = "grey",colNA = "darkgrey",main = paste("Polygon, subspp:", name))
    cols = c( "black", "red", "blue", "green", "cyan", "magenta",
              "pink", "white", "purple", "orange", "yellow", "sienna",
              "thistle", "palegreen", "powderblue", "aquamarine", "violet", "mediumslateblue",
              "lightsalmon", "lightblue")
    for (i in 1:length(densityPolygons_trim1)) {
      name = names(densityPolygons_trim1)[[i]]
      plot(densityPolygons_trim1[[i]],add = T,border = cols[i],lwd = ((3 * i) / 3))
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
  
  polyLocations = labeledLoc
  
  for (slotA in 1:length(subsppNames)) {
    for (slotB in 1:length(subsppNames)) {
      if (subsppNames[[slotA]] != "unknown" && subsppNames[[slotB]] != "unknown" && slotA != slotB) {
        #print(paste(slotA,slotB,sep=" "))
        polyLocations = locatePolygonPoints(test_points = polyLocations,
                                            polygonA = densityPolygons_trim[[slotA]],polygonB = densityPolygons_trim[[slotB]],
                                            nameA = subsppNames[[slotA]],nameB = subsppNames[[slotB]],
                                            setcoord = T)
      }
    }
  }
  
  
  
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
    
  }
  
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
#    quantile=0.95,
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
#    quantile=0.95,
#    xmin=-125,
#    xmax=-60,
#    ymin=10,
#    ymax=50,
#    plotIt=T,
#    bgLayer=bg,
#    outputDir="/Users/kprovost/Documents/Classes/Spatial Bioinformatics/project/")
