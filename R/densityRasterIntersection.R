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
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' 
#' 
#' 
densityRasterIntersection = function(densA,densB,verbose=F) {
  
  library(sf)
  library(raster)
  
  densityStack = raster::stack("/Users/kprovost/Documents/GitHub/subsppLabelR/Phainopepla nitens/DensityRaster_Phainopepla nitens_0.75.tif")
  densA = densityStack[[2]]
  densB = densityStack[[3]]
  plot(densityStack[[2:3]])
  
  densAR = raster("/Users/kprovost/Documents/GitHub/subsppLabelR/Phainopepla nitens/Phainopepla nitens lepida_raw_raster.tif")
  densBR = raster("/Users/kprovost/Documents/GitHub/subsppLabelR/Phainopepla nitens/Phainopepla nitens nitens_raw_raster.tif")
  
  presAbs = function(densA,densB) {
    ## convert to a pres/abs
    densA_pres = densA
    densB_pres = densB
    densA_pres[!(is.na(densA_pres))] = 1 
    densA_pres[(is.na(densA_pres))] = 0 
    densB_pres[!(is.na(densB_pres))] = 1 
    densB_pres[(is.na(densB_pres))] = 0 
    dens_AB_pres_dif = densA_pres - densB_pres
    dens_AB_pres_dif[densA_pres==0 & densB_pres == 0] = NA
    ## dens_AB_pres_dif is a raster where if the two rasters overlap, it has a 
    ## value of 0. if only one raster is present on that cell, it is a -1 or +1
    ## where it is -1 if B is present and +1 if A is present  
    return(dens_AB_pres_dif)
  }
  
  dens_AB_pres_dif = presAbs(densA,densB)
  plot(dens_AB_pres_dif) ## if 0, overlap. if -1 or +1, only on one side. 
  ## note: presAbs doesn't work with the raw rasters, they are all complete
  
  relOverlap = function(densA,densB) {
    densA_rel = densA
    densA_rel[is.na(densA_rel)] = 0  
    densB_rel = densB
    densB_rel[is.na(densB_rel)] = 0
    dens_AB_dif = densA_rel - densB_rel
    dens_AB_dif[densA_rel==0 & densB_rel == 0] = NA
    return(dens_AB_dif)
  }
  
  dens_AB_dif = relOverlap(densA,densB)
  dens_ABR_dif = relOverlap(densAR,densBR)
  par(mfrow=c(1,2))
  plot(dens_AB_dif)
  plot(dens_ABR_dif)
  
  convertRelToPresAbs = function(dens_AB_dif){
    dens_AB_dif[dens_AB_dif>0] = 1
    dens_AB_dif[dens_AB_dif<0] = -1
    return(dens_AB_dif)
  }
  
  dens_AB_dif = convertRelToPresAbs(dens_AB_dif)
  dens_ABR_dif = convertRelToPresAbs(dens_ABR_dif)
  par(mfrow=c(1,2))
  plot(dens_AB_dif)
  plot(dens_ABR_dif)
  
  ## can we take majority rule?
  dens_AB_pres_overlap = dens_AB_pres_dif
  dens_AB_pres_overlap[dens_AB_pres_overlap!=0] = NA
  
  ## calculate with relative densities, which range from 0 to 1 
  plot(dens_AB_pres_overlap)
  ## now get from the other one 
  dens_AB_rel_overlap = dens_AB_dif
  dens_AB_rel_overlap[is.na(dens_AB_pres_overlap)] = NA
  plot(dens_AB_rel_overlap)
  ## can we dissolve the rasters somehow
  
  ## fine lets just loop over the damn raster and find all the adjacent cells
  ## dave helped me do the recursive code so let's do the recursive code he coded out

  ras = dens_AB_rel_overlap
  myCells = which(!is.na(values(ras)))
  target = myCells
  visited = c()
  clusterList = list()
  for(cell in myCells) {
    #print("MY CELL")
    #print(cell)
    if(cell %in% visited) { next }  else {
      returned = recurseNeighbors(cell,visited)
      cluster = c(returned$cluster)
      visited = c(returned$visited)
      #print("MY CLUSTER")
      #print(cluster)
      clusterList = c(clusterList,list(cluster))
      #print("MY VISITED")
      #print(visited)
    }
  }
  print(clusterList)
  
  ## the way this is set up it won't work properly if there are NA cells in myCells
  recurseNeighbors = function(cell,visited) {
    #print("CELLS VISITED")
    #print(visited)
    ## base case: cell is visited
    if (cell %in% visited) {
      #print(paste("VISITED",cell))
      return(list(cluster=c(),visited=c(visited,cell)))
    }
    ## base case: cell is invalid (should not happen)
    ## base case: cell is empty (should not happen)
    ## recursive case: cell is valid and not visited
    ## get neighbors of cell 
    #print(paste("NOT VISITED",cell))
    cellNeighbors = adjacent(ras,cell,directions=8,target,include=F,sorted=T,pairs=F)
    cluster = c()
    visited = sort(unique(c(visited,cell)))
    #print("NEIGHBORS")
    #print(cellNeighbors)
    for (neighbor in cellNeighbors) {
      #print(paste("NEIGHBOR",neighbor))
      returned = recurseNeighbors(neighbor,visited)
      cluster = sort(unique(c(cluster,returned$cluster)))
      visited = sort(unique(c(visited,returned$visited)))
      #print("CLUSTER")
      #print(cluster)
    }
    return(list(cluster=unique(c(cluster,cell)),
                visited=unique(c(visited,cell))))
  }
  
  ## i can use this function to check if the highest density area is connected to the polygon tho which is nice
  
  ras = densB
  densCells = which(!is.na(values(ras)))
  target = densCells
  ## get the max density cell
  densestCell = which(values(ras)==max(values(ras),na.rm=T))
  visited = c()
  clusterList = list()
  print("MY CELL")
  print(densestCell)
  returned = recurseNeighbors(densestCell,visited)
  cluster = c(returned$cluster)
  visited = c(returned$visited)
  
  plot(densB)
  densBtest = densB
  densBtest[cluster] = 10
  plot(densBtest)
  
  
  spatrast = terra::rast(dens_AB_rel_overlap)
  starsrast = st_as_stars(spatrast)
  polygon = st_as_sf(starsrast)
  st_union(polygon)
  pbuffer = st_buffer(polygon,0.1)
  polygon1 = aggregate(polygon,pbuffer,mean,do_union=T,simplify=T,
                       join = function(x, y) st_is_within_distance(x, y, dist = 0.3))
  plot(polygon1)
  
  
  Mode <- function(x, na.rm = TRUE) {
    if(na.rm){
      x = x[!is.na(x)]
    }
    
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  
  polygon2 = aggregate(polygon, by = list(diss = polygon$layer), 
                       FUN = function(x)x[1], do_union = TRUE)
  plot(polygon2)
  pbuffer = st_buffer(polygon2,0.2)
  #polygon3 = st_union(polygon2, by_feature = TRUE)
  polygon3 = aggregate(polygon2,pbuffer,mean,do_union=T,simplify=F,
                       join = function(x, y) st_is_within_distance(x, y, dist = 0.3))
  plot(polygon3)
  
  
  
  
  polygon = st_as_sf(dens_AB_rel_overlap)
  polygon = raster::rasterToPolygons(dens_AB_rel_overlap,fun = NULL,na.rm = T,dissolve = F)
  plot(polygon,col=1:2)
  
  ## this function locates intersections between rasters
  
  #densA <- raster(nrows=40, ncols=40, xmn=0, xmx=2, ymn=0, ymx=2)
  #densA[] <- seq(1, 100, length.out=ncell(densA))
  #densB <- raster(outer(1:20,20:1), xmn=0, xmx=1, ymn=0, ymx=1)
  #densB <- reclassify(densB, c(0,100,NA))
  plot(densA)
  overlapRas = mask(crop(densB, densA), densA)
  plot(overlapRas)
  overlapRas = mask(crop(densA, densB), densB)
  plot(overlapRas)
  
  return(polygon)
}
