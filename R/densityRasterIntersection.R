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
  ## this needs to take in two density rasters
  ## and then it needs to remove the overlaps of the density rasters
  ## the logic for this:
  ## find where they overlap
  library(sf)
  library(raster)

  myFolder = "C:/Users/kaiya/Documents/Work/GitHub/subsppLabelR/"
  setwd(myFolder)

  densityStack = raster::stack("Phainopepla nitens/DensityRaster_Phainopepla nitens_0.75.tif")
  densA = densityStack[[2]]
  densB = densityStack[[3]]
  plot(densityStack[[2:3]])

  densAR = raster("Phainopepla nitens/Phainopepla nitens lepida_raw_raster.tif")
  densBR = raster("Phainopepla nitens/Phainopepla nitens nitens_raw_raster.tif")

  presAbs = function(densA,densB) {
    ## convert to a pres/abs
    densA_pres = densA
    densB_pres = densB
    densA_pres[!(is.na(densA_pres))] = 1
    densA_pres[(is.na(densA_pres))] = 0
    densB_pres[!(is.na(densB_pres))] = 1
    densB_pres[(is.na(densB_pres))] = 0
    dens_pres_abs_overlap = densA_pres - densB_pres
    dens_pres_abs_overlap[densA_pres==0 & densB_pres == 0] = NA
    ## dens_pres_abs_overlap is a raster where if the two rasters overlap, it has a
    ## value of 0. if only one raster is present on that cell, it is a -1 or +1
    ## where it is -1 if B is present and +1 if A is present
    return(dens_pres_abs_overlap)
  }

  dens_pres_abs_overlap = presAbs(densA,densB)
  plot(dens_pres_abs_overlap) ## if 0, overlap. if -1 or +1, only on one side.
  ## note: presAbs doesn't work with the raw rasters, they are all complete

  relOverlap = function(densA,densB) {
    densA_rel = densA
    densA_rel[is.na(densA_rel)] = 0
    densB_rel = densB
    densB_rel[is.na(densB_rel)] = 0
    dens_relative_bool = densA_rel - densB_rel
    dens_relative_bool[densA_rel==0 & densB_rel == 0] = NA
    return(dens_relative_bool)
  }

  dens_relative_values = relOverlap(densA,densB)
  dens_relative_raw_values = relOverlap(densAR,densBR)
  par(mfrow=c(1,2))
  plot(dens_relative_values)
  plot(dens_relative_raw_values)

  convertRelToPresAbs = function(dens_relative_bool){
    dens_relative_bool[dens_relative_bool>0] = 1
    dens_relative_bool[dens_relative_bool<0] = -1
    return(dens_relative_bool)
  }

  dens_relative_bool = convertRelToPresAbs(dens_relative_values)
  dens_relative_raw_bool = convertRelToPresAbs(dens_relative_raw_values)
  par(mfrow=c(1,2))
  plot(dens_relative_bool)
  plot(dens_relative_raw_bool)

  ## can we take majority rule?
  dens_pres_abs_only_overlap = dens_pres_abs_overlap
  dens_pres_abs_only_overlap[dens_pres_abs_only_overlap!=0] = NA

  ## calculate with relative densities, which range from 0 to 1
  plot(dens_pres_abs_only_overlap)
  ## now get from the other one
  dens_AB_rel_overlap = dens_relative_bool
  dens_AB_rel_overlap[is.na(dens_pres_abs_only_overlap)] = NA
  plot(dens_AB_rel_overlap)

  ## assign the overlapping raster area to individual clusters
  ## run the recursion for the non raw rasters
  ras = dens_AB_rel_overlap
  myCells = which(!is.na(values(ras)))
  targetCells = myCells

  loopFindClusters <- function(ras,allCells=which(!is.na(values(ras))),targetCells=allCells,verbose=F) {
    visited = c()
    clusterList = list()
    for(cell in targetCells) {
      if(verbose) {print("MY CELL")}
      if(verbose) {print(cell)}
      if(cell %in% visited) { next }  else {
        if(verbose) {print("TARGET CELLS")}
        if(verbose) {print(targetCells)}
        returned = recurseNeighbors(cell=cell,visited=visited,ras=ras,allCells=allCells,verbose=verbose)
        cluster = c(returned$cluster)
        visited = c(returned$visited)
        if(verbose) {print("MY CLUSTER")}
        if(verbose) {print(cluster)}
        clusterList = c(clusterList,list(cluster))
        if(verbose) {print("MY VISITED")}
        if(verbose) {print(visited)}
      }

    }
    return(clusterList)
  }

  clusterList = loopFindClusters(ras=dens_AB_rel_overlap,allCells=myCells,targetCells=myCells,verbose=F)

  identified_overlap_clusters = dens_AB_rel_overlap
  for(i in 1:length(clusterList)) {
    clusterCells = clusterList[[i]]
    identified_overlap_clusters[clusterCells] = i
  }
  plot(identified_overlap_clusters)






  ## find the cluster with the densest cell for each raster
  ## i can use this function to check if the highest density area is connected to the polygon tho which is nice

  densCellsB = which(!is.na(values(densB)))
  ## get the max density cell
  targetB = which(values(densB)==max(values(densB),na.rm=T))
  densestClusterListB = loopFindClusters(ras=densB,
                                        allCells=densCellsB,
                                        targetCells=targetB,
                                        verbose=F)

  plot(densB)
  densest_B = densB
  densest_B[densCellsB] = 0
  densest_B[densestClusterListB[[1]]] = 1
  plot(densest_B)

  densCellsA = which(!is.na(values(densA)))
  ## get the max density cell
  targetA = which(values(densA)==max(values(densA),na.rm=T))
  densestClusterListA = loopFindClusters(ras=densA,
                                        allCells=densCellsA,
                                        targetCells=targetA,
                                        verbose=F)

  plot(densA)
  densest_A = densA
  densest_A[densCellsA] = 0
  densest_A[densestClusterListA[[1]]] = 1
  plot(densest_A)



  par(mfrow=c(1,3))
  plot(identified_overlap_clusters,col=c("red","blue"))
  plot(densest_A); plot(identified_overlap_clusters,add=T,col="red")
  plot(densest_B); plot(identified_overlap_clusters,add=T,col="red")

  valid_polygon_A = densest_A
  valid_polygon_B = densest_B

  ## use the densest cluster to assign the thing
  ## for each polygon in the assigned raster, check if that polygon overlaps A and B most dense
  identified_overlap_clusters =

  for(overlapCluster in clusterList) {
    ## find which polygon this cluster belongs to on both A and B


    is_A_densest = (1 %in% densest_A[overlapCluster])
    is_B_densest = (1 %in% densest_B[overlapCluster])
    ## A is densest and B is not densest
    if(is_A_densest && !is_B_densest){
      ## identify the polygon for B
      ## get valid cells for B
      valid_B = intersect(densCellsB, overlapCluster)
      b_cluster = loopFindClusters(ras=densest_B,allCells=densCellsB,targetCells=valid_B,verbose=F)
      valid_polygon_B[b_cluster[[1]]]=NA

    }
    ## A is not densest and B is densest
    if(!is_A_densest && is_B_densest){
      ## identify the polygon for A
      ## get valid cells for A
      valid_A = intersect(densCellsA, overlapCluster)
      a_cluster = loopFindClusters(ras=densest_A,allCells=densCellsA,targetCells=valid_A,verbose=F)
      valid_polygon_A[a_cluster[[1]]]=NA

    }

    ## both A and B are not densest
    ## both A and B are densest

  }




  ## if tied, then break the tie with which has more density
  ## return each raster as a boolean after the assignment
















  ## can we dissolve the rasters somehow

  ## fine lets just loop over the damn raster and find all the adjacent cells
  ## dave helped me do the recursive code so let's do the recursive code he coded out

  dens_pres_abs_overlap_reassigned = dens_pres_abs_overlap
  ## ok now we have the cluster list
  ## check each cluster's cells and get the majority rule or the sum
  for(cluster in clusterList) {
    valuesA = densA[cluster]
    valuesB = densB[cluster]
    sum(valuesA)
    sum(valuesB)
    winner = valuesA
    winner[valuesA > valuesB] = 1 ## +1
    winner[valuesB > valuesA] = -1 ## -1
    assigned = names(which(table(winner)==max(table(winner))))
    dens_pres_abs_overlap_reassigned[cluster] = as.numeric(assigned)
  }

  plot(dens_pres_abs_overlap)
  plot(dens_AB_rel_overlap)
  plot(dens_pres_abs_overlap_reassigned)






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
#' Recursively find neighbors
#'
#' This function XXXXXXXX
#'
#' @param cell XXX
#' @param visited xXX
#' @param ras xXX
#' @param allCells xXX
#' @param verbose xXX
#'
#' @export
#' @examples
#'
recurseNeighbors = function(cell,visited,ras,allCells,verbose=F) {
  if(verbose) {print("CELLS VISITED")}
  if(verbose) {print(visited)}
  ## base case: cell is visited
  if (cell %in% visited) {
    if(verbose) {print(paste("VISITED",cell))}
    return(list(cluster=c(),visited=c(visited,cell)))
  }
  ## base case: cell is invalid (should not happen)
  ## base case: cell is empty
  if(is.na(ras[cell])) {
    if(verbose) {print(paste("EMPTY",cell))}
    return(list(cluster=c(),visited=c(visited,cell)))
  }
  ## recursive case: cell is valid and not visited
  ## get neighbors of cell
  if(verbose) {print(paste("NOT VISITED",cell))}
  cellNeighbors = adjacent(x=ras,cells=cell,directions=8,target=allCells,include=F,sorted=T,pairs=F)
  cluster = c()
  visited = sort(unique(c(visited,cell)))
  if(verbose) {print("NEIGHBORS")}
  if(verbose) {print(cellNeighbors)}
  for (neighbor in cellNeighbors) {
    if(verbose) {print(paste("NEIGHBOR",neighbor))}
    returned = recurseNeighbors(neighbor,visited,ras,allCells,verbose)
    cluster = sort(unique(c(cluster,returned$cluster)))
    visited = sort(unique(c(visited,returned$visited)))
    if(verbose) {print("CLUSTER")}
    if(verbose) {print(cluster)}
  }
  return(list(cluster=unique(c(cluster,cell)),
              visited=unique(c(visited,cell))))
}

