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
#' Remove Raster Intersections
#'
#' This function finds intersections between raster and then removes them
#'
#' @param densA The first raster
#' @param densB The second raster
#' @param verbose Whether to print extra statements
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
#' XXXXXXXX
#'
densityRasterRemoveIntersection = function(densA,densB,verbose=F) {
  ## this needs to take in two density rasters
  ## and then it needs to remove the overlaps of the density rasters
  library(sf)
  library(raster)
  
  #myFolder = "C:/Users/kaiya/Documents/Work/GitHub/subsppLabelR/"
  #myFolder = "~/Documents/GitHub/subsppLabelR/"
  #setwd(myFolder)
  
  #densityStack = raster::stack("Phainopepla nitens/DensityRaster_Phainopepla nitens_0.75.tif")
  #densA = densityStack[[2]]
  #densB = densityStack[[3]]
  #plot(densityStack[[2:3]])
  
  #densAR = raster("Phainopepla nitens/Phainopepla nitens lepida_raw_raster.tif")
  #densBR = raster("Phainopepla nitens/Phainopepla nitens nitens_raw_raster.tif")
  
  dens_pres_abs_overlap = presAbs(densA,densB)
  #plot(dens_pres_abs_overlap) ## if 0, overlap. if -1 or +1, only on one side.
  ## note: presAbs doesn't work with the raw rasters, they are all complete
  
  dens_relative_values = relOverlap(densA,densB)
  #dens_relative_raw_values = relOverlap(densAR,densBR)
  #par(mfrow=c(1,2))
  #plot(dens_relative_values)
  #plot(dens_relative_raw_values)
  
  dens_relative_bool = convertRelToPresAbs(dens_relative_values)
  #dens_relative_raw_bool = convertRelToPresAbs(dens_relative_raw_values)
  #par(mfrow=c(1,2))
  #plot(dens_relative_bool)
  #plot(dens_relative_raw_bool)
  
  ## can we take majority rule?
  dens_pres_abs_only_overlap = dens_pres_abs_overlap
  dens_pres_abs_only_overlap[dens_pres_abs_only_overlap!=0] = NA
  
  ## calculate with relative densities, which range from 0 to 1
  #plot(dens_pres_abs_only_overlap)
  ## now get from the other one
  dens_AB_rel_overlap = dens_relative_bool
  dens_AB_rel_overlap[is.na(dens_pres_abs_only_overlap)] = NA
  #plot(dens_AB_rel_overlap)
  
  ## assign the overlapping raster area to individual clusters
  ## run the recursion for the non raw rasters
  ras = dens_AB_rel_overlap
  myCells = which(!is.na(values(ras)))
  targetCells = myCells
  
  clusterList = loopFindClusters(ras=dens_AB_rel_overlap,allCells=myCells,targetCells=myCells,verbose=F)
  
  identified_overlap_clusters = dens_AB_rel_overlap
  for(i in 1:length(clusterList)) {
    clusterCells = clusterList[[i]]
    identified_overlap_clusters[clusterCells] = i
  }
  #plot(identified_overlap_clusters)
  
  ## find the cluster with the densest cell for each raster
  ## i can use this function to check if the highest density area is connected to the polygon tho which is nice
  
  densCellsB = which(!is.na(values(densB)))
  ## get the max density cell
  targetB = which(values(densB)==max(values(densB),na.rm=T))
  densestClusterListB = loopFindClusters(ras=densB,
                                         allCells=densCellsB,
                                         targetCells=targetB,
                                         verbose=F)
  
  #plot(densB)
  densest_B = densB
  densest_B[densCellsB] = 0
  densest_B[densestClusterListB[[1]]] = 1
  #plot(densest_B)
  
  densCellsA = which(!is.na(values(densA)))
  ## get the max density cell
  targetA = which(values(densA)==max(values(densA),na.rm=T))
  densestClusterListA = loopFindClusters(ras=densA,
                                         allCells=densCellsA,
                                         targetCells=targetA,
                                         verbose=F)
  
  #plot(densA)
  densest_A = densA
  densest_A[densCellsA] = 0
  densest_A[densestClusterListA[[1]]] = 1
  #plot(densest_A)
  
  #par(mfrow=c(1,3))
  #plot(identified_overlap_clusters,col=c("red","blue"))
  #plot(densest_A); plot(identified_overlap_clusters,add=T,col="red")
  #plot(densest_B); plot(identified_overlap_clusters,add=T,col="red")
  
  valid_polygon_A = densest_A
  valid_polygon_B = densest_B
  
  ## use the densest cluster to assign the thing
  ## for each polygon in the assigned raster, check if that polygon overlaps A and B most dense
  #identified_overlap_clusters
  
  for(overlapCluster in clusterList) {
    ## overlapCluster = clusterList[[1]]
    ## find which polygon this cluster belongs to on both A and B
    ## get valid cells for A
    valid_A = intersect(densCellsA, overlapCluster)
    a_cluster = loopFindClusters(ras=densest_A,allCells=densCellsA,targetCells=valid_A,verbose=F)
    ## get valid cells for B
    valid_B = intersect(densCellsB, overlapCluster)
    b_cluster = loopFindClusters(ras=densest_B,allCells=densCellsB,targetCells=valid_B,verbose=F)
    ## check if these are the densest polygons for A and B
    is_A_densest = (1 %in% densest_A[overlapCluster])
    is_B_densest = (1 %in% densest_B[overlapCluster])
    ## calculate the density values for each cluster
    density_values_A = densA[a_cluster[[1]]]
    density_values_B = densB[b_cluster[[1]]]
    
    ## A is densest and B is not densest
    if(is_A_densest && !is_B_densest){
      ## remove those cells from B
      valid_polygon_B[b_cluster[[1]]]=NA
      
    }
    ## A is not densest and B is densest
    else if(!is_A_densest && is_B_densest){
      ## remove those cells from A
      valid_polygon_A[a_cluster[[1]]]=NA
    }
    ## both A and B are not densest
    else if(!is_A_densest && !is_B_densest){
      ## if tied, then break the tie with which has more density and remove the 
      ## polygon that is less dense
      ## A is more dense, get rid of B
      if(sum(density_values_A,na.rm=T) > sum(density_values_B,na.rm=T)) {
        valid_polygon_B[b_cluster[[1]]]=NA
      }
      ## B is more dense, get rid of A
      else if(sum(density_values_A,na.rm=T) < sum(density_values_B,na.rm=T)) {
        valid_polygon_A[a_cluster[[1]]]=NA
      }
      ## they are equally dense, DIE
      else if(sum(density_values_A,na.rm=T) == sum(density_values_B,na.rm=T)) {
        stop("ERROR: Density rasters should differ in their sums -- check your inputs!")
      }
      
    }
    ## both A and B are densest
    else if(is_A_densest && is_B_densest){
      ## if tied, then break the tie with which has more density, but keep both
      ## polygons, just remove the cells based on which is more dense 
      ## A is more dense, remove overlap cells from B
      if(sum(density_values_A,na.rm=T) > sum(density_values_B,na.rm=T)) {
        valid_polygon_B[overlapCluster[[1]]]=NA
      }   
      ## B is more dense, remove overlap cells from A
      else if(sum(density_values_A,na.rm=T) < sum(density_values_B,na.rm=T)) {
        valid_polygon_A[overlapCluster[[1]]]=NA
      }  
      ## they are equally dense, DIE
      else if(sum(density_values_A,na.rm=T) == sum(density_values_B,na.rm=T)) {
        stop("ERROR: Density rasters should differ in their sums -- check your inputs!")
      }
      
    }
  }
  
  #par(mfrow=c(2,2))
  #plot(densest_A)
  #plot(densest_B)
  #plot(valid_polygon_A)
  #plot(valid_polygon_B)
  
  ## need to return valid_polygon_A and valid_polygon_B but as the original density
  valid_polygon_A[!(is.na(valid_polygon_A))] = densA[!(is.na(valid_polygon_A))]
  valid_polygon_B[!(is.na(valid_polygon_B))] = densB[!(is.na(valid_polygon_B))]
  
  return(list(valid_polygon_A,valid_polygon_B))
}
#' Recursively find neighbors
#'
#' This function finds neighbors of a cell recursively
#'
#' @param cell The cell you are checking
#' @param visited Cells you have already visited
#' @param ras The raster you are checking
#' @param allCells All possible valid cells
#' @param verbose Whether or not to print extra statements
#'
#' @export
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
#' Calculate presence/absence raster from two density rasters
#'
#' This function turns two density rasters into a presence/absence boolean
#'
#' @param densA The first raster
#' @param densB The second raster
#'
#' @export
#'
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
#' Calculate relative overlap raster from two density rasters
#'
#' This function calculates the relative overlap from two density rasters
#'
#' @param densA The first raster
#' @param densB The second raster
#'
#' @export
#' @examples
#'
relOverlap = function(densA,densB) {
  densA_rel = densA
  densA_rel[is.na(densA_rel)] = 0
  densB_rel = densB
  densB_rel[is.na(densB_rel)] = 0
  dens_relative_bool = densA_rel - densB_rel
  dens_relative_bool[densA_rel==0 & densB_rel == 0] = NA
  return(dens_relative_bool)
}
#' Loop to find recursive clusters and returns list of clusters
#'
#' This function loops over clusters to find them
#'
#' @param ras The raster
#' @param allCells All possible cells
#' @param targetCells The cells to be clustered
#' @param verbose Whether or not to be verbose
#'
#' @export
#'
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
#' Convert relative raster to presence/absence raster
#'
#' This function converts relative values of a raster into a presence/absence
#'
#' @param ras dens_relative_bool
#'
#' @export
#'
convertRelToPresAbs = function(dens_relative_bool){
  dens_relative_bool[dens_relative_bool>0] = 1
  dens_relative_bool[dens_relative_bool<0] = -1
  return(dens_relative_bool)
}