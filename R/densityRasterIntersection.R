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
  plot(densityStack)
  
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
  myCells = which(!is.na(values(dens_AB_rel_overlap)))  
  
  adjCells = lapply(myCells,FUN=function(x){
    adjacent(dens_AB_rel_overlap,x,directions=8,target=myCells,sorted=T,include=T,pairs=F)
  })
  names(adjCells) = myCells
  
  posDirection <- matrix(c(NA, NA, NA, 
                          NA, 0,  1, 
                          1, 1, 1), ncol=3, byrow=TRUE)
  adjCells = adjacent(dens_AB_rel_overlap,myCells,directions=posDirection,target=myCells,sorted=T,include=F)
  adjCells2 = as.data.frame(adjCells)
  adjCells2$keep = adjCells2$from<adjCells2$to
  adjCells2 = adjCells2[adjCells2$keep==T,]
  
  
  
  
  findLabelAdjCell = function(myCell=myCell,ras=dens_AB_rel_overlap,allCells=myCells,label=myCell,cellLabels=NULL,verbose=T) {
    if(verbose==T) { print(myCell) }
    
    adjCells = adjacent(ras,myCell,directions=8,target=allCells[allCells>=myCell],sorted=T,include=T,pairs=F)

    for(newCell in adjCells[adjCells>myCell]){
      if(verbose==T) { print("LOOP") }
      
      newAdjCells = findLabelAdjCell(myCell=newCell,ras=ras,allCells=allCells[!(allCells %in% adjCells)],label=myCell,cellLabels=cellLabels)
      adjCells = sort(intersect(adjCells,newAdjCells))
    }

    return(adjCells)
  }
  
  adjCells = findLabelAdjCell(myCell=myCell,ras=dens_AB_rel_overlap,allCells=myCells,label=myCell,cellLabels=NULL,verbose=T)
  
  ## given a cell, find cells that are adjacent to it or to a neighbor etc the whole way down
  getAllTouchingCells = function(myCell,foundCells=NULL,ras,allCells,checkedCells=NULL) {
    if(is.null(foundCells)) { foundCells = c(myCell) }
    if(verbose==T) { print(myCell) }
    if(verbose==T) { print(foundCells) }
    
    ## get all cells touching this cell
    adjCells = adjacent(ras,myCell,directions=8,target=allCells[!(allCells %in% foundCells)],sorted=T,include=F,pairs=F)
    foundCells = sort(union(adjCells,foundCells))
    
    ## see if all of the found cells were checked 
    outersect <- function(x, y) {
      sort(c(setdiff(x, y),
             setdiff(y, x)))
    }
    
    if(is.null(checkedCells)) { checkedCells = c(myCell) }
    cellsLeft = outersect(checkedCells,foundCells)
    print(cellsLeft)
    if(length(cellsLeft)>0) {
      if(verbose==T) { print("RECURSE") }
      myList = getAllTouchingCells(myCell=cellsLeft[1],
                          foundCells=c(foundCells,cellsLeft[1]),
                          ras=ras,
                          allCells=allCells[!(allCells %in% intersect(foundCells,checkedCells))],
                          checkedCells=checkedCells)
      reFound = myList[[1]]
      reChecked = myList[[2]]
      if(verbose==T) { print("DONE RECURSE") }
      foundCells = sort(union(reFound,foundCells,myCell,cellsLeft[1]))
      checkedCells = sort(union(reChecked,checkedCells,cellsLeft[1]))
    }
    
    return(list(foundCells=foundCells,checkedCells=checkedCells))
  }
  getAllTouchingCells(myCell=myCell,foundCells=NULL,ras=dens_AB_rel_overlap,allCells,checkedCells=NULL)
  
  ## base case: there are no cells in a direction touching the cell that have a cell index larger than that cell
  ## if that is the case, it needs to return itself?
  ## you just need to check if every thing it returns is a cell or not 
  
  touchingCells = function(myCell,neighborhood,ras,allCells) {
    
    ##set up adjacency matrix to look at each direction
    middleRight <- matrix(c(NA, NA, NA, 
                            NA, 0,  1, 
                            NA, NA, NA), ncol=3, byrow=TRUE)
    bottomLeft <- matrix(c(NA, NA, NA, 
                           NA, 0,  NA, 
                           1, NA, NA), ncol=3, byrow=TRUE)
    bottomCenter <- matrix(c(NA, NA, NA, 
                             NA, 0,  NA, 
                             NA, 1, NA), ncol=3, byrow=TRUE)
    bottomRight <- matrix(c(NA, NA, NA, 
                            NA, 0,  NA, 
                            NA, NA, 1), ncol=3, byrow=TRUE)
    
    if(neighborhood=="middleRight") {
      adjCell = adjacent(ras,myCell,directions=middleRight,target=allCells,sorted=T,include=F,pairs=F)
    } else if (neighborhood=="bottomLeft") {
      adjCell = adjacent(ras,myCell,directions=bottomLeft,target=allCells,sorted=T,include=F,pairs=F)
    } else if (neighborhood=="bottomCenter") {
      adjCell = adjacent(ras,myCell,directions=bottomCenter,target=allCells,sorted=T,include=F,pairs=F)
    } else if (neighborhood=="bottomRight") {
      adjCell = adjacent(ras,myCell,directions=bottomRight,target=allCells,sorted=T,include=F,pairs=F)
    } else {
      stop("neighborhood value not accepted")
      adjCell = adjacent(ras,myCell,directions=8,target=allCells,sorted=T,include=F,pairs=F)
    }
    
    if(length(adjCell) == 0 ) {
      ## we have reached the end of this path 
      print("END OF PATH")
      print(myCell)
      return(myCell)
    } else {
      ## there are more cells to check
      print("CHECK EACH DIRECTION")
      print("MR")
      mR = touchingCells(adjCell,neighborhood="middleRight",ras,allCells)
      print("BL")
      bL = touchingCells(adjCell,neighborhood="bottomLeft",ras,allCells)
      print("BC")
      bC = touchingCells(adjCell,neighborhood="bottomCenter",ras,allCells)
      print("BR")
      bR = touchingCells(adjCell,neighborhood="bottomRight",ras,allCells)
      outputCells = sort(unique(c(mR,bL,bC,bR)))
      return(outputCells)
      
    }
  } 
  
  touchingCells(myCell=410,neighborhood="bottomLeft",ras,allCells)
  
  
  ## time for a list of lists
  for(myCell in myCells) {
    touching_cells_from = adjCells2$from[adjCells2$from==myCell | adjCells2$to==myCell]
    touching_cells_to = adjCells2$to[adjCells2$from==myCell | adjCells2$to==myCell]
    touching_cells = intersect(touching_cells_from,touching_cells_to)
  }
  
  
  
  
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
