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
  
  densityStack = raster::stack("/Users/kprovost/Documents/GitHub/subsppLabelR/Phainopepla nitens/DensityRaster_Phainopepla nitens_0.75.tif")
  densA = densityStack[[2]]
  densB = densityStack[[3]]
  plot(densityStack)
  
  ## convert to a pres/abs
  densA_pres = densA
  densB_pres = densB
  densA_pres[!(is.na(densA_pres))] = 1 
  densA_pres[(is.na(densA_pres))] = 0 
  densB_pres[!(is.na(densB_pres))] = 1 
  densB_pres[(is.na(densB_pres))] = 0 
  dens_AB_pres_dif = densA_pres - densB_pres
  dens_AB_pres_dif[densA_pres==0 & densB_pres == 0] = NA
  plot(dens_AB_pres_dif) ## if 0, overlap. if -1 or +1, only on one side. 
  
  ## probably should do this with the relative densities 
  densA_rel = densA
  densA_rel[is.na(densA_rel)] = 0  
  densB_rel = densB
  densB_rel[is.na(densB_rel)] = 0
  dens_AB_dif = densA_rel - densB_rel
  dens_AB_dif[densA_rel==0 & densB_rel == 0] = NA
  plot(dens_AB_dif)
  dens_AB_dif[dens_AB_dif>0] = 1
  dens_AB_dif[dens_AB_dif<0] = -1
  plot(dens_AB_dif)
  plot(dens_AB_pres_dif)
  ## can we take majority rule?
  dens_AB_pres_overlap = dens_AB_pres_dif
  dens_AB_pres_overlap[dens_AB_pres_overlap!=0] = NA
  plot(dens_AB_pres_overlap)
  
  ## this function locates intersections between rasters
  library(sp)
  library(raster)
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
