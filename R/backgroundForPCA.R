#' @import caret
#' @import doFuture
#' @import dplyr
#' @import ggplot2
#' @import grid
#' @import gtools
#' @import gridExtra
#' @import h2o
#' @import MASS
#' @import plyr
#' @import raster
#' @import rebird
#' @import sp
#' @import spocc
#' @import viridis
#' @import AppliedPredictiveModeling
#' @import RColorBrewer
#' @import ENMeval
#' @import ENMTools
#' @import ecospat
#' @import ade4
#' @import adehabitatMA
#' @import adehabitatHR
NULL
#' Generating background points for a pca
#'
#' This function generates background points for making a PCA of environment.
#'
#' @param localities Cleaned and thinned localities
#' @param r The radius in m with which to get background points
#' @param num The number of background points
#' @param e The environmental raster for niche
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#' loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
backgroundForPCA = function(localities=locs_thinned,r=1,num=(100*nrow(localities)),e=Env,verbose=T){
  localities=localities[,c("Longitude","Latitude")]
  localities_pol = terra::vect(as.matrix(localities),"points")
  bg1 = ENMTools::background.buffer(points=localities_pol, buffer.width = r,
                                    n = num, mask = as(e[[1]],"SpatRaster"),
                                    buffer.type="circles",return.type="points")
  bg1 = terra::as.data.frame(bg1,geom="XY")
  extract1 = na.omit(cbind(localities,
                           raster::extract(e, localities), rep(1, nrow(localities))))
  colnames(extract1)[ncol(extract1)] = 'occ'
  extbg1 = na.omit(cbind(bg1, raster::extract(e, bg1), rep(0, nrow(bg1))))
  colnames(extbg1)[ncol(extbg1)] = 'occ'
  colnames(extbg1)[1:2] = colnames(extract1)[1:2]
  dat1 = rbind(extract1, extbg1)
  return(list(bgenv=dat1,bgpoints=bg1,bgext=extract1,bgextbg=extbg1))
}
