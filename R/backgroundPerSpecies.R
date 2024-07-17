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
#' Generating background points for a species
#'
#' This function generates background points for each species.
#'
#' @param localities Thinned cleaned localities.
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#'loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
#' perspecies_bgstuff = backgroundPerSpecies(loc_thin)
backgroundPerSpecies = function(localities=loc_thin,verbose=T,name="assigned",e=Env){
  if(verbose==T){print("starting backgroundPerSpecies")}
  name_col = which(colnames(localities) == name)
  loc_thin_by_subspecies = split(localities, localities[,name_col])
  bgenv_by_subspecies = list()

  bgext_by_subspecies = list()
  bgpoints_by_subspecies = list()
  bgextbg_by_subspecies = list()

  for(i in 1:length(names(loc_thin_by_subspecies))){
    singleSubspp = loc_thin_by_subspecies[[i]]
    subsppName = names(loc_thin_by_subspecies)[[i]]
    single_bgstuff = backgroundForPCA(singleSubspp,e=e)
    bgenv_by_subspecies[[i]] = single_bgstuff$bgenv
    bgext_by_subspecies[[i]] = single_bgstuff$bgext
    bgpoints_by_subspecies[[i]] = single_bgstuff$bgpoints
    bgextbg_by_subspecies[[i]] = single_bgstuff$bgextbg

  }
  names(bgenv_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgext_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgpoints_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgextbg_by_subspecies) = names(loc_thin_by_subspecies)

  return(list(bgenv_by_subspecies=bgenv_by_subspecies,
              bgext_by_subspecies=bgext_by_subspecies,
              bgpoints_by_subspecies=bgpoints_by_subspecies,
              bgextbg_by_subspecies=bgextbg_by_subspecies))
}
