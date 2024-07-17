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
#' @import sf
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
#' Clean points by the environment
#'
#' This function removes any localities that do not have environmental data
#' so that the niche can be calculated. It assumes the second and third columns
#' are longitude and latitude, respectively.
#'
#' @param Env The environmental variables to calculate niche from
#' @param loc The locality data
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
cleanByEnvironment = function(Env,loc,latname="latitude",lonname="longitude"){
  latnum=which(colnames(loc)==latname)
  lonnum=which(colnames(loc)==lonname)
  loc[,latnum] = as.numeric(loc[,latnum])
  loc[,lonnum] = as.numeric(loc[,lonnum])
  extr = raster::extract(Env, loc[,c(lonnum,latnum)]) ## gets values from Env that are at loc
  head(extr)
  loc_clean = loc[!is.na(extr[,1]),]
  print(paste("Removed",nrow(loc)-nrow(loc_clean),"rows with no Env data"))

  return(loc_clean)
}

