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
#' Outlier Detection
#'
#' Calculates probability values for p_x given an epsilon
#'
#' @param localities A list of localities to check for anomalies
#' @param epsilon An value with which to flag anomalies with probability less than epsilon

#' @export
#'
#'
mgdad = function(localities,epsilon=0.0001) {
  p_x = mgdad_px(localities)
  anomalies = which(p_x < epsilon)
  ## anything where the probability is less than epsilon is considered too improbable
  return(anomalies)
}
