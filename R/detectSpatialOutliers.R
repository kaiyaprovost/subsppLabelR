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
#' This does outlier detection on points using Multivariate Gaussian
#' Distribution anomaly detection. Returns anomalies, which are a named
#' set of integers where the names are the original rowname of the point
#' and the integers are the rowname of the dataframe, if they are not the
#' same.
#'
#' @param localities A list of localities to check for anomalies
#' @param fulldataset The full dataset this is taken from
#' @param epsilon An value with which to flag anomalies with probability less than epsilon
#'
#' @export
#'
#'
detectSpatialOutliers = function(localities,
                                 fulldataset = localities,
                                 epsilon = 0.0001) {

  space = localities[,c("longitude","latitude")]
  anomalies = mgdad(localities=space,epsilon=epsilon)
  #space_anomalies=mgdad(localities=space,epsilon=epsilon)
  #space_only_anomalies = space[space_anomalies,]
  #anomalies = which(do.call(paste0, fulldataset[,c("longitude","latitude")]) %in% do.call(paste0, space_only_anomalies))
  return(anomalies)
}

