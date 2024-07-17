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
#' Distribution anomaly detection. Based off of
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
  ## TODO: why is this removing central anomalies?
  ## TODO: fix the code, taking the unique data is what is causing the issue, need to match

  mgdad = function(localities,epsilon=0.0001) {
    ## Multivariate Gaussian Distribution anomaly detection
    m=nrow(localities) ## the number of rows
    n = ncol(localities) ## the number of columns, aka dimensions
    mu = colMeans(localities,na.rm=T) ## takes the mean of each column
    Sigma = cov(localities) ## get the covariance matrix
    centered <- caret::preProcess(localities, method = "center") ## generates a scaling with which to center the values
    space_mu <- as.matrix(predict(centered, localities)) ## converts the data to the centered values
    A = (2 * pi) ^ (-n / 2) * det(Sigma) ^ (-0.5)
    B = exp(-0.5 * rowSums((space_mu %*% MASS::ginv(Sigma)) * space_mu))
    p_x = A * B ## this is the probability of each value coming from the same normal distribution
    anomalies = which(p_x < epsilon) ## anything where the probability is less than epsilon is considered too improbable
    return(anomalies)
  }

  space = localities[,c("longitude","latitude")]
  anomalies = mgdad(localities=space,epsilon=epsilon)
  #space_anomalies=mgdad(localities=space,epsilon=epsilon)
  #space_only_anomalies = space[space_anomalies,]
  #anomalies = which(do.call(paste0, fulldataset[,c("longitude","latitude")]) %in% do.call(paste0, space_only_anomalies))
  return(anomalies)
}

