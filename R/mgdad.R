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
#' Outlier Probability
#'
#' Calculates Multivariate Gaussian Distribution values for p_x
#'
#' @param localities A list of localities to check for anomalies
#'
#' @export
#'
#'
mgdad_px = function(localities) {
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
  return(p_x)
}

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
