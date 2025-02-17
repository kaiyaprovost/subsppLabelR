% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/densityRasterRemoveIntersection.R
\name{loopFindClusters}
\alias{loopFindClusters}
\title{Loop to find recursive clusters and returns list of clusters}
\usage{
loopFindClusters(
  ras,
  allCells = which(!is.na(values(ras))),
  targetCells = allCells,
  verbose = F
)
}
\arguments{
\item{ras}{The raster}

\item{allCells}{All possible cells}

\item{targetCells}{The cells to be clustered}

\item{verbose}{Whether or not to be verbose}
}
\description{
This function loops over clusters to find them
}
