#' @import raster
#' @import MASS
#' @import spocc
#' @import rgeos
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
NULL
#' Pull Subspecies Occurrences
#'
#' This function uses spocc::occ to query GBIF and other
#' repositories for one species and N number of
#' subspecies, specified by the user. Returns a list of
#' occurrence record objects.
#' Note: currently only tested for N=2 and N=3 subspecies.
#'
#' @param spp Genus and species to query, as string
#' @param subsppList Strings of subspecies to query
#' @param pointLimit Maximum point limit to return for each database -- see spocc::occ
#' @param dbToQuery List of databases to search through -- see spocc::occ
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
subspeciesOccQuery = function(spp,
                              subsppList = NULL,
                              pointLimit = 500,
                              dbToQuery = c("gbif", "bison", "inat", "ecoengine", "vertnet"),
                              ...) {
  ## TODO: support for no subspecies given
  ## this function uses spocc to query for one species and multiple subspecies
  #library(spocc)
  print(paste("Getting Species: ", spp))
  sppOcc = spocc::occ(query = spp,limit = pointLimit,has_coords = T,from = dbToQuery,...)
  if(is.null(subsppList)) {
    print("WARNING: No subspecies given -- not doing subspecies labeling")
    subSppListOcc = NULL
  } else {
    subSppListOcc = lapply(subsppList, function(x) {
      print(paste("     Getting Subspecies: ", x))
      to_return = spocc::occ(query = paste(spp, x, sep = " "),
                             limit = pointLimit,has_coords = T,from = dbToQuery)
      return(to_return)
    })
    names(subSppListOcc) = subsppList
    print(sppOcc)
    print(subSppListOcc)
  }
  toReturn = list(sppOcc, subSppListOcc)
  names(toReturn) = c("unknown", "labeled")
  return(toReturn)
}

