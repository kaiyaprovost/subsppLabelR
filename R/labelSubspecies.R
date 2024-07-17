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
#' Convert Occurence List to Labeled Dataframe
#'
#' This function converts a labeled list of subspecies occurences
#' into a dataframe of occurences, with a column for subspecies, using occ2dfSubspeciesLabels().
#' The latter function works only on a single subspecies whereas this function
#' performs the task for all subspecies in the list.
#'
#' @param subsppOccList A list of subspecies occurences from subspeciesOccQuery() or similar
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' subsppNames = unique(labeledLoc$subspecies)
labelSubspecies = function(subsppOccList,spp,subsppList) {
  ## this function takes a list of three taxa and labels them with subspecies information
  ## TODO: turn this into a function to give multiple subspecies and return it
  ## TODO: doesn't work if one subspp has zero records
  sppOcc = subsppOccList[[1]]
  name_sppOcc = names(subsppOccList)[1]
  #print("Giving occ2df labels")
  sppLocLab = occ2dfSubspeciesLabels(subsppOccList_object = sppOcc,
                                      subsppOccList_name = name_sppOcc)
  labeledOcc = subsppOccList[[2]]
  if(is.null(labeledOcc)) {
    print("No subspecies present -- returning for species only")
  } else {
    #print(paste("Length labeledOcc:",length(labeledOcc)))
    for (occ in 1:length(labeledOcc)) {
      print(names(labeledOcc)[[occ]])
      subsppLoc = occ2dfSubspeciesLabels(subsppOccList_object = labeledOcc[[occ]],subsppOccList_name = names(labeledOcc)[[occ]])
      #print("check1")
      sppLocLab = rbind(sppLocLab, subsppLoc)
      #print("check2")
    }
  }
  return(sppLocLab)
}
