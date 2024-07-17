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
occ2dfSubspeciesLabels = function(subsppOccList_object,
                                   subsppOccList_name) {
  ## thus function turns an occ object into a dataframe with a column for subspecies
  ## TODO: make it optional to do the "unique" thing for future processing
  ## TODO: does not work if zero records
  sppDf = data.frame(spocc::occ2df(subsppOccList_object))
  if (nrow(sppDf) <= 0) {
    sppLocLab = sppDf
    print("THIS SUBSPECIES HAS ZERO RECORDS")
  } else {
    sppLoc = unique(na.omit(sppDf[, 1:3]))
    sppLocLab = sppDf
    sppLocLab$subspecies = subsppOccList_name
  }
  ## TODO: add a way to make this output the dataframe
  #print("return")
  return(sppLocLab)
}
