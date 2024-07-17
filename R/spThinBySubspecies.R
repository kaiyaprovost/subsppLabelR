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
#' @import rgeos
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
#' Thinning points by subspecies
#'
#' This function takes cleaned localities and applies a thinning algorithm so
#' that points are no closer than the thin.par km. Does this per species, or
#' whatever level of taxonomy is represented in the species column
#'
#' @param loc_good_clean Cleaned localities with associated Env data
#' @param thin.par The amount of km points should not be closer than
#' @param reps How many times to repeat the thinning
#' @param lat.col the column name for latitude
#' @param long.col The column name for longitude
#' @param spec.col The column name for species
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
spThinBySubspecies = function(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",
                              long.col="longitude",spec.col="assigned",verbose=T,write.files=T,
                              overwrite=F,species){
  if(verbose==T){
    print("starting spThinBySubspecies")
    print(unique(loc_good_clean$assigned))
  }

  locs_thinned=lapply(unique(loc_good_clean$assigned),FUN=function(subspp){
    print(subspp)
    loc_temp = loc_good_clean[loc_good_clean$assigned==subspp,]
    this_base= paste("thinned_data_",species,"_",subspp,sep="")
    this_file = paste(this_base,"_thin1.csv",sep="")

    if(overwrite==T){
      loc_thin = spThin::thin(loc.data = loc_temp,
                              lat.col = lat.col,
                              long.col = long.col,
                              spec.col = spec.col,
                              thin.par = thin.par, ## km distance that records need to be separated by
                              reps = reps, ## number of times to repeat thinning process
                              locs.thinned.list.return = T,
                              write.files = write.files,
                              max.files = 1,
                              write.log.file = F,
                              out.dir=getwd(),
                              out.base = this_base)[[1]]
      loc_thin$assigned = subspp
    } else {
      if(file.exists(this_file)){
        print("SKIPPING THINNING, FILE EXISTS")
        loc_thin = read.table(this_file,header=T,sep=",")
      } else {
        ## note: if you do not  change the  "name" column, it will error out and only use  the first subspecies.
        loc_thin = spThin::thin(loc.data = loc_temp,
                                lat.col = lat.col,
                                long.col = long.col,
                                spec.col = spec.col,
                                thin.par = thin.par, ## km distance that records need to be separated by
                                reps = reps, ## number of times to repeat thinning process
                                locs.thinned.list.return = T,
                                write.files = write.files,
                                max.files = 1,
                                write.log.file = F,
                                out.dir=getwd(),
                                out.base = this_base)[[1]]
        loc_thin$assigned = subspp
      }
    }
    colnames(loc_thin)[colnames(loc_thin)=="longitude"] = "Longitude"
    colnames(loc_thin)[colnames(loc_thin)=="latitude"] = "Latitude"
    return(loc_thin)
  })
  locs_thinned_df = do.call(gtools::smartbind,locs_thinned)
  return(locs_thinned_df)
}

