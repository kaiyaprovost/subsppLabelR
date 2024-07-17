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
#' Make pairwise niche overlaps
#'
#' This function uses the background PCA grid for climate that was generated to
#' calculate pairwise niche overlaps for the species.
#'
#' @param pca_grid_clim The grid for the climate PCA
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#' loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
#' perspecies_bgstuff = backgroundPerSpecies(loc_thin)
#' pcaOutput = createPcaToCompare(loc_thin_bgstuff,perspecies_bgstuffspecies)
#' pca_grid_clim = pcaOutput$grid_clim
#' overlap_df  = pairwiseNicheOverlap(pca_grid_clim)
pairwiseNicheOverlap = function(pca_grid_clim,verbose=T){
  if(verbose==T){print("starting pairwiseNicheOverlap")}
  overlap_df = data.frame(spp1=character(),
                          spp2=character(),
                          SchoenersD=numeric(),
                          modifiedHellingersI=numeric())

  for(i in 1:length(pca_grid_clim)){
    for(j in 1:length(pca_grid_clim)){
      if(i<j){
        if(verbose==T){print(paste(i,j))}
        spp1_name = names(pca_grid_clim)[[i]]
        spp1 = pca_grid_clim[[i]]
        spp2_name = names(pca_grid_clim)[[j]]
        spp2 = pca_grid_clim[[j]]
        overlap <- ecospat::ecospat.niche.overlap(spp1, spp2, cor=T)
        rowToAdd = cbind(as.character(spp1_name),
                         as.character(spp2_name),
                         as.numeric(overlap$D),
                         as.numeric(overlap$I))
        colnames(rowToAdd) = colnames(overlap_df)
        overlap_df = rbind(overlap_df,rowToAdd)
      }
    }
  }
  return(overlap_df)
}

