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
#' @import sp
#' @import sf
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
#' Generating a pca for comparison across niches.
#'
#' This function uses ecospat to create PCA across Environments to compare
#' niches. Requires some custom ecospat functions which are given below.
#'
#' @param loc_thin_bgstuff Background points for all species.
#' @param perspecies_bgstuff Background points for each species.
#' @param species Species name, for plotting
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#' loc_thin_bgstuff = generateBackgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
#' perspecies_bgstuff = generateBackgroundPerSpecies(loc_thin)
#' pcaOutput = createPcaToCompare(loc_thin_bgstuff,perspecies_bgstuffspecies)
#' pca_grid_clim = pcaOutput$grid_clim
createPcaToCompare = function(loc_thin_bgstuff,perspecies_bgstuff,species,verbose=T) {
  if(verbose==T){print("starting createPcaToCompare")}
  bg_dat = loc_thin_bgstuff$bgenv
  bg_bg = loc_thin_bgstuff$bgpoints

  bgext_by_subspecies = perspecies_bgstuff$bgext_by_subspecies
  bgenv_by_subspecies = perspecies_bgstuff$bgenv_by_subspecies

  ## pca bg points
  pca.env <- ade4::dudi.pca(bg_dat[,3:(ncol(bg_dat)-1)],scannf=F,nf=2)
  png(paste("PCAcorrelationCircle_",species,".png",sep=""))
  ecospat::ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
  dev.off()

  ## pca scores whole study area, all points, all subspecies
  scores_globclim<-pca.env$li # PCA scores for the whole study area (all points)

  ## now get pca scores per species instead

  scores = list()
  scores_clim = list()
  grid_clim = list()

  for(i in 1:length(names(bgext_by_subspecies))){
    singleSubspp_bgext = bgext_by_subspecies[[i]]
    singleSubspp_bgenv = bgenv_by_subspecies[[i]]
    subsppName = names(bgext_by_subspecies)[[i]]

    scores_subspp = ade4::suprow(pca.env,
                                 singleSubspp_bgext[which(
                                   singleSubspp_bgext[,ncol(singleSubspp_bgext)]==1)
                                   ,3:(ncol(singleSubspp_bgext)-1)])$li # PCA scores for the species 1 distribution


    scores[[i]] = scores_subspp

    scores_clim_subspp = ade4::suprow(pca.env,
                                      singleSubspp_bgenv[,3:(ncol(singleSubspp_bgenv)-1)])$li # PCA scores for the whole native study area species 1 ## bgenv

    scores_clim[[i]] = scores_clim_subspp

    ## make a dynamic occurrence densities grid
    ## grid of occ densities along one or two environmental gradients
    ## glob = env variables in background
    ## glob1 = env variables for species
    ## sp = occurrences of species
    ## R = resolution
    ## th.sp = a threshhold to elimite low density values of species occurrences
    # grid_clim_subspp <- ecospat.grid.clim.dyn_custom(glob = scores_globclim,
    #                                                  glob1 = scores_clim_subspp,
    #                                                  sp = scores_subspp,
    #                                                  R = 100,
    #                                                  th.sp = 0,
    #                                                  th.env = 0,
    #                                                  removeNA=T)
    grid_clim_subspp <- ecospat::ecospat.grid.clim.dyn(glob = scores_globclim,
                                                       glob1 = scores_clim_subspp,
                                                       sp = scores_subspp,
                                                       R = 100,
                                                       th.sp = 0,
                                                       th.env = 0)


    grid_clim[[i]] = grid_clim_subspp

  }

  names(scores) = names(bgext_by_subspecies)
  names(scores_clim) = names(bgext_by_subspecies)
  names(grid_clim) = names(bgext_by_subspecies)

  pdf(paste("NicheSpaceComparison_",species,".pdf",sep=""))
  par(mfrow=n2mfrow(length(grid_clim)),
      ask=F)
  for(i in 1:length(grid_clim)){
    plot(grid_clim[[i]]$w,main=names(grid_clim)[[i]])
  }
  dev.off()

  return(list(scores_globclim=scores_globclim,
              scores=scores,
              scores_clim=scores_clim,
              grid_clim=grid_clim))

}

