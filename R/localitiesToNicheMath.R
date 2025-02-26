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
#' Wrapper for calculating niche math
#'
#' This function serves as a wrapper for the rest of the functions.
#'
#' @param species The species whose name to print.
#' @param Env The environmental variables to compare.
#' @param locs The localities to be processed.
#'
#' @export
#' @examples
#'
#' localitiesToNicheMath(Env,loc,species)
localitiesToNicheMath = function(Env,loc,species,rep1=10,rep2=1000,
                                 RMvalues=seq(0.5,4,0.5),
                                 fc=c("L", "LQ", "H"),numCores=1,
                                 method='block',verbose=T,
                                 runNicheModels=T,overwrite=F){

  loc_good_clean = cleanByEnvironment(Env,loc)
  loc_thin = spThinBySubspecies(loc_good_clean,species=occ_name,overwrite=overwrite)
  all_subspp = unique(loc_thin$assigned)
  #if(verbose==T){View(loc_thin)}

  bg_dat_filename = paste(species,"_bg_dat.txt",sep="")
  bg_bg_filename = paste(species,"_bg_bg.txt",sep="")

  if(!file.exists(bg_dat_filename) | !file.exists(bg_bg_filename) | overwrite==T) {
    print("Making bg_dat and bg_bg")
    loc_thin_bgstuff = generateBackgroundForPCA(localities = loc_thin,e=Env)
    bg_dat = loc_thin_bgstuff$bgenv
    bg_bg = loc_thin_bgstuff$bgpoints
    write.table(bg_dat,bg_dat_filename)
    write.table(bg_bg,bg_bg_filename)
  } else {
    print("Importing bg_dat and bg_bg")
    bg_dat = read.table(bg_dat_filename)
    bg_bg = read.table(bg_bg_filename)
  }

  perspecies_bgstuff = generateBackgroundPerSpecies(localities = loc_thin,e=Env)

  bgenv_by_subspecies = perspecies_bgstuff$bgenv_by_subspecies
  ## bgenv_by_subspecies has list of one df per subspecies
  lapply(1:length(bgenv_by_subspecies), FUN=function(i) {
    ## iterates over the df
    name = names(bgenv_by_subspecies)[i]
    df = bgenv_by_subspecies[[i]]
    bgenv_sub_filename = paste(species,"_",name,"_bgenv_by_subspecies.txt",sep="")
    write.table(df,bgenv_sub_filename,row.names = F)
  })

  bgext_by_subspecies = perspecies_bgstuff$bgext_by_subspecies
  ## bgext_by_subspecies has list of one df per subspecies
  lapply(1:length(bgext_by_subspecies), FUN=function(i) {
    ## iterates over the df
    name = names(bgext_by_subspecies)[i]
    df = bgext_by_subspecies[[i]]
    bgext_sub_filename = paste(species,"_",name,"_bgext_by_subspecies.txt",sep="")
    write.table(df,bgext_sub_filename,row.names = F)
  })

  bgpoints_by_subspecies = perspecies_bgstuff$bgpoints_by_subspecies
  ## bgpoints_by_subspecies has list of one df per subspecies
  lapply(1:length(bgpoints_by_subspecies), FUN=function(i) {
    ## iterates over the df
    name = names(bgpoints_by_subspecies)[i]
    df = bgpoints_by_subspecies[[i]]
    bgpoints_sub_filename = paste(species,"_",name,"_bgpoints_by_subspecies.txt",sep="")
    write.table(df,bgpoints_sub_filename,row.names = F)
  })

  bgextbg_by_subspecies = perspecies_bgstuff$bgextbg_by_subspecies
  ## bgextbg_by_subspecies has list of one df per subspecies
  lapply(1:length(bgextbg_by_subspecies), FUN=function(i) {
    ## iterates over the df
    name = names(bgextbg_by_subspecies)[i]
    df = bgextbg_by_subspecies[[i]]
    bgextbg_sub_filename = paste(species,"_",name,"_bgextbg_by_subspecies.txt",sep="")
    write.table(df,bgextbg_sub_filename,row.names = F)
  })

  pcaOutput = createPcaToCompare(loc_thin_bgstuff,perspecies_bgstuff,species) ## apparently not fixed yet
  pca_grid_clim = pcaOutput$grid_clim
  if(verbose==T){print("finished createPcaToCompare")}
  overlap_filename = paste(species,"_overlap.txt",sep="")
  if(!file.exists(overlap_filename) | overwrite==T) {
    pairwiseNicheEquivalence(pca_grid_clim,rep1=rep1,rep2=rep2,species=species)
    if(verbose==T){print("finished pairwiseNicheEquivalence")}
    overlap_df = pairwiseNicheOverlap(pca_grid_clim)
    write.table(overlap_df,overlap_filename)
    if(verbose==T){print("finished pairwiseNicheOverlap")}
  } else {
    print("Overlap file already exists! Not running")
  }
  if(runNicheModels==T){
    listENMresults = lapply(1:length(perspecies_bgstuff$bgpoints_by_subspecies),function(i){
      ##TODO: add in bg points from above
      subspp = names(perspecies_bgstuff$bgpoints_by_subspecies)[[i]]
      print(paste("Running",subspp))
      res = ENMevaluate(occs=perspecies_bgstuff$bgpoints_by_subspecies[[i]], envs = Env, partitions=method,
                        algorithm="maxnet",tune.args=list(fc=fc,rm=RMvalues),
                        parallel=T, numCores=numCores)
      #names(res) = names(nitens_by_subspp)[[i]]
      return(res)
    })
    names(listENMresults) = names(perspecies_bgstuff$bgpoints_by_subspecies)
    #if(verbose==T){View(listENMresults)}
    return(listENMresults)
  }
  else {
    print("SKIPPING NICHE MODELING")
    return(NULL)
  }

}

