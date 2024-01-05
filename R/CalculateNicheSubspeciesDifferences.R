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
#' Clean points by the environment
#'
#' This function removes any localities that do not have environmental data
#' so that the niche can be calculated. It assumes the second and third columns
#' are longitude and latitude, respectively.
#'
#' @param Env The environmental variables to calculate niche from
#' @param loc The locality data
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
cleanByEnvironment = function(Env,loc,latname="latitude",lonname="longitude"){
  latnum=which(colnames(loc)==latname)
  lonnum=which(colnames(loc)==lonname)
  loc[,latnum] = as.numeric(loc[,latnum])
  loc[,lonnum] = as.numeric(loc[,lonnum])
  extr = raster::extract(Env, loc[,c(lonnum,latnum)]) ## gets values from Env that are at loc
  head(extr)
  loc_clean = loc[!is.na(extr[,1]),]
  print(paste("Removed",nrow(loc)-nrow(loc_clean),"rows with no Env data"))

  return(loc_clean)
}
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
#' Generating background points for a pca
#'
#' This function generates background points for making a PCA of environment.
#'
#' @param localities Cleaned and thinned localities
#' @param r The radius in m with which to get background points
#' @param num The number of background points
#' @param e The environmental raster for niche
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#' loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
backgroundForPCA = function(localities=locs_thinned,r=1,num=(100*nrow(localities)),e=Env,verbose=T){
  localities=localities[,c("Longitude","Latitude")]
  localities_pol = terra::vect(as.matrix(localities),"points")
  bg1 = ENMTools::background.buffer(points=localities_pol, buffer.width = r,
                                    n = num, mask = as(e[[1]],"SpatRaster"),
                                    buffer.type="circles",return.type="points")
  bg1 = terra::as.data.frame(bg1,geom="XY")
  extract1 = na.omit(cbind(localities,
                           raster::extract(e, localities), rep(1, nrow(localities))))
  colnames(extract1)[ncol(extract1)] = 'occ'
  extbg1 = na.omit(cbind(bg1, raster::extract(e, bg1), rep(0, nrow(bg1))))
  colnames(extbg1)[ncol(extbg1)] = 'occ'
  colnames(extbg1)[1:2] = colnames(extract1)[1:2]
  dat1 = rbind(extract1, extbg1)
  return(list(bgenv=dat1,bgpoints=bg1,bgext=extract1,bgextbg=extbg1))
}
#' Generating background points for a species
#'
#' This function generates background points for each species.
#'
#' @param localities Thinned cleaned localities.
#'
#' @export
#' @examples
#'
#' loc_good_clean = cleanByEnvironment(Env, loc)
#' locs_thinned = spThinBySubspecies(loc_good_clean,thin.par=10,reps=1,lat.col="latitude",long.col="longitude",spec.col="assigned")
#'loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
#' perspecies_bgstuff = backgroundPerSpecies(loc_thin)
backgroundPerSpecies = function(localities=loc_thin,verbose=T,name="assigned",e=Env){
  if(verbose==T){print("starting backgroundPerSpecies")}
  name_col = which(colnames(localities) == name)
  loc_thin_by_subspecies = split(localities, localities[,name_col])
  bgenv_by_subspecies = list()

  bgext_by_subspecies = list()
  bgpoints_by_subspecies = list()
  bgextbg_by_subspecies = list()

  for(i in 1:length(names(loc_thin_by_subspecies))){
    singleSubspp = loc_thin_by_subspecies[[i]]
    subsppName = names(loc_thin_by_subspecies)[[i]]
    single_bgstuff = backgroundForPCA(singleSubspp,e=e)
    bgenv_by_subspecies[[i]] = single_bgstuff$bgenv
    bgext_by_subspecies[[i]] = single_bgstuff$bgext
    bgpoints_by_subspecies[[i]] = single_bgstuff$bgpoints
    bgextbg_by_subspecies[[i]] = single_bgstuff$bgextbg

  }
  names(bgenv_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgext_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgpoints_by_subspecies) = names(loc_thin_by_subspecies)
  names(bgextbg_by_subspecies) = names(loc_thin_by_subspecies)

  return(list(bgenv_by_subspecies=bgenv_by_subspecies,
              bgext_by_subspecies=bgext_by_subspecies,
              bgpoints_by_subspecies=bgpoints_by_subspecies,
              bgextbg_by_subspecies=bgextbg_by_subspecies))
}
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
#' loc_thin_bgstuff = backgroundForPCA(localities=loc_good[,c("Longitude","Latitude")],r=200000,num=(100*nrow(localities)),e=Env)
#' perspecies_bgstuff = backgroundPerSpecies(loc_thin)
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
#' Test for pairwise niche equivalence
#'
#' This function uses the background PCA grid for climate that was generated to
#' calculate pairwise niche equivalences for the species. This uses
#' randomization to pairwise equivalence so you must set the parameters for
#' randomizing.
#'
#' @param pca_grid_clim The grid for the climate PCA
#' @param rep1 The number of repetitions for niche equivalency
#' @param rep2 The number of repeititons for niche similarity
#'
#' @export
#' @examples
#'
#' printPointsPdfSuspect(species,subspecies,bg,loc_suspect)
pairwiseNicheEquivalence = function(pca_grid_clim,rep1=10,rep2=1000,species,verbose=T){
  if(verbose==T){print("starting pairwiseNicheEquivalence")}
  for(i in 1:length(pca_grid_clim)){
    for(j in 1:length(pca_grid_clim)){
      if(i<j){
        print(paste(i,j))
        spp1_name = names(pca_grid_clim)[[i]]
        spp1 = pca_grid_clim[[i]]
        spp2_name = names(pca_grid_clim)[[j]]
        spp2 = pca_grid_clim[[j]]

        #eq.test <- ecospat.niche.equivalency.test_custom(z1=spp1, z2=spp2,
        #                                                 rep=rep1, alternative = "higher"
        #                                                 )
        eq.test = ecospat::ecospat.niche.equivalency.test(z1=spp1, z2=spp2,rep=rep1,
                                                          overlap.alternative = "higher", ## testing for niche conservatism
                                                          expansion.alternative = "lower",
                                                          stability.alternative = "higher",
                                                          unfilling.alternative = "lower"
        )
        pdf(paste("EquivalencyOverlapTests_",species,"_",spp1_name,"_",spp2_name,".pdf",sep=""))
        par(mfrow=c(2,1))
        ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
        ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
        dev.off()
        print(paste("Running niche similarity test for",spp1_name,"-",spp2_name))
        sim.test <- ecospat::ecospat.niche.similarity.test(z1=spp1, z2=spp2,
                                                           rep=rep2, overlap.alternative = "higher", ## testing for niche conservatism
                                                           expansion.alternative = "lower",
                                                           stability.alternative = "higher",
                                                           unfilling.alternative = "lower",
                                                           rand.type=2)
        pdf(paste("EquivalencyOverlapTests_",species,"_",spp1_name,"_",spp2_name,".pdf",sep=""))
        par(mfrow=c(2,2))
        ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "D", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "I", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        dev.off()
        print(paste("Running niche similarity test for",spp2_name,"-",spp1_name))
        sim.test2 <- ecospat::ecospat.niche.similarity.test(z1=spp2, z2=spp1,
                                                            rep=rep2, overlap.alternative = "higher", ## testing for niche conservatism
                                                            expansion.alternative = "lower",
                                                            stability.alternative = "higher",
                                                            unfilling.alternative = "lower",
                                                            rand.type=2)

        pdf(paste("EquivalencyOverlapTests_",species,"_",spp1_name,"_",spp2_name,".pdf",sep=""))
        par(mfrow=c(2,3))
        ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "D", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        ecospat.plot.overlap.test(sim.test2, "D", paste("Similarity ",spp2_name,"->",spp1_name,sep=""))
        ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "I", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        ecospat.plot.overlap.test(sim.test2, "I", paste("Similarity ",spp2_name,"->",spp1_name,sep=""))
        dev.off()

      }
    }
  }
}
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
  #if(verbose==T){View(loc_thin)}
  loc_thin_bgstuff = backgroundForPCA(localities = loc_thin,e=Env)
  bg_dat = loc_thin_bgstuff$bgenv
  bg_bg = loc_thin_bgstuff$bgpoints
  perspecies_bgstuff = backgroundPerSpecies(localities = loc_thin,e=Env)
  pcaOutput = createPcaToCompare(loc_thin_bgstuff,perspecies_bgstuff,species)
  if(verbose==T){print("finished createPcaToCompare")}
  overlap_filename = paste(species,"_overlap.txt",sep="")
  if(!file.exists(overlap_filename) | overwrite==T) {
    pca_grid_clim = pcaOutput$grid_clim
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
#'
