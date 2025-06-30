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
#' Test for pairwise niche equivalence
#'
#' This function uses the background PCA grid for climate that was generated to
#' calculate pairwise niche equivalences for the species. This uses
#' randomization to pairwise equivalence so you must set the parameters for
#' randomizing.
#' 
#' The default alternative settings will test for niche conservatism.
#'
#' @param pca_grid_clim The grid for the climate PCA
#' @param rep1 The number of repetitions for niche equivalency
#' @param rep2 The number of repeititons for niche similarity
#'
#' @export
#' @examples
#'
#' printPointsPdfSuspect(species,subspecies,bg,loc_suspect)
pairwiseNicheEquivalence = function(pca_grid_clim,rep1=1000,rep2=1000,species,verbose=T,
                                    overlap.alternative = "different",
                                    expansion.alternative = "different",
                                    stability.alternative = "different",
                                    unfilling.alternative = "different"){
  
  if(verbose==T){print("starting pairwiseNicheEquivalence")}
  for(i in 1:length(pca_grid_clim)){
    for(j in 1:length(pca_grid_clim)){
      if(i<j){
        print(paste(i,j))
        spp1_name = names(pca_grid_clim)[[i]] ## list?
        spp1 = pca_grid_clim[[i]]
        spp2_name = names(pca_grid_clim)[[j]] ## list
        spp2 = pca_grid_clim[[j]]

        eq.test = ecospat.niche.equivalency.test_custom(z1=spp1, z2=spp2,rep=rep1,
                                                          overlap.alternative = overlap.alternative, 
                                                          expansion.alternative = expansion.alternative,
                                                          stability.alternative =  stability.alternative,
                                                          unfilling.alternative = unfilling.alternative
        )
        pdf(paste("EquivalencyOverlapTests_",species,"_",spp1_name,"_",spp2_name,".pdf",sep=""))
        par(mfrow=c(2,1))
        ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
        ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
        dev.off()
        print(paste("Running niche similarity test for",spp1_name,"-",spp2_name))
        sim.test <- ecospat.niche.similarity.test_custom(z1=spp1, z2=spp2,
                                                           rep=rep2, overlap.alternative = overlap.alternative,
                                                           expansion.alternative = expansion.alternative,
                                                           stability.alternative = stability.alternative,
                                                           unfilling.alternative = unfilling.alternative,
                                                           rand.type=2)
        pdf(paste("EquivalencyOverlapTests_",species,"_",spp1_name,"_",spp2_name,".pdf",sep=""))
        par(mfrow=c(2,2))
        ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "D", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
        ecospat.plot.overlap.test(sim.test, "I", paste("Similarity ",spp1_name,"->",spp2_name,sep=""))
        dev.off()
        print(paste("Running niche similarity test for",spp2_name,"-",spp1_name))
        sim.test2 <- ecospat.niche.similarity.test_custom(z1=spp2, z2=spp1,
                                                            rep=rep2, overlap.alternative = overlap.alternative,
                                                            expansion.alternative = expansion.alternative,
                                                            stability.alternative = stability.alternative,
                                                            unfilling.alternative = unfilling.alternative,
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

