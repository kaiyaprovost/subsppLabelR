##### NOTES FROM LAB MEETING #####
# ideas for subsppR
# -- mccormack 2010 jniche paper -- niche tests are not very robust
# -- scale issue -- why not just do by longitude or whatever?
#   -- what about the rotation of points or something?
#   -- GIGO though, the data points aren't good.
# -- are the son/chi desert different or the same? scale?
# -- are the climatic envelopes different at a biome scale?
# -- also phaino is a weird one to test
# -- have you tried using other layers?
# -- eclogy keeps them separated or not?
# -- first one -- are the deserts different -- does env drive differentiation, what about it? local adaptation? flipside of CFB -- locally adapted is same as saying cant cross CFB. but can be locally adapted and have the same niche so what's important s why can't cross barrier.
# -- second one is getting at the latter. try to understand why the barrier affects your species. then can set up to test ability of them to cross using niche models in comparative framework.
# -- consistent biases allows for answers
# -- to test the strength of the barrier - does resistance explain differentiation across the barrier.
# -- ideally with spatially explicit simulation. going to be difficult to get aorund points and scale issues of niche model.
# -- barrier is small geographic area with large scale data. but can refine data down.
# -- ditch ebird stuff totally from the gbif data. probably just stationary counts less than 10 minutes and will get tons of recrods.
# -- specimens may have too much lat long error.
# -- remember the error estimates.
# -- stationary counts, like <1 hour or something? see how much data that is. past 5 years? people using ebird app which does very good lat long. high res high quality points. factor in cell phone shittiness.
# -- most of these birds not a climatic issue. but is correlated with broad scale habitat?
# -- start with large scale stuff to see if it works? then it probalby wont and go for small stuff.
# -- how habitat associations correlate with climate?
# -- could do the ML though. LandSat? landscape classification?
# -- mary has a remose sensing package?

## -- hypervolumes are more accurate than minimum convex polygons! r package hypervolume. see the woodpecker color paper https://www.biorxiv.org/content/biorxiv/early/2018/07/23/375261.full.pdf




## input object: results of PairsOfBirdSpeciesNicheOverlap-3orMoreSpecies.R
## list of two objects
## first of which is a dataframe with columns: name (spp name), lat, long, subspecies (a priori),
## then some number of columns of boolean values as to whether a point
## is assigned to a subspecies or not
## second of which are polygons associated with each subspecies

## problem: multiple tests between 19 subspecies? bonferoni corrections
## maybe only do if adjacent? -- probably not, says the instructors
## may want to figure out another way to measure similarity to walk through subspp space
##TODO: consider lumping?


## for now testing with phainopeplaNitens
##TODO: change all the lapply to function if possible
##TODO: update plotting here (and in other file) so that everything plots iteratively, esp if pairwise
##TODO: add in species name not just subspp name

devtools::install_github('kaiyaprovost/subsppLabelR')
library(subsppLabelR)
## need to add a check in here -- remove unknown points with exact same lat/long as a labeled point

setwd("~/Documents/Classes/Finished Classes/Spatial Bioinformatics/project/")

THIN_BEYOND = FALSE

## also split these into subspecies assignments with $assigned
## this is a list, can access with $x or [[i]]

## if you want you can check through the points that suck to add in more localities
## but for now not doing that

## now we have groups that are assigned to single areas and don't mismatch
## let's split into two different groupings and make some niche models!
## then test overlap

# step 1 -- import environmental variables
##TODO: add in other data than Worldclim
##TODO: water layer, how often water there is at that spot? 30m layer!!!
Env = raster::stack(list.files(
  path='/Users/kprovost/Documents/Classes/Finished Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  pattern="\\.bil$",
  full.names=T))
ext = raster::extent(c(-125,-60,10,50)) ## make sure this will play nice with your points
Env = raster::crop(Env, ext)
bg = Env[[1]] ## just for plotting

## get locs
#species = "Phainopepla nitens"
#subspecies = c("nitens","lepida")

species = "Cardinalis sinuatus"
subspecies = c("sinuatus","fulvescens","peninsulae")
#species = "Cardinalis cardinalis"
#subspecies = c("affinis","canicaudus","cardinalis","carneus",
#                "clintoni","coccineus","flammiger","floridanus",
#                "igneus","littoralis","magnirostris","mariae",
#                "phillipsi","saturatus","seftoni","sinaloensis",
#                "superbus","townsendi","yucatanicus")

#subspecies_igne = c("affinis","clintoni","igneus","seftoni","sinaloensis","superbus","townsendi")
#subspecies_card = c("canicaudus","cardinalis","floridanus","magnirostris")
#subspecies_cocc = c("coccineus","flammiger","littoralis","phillipsi","yucatanicus")
#subspecies_rest = c("carneus","mariae","saturatus")
#test = c("clintoni","affinis")

## there is a bug -- if one subspp range is entirely subsumed within another polygon,
## will delete that subspecies. no bueno

alllocs = "/Users/kprovost/Documents/Classes/Finished Classes/Spatial Bioinformatics/project/big_sinuatus_testrun_NOTWORKING/AllLoci_Cardinalis sinuatus_sinuatus fulvescens peninsulae.txt"

library(subsppLabelR)

processedSpecies = databaseToAssignedSubspecies(spp=species,
                                                subsppList=subspecies,
                                                pointLimit=10,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
                                                quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
                                                outputDir="~/Documents/Classes/Finished Classes/Spatial Bioinformatics/project/",
                                                datafile=alllocs)

## THIS STILL ISN'T TRIMMING PROPERLY -- maybe do the loop twice?

















#processedSpecies_rest = databaseToAssignedSubspecies(spp=species,
#                                                 subsppList=subspecies_rest,
#                                                 pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#                                                 quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#                                                 outputDir="~/Documents/Classes/Spatial Bioinformatics/project/")
#processedSpecies_cocc = databaseToAssignedSubspecies(spp=species,
#                                                     subsppList=subspecies_cocc,
#                                                     pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#                                                     quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#                                                     outputDir="~/Documents/Classes/Spatial Bioinformatics/project/")
#processedSpecies_card = databaseToAssignedSubspecies(spp=species,
#                                                     subsppList=subspecies_card,
#                                                     pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#                                                     quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#                                                     outputDir="~/Documents/Classes/Spatial Bioinformatics/project/")
#processedSpecies_igne = databaseToAssignedSubspecies(spp=species,
#                                                     subsppList=subspecies_igne,
#                                                     pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
#                                                     quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
#                                                     outputDir="~/Documents/Classes/Spatial Bioinformatics/project/")


loc_suspect = processedSpecies$loc_suspect
loc_good = processedSpecies$loc_good
pol = processedSpecies$pol

write.table(loc_suspect,paste(paste("SuspectLoci",species,paste(subspecies,collapse=" "),sep="_"),".txt",sep=""),
            quote=FALSE,sep="\t",row.names=FALSE)
write.table(loc_good,paste(paste("GoodLoci",species,paste(subspecies,collapse=" "),sep="_"),".txt",sep=""),
            quote=FALSE,sep="\t",row.names=FALSE)

library(rgdal)
library(sp)
for (i in 1:length(pol)) {
  obj = pol[[i]]
  name = names(pol)[i]
  #print(name)
  obj2 = SpatialPolygonsDataFrame(obj,data=as.data.frame(rep(1,length(obj))),
                                  match.ID = FALSE)
  #print(obj)
  #print("---")
  writeOGR(obj=obj2,dsn=paste(paste("Polygon",species,name,sep="_"),".shp",sep=""),
           layer=name,driver="ESRI Shapefile")
  #print("end")
}

## BELOW: WRITE TO FIX FUNCTION DECISION TREE

flagPolygonOverlap2 = function(subsppPoly1=polA,subsppPoly2=polB){
  ##function(subsppPoly1,subsppPoly2)
  ## this function checks for overlaps between polygons
  ## TODO: remove polygon if not touching another of same spp
  ## but also closer to polygon of other species
  #library(rgeos)
  #library(raster)

  ## there is a bug -- if one subspp range is entirely subsumed within another polygon,
  ## will delete that subspecies. no bueno

  badList_subsppA_features = c()
  badList_subsppB_features = c()
  overlapsToRemove_subsppA = c()
  overlapsToRemove_subsppB = c()

  for (feature_subsppA in (1:length(subsppPoly1))){ ## get the features within the subspecies polygon
    for(feature_subsppB in (1:length(subsppPoly2))) { ## get the features within the subspecies polygon
      ## check areas
      totArea1 = rgeos::gArea(subsppPoly1) ## the whole area of the subspecies
      totArea2 = rgeos::gArea(subsppPoly2) ## the whole area of the subspecies
      area1 = rgeos::gArea(subsppPoly1[feature_subsppA,]) ## the area of the single feature
      area2 = rgeos::gArea(subsppPoly2[feature_subsppB,]) ## the area of the single feature

      # subsppPoly1$totalArea = totArea1
      # subsppPoly2$totalArea = totArea2
      # subsppPoly1[feature_subsppA,]$featureArea = area1
      # subsppPoly2[feature_subsppB]$featureArea = area2

      if(rgeos::gIntersects(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,])){
        #print("INTERSECTS")
        testInt = rgeos::gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,])
        #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
        #print(length(intersect))


        ## if they overlap

        #plot(subsppPoly1[feature_subsppA,],border="red",add=T)
        #plot(subsppPoly2[feature_subsppB,],border="cyan",add=T)
        testArea = rgeos::gArea(testInt)

        if(testArea > 0) {
          #print("OVERLAP AREA NOT ZERO")
          intersect = raster::intersect(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,])
          #intersect = gIntersection(subsppPoly1[feature_subsppA,],subsppPoly2[feature_subsppB,],byid=T)
          areaInt=rgeos::gArea(intersect)

          ## if they overlap

          area1percent = areaInt/totArea1 ## the area of the feature as a percent of the area of the subspecies
          area2percent = areaInt/totArea2

          ## NEW IF STATEMENTS
          if (area1 > area2) { ## if A big B small -- area1 vs area2
            if (testArea < area2) { ## if B not subsumed, remove from A -- testArea
              badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
            }
            else if (testArea >= area2) { ## else if B subsumed
              if (area2percent == 1) { ## if B = 100% of subspecies, remove from A -- area2percent
                badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              }
              else if (area2percent != 1) { ## else if B not 100% of subspecies, remove from B
                badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              }
            }
          }
          else if (area1 < area2) { ## else if A small B big
            if (testArea < area1) {## if A not subsumed, remove from B -- area1percent
              badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
            }
            else if (testArea >= area1) { ## else if A subsumed
              if (area1percent == 1) { ## if A = 100% of subspecies, remove from B
                badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
              }
              else if (area1percent != 1) { ## else if A not 100% of subspecies, remove from A
                badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
              }
            }
          }
          else if (area1 == area2) { ## else if A = B
            ## check total areas and remove larger
            if (totArea1 > totArea2) {
              ## if total A is greater than total B, remove A
              badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
            }
            else if (totArea1 < totArea2) { ## if total B is greater than total A, remove B
              badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
            }
            else if (totArea1 == totArea2) { ## if they are the exact same size
              if (testArea < area1 && testArea < area2) { ## if neither subsumed, remove arbitrary
                flip = sample(1:2,1)
                if (flip == 1) { ## flip coin remove A
                  badList_subsppA_features = c(badList_subsppA_features,feature_subsppA)
                }
                else if (flip == 2) { ## flip coin remove B
                  badList_subsppB_features = c(badList_subsppB_features,feature_subsppB)
                }
              }
              else if (testArea == area1 && testArea == area2) { ## if both are subsumed, keep both
                print("EXACT MATCH, KEEPING BOTH ")
              }
            }
          }
        }
      }
    }




    toReturn = list(subsppApoly_toremove = badList_subsppA_features,
                    subsppBpoly_toremove = badList_subsppB_features,
                    subsppA_intToRemove = overlapsToRemove_subsppA,
                    subsppB_intToRemove = overlapsToRemove_subsppB)
    return(toReturn)
  }
}

## END





# loc_suspect = rbind(processedSpecies_card$loc_suspect,
#                     processedSpecies_cocc$loc_suspect,
#                     processedSpecies_rest$loc_suspect,
#                     processedSpecies_igne$loc_suspect)
# loc_good = rbind(processedSpecies_card$loc_good,
#                  processedSpecies_cocc$loc_good,
#                  processedSpecies_rest$loc_good,
#                  processedSpecies_igne$loc_good)
# loc_good = rbind(processedSpecies_card$loc_good,
#                  processedSpecies_cocc$loc_good,
#                  processedSpecies_rest$loc_good,
#                  processedSpecies_igne$loc_good)

## build png of the points used

printPointsPng = function(species,subspecies,bg,loc_good){
  png(paste("Subspecies_assignment_",species,".png",sep=""))
  #print(length(subspecies))
  #col = floor(sqrt(length(subspecies)))
  #row = ceiling((sqrt(length(subspecies))))
  mf = n2mfrow(length(subspecies))
  if(mf[[1]] > mf[[2]]){mf = rev(mf)}
  par(mfrow=mf)

  for(sub in subspecies){
    print(sub)
    raster::plot(bg, col="grey",colNA="darkgrey",main=paste("assigned",sub,sep=" "),
                 legend=F)
    temp = loc_good[loc_good$assigned==sub,]
    points(temp[temp$assigned==sub,2:3],
           col=as.factor(temp$subspecies),
           pch=as.numeric(as.factor(temp$subspecies)))
    legend("top", legend=as.factor(unique(temp$subspecies)),
           pch=unique(as.numeric(as.factor(temp$subspecies))),
           bty="n",
           col=as.factor(unique(temp$subspecies)))
  }
  dev.off()
}

printPointsPng(species=species,subspecies=subspecies,bg=bg,loc_good=loc_good)

printPointsPdfGood = function(species,subspecies,bg,loc_good){
  pdf(paste("Subspecies_assignment_goodLoci_",species,".pdf",sep=""))
  #print(length(subspecies))
  #col = floor(sqrt(length(subspecies)))
  #row = ceiling((sqrt(length(subspecies))))
  #mf = n2mfrow(length(subspecies))
  #if(mf[[1]] > mf[[2]]){mf = rev(mf)}
  #par(mfrow=mf)

  for(sub in subspecies){
    print(sub)

    raster::plot(bg, col="grey",colNA="darkgrey",main=paste("assigned",sub,sep=" "),
                 legend=F)
    temp = loc_good[loc_good$assigned==sub,]
    points(temp[temp$assigned==sub,2:3],
           col=as.factor(temp$subspecies),
           pch=as.numeric(as.factor(temp$subspecies)))
    legend("top", legend=as.factor(unique(temp$subspecies)),
           pch=unique(as.numeric(as.factor(temp$subspecies))),
           bty="n",
           col=as.factor(unique(temp$subspecies)))
  }
  dev.off()
}

printPointsPdfSuspect = function(species,subspecies,bg,loc_suspect){
  pdf(paste("Subspecies_assignment_suspectLoci_",species,".pdf",sep=""))
  #print(length(subspecies))
  #col = floor(sqrt(length(subspecies)))
  #row = ceiling((sqrt(length(subspecies))))
  #mf = n2mfrow(length(subspecies))
  #if(mf[[1]] > mf[[2]]){mf = rev(mf)}
  #par(mfrow=mf)

  for(sub in subspecies){
    print(sub)
    raster::plot(bg, col="grey",colNA="darkgrey",main=paste("assigned",sub,sep=" "),
                 legend=F)
    temp = loc_suspect[loc_suspect$assigned==sub,]
    points(temp[temp$assigned==sub,2:3],
           col=as.factor(temp$subspecies),
           pch=as.numeric(as.factor(temp$subspecies)))
    legend("top", legend=as.factor(unique(temp$subspecies)),
           pch=unique(as.numeric(as.factor(temp$subspecies))),
           bty="n",
           col=as.factor(unique(temp$subspecies)))
  }
  dev.off()
}

printPointsPdfGood(species=species,subspecies=subspecies,bg=bg,loc_good=loc_good)
printPointsPdfSuspect(species=species,subspecies=subspecies,bg=bg,loc_suspect=loc_suspect)


# png("Phainopepla nitens assignment.png")
# par(mfrow=c(1,2))
# plot(bg, col="grey",colNA="darkgrey",main="assigned lepida")
# points(loc_good[loc_good$assigned=="lepida",2:3],col=as.factor(loc_good$subspecies),
#        pch=as.numeric(as.factor(loc_good$subspecies)))
# legend("top", legend=as.factor(unique(loc_good$subspecies)),
#        pch=unique(as.numeric(as.factor(loc_good$subspecies))),
#        bty="n",
#        col=as.factor(unique(loc_good$subspecies)))
# plot(bg, col="grey",colNA="darkgrey",main="assigned nitens")
# points(loc_good[loc_good$assigned=="nitens",2:3],col=as.factor(loc_good$subspecies),
#        pch=as.numeric(as.factor(loc_good$subspecies)))
# legend("top", legend=as.factor(unique(loc_good$subspecies)),
#        pch=unique(as.numeric(as.factor(loc_good$subspecies))),
#        bty="n",
#        col=as.factor(unique(loc_good$subspecies)))
# dev.off()

# raster::plot(bg, col="grey",colNA="darkgrey")
# points(loc_good[loc_good$subspecies=="nitens",2:3],col=as.factor(loc_good$assigned),
#        pch=as.numeric(as.factor(loc_good$subspecies)),main="nitens prior")
# raster::plot(bg, col="grey",colNA="darkgrey")
# points(loc_good[loc_good$subspecies=="lepida",2:3],col=as.factor(loc_good$assigned),
#        pch=as.numeric(as.factor(loc_good$subspecies)),main="lepida prior")

## step 2 -- remove any locations with no Env data

if THIN_BEYOND == TRUE {

  cleanByEnvironment = function(Env,loc){
    loc[,2] = as.numeric(loc[,2])
    loc[,3] = as.numeric(loc[,3])
    extr = raster::extract(Env, loc[,2:3]) ## gets values from Env that are at loc
    head(extr)
    loc_clean = loc[!is.na(extr[,1]),]
    print(paste("Removed",nrow(loc)-nrow(loc_clean),"rows with no Env data"))

    return(loc_clean)
  }

  loc_good_clean = cleanByEnvironment(Env=Env,loc=loc_good)

  spThinBySubspecies = function(loc_good_clean,species,overwrite=T){
    ## step 3 -- subset data based on subspecies, then run spthin
    ## this only gives lat longs, results in a named list with subspp names as names
    loc_good_by_subspp = split(loc_good_clean[,2:3], loc_good_clean$assigned)

    pathlist <- sapply(names(loc_good_by_subspp),function(x) NULL)

    for(i in 1:length(names(loc_good_by_subspp))){
      sub = names(loc_good_by_subspp)[[i]]
      dir = paste(getwd(),"/Thinned/",sep="")
      if (!(file.exists(dir))){
        dir.create(file.path(dir))
        #setwd(file.path(subDir))
      }
      subDir = paste(getwd(),"/Thinned/",species," ",sub,sep="")
      if (!(file.exists(subDir))){
        dir.create(file.path(subDir))
        #setwd(file.path(subDir))
      }

      loc_sub = loc_good_by_subspp[[i]]
      loc_sub$name = sub

      print(sub)

      ## select a single subspecies and thin it

      if(overwrite==T){
        if ((file.exists(paste(subDir,"/thinned_data_thin1.csv",sep="")))){
          file.remove(paste(subDir,"/thinned_data_thin1.csv",sep=""))
          #setwd(file.path(subDir))
        }
      }

      path = paste(subDir,"/thinned_data_thin1.csv",sep="")
      pathlist[[i]] = path

      thin<-spThin::thin(loc.data = loc_sub,
                         lat.col = "latitude",
                         long.col = "longitude",
                         spec.col = "name",
                         thin.par = 10, ## km distance that records need to be separated by
                         reps = 10, ## number of times to repeat thinning process
                         locs.thinned.list.return = T,
                         write.files = T,
                         max.files = 1,
                         out.dir = subDir,
                         write.log.file = T)
    }

    ## convert a list of the paths to the thinned stuff, named by subspecies
    newThin = data.frame()
    for(path in pathList){
      addThin = read.csv(path)
      newThin = rbind(newThin,addThin)
    }
    return(newThin)

  }

  loc_thin = spThinBySubspecies(loc_good_clean = loc_good_clean,species=species,overwrite = T)

  ## import the thinned data one at a time to run enmeval on
  ## step 4 -- run EMNeval
  ##TODO: bg points, buffer or something
  ##TODO: run all subspp at same time?
  ##TODO: parallel it?

  ## buffer the points!


  ##TODO: get this running on a list rather than doing manually
  ## and do all the pairwise stuff
  backgroundForPCA = function(localities=loc_good[,2:3],
                              r=200000,
                              num=(100*nrow(localities)),
                              e=Env){
    library(ENMTools)
    bg1 = ENMTools::background.points.buffer(localities, radius = r,
                                             n = num, mask = e[[1]])
    extract1 = na.omit(cbind(localities,
                             extract(e, localities), rep(1, nrow(localities))))
    colnames(extract1)[ncol(extract1)] = 'occ'
    extbg1 = na.omit(cbind(bg1, extract(e, bg1), rep(0, nrow(bg1))))
    colnames(extbg1)[ncol(extbg1)] = 'occ'
    dat1 = rbind(extract1, extbg1)
    return(list(bgenv=dat1,bgpoints=bg1,bgext=extract1,bgextbg=extbg1))
  }


  loc_thin_bgstuff = backgroundForPCA(localities=loc_thin[,2:3])
  bg_dat = loc_thin_bgstuff$bgenv
  bg_bg = loc_thin_bgstuff$bgpoints


  backgroundPerSpecies = function(localities=loc_thin){
    loc_thin_by_subspecies = split(loc_thin, loc_thin$name)
    bgenv_by_subspecies = list()

    bgext_by_subspecies = list()
    bgpoints_by_subspecies = list()
    bgextbg_by_subspecies = list()

    for(i in 1:length(names(loc_thin_by_subspecies))){
      singleSubspp = loc_thin_by_subspecies[[i]]
      subsppName = names(loc_thin_by_subspecies)[[i]]
      single_bgstuff = backgroundForPCA(singleSubspp[,2:3])
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

  perspecies_bgstuff = backgroundPerSpecies(localities=loc_thin)

  ## the base ecospat function didn't check for na.rm so I did so
  ecospat.grid.clim.dyn_custom <- function(glob, glob1, sp, R, th.sp = 0, th.env = 0,
                                           geomask = NULL,removeNA=T) {

    glob <- as.matrix(glob)
    glob1 <- as.matrix(glob1)
    sp <- as.matrix(sp)
    l <- list()

    if (ncol(glob) > 2)
      stop("cannot calculate overlap with more than two axes")

    if (ncol(glob) == 1) {
      # if scores in one dimension (e.g. LDA,SDM predictions,...)
      xmax <- max(glob[, 1])
      xmin <- min(glob[, 1])
      x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), length.out = R)  # breaks on score gradient 1
      sp.dens <- density(sp[, 1], kernel = "gaussian", from = xmin, to = xmax,
                         n = R, cut = 0)  # calculate the density of occurrences in a vector of R pixels along the score gradient
      # using a gaussian kernel density function, with R bins.
      glob1.dens <- density(glob1[, 1], kernel = "gaussian", from = xmin,
                            to = xmax, n = R, cut = 0)  # calculate the density of environments in glob1
      z <- sp.dens$y * nrow(sp)/sum(sp.dens$y)  # rescale density to the number of occurrences in sp
      # number of occurrence/pixel
      Z <- glob1.dens$y * nrow(glob)/sum(glob1.dens$y)  # rescale density to the number of sites in glob1
      glob1r <- sapply(glob1, findInterval, glob1.dens$x)
      th.env <- quantile(glob1.dens$y[glob1r], th.env)
      glob1rm <- which(Z < th.env)
      spr <- sapply(sp, findInterval, sp.dens$x)
      th.sp <- quantile(sp.dens$y[spr], th.sp)
      sprm <- which(z < th.sp)
      z[sprm] <- 0  # remove infinitesimally small number generated by kernel density function
      Z[glob1rm] <- 0  # remove infinitesimally small number generated by kernel density function

      z.uncor <- z/max(z)  # rescale between [0:1] for comparison with other species
      z.cor <- z/Z  # correct for environment prevalence
      z.cor[is.na(z.cor)] <- 0  # remove n/0 situations
      z.cor[z.cor == "Inf"] <- 0  # remove 0/0 situations
      z.cor <- z.cor/max(z.cor)  # rescale between [0:1] for comparison with other species
      w <- z.uncor
      w[w > 0] <- 1
      l$x <- x
      l$z <- z
      l$z.uncor <- z.uncor
      l$z.cor <- z.cor
      l$Z <- Z
      l$glob <- glob
      l$glob1 <- glob1
      l$sp <- sp
      l$w <- w
    }

    if (ncol(glob) == 2) {
      # if scores in two dimensions (e.g. PCA)

      xmin <- min(glob[, 1])
      xmax <- max(glob[, 1])
      ymin <- min(glob[, 2])
      ymax <- max(glob[, 2])  # data preparation
      glob1r <- data.frame(cbind((glob1[, 1] - xmin)/abs(xmax - xmin), (glob1[,
                                                                              2] - ymin)/abs(ymax - ymin)))  # data preparation
      spr <- data.frame(cbind((sp[, 1] - xmin)/abs(xmax - xmin), (sp[, 2] -
                                                                    ymin)/abs(ymax - ymin)))  # data preparation
      mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))),
                     nrcol = R-2, count = FALSE) # data preparation
      sp.dens <- kernelUD(SpatialPoints(spr[, 1:2]), h = "href", grid = mask,
                          kern = "bivnorm")  # calculate the density of occurrences in a grid of RxR pixels along the score gradients
      sp.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax, matrix(sp.dens$ud,
                                                                               nrow = R))
      # using a gaussian kernel density function, with RxR bins.
      # sp.dens$var[sp.dens$var>0 & sp.dens$var<1]<-0
      glob1.dens <- kernelUD(SpatialPoints(glob1r[, 1:2]), grid = mask, kern = "bivnorm")
      glob1.dens <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax,
                           matrix(glob1.dens$ud, nrow = R))
      # glob1.dens$var[glob1.dens$var<1 & glob1.dens$var>0]<-0

      x <- seq(from = min(glob[, 1]), to = max(glob[, 1]), length.out = R)  # breaks on score gradient 1
      y <- seq(from = min(glob[, 2]), to = max(glob[, 2]), length.out = R)  # breaks on score gradient 2
      glob1r <- extract(glob1.dens, glob1)
      Z.th <- quantile(glob1r, th.env,na.rm=removeNA)
      glob1.dens[glob1.dens < Z.th] <- 0
      if (!is.null(geomask)) {
        proj4string(geomask) <- NA
        glob1.dens <- mask(glob1.dens, geomask, updatevalue = 0)  # Geographical mask in the case if the analysis takes place in the geographical space
      }
      Z <- glob1.dens * nrow(glob1)/cellStats(glob1.dens, "sum")

      spr <- extract(sp.dens, sp)
      z.th <- quantile(spr, th.sp)
      sp.dens[Z == 0] <- 0
      sp.dens[sp.dens < z.th] <- 0
      if (!is.null(geomask)) {
        sp.dens <- mask(sp.dens, geomask, updatevalue = 0)  # Geographical mask in the case if the analysis takes place in the geographical space
      }
      z <- sp.dens * nrow(sp)/cellStats(sp.dens, "sum")
      z.uncor <- z/cellStats(z, "max")
      w <- z.uncor  # remove infinitesimally small number generated by kernel density function
      w[w > 0] <- 1
      z.cor <- z/Z  # correct for environment prevalence
      z.cor[is.na(z.cor)] <- 0  # remove n/0 situations
      z.cor <- z.cor/cellStats(z.cor, "max")
      l$x <- x
      l$y <- y
      l$z <- z
      l$z.uncor <- z.uncor
      l$z.cor <- z.cor
      l$Z <- Z
      l$glob <- glob
      l$glob1 <- glob1
      l$sp <- sp
      l$w <- w

    }

    return(l)
  }


  createPcaToCompare = function(loc_thin_bgstuff,perspecies_bgstuff,species) {
    bg_dat = loc_thin_bgstuff$bgenv
    bg_bg = loc_thin_bgstuff$bgpoints

    bgext_by_subspecies = perspecies_bgstuff$bgext_by_subspecies
    bgenv_by_subspecies = perspecies_bgstuff$bgenv_by_subspecies

    ## pca bg points
    pca.env <- dudi.pca(bg_dat[,3:(ncol(bg_dat)-1)],scannf=F,nf=2)
    png(paste("PCAcorrelationCircle_",species,".png",sep=""))
    ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
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

      scores_subspp = suprow(pca.env,
                             singleSubspp_bgext[which(
                               singleSubspp_bgext[,ncol(singleSubspp_bgext)]==1)
                               ,3:(ncol(singleSubspp_bgext)-1)])$li # PCA scores for the species 1 distribution


      scores[[i]] = scores_subspp

      scores_clim_subspp = suprow(pca.env,
                                  singleSubspp_bgenv[,3:(ncol(singleSubspp_bgenv)-1)])$li # PCA scores for the whole native study area species 1 ## bgenv

      scores_clim[[i]] = scores_clim_subspp

      ## make a dynamic occurrence densities grid
      ## grid of occ densities along one or two environmental gradients
      ## glob = env variables in background
      ## glob1 = env variables for species
      ## sp = occurrences of species
      ## R = resolution
      ## th.sp = a threshhold to elimite low density values of species occurrences
      grid_clim_subspp <- ecospat.grid.clim.dyn_custom(glob = scores_globclim,
                                                       glob1 = scores_clim_subspp,
                                                       sp = scores_subspp,
                                                       R = 100,
                                                       th.sp = 0,
                                                       th.env = 0,
                                                       removeNA=T)
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

  pcaOutput = createPcaToCompare(loc_thin_bgstuff=loc_thin_bgstuff,
                                 perspecies_bgstuff=perspecies_bgstuff,
                                 species=species)


  pca_grid_clim = pcaOutput$grid_clim

  ## PCA the background points
  # pca.env <- dudi.pca(bg_dat[,3:(ncol(bg_dat)-1)],scannf=F,nf=2)
  # png(paste("PCAcorrelationCircle_",species,".png",sep=""))
  # ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
  # dev.off()
  # scores.sp1 <- suprow(pca.env,
  #                      bgext_lepida[which(bgext_lepida[,22]==1),3:21])$li # PCA scores for the species 1 distribution
  # scores.sp2 <- suprow(pca.env,
  #                      bgext_nitens[which(bgext_nitens[,22]==1),3:21])$li # PCA scores for the species 2 distribution
  # scores.clim1 <- suprow(pca.env,bg_lepida[,3:21])$li # PCA scores for the whole native study area species 1 ## bgenv
  # scores.clim2 <- suprow(pca.env,bg_nitens[,3:21])$li # PCA scores for the whole native study area species 2
  #
  # ## make a dynamic occurrence densities grid
  # ## grid of occ densities along one or two environmental gradients
  # ## glob = env variables in background
  # ## glob1 = env variables for species
  # ## sp = occurrences of species
  # ## R = resolution
  # ## th.sp = a threshhold to elimite low density values of species occurrences
  # grid.clim1 <- ecospat.grid.clim.dyn(
  #   glob = scores.globclim,
  #   glob1 = scores.clim1,
  #   sp = scores.sp1,
  #   R = 100,
  #   th.sp = 0
  # )
  # grid.clim2 <- ecospat.grid.clim.dyn(
  #   glob = scores.globclim,
  #   glob1 = scores.clim2,
  #   sp = scores.sp2,
  #   R = 100,
  #   th.sp = 0
  # )

  ## get overlaps and test for niche equivalence pairwise

  pairwiseNicheOverlap = function(pca_grid_clim=pca_grid_clim){

    overlap_df = data.frame(spp1=character(),
                            spp2=character(),
                            SchoenersD=numeric(),
                            modifiedHellingersI=numeric())

    for(i in 1:length(pca_grid_clim)){
      for(j in 1:length(pca_grid_clim)){
        if(i<j){
          #print(paste(i,j))
          spp1_name = names(pca_grid_clim)[[i]]
          spp1 = pca_grid_clim[[i]]
          spp2_name = names(pca_grid_clim)[[j]]
          spp2 = pca_grid_clim[[j]]
          overlap <- ecospat.niche.overlap(spp1, spp2, cor=T)
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

  overlap_df = pairwiseNicheOverlap(pca_grid_clim=pca_grid_clim)

  # ## calculate overlap between the species niches
  # ## Schoener's D
  # D.overlap <- ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$D
  # D.overlap # 0.3298541 for Toxostoma cris vs leco
  # ## modified Hellinger's I
  # I.overlap = ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$I
  # I.overlap # 0.5620279 for Toxostoma cris vs leco
  # ## this is looking for overlap in PCA-environmental space of the two species


  ## is this a projection of the niche?
  # ## its the pca space occupied by the two species
  # png(paste("NicheSpaceComparison_",species,".png",sep=""))
  # par(mfrow=c(1,2))
  # plot(grid.clim1$w,main="lepida",sub=paste("Sch. D =",round(D.overlap,digits=6)))
  # plot(grid.clim2$w,main="nitens",sub=paste("Hel. I =",round(I.overlap,digits=6)))
  # dev.off()

  ## test for niche equvalence pairwise
  ## first need to modify function again to remove NA
  ecospat.niche.equivalency.test_custom <- function(z1, z2, rep, alternative = "greater", ncores=1) {

    R <- length(z1$x)
    l <- list()

    obs.o <- ecospat.niche.overlap(z1, z2, cor = TRUE)  #observed niche overlap

    if (ncores == 1){
      sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.eq.gen_custom, z1, z2)), byrow = TRUE,
                                    ncol = 2))  #simulate random overlap
    }else{
      #number of cores attributed for the permutation test
      cl <- makeCluster(ncores)  #open a cluster for parallelization
      invisible(clusterEvalQ(cl))  #import the internal function into the cluster
      sim.o <- as.data.frame(matrix(unlist(parLapply(cl, 1:rep, overlap.eq.gen_custom, z1, z2)), byrow = TRUE,
                                    ncol = 2))  #simulate random overlap
      stopCluster(cl)  #shutdown the cluster
    }
    colnames(sim.o) <- c("D", "I")
    l$sim <- sim.o  # storage
    l$obs <- obs.o  # storage

    if (alternative == "greater") {
      l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
      l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = greater -> test for niche conservatism/convergence
    }
    if (alternative == "lower") {
      l$p.D <- (sum(sim.o$D <= obs.o$D) + 1)/(length(sim.o$D) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
      l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)  # storage of p-values alternative hypothesis = lower -> test for niche divergence
    }

    return(l)
  }

  ##################################################################################################

  #### internal functions from ENMtools that are needed for custom ones
  overlap.sim.gen <- function(repi, z1, z2, rand.type = rand.type) {
    R1 <- length(z1$x)
    R2 <- length(z2$x)
    if (is.null(z1$y) & is.null(z2$y)) {
      if (rand.type == 1) {
        # if rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly
        # shifted
        center.z1 <- which(z1$z.uncor == 1)  # define the centroid of the observed niche
        Z1 <- z1$Z/max(z1$Z)
        rand.center.z1 <- sample(1:R1, size = 1, replace = FALSE, prob = Z1)  # randomly (weighted by environment prevalence) define the new centroid for the niche
        xshift.z1 <- rand.center.z1 - center.z1  # shift on x axis
        z1.sim <- z1
        z1.sim$z <- rep(0, R1)  # set intial densities to 0
        for (i in 1:length(z1$x)) {
          i.trans.z1 <- i + xshift.z1
          if (i.trans.z1 > R1 | i.trans.z1 < 0)
            (next)()  # densities falling out of the env space are not considered
          z1.sim$z[i.trans.z1] <- z1$z[i]  # shift of pixels
        }
        z1.sim$z <- (z1$Z != 0) * 1 * z1.sim$z  # remove densities out of existing environments
        z1.sim$z.cor <- (z1.sim$z/z1$Z)/max((z1.sim$z/z1$Z), na.rm = TRUE)  #transform densities into occupancies
        z1.sim$z.cor[which(is.na(z1.sim$z.cor))] <- 0
        z1.sim$z.uncor <- z1.sim$z/max(z1.sim$z, na.rm = TRUE)
        z1.sim$z.uncor[which(is.na(z1.sim$z.uncor))] <- 0
      }

      center.z2 <- which(z2$z.uncor == 1)  # define the centroid of the observed niche
      Z2 <- z2$Z/max(z2$Z)
      rand.center.z2 <- sample(1:R2, size = 1, replace = FALSE, prob = Z2)  # randomly (weighted by environment prevalence) define the new centroid for the niche

      xshift.z2 <- rand.center.z2 - center.z2  # shift on x axis
      z2.sim <- z2
      z2.sim$z <- rep(0, R2)  # set intial densities to 0
      for (i in 1:length(z2$x)) {
        i.trans.z2 <- i + xshift.z2
        if (i.trans.z2 > R2 | i.trans.z2 < 0)
          (next)()  # densities falling out of the env space are not considered
        z2.sim$z[i.trans.z2] <- z2$z[i]  # shift of pixels
      }
      z2.sim$z <- (z2$Z != 0) * 1 * z2.sim$z  # remove densities out of existing environments
      z2.sim$z.cor <- (z2.sim$z/z2$Z)/max((z2.sim$z/z2$Z), na.rm = TRUE)  #transform densities into occupancies
      z2.sim$z.cor[which(is.na(z2.sim$z.cor))] <- 0
      z2.sim$z.uncor <- z2.sim$z/max(z2.sim$z, na.rm = TRUE)
      z2.sim$z.uncor[which(is.na(z2.sim$z.uncor))] <- 0
    }

    if (!is.null(z2$y) & !is.null(z1$y)) {
      if (rand.type == 1) {
        # if rand.type = 1, both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly
        # shifted
        centroid.z1 <- which(z1$z.uncor == 1, arr.ind = TRUE)[1, ]  # define the centroid of the observed niche
        Z1 <- z1$Z/max(z1$Z)
        rand.centroids.z1 <- which(Z1 > 0, arr.ind = TRUE)  # all pixels with existing environments in the study area
        weight.z1 <- Z1[Z1 > 0]
        rand.centroid.z1 <- rand.centroids.z1[sample(1:nrow(rand.centroids.z1), size = 1, replace = FALSE,
                                                     prob = weight.z1), ]  # randomly (weighted by environment prevalence) define the new centroid for the niche
        xshift.z1 <- rand.centroid.z1[1] - centroid.z1[1]  # shift on x axis
        yshift.z1 <- rand.centroid.z1[2] - centroid.z1[2]  # shift on y axis
        z1.sim <- z1
        z1.sim$z <- matrix(rep(0, R1 * R1), ncol = R1, nrow = R1)  # set intial densities to 0
        for (i in 1:R1) {
          for (j in 1:R1) {
            i.trans.z1 <- i + xshift.z1
            j.trans.z1 <- j + yshift.z1
            if (i.trans.z1 > R1 | i.trans.z1 < 0)
              (next)()  # densities falling out of the env space are not considered
            if (j.trans.z1 > R1 | j.trans.z1 < 0)
              (next)()
            z1.sim$z[i.trans.z1, j.trans.z1] <- z1$z[i, j]  # shift of pixels
          }
        }
        z1.sim$z <- (z1$Z != 0) * 1 * z1.sim$z  # remove densities out of existing environments
        z1.sim$z.cor <- (z1.sim$z/z1$Z)/max((z1.sim$z/z1$Z), na.rm = TRUE)  #transform densities into occupancies
        z1.sim$z.cor[which(is.na(z1.sim$z.cor))] <- 0
        z1.sim$z.uncor <- z1.sim$z/max(z1.sim$z, na.rm = TRUE)
        z1.sim$z.uncor[which(is.na(z1.sim$z.uncor))] <- 0
      }
      centroid.z2 <- which(z2$z.uncor == 1, arr.ind = TRUE)[1, ]  # define the centroid of the observed niche
      Z2 <- z2$Z/max(z2$Z)
      rand.centroids.z2 <- which(Z2 > 0, arr.ind = TRUE)  # all pixels with existing environments in the study area
      weight.z2 <- Z2[Z2 > 0]
      rand.centroid.z2 <- rand.centroids.z2[sample(1:nrow(rand.centroids.z2), size = 1, replace = FALSE,
                                                   prob = weight.z2), ]  # randomly (weighted by environment prevalence) define the new centroid for the niche
      xshift.z2 <- rand.centroid.z2[1] - centroid.z2[1]  # shift on x axis
      yshift.z2 <- rand.centroid.z2[2] - centroid.z2[2]  # shift on y axis
      z2.sim <- z2
      z2.sim$z <- matrix(rep(0, R2 * R2), ncol = R2, nrow = R2)  # set intial densities to 0
      for (i in 1:R2) {
        for (j in 1:R2) {
          i.trans.z2 <- i + xshift.z2
          j.trans.z2 <- j + yshift.z2
          if (i.trans.z2 > R2 | i.trans.z2 < 0)
            (next)()  # densities falling out of the env space are not considered
          if (j.trans.z2 > R2 | j.trans.z2 < 0)
            (next)()
          z2.sim$z[i.trans.z2, j.trans.z2] <- z2$z[i, j]  # shift of pixels
        }
      }
      z2.sim$z <- (z2$Z != 0) * 1 * z2.sim$z  # remove densities out of existing environments
      z2.sim$z.cor <- (z2.sim$z/z2$Z)/max((z2.sim$z/z2$Z), na.rm = TRUE)  #transform densities into occupancies
      z2.sim$z.cor[which(is.na(z2.sim$z.cor))] <- 0
      z2.sim$z.uncor <- z2.sim$z/max(z2.sim$z, na.rm = TRUE)
      z2.sim$z.uncor[which(is.na(z2.sim$z.uncor))] <- 0
    }

    if (rand.type == 1) {
      o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)
    }
    if (rand.type == 2)
    {
      o.i <- ecospat.niche.overlap(z1, z2.sim, cor = TRUE)
    }  # overlap between random and observed niches
    sim.o.D <- o.i$D  # storage of overlaps
    sim.o.I <- o.i$I
    return(c(sim.o.D, sim.o.I))
  }
  overlap.eq.gen_custom <- function(repi, z1, z2) {
    if (is.null(z1$y)) {
      # overlap on one axis

      occ.pool <- c(z1$sp, z2$sp)  # pool of random occurrences
      rand.row <- sample(1:length(occ.pool), length(z1$sp))  # random reallocation of occurrences to datasets
      sp1.sim <- occ.pool[rand.row]
      sp2.sim <- occ.pool[-rand.row]
    }

    if (!is.null(z1$y)) {
      # overlap on two axes

      occ.pool <- rbind(z1$sp, z2$sp)  # pool of random occurrences
      row.names(occ.pool)<-c()  # remove the row names
      rand.row <- sample(1:nrow(occ.pool), nrow(z1$sp))  # random reallocation of occurrences to datasets
      sp1.sim <- occ.pool[rand.row, ]
      sp2.sim <- occ.pool[-rand.row, ]
    }

    z1.sim <- ecospat.grid.clim.dyn_custom(z1$glob, z1$glob1, data.frame(sp1.sim), R = length(z1$x))  # gridding
    z2.sim <- ecospat.grid.clim.dyn_custom(z2$glob, z2$glob1, data.frame(sp2.sim), R = length(z2$x))

    o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)  # overlap between random and observed niches
    sim.o.D <- o.i$D  # storage of overlaps
    sim.o.I <- o.i$I
    return(c(sim.o.D, sim.o.I))
  }



  pairwiseNicheEquivalence = function(pca_grid_clim=pca_grid_clim,rep1=10,rep2=1000){
    for(i in 1:length(pca_grid_clim)){
      for(j in 1:length(pca_grid_clim)){
        if(i<j){
          print(paste(i,j))
          spp1_name = names(pca_grid_clim)[[i]]
          spp1 = pca_grid_clim[[i]]
          spp2_name = names(pca_grid_clim)[[j]]
          spp2 = pca_grid_clim[[j]]

          eq.test <- ecospat.niche.equivalency.test_custom(z1=spp1, z2=spp2,
                                                           rep=rep1, alternative = "greater")
          sim.test <- ecospat.niche.similarity.test(z1=spp1, z2=spp2,
                                                    rep=rep2, alternative = "greater",
                                                    rand.type=2)
          sim.test2 <- ecospat.niche.similarity.test(z1=spp1, z2=spp2,
                                                     rep=rep2, alternative = "greater",
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

  pairwiseNicheEquivalence(pca_grid_clim=pca_grid_clim,rep1=10,rep2=1000)

  #
  # ## test for niche equivalence
  # ## this is the observed niche overlap when you randomize species id
  # ## greater means you test for niche conservatism
  # ## lower means you test for niche divergence
  # eq.test <- ecospat.niche.equivalency.test(grid.clim1, grid.clim2,
  #                                           rep=10, alternative = "greater") ##rep = 1000 recommended for operational runs
  # ## for 10 takes 1 minute
  #
  # ## then test for niche similarity --
  # ## overlap between spp 1 and overlaps between random spp 2 bg niches
  # ## rand.type 2 means only z2 is randomly shifted
  # sim.test <- ecospat.niche.similarity.test(grid.clim1, grid.clim2,
  #                                           rep=1000, alternative = "greater",
  #                                           rand.type=2)
  # sim.test2 <- ecospat.niche.similarity.test(grid.clim2,grid.clim1,
  #                                            rep=1000, alternative = "greater",
  #                                            rand.type=2)
  # png(paste("EquivalencyOverlapTests_",species,"_nitens vs lepida",".png",sep=""))
  # par(mfrow=c(2,3))
  # ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
  # ecospat.plot.overlap.test(sim.test, "D", "Similarity 1->2")
  # ecospat.plot.overlap.test(sim.test2, "D", "Similarity 2->1")
  # ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
  # ecospat.plot.overlap.test(sim.test, "I", "Similarity 1->2")
  # ecospat.plot.overlap.test(sim.test2, "I", "Similarity 2->1")
  # dev.off()


  ## generate the actual models

  library(ENMeval)

  listENMresults = lapply(1:length(nitens_by_subspp),function(i){
    ##TODO: add in bg points from above
    subspp = names(nitens_by_subspp)[[i]]
    print(paste("Running",subspp))
    res = ENMevaluate(occ=nitens_by_subspp[[i]], env = Env, method='block',
                      parallel=T, numCores=4, fc=c("L", "LQ", "H"), #c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")
                      RMvalues=seq(0.5,4,0.5), rasterPreds=F)
    #names(res) = names(nitens_by_subspp)[[i]]
    return(res)
  })
  names(listENMresults) = names(nitens_by_subspp)

  # Mean.AUC is the one to use,
  ## full.aUC it's the AUC of the model used for AICc comparisons
  ## using all the points and no withheld data
  ## for omission rate, if small sample size and high confidence
  ## do mean.ORmin. if don't know, do mean.OR10

  ## now find the best models
  ## minimize omission rate and optimize AUC values
  ## remember to change ORmin to OR10 depending on sample size

  minORmaxAUCmodel = function(res,or="Mean.OR10",auc="Mean.AUC"){
    setsort = res@results[order(res@results[,or]),]
    setsort2 = setsort[order(setsort[,auc], decreasing=TRUE),]
    top = setsort2[1,]
    ## extract out and put in maxent
    best = which(as.character(res@results[,1]) == as.character(setsort2[1,1]))
    return(best)
  }

  bestModels = lapply(1:length(listENMresults),function(i){
    res = listENMresults[[i]]
    name = names(listENMresults)[[i]]
    print(paste("Choosing best model for",name))
    best = minORmaxAUCmodel(res=listENMresults[[i]])
    return(best)
  })
  names(bestModels) = names(listENMresults)


  # getBestRaster = function(best,res,wdToPrint=paste(getwd(),"/",sep=""),prefix,suffix,species,Env){
  #   setwd(wdToPrint)
  #   pred.raw = predict(Env, res@models[[best]]) ## need to run a new model in maxent with those settings
  #   ## need to extract the features and the reg multiplier for the model and put in the arguments into the maxent call
  #   writeRaster(pred.raw,filename=paste(wdToPrint,prefix,"_BestModel_",species," ",suffix,".asc",sep=""),
  #               format="ascii",overwrite=T)
  #   return(pred.raw)
  # }

  #####

  getBestRaster = function(best,res,wdToPrint=paste(getwd(),"/",sep=""),prefix,suffix,species,Env,locs){
    setwd(wdToPrint)
    mod.table<-res@results
    no.zero.param <- mod.table[mod.table$nparam != 0,]
    ordered<-no.zero.param[with(no.zero.param, order(delta.AICc)), ]
    opt.mod<-ordered[1,]

    ## beta multiplier
    b.m<-opt.mod$rm
    beta.mulr<- paste('betamultiplier=',b.m,sep='')
    ## java maxent needs everything to be set false
    ## anything you want then needs to be true
    false.args<-c('noautofeature','noproduct','nothreshold','noquadratic','nohinge','nolinear')
    ## then you set anything as true
    feat<-strsplit(as.vector(opt.mod[,2]), ',')[[1]]
    if(feat == 'L'){
      feats = 'linear'
    } else if(feat == 'LQ'){
      feats = c('quadratic', 'linear')
    } else if(feat == 'H'){
      feats = 'hinge'
    } else if(feat == 'P'){
      feats = 'product'
    } else if(feat == 'T'){
      feats = 'threshold'
    } else if(feat == 'LQH'){
      feats = c('linear', 'quadratic', 'hinge')
    } else if(feat == 'LQHP'){
      feats = c('linear', 'quadratic', 'hinge', 'product')
    } else if(feat == 'LQHPT'){
      feats = c('linear', 'quadratic', 'hinge', 'product', 'threshold')
    }
    for (j in 1:length(feats)){false.args[which(sub('no','',false.args)==feats[j])] = feats[j]}
    m <-maxent(Env, locs, args=c(false.args, beta.mulr, 'noremoveduplicates', 'noaddsamplestobackground'))
    pred.raw  <- predict(object= m, x=Env, na.rm=TRUE, format='GTiff',overwrite=TRUE, progress='text',args='logistic')
    ## need to extract the features and the reg multiplier for the model and put in the arguments into the maxent call
    writeRaster(pred.raw,filename=paste(wdToPrint,prefix,"_BestModel_",species," ",suffix,".asc",sep=""),
                format="ascii",overwrite=T)
  }



  bestRasters = lapply(1:length(bestModels),function(i){
    locs = nitens_by_subspp[[i]]
    name = names(listENMresults)[[i]]
    print(paste("Making rasters for",name))
    test = getBestRaster(best=bestModels[[i]],res=listENMresults[[i]],prefix="PredRaw",suffix=name,species=species,Env=Env,locs=locs)
    plot(test, col=viridis::viridis(99))
    return(test)
  })
  names(bestRasters) = names(listENMresults)

  png(paste("BestRastersTest_",species," ","nitens_vs_lepida.png",sep=""))
  par(mfrow=c(1,2))
  plot(bestRasters[[1]],col=viridis::viridis(99),main=names(bestRasters)[[1]])
  plot(bestRasters[[2]],col=viridis::viridis(99),main=names(bestRasters)[[2]])
  dev.off()

  ## now do thresholding
  ## get the enmeval output

  threshRasters = lapply(1:length(bestRasters),function(i){
    pred = bestRasters[[i]]
    name = names(bestRasters)[[i]]
    best = bestModels[[i]]
    ev.set <- evaluate(nitens_by_subspp[[i]], listENMresults[[i]]@bg.pts, listENMresults[[i]]@models[[best]], Env)
    th1 = threshold(ev.set) ## omission options
    p1.kappa = pred >= th1$kappa ## equal sensitivity and specificity according to ROC curve
    p1.spsen = pred >= th1$spec_sens
    p1.nomit = pred >= th1$no_omission
    p1.preva = pred >= th1$prevalence
    p1.equal = pred >= th1$equal_sens_spec
    p1.sns09 = pred >= th1$sensitivity

    png(paste("Thresholds_",species," ",name,".png",sep=""))
    par(mfrow=c(2,3))
    plot(p1.kappa, col=viridis::viridis(99),main=paste(name,"kappa"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    plot(p1.spsen, col=viridis::viridis(99),main=paste(name,"spec_sens"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    plot(p1.nomit, col=viridis::viridis(99),main=paste(name,"no_omisison"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    plot(p1.preva, col=viridis::viridis(99),main=paste(name,"prevalence"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    plot(p1.equal, col=viridis::viridis(99),main=paste(name,"equal_sens_spec"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    plot(p1.sns09, col=viridis::viridis(99),main=paste(name,"sensitivity0.9"))#; points(nitens_by_subspp[[i]],col=rgb(0,0,0,0.01))
    dev.off()

    writeRaster(p1.kappa,filename=paste("Threshold_",species," ",name,"_kappa.asc",sep=""),format="ascii",overwrite=T)
    writeRaster(p1.spsen,filename=paste("Threshold_",species," ",name,"_spec_sens.asc",sep=""),format="ascii",overwrite=T)
    writeRaster(p1.nomit,filename=paste("Threshold_",species," ",name,"_no_omisison.asc",sep=""),format="ascii",overwrite=T)
    writeRaster(p1.preva,filename=paste("Threshold_",species," ",name,"_prevalence.asc",sep=""),format="ascii",overwrite=T)
    writeRaster(p1.equal,filename=paste("Threshold_",species," ",name,"_equal_sens_spec.asc",sep=""),format="ascii",overwrite=T)
    writeRaster(p1.sns09,filename=paste("Threshold_",species," ",name,"_sensitivity0.9.asc",sep=""),format="ascii",overwrite=T)

    th1ras = list(kappa=p1.kappa,
                  specSens = p1.spsen,
                  noOmit = p1.nomit,
                  prevalence = p1.preva,
                  equalSensSpec = p1.equal,
                  sensitivity_0.9 = p1.sns09)
    return(th1ras)

  })
  names(threshRasters) = names(bestRasters)

  print("Are these really subspecies, folks?")




  ### plotting stuff
  plot(bg, col="grey",colNA="darkgrey")
  points(loc_good$longitude,loc_good$latitude,col=as.factor(loc_good$assigned))

  plot(bg, col="grey",colNA="darkgrey")
  points(loc_suspect$longitude,loc_suspect$latitude,col=as.factor(loc_suspect$assigned),
         pch="x")

}
