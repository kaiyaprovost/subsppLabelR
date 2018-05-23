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

install_github('kaiyaprovost/subsppLabelR')
library(subsppLabelR)
setwd("~/Documents/Classes/Spatial Bioinformatics/project/")

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
  path='/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  pattern="\\.bil$",
  full.names=T))
ext = raster::extent(c(-125,-60,10,50)) ## make sure this will play nice with your points
Env = raster::crop(Env, ext)
bg = Env[[1]] ## just for plotting 

## get locs 
species = "Phainopepla nitens"
subspecies = c("nitens","lepida")

phainopeplaNitens = databaseToAssignedSubspecies(spp=species,
                                                 subsppList=subspecies,
                                                 pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
                                                 quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg,
                                                 outputDir="~/Documents/Classes/Spatial Bioinformatics/project/")
nitens_loc_suspect = phainopeplaNitens$loc_suspect
nitens_loc_good = phainopeplaNitens$loc_good
nitens_pol = phainopeplaNitens$pol

png("Phainopepla nitens assignment.png")
par(mfrow=c(1,2))
plot(bg, col="grey",colNA="darkgrey",main="assigned lepida")
points(nitens_loc_good[nitens_loc_good$assigned=="lepida",2:3],col=as.factor(nitens_loc_good$subspecies),
       pch=as.numeric(as.factor(nitens_loc_good$subspecies)))
legend("top", legend=as.factor(unique(nitens_loc_good$subspecies)),
       pch=unique(as.numeric(as.factor(nitens_loc_good$subspecies))),
       bty="n", 
       col=as.factor(unique(nitens_loc_good$subspecies)))
plot(bg, col="grey",colNA="darkgrey",main="assigned nitens")
points(nitens_loc_good[nitens_loc_good$assigned=="nitens",2:3],col=as.factor(nitens_loc_good$subspecies),
       pch=as.numeric(as.factor(nitens_loc_good$subspecies)))
legend("top", legend=as.factor(unique(nitens_loc_good$subspecies)),
       pch=unique(as.numeric(as.factor(nitens_loc_good$subspecies))),
       bty="n", 
       col=as.factor(unique(nitens_loc_good$subspecies)))
dev.off()

plot(bg, col="grey",colNA="darkgrey")
points(nitens_loc_good[nitens_loc_good$subspecies=="nitens",2:3],col=as.factor(nitens_loc_good$assigned),
       pch=as.numeric(as.factor(nitens_loc_good$subspecies)),main="nitens prior")
plot(bg, col="grey",colNA="darkgrey")
points(nitens_loc_good[nitens_loc_good$subspecies=="lepida",2:3],col=as.factor(nitens_loc_good$assigned),
       pch=as.numeric(as.factor(nitens_loc_good$subspecies)),main="lepida prior")

## step 2 -- remove any locations with no Env data
extr = raster::extract(Env, nitens_loc_good[,2:3]) ## gets values from Env that are at loc 
head(extr)
nitens_loc_good = nitens_loc_good[!is.na(extr[,1]),]

## step 3 -- subset data based on subspecies, then run spthin
## this only gives lat longs, results in a named list with subspp names as names
nitens_by_subspp = split(nitens_loc_good[,2:3], nitens_loc_good$assigned) 

## TODO: RUN SPTHIN!!!!!!!!!!!

## step 4 -- run EMNeval
##TODO: bg points, buffer or something
##TODO: run all subspp at same time?
##TODO: parallel it?

## buffer the points!
library(ENMTools)

##TODO: get this running on a list rather than doing manually 
## and do all the pairwise stuff
backgroundForPCA = function(localities=nitens_loc_good[,2:3],r=200000,
                            num=(100*nrow(localities)),e=Env){
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

  bgstuff = backgroundForPCA(localities=nitens_loc_good[,2:3])
  bg_dat = bgstuff$bgenv
  bg_bg = bgstuff$bgpoints
  
  ## get bg_dat for each species
  bg_lepida = backgroundForPCA(nitens_by_subspp$lepida)$bgenv
  bg_nitens = backgroundForPCA(nitens_by_subspp$nitens)$bgenv
  bgext_lepida = backgroundForPCA(nitens_by_subspp$lepida)$bgext
  bgext_nitens = backgroundForPCA(nitens_by_subspp$lepida)$bgext
  
## PCA the background points 
pca.env <- dudi.pca(bg_dat[,3:21],scannf=F,nf=2)
png(paste("PCAcorrelationCircle_",species,".png",sep=""))
dev.off()
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
scores.globclim<-pca.env$li # PCA scores for the whole study area (all points)
scores.sp1 <- suprow(pca.env,
                     bgext_lepida[which(bgext_lepida[,22]==1),3:21])$li # PCA scores for the species 1 distribution
scores.sp2 <- suprow(pca.env,
                     bgext_nitens[which(bgext_nitens[,22]==1),3:21])$li # PCA scores for the species 2 distribution
scores.clim1 <- suprow(pca.env,bg_lepida[,3:21])$li # PCA scores for the whole native study area species 1
scores.clim2 <- suprow(pca.env,bg_nitens[,3:21])$li # PCA scores for the whole native study area species 2

## make a dynamic occurrence densities grid
## grid of occ densities along one or two environmental gradients
## glob = env variables in background 
## glob1 = env variables for species
## sp = occurrences of species
## R = resolution
## th.sp = a threshhold to elimite low density values of species occurrences
grid.clim1 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim1,
  sp = scores.sp1,
  R = 100,
  th.sp = 0
)
grid.clim2 <- ecospat.grid.clim.dyn(
  glob = scores.globclim,
  glob1 = scores.clim2,
  sp = scores.sp2,
  R = 100,
  th.sp = 0
)

## calculate overlap between the species niches
## Schoener's D
D.overlap <- ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$D 
D.overlap # 0.3298541 for Toxostoma cris vs leco
## modified Hellinger's I
I.overlap = ecospat.niche.overlap (grid.clim1, grid.clim2, cor=T)$I
I.overlap # 0.5620279 for Toxostoma cris vs leco
## this is looking for overlap in PCA-environmental space of the two species


## is this a projection of the niche?
## its the pca space occupied by the two species
png(paste("NicheSpaceComparison_",species,".png",sep=""))
par(mfrow=c(1,2))
plot(grid.clim1$w,main="lepida",sub=paste("Sch. D =",round(D.overlap,digits=6)))
plot(grid.clim2$w,main="nitens",sub=paste("Hel. I =",round(I.overlap,digits=6)))
dev.off()

## test for niche equivalence
## this is the observed niche overlap when you randomize species id
## greater means you test for niche conservatism
## lower means you test for niche divergence 
eq.test <- ecospat.niche.equivalency.test(grid.clim1, grid.clim2,
                                          rep=10, alternative = "greater") ##rep = 1000 recommended for operational runs
## for 10 takes 1 minute

## then test for niche similarity -- 
## overlap between spp 1 and overlaps between random spp 2 bg niches
## rand.type 2 means only z2 is randomly shifted
sim.test <- ecospat.niche.similarity.test(grid.clim1, grid.clim2,
                                          rep=1000, alternative = "greater",
                                          rand.type=2) 
sim.test2 <- ecospat.niche.similarity.test(grid.clim2,grid.clim1,
                                           rep=1000, alternative = "greater",
                                           rand.type=2) 
png(paste("EquivalencyOverlapTests_",species,"_nitens vs lepida",".png",sep=""))
par(mfrow=c(2,3))
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
ecospat.plot.overlap.test(sim.test, "D", "Similarity 1->2")
ecospat.plot.overlap.test(sim.test2, "D", "Similarity 2->1")
ecospat.plot.overlap.test(eq.test, "I", "Equivalency")
ecospat.plot.overlap.test(sim.test, "I", "Similarity 1->2")
ecospat.plot.overlap.test(sim.test2, "I", "Similarity 2->1")
dev.off()


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
points(nitens_loc_good$longitude,nitens_loc_good$latitude,col=as.factor(nitens_loc_good$assigned))

plot(bg, col="grey",colNA="darkgrey")
points(nitens_loc_suspect$longitude,nitens_loc_suspect$latitude,col=as.factor(nitens_loc_suspect$assigned),
       pch="x")

