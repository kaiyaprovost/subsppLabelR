twoSubspeciesOccGbif = function(spp="Phainopepla nitens",subspp1="nitens",subspp2="lepida",
                             pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet")) {
  ## this function uses spocc to query for one species and two subspecies
  ## TODO: turn this into a function to give multiple subspecies and return it
  library(spocc)
  
  print(paste("Getting Species: ",spp))
  sppOcc = occ(query=spp,limit=pointLimit,has_coords=T,from=dbToQuery)
  print(paste("     Getting Subspecies: ",subspp1))
  subSppOcc1 = occ(query=paste(spp,subspp1,sep=" "),
                   limit=pointLimit,has_coords=T,from=dbToQuery)
  print(paste("     Getting Subspecies: ",subspp2))
  subSppOcc2 = occ(query=paste(spp,subspp2,sep=" "),
                   limit=pointLimit,has_coords=T,from=dbToQuery)

  print(sppOcc)
  print(subSppOcc1)
  print(subSppOcc2)
  
  toReturn = list(sppOcc,subSppOcc1,subSppOcc2)
  names(toReturn) = c("unknown",subspp1,subspp2)
  
  return(toReturn)

  }

occ2df_subspeciesLabels = function(subsppOccList_object,subsppOccList_name){
  ## thus function turns an occ object into a dataframe with a column for subspecies
  ## TODO: make it optional to do the "unique" thing for future processing
  
  sppDf = data.frame(occ2df(subsppOccList_object))
  sppLoc = unique(na.omit(sppDf[,1:3]))
  sppLocLab = sppLoc
  sppLocLab$subspecies = subsppOccList_name
  
  return(sppLocLab)
}

labelSubspecies = function(subsppOccList) {
  ## this function takes a list of three taxa and labels them with subspecies information
  ## TODO: turn this into a function to give multiple subspecies and return it
  
  sppOcc = subsppOccList[[1]]
  name_sppOcc = names(subsppOccList)[1]
  subsppOcc1 = subsppOccList[[2]]
  name_subsppOcc1 = names(subsppOccList)[2]
  subsppOcc2 = subsppOccList[[3]]
  name_subsppOcc2 = names(subsppOccList)[3]
  
  #print("Giving occ2df labels")
  sppLocLab = occ2df_subspeciesLabels(subsppOccList_object=sppOcc,subsppOccList_name=name_sppOcc)
  #head(sppLocLab)
  subspp1LocLab = occ2df_subspeciesLabels(subsppOccList_object=subsppOcc1,subsppOccList_name=name_subsppOcc1)
  #head(subspp1LocLab)
  subspp2LocLab = occ2df_subspeciesLabels(subsppOccList_object=subsppOcc2,subsppOccList_name=name_subsppOcc2)
  #head(subspp2LocLab)

  plot(sppLocLab$longitude,sppLocLab$latitude)
  plot(subspp1LocLab$longitude,subspp1LocLab$latitude)
  plot(subspp2LocLab$longitude,subspp2LocLab$latitude)
  
    
  fullPointLocLab = unique(rbind(sppLocLab,subspp1LocLab,subspp2LocLab))
  toReturn = list(fullPointLocLab=fullPointLocLab,
       sppOcc=sppOcc,name_sppOcc=name_sppOcc,
       subsppOcc1=subsppOcc1,name_subsppOcc1=name_subsppOcc1,
       subsppOcc2=subsppOcc2,name_subsppOcc2=name_subsppOcc2)
  return(toReturn)
  
}

subspeciesDensityMap = function(localities,quantile=0.95,xmin=-125,
                                xmax=-60,ymin=10,ymax=50) {
  ## this function uses kernel density to make a raster that will then be used to filter
  ## the data to remove all but the 5% (or 1-quantile) most dense cells
  
  library(MASS)
  library(raster)
  range = c( xmin, xmax, ymin, ymax) 
  ext = raster::extent(range)
  
  w1 = matrix(1,3,3)
  ## generate the two dimensional kernel density estimation
  density = kde2d(localities$longitude, localities$latitude, lims=range)
  ## convert to raster
  densRas = raster(density)
  ## take the top percentile of the points, only the densest areas
  quan = quantile(densRas[densRas],quantile)
  densRas_trim = densRas
  densRas_trim[densRas_trim <= quan] = NA
  plot(densRas_trim,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  
  return(densRas_trim)
  
}

densityMapToPolygons = function(densityMap) {
  ## this function converts density maps to polygons
  ## will work on other kinds of polygons as well
  polygon = densityMap
  polygon[!is.na(polygon)] = 1
  polygon <- disaggregate(rasterToPolygons(polygon, fun=NULL, na.rm=T,dissolve=T))
  plot(polygon)
  return(polygon)
}

whichPolygonsCompletelyOverlap = function(subsppPoly1,subsppPoly2){
  ## this function checks for overlaps between polygons
  ## it only returns a list of polygons that are completely within other polygons
  ## TODO: what about things that are in neither polygon?
  ## TODO: what about things that are in both?
  
  ## TODO: if overlapping, remove from larger polygon size! 
  library(rgeos)
  library(raster)
  badList_subspp1_features = c()
  badList_subspp2_features = c()
  overlapsToRemove_subspp1 = c()
  overlapsToRemove_subspp2 = c()
  
  for (feature_subspp1 in range(1,length(subsppPoly1))){
    for(feature_subspp2 in range(1,length(subsppPoly2))) {
      #intersect = gIntersection(subsppPoly1[feature_subspp1,],subsppPoly2[feature_subspp2,],byid=T)
      intersect = raster::intersect(subsppPoly1[feature_subspp1,],subsppPoly2[feature_subspp2,])
      #print(length(intersect))
      if (length(intersect) != 0){ 
        
        ## if they overlap
        
        #plot(subsppPoly1[feature_subspp1,],border="red",add=T)
        #plot(subsppPoly2[feature_subspp2,],border="cyan",add=T)
        areaInt = gArea(intersect)
        
        if(areaInt > 0) {
          
          ## if they overlap
          ## check areas 
          area1 = gArea(subsppPoly1[feature_subspp1,])
          area2 = gArea(subsppPoly2[feature_subspp2,])
          
          area1percent = areaInt/area1
          area2percent = areaInt/area1 
          
          if (areaInt >= area1 || areaInt >= area2) {
            if (areaInt >= area1) {
              ## if the overlap entirely subsumes an area 
              ## remove the area 
              ## TODO: change to check for density?
              
              #print("remove area 2")
              badList_subspp1_features = c(badList_subspp1_features,feature_subspp1)
            } 
            if (areaInt >= area2) {
              ## if the overlap entirely subsumes an area 
              ## remove the area 
              ## TODO: change to check for density?
              
              badList_subspp2_features = c(badList_subspp2_features,feature_subspp2)
            }
          }
          else {
            if (area1percent <= area2percent) {
              #print("remove area1")
              overlapsToRemove_subspp1 = c(overlapsToRemove_subspp1,intersect)
            }
            else if (area1percent >= area2percent) {
              #print("remove area1")
              overlapsToRemove_subspp2 = c(overlapsToRemove_subspp2,intersect)
            }
          }
        }
        
        
      }
    }
  }
  toReturn = list(subspp1poly_toremove = badList_subspp1_features,
                  subspp2poly_toremove = badList_subspp2_features,
                  subspp1_intToRemove = overlapsToRemove_subspp1,
                  subspp2_intToRemove = overlapsToRemove_subspp2)
  return(toReturn)
}

trimPolygonsByOverlap = function(polygon,idList=NULL,intList=NULL) {
  ## it then removes polygons that are completely within another polygon as ID'd by whichPolygonsCompletelyOverlap()
  ## TODO: what about things that are in neither polygon?
  ## TODO: what about things that are in both?
  polytrim = polygon
  if(is.null(idList)){
    if(is.null(intList)){
      print("no overlaps to remove")
    }
    else {
      for(int in 1:length(intList)){
        polytrim = gDifference(polytrim,intList[[int]])
      }
      print("removing differences from larger shapefile")
    }
  }
  else {
    polytrim = polygon[-idList,]
    if(is.null(intList)){
      print("removing subsumed polygons")
    }
    else {
      for(int in 1:length(intList)){
        polytrim = gDifference(polytrim,intList[[int]])
      }
      print("removing subsumed polygons and differences from larger shape")
    }
    
  }
  return(polytrim)
}

locatePolygonPoints = function(test_points,polygonA,polygonB,crs="+proj=longlat +ellps=WGS84",setcoord=T,
                               nameA="A",nameB="B") {
  ## this function determines whether points overlap with polygon A, polygon B, neither, or both
  ## and then adds columns to label them
  library(dplyr)
  pts = test_points
  if(setcoord==T){
    coordinates(pts) = ~longitude+latitude
    proj4string(pts) = CRS(crs)
    polygonA = spTransform(polygonA,crs)
    polygonB = spTransform(polygonB,crs)
  }
  
  inpolygonA = test_points[which(!is.na(over(pts,polygonA))),]
  inpolygonB = test_points[which(!is.na(over(pts,polygonB))),]
  notInpolygonA = test_points[which(is.na(over(pts,polygonA))),]
  notInpolygonB = test_points[which(is.na(over(pts,polygonB))),]
  
  inBothPolygons = intersect(inpolygonA,inpolygonB)
  inNeitherPolygon = intersect(notInpolygonA,notInpolygonB)
  onlypolygonA = intersect(inpolygonA,notInpolygonB)
  onlypolygonB = intersect(inpolygonB,notInpolygonA)
  
  inBothPolygons_1 = cbind(inBothPolygons,
                           rep(1,length(inBothPolygons[,1])),
                           rep(1,length(inBothPolygons[,1])),
                           rep("both",length(inBothPolygons[,1])))
  inNeitherPolygon_1 = cbind(inNeitherPolygon,
                           rep(0,length(inNeitherPolygon[,1])),
                           rep(0,length(inNeitherPolygon[,1])),
                           rep("neither",length(inNeitherPolygon[,1])))  
  onlypolygonA_1 = cbind(onlypolygonA,
                             rep(1,length(onlypolygonA[,1])),
                             rep(0,length(onlypolygonA[,1])),
                             rep(nameA,length(onlypolygonA[,1])))  
  onlypolygonB_1 = cbind(onlypolygonB,
                         rep(0,length(onlypolygonB[,1])),
                         rep(1,length(onlypolygonB[,1])),
                         rep(nameB,length(onlypolygonB[,1])))  
  colnames(inBothPolygons_1) = c(colnames(inBothPolygons),nameA,nameB,"assigned_subspecies")
  colnames(inNeitherPolygon_1) = c(colnames(inNeitherPolygon),nameA,nameB,"assigned_subspecies")
  colnames(onlypolygonA_1) = c(colnames(onlypolygonA),nameA,nameB,"assigned_subspecies")
  colnames(onlypolygonB_1) = c(colnames(onlypolygonB),nameA,nameB,"assigned_subspecies")
  
  toReturn = rbind(inBothPolygons_1,inNeitherPolygon_1,
                   onlypolygonA_1,onlypolygonB_1)
  
  return(toReturn)
}

databaseToAssignedSubspecies = function(spp,subspp1,subspp2,pointLimit,dbToQuery,quantile=0.95,xmin=-125,
                               xmax=-60,ymin=10,ymax=50,plotIt=F,bgLayer,...) {
  
  library(dplyr)
  library(spocc)
  library(MASS)
  library(rgeos)
  library(raster)
  ## step 1: get the species
  ## in this case, subspecies
  print("Downloading species occurrences")
  listFromSubspeciesOcc = twoSubspeciesOccGbif(spp=spp,subspp1=subspp1,subspp2=subspp2,
                                               pointLimit=pointLimit,dbToQuery=dbToQuery)
  
  ## label the data by subspeces
  print("Labeling data by subspecies")
  output = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
  labeledLoc = output$fullPointLocLab
  sppOcc = output$sppOcc
  name_sppOcc = output$name_sppOcc
  subsppOcc1 = output$subsppOcc1
  name_subsppOcc1 = output$name_subsppOcc1
  subsppOcc2 = output$subsppOcc2
  name_subsppOcc2 = output$name_subsppOcc2
  
  if(plotIt==T){
    print("Plotting")
    plot(bgLayer, col="grey",colNA="darkgrey")
    points(labeledLoc$longitude,labeledLoc$latitude,
           col=as.factor(labeledLoc$subspecies))
  }
  
  ## to retuce error take only subspecies within main density
  ## clean up the polygons so that if grouping way out in middle of nowhere, get rid of it
  ## remove points that fall within the other subspecies' polygon 
  ## and account for data being poor
  ## build the density of the points
  ## remove all but 95% (or quantile) most dense cells
  
  print("Buidling species kernel density maps")
  spp_densRas = subspeciesDensityMap(localities=labeledLoc[labeledLoc$subspecies==name_sppOcc,],
                                     quantile=quantile,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  subspp1_densRas = subspeciesDensityMap(localities=labeledLoc[labeledLoc$subspecies==name_subsppOcc1,],
                                         quantile=quantile,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  subspp2_densRas = subspeciesDensityMap(localities=labeledLoc[labeledLoc$subspecies==name_subsppOcc2,],
                                         quantile=quantile,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
  
  if (plotIt==T) {
    par(mfrow=c(2,2))
    plot(bgLayer, col="grey",colNA="darkgrey",main="all")
    plot(spp_densRas,add=T,col=viridis::viridis(99))
    plot(subspp1_densRas,add=T,col=viridis::viridis(99))
    plot(subspp2_densRas,add=T,col=viridis::viridis(99))
    
    plot(bgLayer, col="grey",colNA="darkgrey",main=name_sppOcc)
    plot(spp_densRas,add=T,col=viridis::viridis(99))
    plot(bgLayer, col="grey",colNA="darkgrey",main=name_subsppOcc1)
    plot(subspp1_densRas,add=T,col=viridis::viridis(99))
    plot(bgLayer, col="grey",colNA="darkgrey",main=name_subsppOcc2)
    plot(subspp2_densRas,add=T,col=viridis::viridis(99))
    par(mfrow=c(1,1))
  }
  
  ## convert to polygons
  print("Converting density maps to polygons")
  spp_densityPolygon = densityMapToPolygons(densityMap=spp_densRas)
  subspp1_densityPolygon = densityMapToPolygons(densityMap=subspp1_densRas)
  subspp2_densityPolygon = densityMapToPolygons(densityMap=subspp2_densRas)
  
  if (plotIt==T) {
    plot(bgLayer, col="grey",colNA="darkgrey",main="all")
    plot(spp_densityPolygon,add=T,border="black",lwd=7)
    plot(subspp1_densityPolygon,add=T,border="red",lwd=4)
    plot(subspp2_densityPolygon,add=T,border="cyan",lwd=1)
  }
  
  ## check overlaps between polygons
  
  print("Checking Overlaps of Polygons")
  polygonsToRemove = whichPolygonsCompletelyOverlap(subsppPoly1 = subspp1_densityPolygon,
                                                    subsppPoly2 = subspp2_densityPolygon)
  
  subspp1_polygonsToRemove = polygonsToRemove$subspp1poly_toremove
  subspp2_polygonsToRemove = polygonsToRemove$subspp2poly_toremove
  overlapToRemove_subspp1 = polygonsToRemove$subspp1_intToRemove
  overlapToRemove_subspp2 = polygonsToRemove$subspp2_intToRemove
  
  ## remove polygons that are completely within other polygon
  ## TODO: what about things that are in neither polgon?
  ## TODO: what about things that are in both?
  
  print("Removing overlapping polygons")
  subspp1_densityPolygon_trim = trimPolygonsByOverlap(polygon=subspp1_densityPolygon,
                                                      idList = subspp1_polygonsToRemove,
                                                      intList=overlapToRemove_subspp1)
  subspp2_densityPolygon_trim = trimPolygonsByOverlap(polygon=subspp2_densityPolygon,
                                                      idList = subspp2_polygonsToRemove,
                                                      intList=overlapToRemove_subspp2)
  
  if (plotIt==T) {
    plot(bgLayer, col="grey",colNA="darkgrey",main="all")
    plot(spp_densityPolygon,add=T,col="black",lwd=7)
    plot(subspp1_densityPolygon_trim,add=T,border="red",lwd=4)
    plot(subspp2_densityPolygon_trim,add=T,border="cyan",lwd=1)
  }
  
  
  ## remove points that are in wrong polygon
  ## if labeled and in wrong polygon, unlabel
  ## label unlabeled points based on polygons
  
  print("Locating points relative to polygons")
  polyLocations = locatePolygonPoints(test_points = labeledLoc,
                                      polygonA=subspp1_densityPolygon_trim,
                                      polygonB=subspp2_densityPolygon_trim,
                                      nameA=name_subsppOcc1,
                                      nameB=name_subsppOcc2)
  
  ## TODO: see if you can label unlabeled points based on polygons or nearest neighbor
  ## or nearest neighbor 
  ## output the new data
  
  ## plot (optional)
  
  if (plotIt == T) { 
    print("Plotting")
    plot(bgLayer, col="grey",colNA="darkgrey")
    plot(subspp1_densityPolygon_trim,add=T,col=rgb(1,0,0,0),border="darkred",lwd=4)
    plot(subspp2_densityPolygon_trim,add=T,col=rgb(0,1,1,0),border="darkgreen",lwd=1)
    points(polyLocations$longitude,polyLocations$latitude,
           col=polyLocations$assigned_subspecies)
    
    }

  
  ## return nice clean data
  return(polyLocations)
  
}

Env = raster::stack(list.files(
  path='/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  pattern="\\.bil$",
  full.names=T))
ext = raster::extent(c(-125,-60,10,50))
Env = raster::crop(Env, ext)
bg = Env[[1]]

phainopeplaNitens = databaseToAssignedSubspecies(spp="Phainopepla nitens",subspp1="nitens",subspp2="lepida",
                                                 pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
                                                 quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg)
  
  
cardinalisSinuatus = databaseToAssignedSubspecies(spp="Cardinalis sinuatus",subspp1="sinuatus",subspp2="fulvescens",
                                                 pointLimit=500,dbToQuery=c("gbif","bison","inat","ebird","ecoengine","vertnet"),
                                                 quantile=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50,plotIt=T,bgLayer=bg)
