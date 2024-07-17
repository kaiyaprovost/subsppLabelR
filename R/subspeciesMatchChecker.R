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
#' Check Subspecies Matches and Labels
#'
#' This function takes the locations of points relative to polygons from locatePolygonPoints() or similar
#' and checks where points are assigned. Points that are in 1) no polygons, 2) multiple polygons, or
#' 3) have a subspecies label that does not match their assignment are flagged as "suspect".
#' Any points that are in single polygons and either 1) were originally unlabeled or 2) have a subspecies
#' label that matches their assignment are flagged as "good".
#'
#' @param locfile Points with subspecies assignments and boolean locations with respect to polygons.
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_sin = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' locs_ful = labeledLoc[labeledLoc$subspecies=="fulvescens",]
#' dens_sin = subspeciesDensityMap(localities=locs_sin,quant=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50)
#' dens_ful = subspeciesDensityMap(localities=locs_ful,quant=0.95,xmin=-125,xmax=-60,ymin=10,ymax=50)
#' densPol_sin = densityMapToPolygons(densityMap=dens_sin)
#' densPol_ful = densityMapToPolygons(densityMap=dens_ful)
#' polyLocations = labeledLoc
#' polyLocations = locatePolygonPoints(test_points=polyLocations,polygonA=densPol_sin,polygonB=densPol_ful,crs="+proj=longlat +ellps=WGS84",setcoord = TRUE,nameA="sinuatus",nameB="fulvescens")
#' checked = subspeciesMatchChecker(locfile = polyLocations)
#' checked_suspect = checked$suspect
#' checked_good = checked$good
subspeciesMatchChecker = function(locfile, subsppNames) {
  #print("a")
  locWithSubspecies = locfile
  #View(locfile)
  #print("b")
  #subsppNames = names(locWithSubspecies[5:length(names(locWithSubspecies))])
  #print("c")
  #print(subsppNames)
  if("unknown" %in% subsppNames){
    numSub = length(subsppNames)-1
  } else {
    numSub = length(subsppNames)
  }
  if(numSub>1){
    #print("d")
    numPoints = length(rownames(locWithSubspecies))
    #print("e")
    lastSubsppCol = length(colnames(locWithSubspecies))
    #View(locWithSubspecies)
    #print("f")
    ## FAILING HERE BECAUSE HAVE REMOVED SUBSPECIES AS WE WENT ALONG
    #if(lastSubsppCol<numSub) {
    subsppAssignCol = locWithSubspecies[, (1+ which(colnames(locWithSubspecies)=="subspecies")):lastSubsppCol]
    #} else {
    #  subsppAssignCol = locWithSubspecies[, (1+lastSubsppCol-numSub):lastSubsppCol]
    #}
    #print(head(locWithSubspecies))
    # for (colnum in 5:length(colnames(locWithSubspecies))){
    #   print(paste("colnum",colnum))
    #   print(head(locWithSubspecies[,colnum]))
    #   num_for_name = as.integer(colnames(locWithSubspecies)[colnum])
    #   print(class(num_for_name))
    #   print(paste("numforname",num_for_name))
    #   print(subsppNames)
    #   name_to_replace = subsppNames[num_for_name]
    #   print(paste("nametoreplace",name_to_replace))
    #   names(locWithSubspecies)[colnum] = name_to_replace
    #   print(names(locWithSubspecies))
    # }
    #print("g")
    subsppPriorCol = locWithSubspecies$subspecies
    #print("h")
    subsppAssignCol = as.data.frame(subsppAssignCol)
    locWithSubspecies$numSubsppGroups <-
      rowSums(subsppAssignCol, na.rm = TRUE)
    #print("i")
    locWithSubspecies$assigned = NA
    #print("j")
    notassigned = locWithSubspecies[which(locWithSubspecies$numSubsppGroups ==
                                            0),]
    #print("k")
    if (nrow(notassigned) != 0) {
      notassigned$assigned = "none"
    }
    #print("l")
    multigroup = locWithSubspecies[which(locWithSubspecies$numSubsppGroups >
                                           1),]
    #print("m")
    if (nrow(multigroup) != 0) {
      multigroup$assigned = "multiple"
    }
    #print("n")
    singlegroup = locWithSubspecies[which(locWithSubspecies$numSubsppGroup ==
                                            1),]
    #print("o")
    ## check whether mismatch between apriori and not
    suspectpoints = rbind(multigroup, notassigned)
    goodpoints = data.frame()
    #print("p")
    for (column in 5:lastSubsppCol) {
      print(column)
      name = colnames(singlegroup)[column]
      #print(name)
      assignedThat = singlegroup[singlegroup[, column] == 1,]
      ## breaking if none assigned
      if (nrow(assignedThat) == 0) {
        print("NONE ASSIGNED:")
        print(name)
      } else {
        ## subspp should be unk or the subspp
        assignedThat$assigned = name
      }
      okaynames = c("unknown", name)
      wrong1 = assignedThat[which(!(assignedThat$subspecies %in% okaynames)),]
      right1 = assignedThat[which(assignedThat$subspecies %in% okaynames),]
      notassignedThat = singlegroup[singlegroup[, column] == 0,]
      wrong2 = notassignedThat[which(notassignedThat$subspecies == name),]
      right2 = notassignedThat[which(notassignedThat$subspecies != name),]
      ## DO NOT ADD RIGHT/WRONG 2, NOT ENOUGH INFORMATION YET
      suspectpoints = rbind(suspectpoints, wrong1)
      goodpoints = rbind(goodpoints, right1)
    }
    #print("q")
    suspectpoints = unique(suspectpoints)
    goodpoints = unique(goodpoints)
    #print("r")
  } else {
    suspectpoints = data.frame()
    goodpoints = locWithSubspecies
  }
  return(list(suspect = suspectpoints,
              good = goodpoints))
}

