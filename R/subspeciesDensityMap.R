#' @import raster
#' @import MASS
#' @import spocc
#' @import dplyr
#' @import sp
#' @import viridis
#' @import caret
#' @import sf
#' @import rebird
NULL
#' Generate Subspecies Density Maps
#'
#' This function uses 2D kernel density estimation to make a density raster
#' of points for each subspecies. It then removes all but the most
#' dense cells as definied by quant.
#'
#' @param localities Labeled localities as generated from labelSubspecies() or subspeciesOccQuery()
#' @param quant quant for density, below which points are removed. E.g., if set to 0.95, removes 95% least dense squares.
#' @param xmin Minimum longitude extent to clip to
#' @param xmax Maximum longitude extent to clip to
#' @param ymin Minimum latitute extent to clip to
#' @param ymax Maximum latitude extent to clip to
#'
#' @export
#' @examples
#'
#' listFromSubspeciesOcc = subspeciesOccQuery(spp="Cardinalis sinuatus",
#'    subsppList=c("sinuatus","peninsulae","fulvescens"),pointLimit=100,
#'    dbToQuery="gbif")
#' labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc)
#' locs = labeledLoc[labeledLoc$subspecies=="sinuatus",]
#' dens = subspeciesDensityMap(localities=locs,quant=0.95,
#'    xmin=-125,xmax=-60,ymin=10,ymax=50)
subspeciesDensityMap = function(localities,
                                quant = 0.95,
                                xmin = NULL,
                                xmax = NULL,
                                ymin = NULL,
                                ymax = NULL,
                                total_range,
                                relative=T,
                                raw_raster=T,
                                subspp,
                                spp,
                                outputDir) {
  ## this function uses kernel density to make a raster that will then be used to filter
  ## the data to remove all but the 5% (or 1-quant) most dense cells
  ## TODO: allow for subspecies-specific quants
  #library(MASS)
  #library(raster)
  if(is.null(xmin)) { xmin = min(localities$longitude,na.rm=T) }
  if(is.null(xmax)) { xmax = max(localities$longitude,na.rm=T) }
  if(is.null(ymin)) { ymin = min(localities$latitude,na.rm=T) }
  if(is.null(ymax)) { ymax = max(localities$latitude,na.rm=T) }
  range = c(xmin, xmax, ymin, ymax)
  ext = raster::extent(range)
  w1 = matrix(1, 3, 3)
  ## generate the two dimensional kernel density estimation
  if (nrow(localities) == 1) {
    print("NOT ENOUGH LOCALITIES")
    return(NULL)
  }
  density = NULL
  try({density = MASS::kde2d(as.numeric(localities$longitude),as.numeric(localities$latitude),
                             lims = range,n = 50)})
  if(is.null(density)){
    print("NULL")
    return(NULL)
  } else {
    ## convert to raster
    densRas = raster::raster(density)
    if(raw_raster==T){
      raw_file = paste(outputDir,spp," ",subspp,"_raw_raster.tif",sep="")
      if(!file.exists(raw_file)){
        writeRaster(densRas,paste(outputDir,spp," ",subspp,"_raw_raster.tif",sep=""),
                    format="GTiff",overwrite=F)
      }

    }

    if(relative==T){
      print("rescaling")
      values(densRas) = scales::rescale(values(densRas),c(0,1))
      values(densRas) = round(values(densRas),digits=10)
    }
    ## we want to print these out before they are clipped

    ## take the top percentile of the points, only the densest areas
    quan = quantile(densRas[densRas], quant)
    densRas_trim = densRas
    densRas_trim[densRas_trim <= quan] = NA
    plot(densRas_trim,colNA="black")    #,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
    total_crs = raster::crs(total_range)
    total_ext = raster::extent(total_range)
    total_res = raster::res(total_range)
    ## crop to the existing layers
    densRas_trim = raster::crop(densRas_trim,total_range)
    densRas_trim = raster::projectRaster(densRas_trim,total_range)
    raster::values(densRas_trim)[is.na(raster::values(total_range))] = NA

    return(densRas_trim)
  }
}
