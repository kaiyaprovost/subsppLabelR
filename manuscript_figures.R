library(raster)
#species = "Cardinalis sinuatus"
#cutoff = "0.9"

cutoff_list = c("0.5","0.6","0.7","0.8","0.9","0.95")
species_list = c("Cardinalis sinuatus","Phainopepla nitens","Geococcyx","Toxostoma curvirostre","Zonotrichia leucophrys")

for(species in species_list) {
  for(cutoff in cutoff_list) {
    print(paste(species,cutoff))
    
    my_bad_file = paste("/Users/kprovost/Documents/Research/subsppLabelR/",species,"/",species,"_",cutoff,"_subspplabelR_loc_suspect.txt",sep="")
    my_good_file = paste("/Users/kprovost/Documents/Research/subsppLabelR/",species,"/",species,"_",cutoff,"_subspplabelR_loc_good.txt",sep="")
    my_good_rasters = paste("/Users/kprovost/Documents/Research/subsppLabelR/",species,"/DensityRaster_",species,"_",cutoff,".tif",sep="")
    ras = stack(my_good_rasters)
    ras = ras[[-1]] ## remove unknown
    df = read.table(my_good_file,header=T,sep="\t")
    df_bad = read.table(my_bad_file,header=T,sep="\t")
    
    png(paste("/Users/kprovost/Documents/Research/subsppLabelR/",species,"/",species,"_",cutoff,".plot.png",sep=""))
    for(i in 1:length(names(ras))) {
      if(i == 1) {
        plot(ras[[i]],main=paste(species,cutoff),col="grey")
      } else {
        plot(ras[[i]],add=T,col="grey")
      }
    }
    points(df_bad$longitude,df_bad$latitude,col="lightgrey",
           pch=0)
    points(df$longitude,df$latitude,col=as.numeric(as.factor(df$assigned)),
           pch=as.numeric(as.factor(df$assigned)))
    legend("topleft",legend=unique(df$assigned),
           col=1:length(ras),pch=1:length(ras))
    dev.off()
  }
}

