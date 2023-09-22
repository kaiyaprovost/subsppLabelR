redo=T; overwrite=T
detach("package:subsppLabelR", unload = TRUE)
#devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
install.packages("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",repos=NULL,type="source",force=T)
library(subsppLabelR)
## parameters
pointLimit=10000
quant=0.95
dbToQuery=c("gbif","inat","bison","vertnet") #EBIRD_KEY = "f49839r87f7g"
xmin = -180; xmax = -50; ymin = 0; ymax = 75;
if(redo==T) {
  species_list = c("Phainopepla nitens",
                   #"Cardinalis sinuatus",
                   #"Campylorhynchus brunneicapillus","Amphispiza bilineata",
                   #"Vireo bellii","Toxostoma crissale","Toxostoma curvirostre","Polioptila melanura",
                   #"Melospiza melodia","Zonotrichia leucophrys","Geococcyx californianus","Geococcyx",
                   #"Melozone fusca","Auriparus flaviceps",
                   NULL
  )

  subspp_list = list(c("lepida","nitens"),
                     #c("sinuatus","fulvescens","peninsulae"),
                     #c("affinis","brunneicapillus","anthonyi","bryanti","guttatus","sandiegensis","seri"),
                     #c("bangsi","bilineata","belvederei","cana","deserticola","carmenae","grisea","opuntia","pacifica","tortugae"),
                     #c("medius","bellii","arizonae","pusillus"),
                     #c("coloradense","crissale","dumosum","trinitatis"),
                     #c("celsum","curvirostre","insularum","maculatum","oberholseri","occidentale","palmeri"),
                     #c("curtata","melanura","lucida"),
                     #c("adusta","amaka","atlantica","beata","caurina","clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella","goldmani","gouldii","graminea","heermanni","inexspectata","insignis","juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia","merrilli","mexicana","micronyx","montana","morphna","pectoralis","pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka","santaecrucis","villai","yuriria"),
                     #c("leucophrys","gambelii","nuttalli","pugetensis","oriantha"),
                     #NULL,
                     #c("californianus","velox"),
                     #c("campoi","fusca","intermedia","jamesi","mesata","mesoleuca","perpallida","potosina","texana","toroi"),
                     #c("acaciarum","flaviceps","hidalgensis","lamprocephalus","ornatus","sinaloae"),
                     NULL
  )

  for(i in 1:length(species_list)) {
    print(i)
    spp = species_list[i]
    subspp = subspp_list[[i]]
    print(spp)
    print(subspp)
    dir.create(paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",sep=""))
    if(overwrite==T | !(file.exists(paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_RAW.txt",sep="")))){
      listFromSubspeciesOcc = subspeciesOccQuery(spp=spp,subsppList=subspp,pointLimit=pointLimit,dbToQuery=dbToQuery)
      labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc,
                                   spp=spp,subsppList=subspp,cleanup_nominate=T)
      head(labeledLoc)
      write.table(labeledLoc,paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_RAW.txt",sep=""),sep="\t",row.names = F)
    } else {
      labeledLoc = read.table(paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_RAW.txt",sep=""),sep="\t",header=T)
    }
    subspp_new = unique(labeledLoc$subspecies)
    subspp_new = subspp_new[subspp_new!="unknown"]
    if(overwrite==T | !(file.exists(paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep="")))){
      labeledLoc$longitude = as.numeric(labeledLoc$longitude)
      labeledLoc$latitude = as.numeric(labeledLoc$latitude)
      labeledLoc$subspecies = as.factor(labeledLoc$subspecies)

      final = subsppLabelR::databaseToAssignedSubspecies(spp=spp,
                                                         subsppList = subspp_new,
                                                         pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                         quantile=quant,
                                                         plotIt=T,
                                                         datafile = labeledLoc,
                                                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                                                         outputDir="C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
      final$loc_suspect
      final$loc_good
      final$pol
      write.table(final$loc_suspect,paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_loc_suspect.txt",sep=""),sep="\t",row.names = F)
      write.table(final$loc_good,paste("C:/Users/kaiya/Documents/GitHub/subsppLabelR/",spp,"_",quant,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep=""),sep="\t",row.names = F)
    }
    print(i)
  }
}
## now get the files
setwd("C:/Users/kaiya/Documents/GitHub/subsppLabelR/")
occfiles=list.files(path="C:/Users/kaiya/Documents/GitHub/subsppLabelR/",
                    pattern="_subspplabelR_loc_good",recursive=T,full.names = T)
wcdata = geodata::worldclim_global("bio",res=10,path="~/",download=F) #wcdata = raster::getData(name="worldclim",download=F,path="~/",var="bio",res=10)
wcdata = raster::stack(wcdata)

for(i in sample(1:length(occfiles))){
  occ_name = basename(occfiles[i])
  occ=read.table(occfiles[i],header=T)
  print(occ_name)
  if(nrow(occ)>0){
    #occ = occ[sample(1:nrow(occ),100),]
    try({
      occ_clean = subsppLabelR::localitiesToNicheMath(Env=wcdata,loc=occ,species=occ_name)
    })
  } else {
    print("NO LOCALITIES!")
  }

}

