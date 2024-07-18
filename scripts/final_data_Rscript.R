## TODO: need to make sure not doing any comparisons with unknowns
## when it is actually determining the values
## TODO: unit testing


#folder = "C:/Users/kaiya/Documents/GitHub/subsppLabelR"
folder = "/Users/kprovost/Documents/GitHub/subsppLabelR/"
setwd(folder)

detach("package:subsppLabelR", unload = TRUE)
devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
#install.packages(folder,repos=NULL,type="source",force=T)
library(subsppLabelR)

## parameters
redo=T; overwrite=F; downloadOnly=T
pointLimit=10000
#quant_list=sort(c(seq(0.05,0.95,0.05),0.99,0.33,0.67))
quant_list=c(0.5)
dbToQuery=c("gbif","inat","bison","vertnet") #EBIRD_KEY = "f49839r87f7g"
xmin = -180; xmax = -50; ymin = 0; ymax = 75;
if(redo==T) {
  species_list = c(
    #"Cardinalis cardinalis",
    #"Phainopepla nitens",
    #"Cardinalis sinuatus",
    #"Campylorhynchus brunneicapillus",
    #"Amphispiza bilineata",
    #"Vireo bellii",
    #"Toxostoma crissale",
    #"Toxostoma curvirostre",
    #"Polioptila melanura",
    #"Melospiza melodia",
    #"Zonotrichia leucophrys",
    #"Geococcyx californianus",
    "Geococcyx",
    #"Melozone fusca",
    #"Auriparus flaviceps",
    NULL
  )
  
  subspp_list = list(
    #c("affinis","canicaudus","cardinalis","carneus","clintoni", "coccineus","flammiger","floridanus","igneus","littoralis", "magnirostris","mariae", "phillipsi","saturatus","seftoni","sinaloensis","superbus","townsendi","yucatanicus"),
    #c("lepida","nitens"),
    #c("fulvescens","peninsulae","sinuatus"),
    #c("affinis","brunneicapillus","anthonyi","bryanti","guttatus","sandiegensis","seri"),
    #c("bangsi","bilineata","belvederei","cana","deserticola","carmenae","grisea","opuntia","pacifica","tortugae"),
    #c("arizonae","bellii","medius","pusillus"),
    #c("coloradense","crissale","dumosum","trinitatis"),
    #c("celsum","curvirostre","insularum","maculatum","oberholseri","occidentale","palmeri"),
    #c("curtata","melanura","lucida"),
    #c("adusta","amaka","atlantica","beata","caurina","clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella","goldmani","gouldii","graminea","heermanni","inexspectata","insignis","juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia","merrilli","mexicana","micronyx","montana","morphna","pectoralis","pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka","santaecrucis","villai","yuriria"),
    # c("leucophrys","gambelii","nuttalli","pugetensis","oriantha"),
    #NULL,
    c("californianus","velox"),
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
    dir.create(paste(folder,spp,"/",sep=""))
    ## if the raw file does not exist, make it. otherwise load it 
    if(!(file.exists(paste(folder,spp,"/",spp,"_subspplabelR_RAW.txt",sep="")))){
      listFromSubspeciesOcc = subspeciesOccQuery(spp=spp,subsppList=subspp,pointLimit=pointLimit,dbToQuery=dbToQuery)
      labeledLoc = labelSubspecies(subsppOccList=listFromSubspeciesOcc,
                                   spp=spp,subsppList=subspp)
      head(labeledLoc)
      write.table(labeledLoc,paste(folder,spp,"/",spp,"_subspplabelR_RAW.txt",sep=""),sep="\t",row.names = F)
    } else {
      labeledLoc = read.table(paste(folder,spp,"/",spp,"_subspplabelR_RAW.txt",sep=""),sep="\t",header=T)
    }
    
    if(downloadOnly==F){
      subspp_new = unique(labeledLoc$subspecies)
      subspp_new = subspp_new[subspp_new!="unknown"]
      
      for(quant in quant_list){
        print(quant)
        if(overwrite==T | !(file.exists(paste(folder,spp,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep="")))){
          labeledLoc$longitude = as.numeric(labeledLoc$longitude)
          labeledLoc$latitude = as.numeric(labeledLoc$latitude)
          labeledLoc$subspecies = as.factor(labeledLoc$subspecies)
          final = subsppLabelR::databaseToAssignedSubspecies(spp=spp,
                                                             subsppList = subspp_new,
                                                             pointLimit=pointLimit,dbToQuery=dbToQuery,
                                                             quant=quant,
                                                             plotIt=T,
                                                             datafile = labeledLoc,
                                                             xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,
                                                             outputDir=paste(folder,spp,"/",sep=""),
                                                             downloadOnly=F)
          
          
          final$loc_suspect
          final$loc_good
          final$pol
          write.table(final$loc_suspect,paste(folder,spp,"/",spp,"_",quant,"_subspplabelR_loc_suspect.txt",sep=""),sep="\t",row.names = F)
          write.table(final$loc_good,paste(folder,spp,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep=""),sep="\t",row.names = F)
          loc_all = rbind(final$loc_suspect,final$loc_good)
          write.table(loc_all,paste(folder,spp,"/",spp,"_",quant,"_subspplabelR_loc_all.txt",sep=""),sep="\t",row.names = F)
          
          
          
        }
        print(i)
      }
    }
  }
}
## now get the files
setwd(folder)
occfiles=list.files(path=folder,
                    pattern="_subspplabelR_loc_good.txt$",recursive=T,full.names = T)
occ_suspect_files = list.files(path=folder,
                               pattern="_subspplabelR_loc_suspect.txt$",recursive=T,full.names = T)
occfiles = c(occ_suspect_files,occfiles)
occfiles = sort(unique(occfiles))
wcdata = geodata::worldclim_global("bio",res=10,path="~/",download=F) #wcdata = raster::getData(name="worldclim",download=F,path="~/",var="bio",res=10)
wcdata = raster::stack(wcdata)
x = file.info(occfiles)
occfiles = occfiles[order(x$size)]
#occfiles = occfiles[grepl("Toxostoma curvirostre",occfiles)]
occfiles = occfiles[grepl("LOCALITIES",occfiles)]

for(i in (1:length(occfiles))){
  occ_name = basename(occfiles[i])
  occ_dir = dirname(occfiles[i])
  setwd(occ_dir)
  try({
    occ=read.table(occfiles[i],header=T)
    #occ = occ[,c("longitude","latitude","assigned")]
    occ = unique(occ)
    print(occ_name)
    unique(occ$assigned)
    if(nrow(occ)>0){
      #occ = occ[sample(1:nrow(occ),100),]
      
      if(length(unique(occ$assigned))>1){
        print(paste("NUMBER SUBSPECIES:",length(unique(occ$assigned))))
        try({
          Env=wcdata
          loc=occ
          species=occ_name
          runNicheModels=F ## turning this on makes it lengthy but runs the other stuff
          occ_clean = subsppLabelR::localitiesToNicheMath(Env=Env,loc=loc,species=species,runNicheModels=runNicheModels,overwrite=overwrite,numCores=4)
          print(occ_clean)
        })
      } else {
        print("ONLY ONE SUBSPECIES!")
      }
    } else {
      print("NO LOCALITIES!")
    }
  })
}

occ_suspect_files = list.files(path=folder,
                               pattern="_subspplabelR_loc_suspect.txt$",recursive=T,full.names = T)
x = file.info(occ_suspect_files)
occ_suspect_files = occ_suspect_files[order(x$size)]
occ_suspect_files = occ_suspect_files[grepl("0.5",occ_suspect_files)]
occ_suspect_files = occ_suspect_files[grepl("Zonotrichia",occ_suspect_files)]
#occ_suspect_files = occ_suspect_files[!grepl("californianus",occ_suspect_files)]
for(j in (1:length(occ_suspect_files))){
  occ_name = basename(occ_suspect_files[j])
  sus_dir = dirname(occ_suspect_files[j])
  setwd(sus_dir)
  try({
    sus_occ=read.table(occ_suspect_files[j],header=T)
    
    ## get rid of any that are subspecies unknown and assigned none
    sus_occ = sus_occ[!(sus_occ$subspecies=="unknown" & sus_occ$assigned=="none"),]
    sus_occ$assigned = sus_occ$subspecies
    
    ## need to get rid of subspecies that have fewer than 10 points so the kernel can be estimated
    sus_occ = sus_occ[,c("longitude","latitude","assigned")]
    sus_occ = unique(sus_occ)
    sus_tab = table(sus_occ$assigned)
    bad_subspp = names(sus_tab)[which(sus_tab<10)]
    sus_occ = sus_occ[!(sus_occ$assigned %in% bad_subspp),]
    
    print(occ_name) ## note: occ_name might accidentally be global
    unique(sus_occ$assigned)
    if(nrow(sus_occ)>0){
      if(length(unique(sus_occ$assigned))>1){
        print(paste("NUMBER SUBSPECIES:",length(unique(sus_occ$assigned))))
        try({
          Env=wcdata
          loc=sus_occ
          species=occ_name
          runNicheModels=F
          sus_clean = subsppLabelR::localitiesToNicheMath(Env=Env,loc=loc,species=species,runNicheModels=runNicheModels)
        })
      } else {
        print("ONLY ONE SUBSPECIES!")
      }
    } else {
      print("NO LOCALITIES!")
    }
  })
  
}


## combine loc_good and loc_suspect into loc_all for each if not already exist
occ_suspect_files = list.files(path=folder,
                               pattern="_subspplabelR_loc_suspect.txt$",recursive=T,full.names = T)
occfiles=list.files(path=folder,
                    pattern="_subspplabelR_loc_good.txt$",recursive=T,full.names = T)

