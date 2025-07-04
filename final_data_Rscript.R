detach("package:subsppLabelR",unload = TRUE)
#pkgDirectory = "C:/Users/kaiya/Documents/Work/GitHub/subsppLabelR/" ## WINDOWS
pkgDirectory = "/Users/kprovost/Documents/GitHub/subsppLabelR/" ## MAC
#devtools::install_github('kaiyaprovost/subsppLabelR',force=T)
roxygen2::roxygenize(package.dir = pkgDirectory)
install.packages(pkgDirectory,repos = NULL,type = "source",force = T)
pkgDirectory = "/Users/kprovost/Documents/Research/subsppLabelR/" ## MAC
library(subsppLabelR)
## parameters
redo = F; overwrite = T
pointLimit = 2000
#quant_list=rev(c(seq(0.5,0.9,0.1),0.95))
quant_list=rev(0.6)
dbToQuery = c("gbif","inat","bison","vertnet") #EBIRD_KEY = "f49839r87f7g"
xmin = -180; xmax = 180; ymin = -90; ymax = 90
spp_epsilon = 10^-6; subspp_epsilon = 10^-3.5
if (redo == T) {
  species_list = c(
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
    #"Chondestes grammacus",
    #"Calidris pusilla",
    #"Calidris mauri",
    #"Calidris alpina",
    #"Calidris alba",
    #"Myiarchus crinitus",
    #"Sayornis phoebe",
    #"Contopus cooperi",
    #"Contopus virens",
    #"Contopus sordidulus",
    NULL
  )
  
  subspp_list = list(#c("lepida","nitens"),
                     #c("fulvescens","peninsulae","sinuatus"),
                     #c("affinis","brunneicapillus","anthonyi","bryanti","guttatus","sandiegensis","seri"),
                     #c("bangsi","bilineata","belvederei","cana","deserticola","carmenae","grisea","opuntia","pacifica","tortugae"),
                     #c("arizonae","bellii","medius","pusillus"),
                     #c("coloradense","crissale","dumosum","trinitatis"),
                     #c("celsum","curvirostre","insularum","maculatum","oberholseri","occidentale","palmeri"),
                     #c("curtata","melanura","lucida"),
                     #c("adusta","amaka","atlantica","beata","caurina","clementae","cleonensis","cooperi","coronatorum","euphonia","fallax","fisherella","goldmani","gouldii","graminea","heermanni","inexspectata","insignis","juddi","kenaiensis","mailliardi","maxillaris","maxima","melodia","merrilli","mexicana","micronyx","montana","morphna","pectoralis","pusillula","rivularis","rufina","saltonis","samuelis","samuelsis","sanaka","santaecrucis","villai","yuriria"),
                     #c("leucophrys","gambelii","nuttalli","pugetensis","oriantha"),
                     #NULL,
                     c("californianus","velox"),
                     #c("campoi","fusca","intermedia","jamesi","mesata","mesoleuca","perpallida","potosina","texana","toroi"),
                     #c("acaciarum","flaviceps","hidalgensis","lamprocephalus","ornatus","sinaloae"),
                     #c("grammacus","strigatus"),
                     #NULL,
                     #NULL,
                     #c("articola","pacifica","hudsonia","arctica","schinzii","alpina","centralis","sakhalina","actites","kistchinski"),
                     #c("alba","rubida"),
                     #NULL,
                     #NULL,
                     #c("cooperi","marjorinus"),
                     #NULL,
                     #c("saturatus","veliei","peninsulae","sordidulus"),
                     NULL)
  
  for (i in sort(1:length(species_list))) {
    print(i)
    spp = species_list[i]
    subspp = subspp_list[[i]]
    print(spp)
    print(subspp)
    dir.create(paste(pkgDirectory,spp,"/",sep = ""))
    if (!(file.exists(paste(pkgDirectory,spp,"/",spp,"_subspplabelR_RAW.txt",sep = "")))) {
      listFromSubspeciesOcc = subspeciesOccQuery(spp = spp,subsppList = subspp,pointLimit = pointLimit,dbToQuery = dbToQuery)
      labeledLoc = labelSubspecies(subsppOccList = listFromSubspeciesOcc,spp = spp,subsppList = subspp) 
      head(labeledLoc)
      write.table(labeledLoc,paste(pkgDirectory,spp,"/",spp,"_subspplabelR_RAW.txt",sep = ""),sep = "\t",row.names = F)
    } else { labeledLoc = read.table(paste(pkgDirectory,spp,"/",spp,"_subspplabelR_RAW.txt",sep = ""),sep = "\t",header = T) }
    subspp_new = unique(labeledLoc$subspecies)
    subspp_new = subspp_new[subspp_new != "unknown"]
    for (quant in quant_list) {
      print(quant)
      if (overwrite == T | !(file.exists(paste(pkgDirectory,spp,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep = "")))) {
        labeledLoc$longitude = as.numeric(labeledLoc$longitude)
        labeledLoc$latitude = as.numeric(labeledLoc$latitude)
        labeledLoc$subspecies = as.factor(labeledLoc$subspecies)
        final = subsppLabelR::createAssignedSubspecies(spp = spp,subsppList = subspp_new,pointLimit = pointLimit,dbToQuery = dbToQuery,quant = quant,plotIt = T,datafile = labeledLoc,xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax,outputDir = paste(pkgDirectory,spp,"/",sep = ""),verbose = T,spp_epsilon=spp_epsilon,subspp_epsilon=subspp_epsilon)
        write.table(final$loc_suspect,paste(pkgDirectory,spp,"/",spp,"_",quant,"_subspplabelR_loc_suspect.txt",sep = ""),sep = "\t",row.names = F)
        write.table(final$loc_good,paste(pkgDirectory,spp,"/",spp,"_",quant,"_subspplabelR_loc_good.txt",sep = ""),sep = "\t",row.names = F)
      }
    }
    print(i)
  }
}
## now get the files
setwd(pkgDirectory)
occfiles = list.files(path = pkgDirectory,pattern = "_subspplabelR_loc_good.txt$",recursive = T,full.names = T)
wcdata = raster::stack(geodata::worldclim_global("bio",res = 10,path = "~/",download = F)) 
x = file.info(occfiles)
#occfiles = occfiles[order(x$size)]
for (i in c(1:length(occfiles))) {
  occ_name = basename(occfiles[i])
  occ_dir = dirname(occfiles[i])
  setwd(occ_dir)
  occ = read.table(occfiles[i],header = T)
  occ = unique(occ)
  print(occ_name)
  unique(occ$assigned)
  if (nrow(occ) > 0) {
    if (length(unique(occ$assigned)) > 1) {
      print(paste("NUMBER SUBSPECIES:",length(unique(occ$assigned))))
      try({
        Env = wcdata
        loc = occ
        species = occ_name
        runNicheModels = F
        occ_clean = subsppLabelR::localitiesToNicheMath(Env = Env,loc = loc,species = species,runNicheModels = runNicheModels,
                                                        overwrite=F)
      })
    } else { print("ONLY ONE SUBSPECIES!") }
  } else { print("NO LOCALITIES!") }
}
occ_suspect_files = list.files(path = pkgDirectory,pattern = "_subspplabelR_loc_suspect.txt$",recursive = T,full.names = T)
for (j in c(1:length(occ_suspect_files))) {
  sus_name = basename(occ_suspect_files[j])
  sus_dir = dirname(occ_suspect_files[j])
  setwd(sus_dir)
  sus_occ = read.table(occ_suspect_files[j],header = T)
  sus_occ = sus_occ[!(sus_occ$subspecies == "unknown" & sus_occ$assigned == "none"),]
  sus_occ$assigned = sus_occ$subspecies
  sus_occ = sus_occ[,c("longitude","latitude","assigned")]
  sus_occ = unique(sus_occ)
  sus_tab = table(sus_occ$assigned)
  bad_subspp = names(sus_tab)[which(sus_tab < 10)]
  sus_occ = sus_occ[!(sus_occ$assigned %in% bad_subspp),]
  print(sus_name)
  unique(sus_occ$assigned)
  if (nrow(sus_occ) > 0) {
    if (length(unique(sus_occ$assigned)) > 1) {
      print(paste("NUMBER SUBSPECIES:",length(unique(sus_occ$assigned))))
      try({
        Env = wcdata
        loc = sus_occ
        species = sus_name
        runNicheModels = T
        sus_clean = subsppLabelR::localitiesToNicheMath(Env = Env,loc = loc,species = species,runNicheModels = runNicheModels)
      })
    } else { print("ONLY ONE SUBSPECIES!") }
  } else { print("NO LOCALITIES!") }
}
