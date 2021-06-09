installLibraries <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = libDir)
  if (!requireNamespace("biomaRt", quietly = TRUE))
    BiocManager::install("biomaRt", lib = libDir)
  if (!requireNamespace("GOSim", quietly = TRUE))
    BiocManager::install("GOSim", lib = libDir)
  if (!requireNamespace("DOSE", quietly = TRUE))
    BiocManager::install("DOSE", lib = libDir)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db", lib = libDir)
  if (!requireNamespace("pathfindR", quietly = TRUE))
    BiocManager::install("pathfindR", lib = libDir)
  if (!requireNamespace("RCy3", quietly = TRUE))
    BiocManager::install("RCy3", lib = libDir)

  if (!requireNamespace("shiny", quietly = TRUE))
    install.packages("shiny", lib = libDir)
  
  if (!requireNamespace("shinyFiles", quietly = TRUE))
    install.packages("shinyFiles", lib = libDir)
  
  if (!requireNamespace("plotly", quietly = TRUE))
    install.packages("plotly", lib = libDir)
  if (!requireNamespace("DT", quietly = TRUE))
    install.packages("DT", lib = libDir)
  if (!requireNamespace("meffil", quietly = TRUE)) {
    install.packages("devtools") # if the devtools package is not installed
    library(devtools)
    install_github("perishky/meffil", lib = libDir)
    #install.packages("meffil")
  }
  #install_version("foobarbaz", "0.1.2")
  if (!requireNamespace("ggforce", quietly = TRUE))
    install.packages("ggforce", lib = libDir)
  if (!requireNamespace("foreach", quietly = TRUE))
    install.packages("foreach", lib = libDir)
  if (!requireNamespace("psych", quietly = TRUE))
    install.packages("psych", lib = libDir)
  if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse", lib = libDir)
  if (!requireNamespace("crosstalk", quietly = TRUE))
    install.packages("crosstalk", lib = libDir)
  if (!requireNamespace("reactable", quietly = TRUE))
    install.packages("reactable", lib = libDir)
  if (!requireNamespace("ggExtra", quietly = TRUE))
    install.packages("ggExtra", lib = libDir)
  # if (!requireNamespace("heatmaply", quietly = TRUE))
  #   install.packages("heatmaply", lib = libDir)
  library(devtools)
  # if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  #   install_github("jokergoo/ComplexHeatmap", lib = libDir)
  # if (!requireNamespace("InteractiveComplexHeatmap", quietly = TRUE))
  #   install_github("jokergoo/InteractiveComplexHeatmap", lib = libDir)
  # 
  # if (!requireNamespace("magick", quietly = TRUE))
  #   install.packages("magick", lib = libDir)
  # if (!requireNamespace("fastcluster", quietly = TRUE))
  #   install.packages("fastcluster", lib = libDir)
  if (!requireNamespace("shinyBS", quietly = TRUE))
    install.packages("shinyBS", lib = libDir)
  if (!requireNamespace("reactlog", quietly = TRUE))
    install.packages("reactlog", lib = libDir)
  if (!requireNamespace("config", quietly = TRUE))
    install.packages("config", lib = libDir)
}

loadLibraries <- function(){
  # library(remotes)
  # library(GOSim)
  # detach(package:GOSim,unload = TRUE, force = TRUE)
  # library(biomaRt)
  # detach(package:biomaRt,unload = TRUE, force = TRUE)
  
  library(shiny)
  #library(shinyFiles)
  library(DT)
  library(ggplot2)
  library(ggExtra) # for ggMarginal
  #library(Cairo)   # For nicer ggplot2 output when deployed on Linux
  library(data.table)
  library(meffil)
  #install.packages("ggforce")
  library(ggforce)
  library(foreach)
  library(psych) # for descriptives and winsor
  library(tidyverse)
  library(dplyr)
  library(plotly)
#  library(crosstalk)
#  library(reactable)
  library(RColorBrewer)
  #library(heatmaply)
#  library(ComplexHeatmap)
#  library(InteractiveComplexHeatmap)
  library(magick)
#  library(fastcluster)
  #library(reactlog)
#  library(org.Hs.eg.db)
  #library(pathfindR)
  library(config)
}

checkConfigVariables <- function() {
#  browser()
  # if (is_null(config$errorTest)) {
  #   stop("errorTest")
  # }
  if (is_null(config$PHENOFileName)) {
    stop("PHENOFileName missing")
  }
  if (is_null(config$firstPHENOVar)) {
    stop("firstPHENOVar missing")
  }
  if (is_null(config$lastPHENOVar)) {
    stop("lastPHENOVar missing")
  }
  if (is_null(config$AdjustFileName)) {
    stop("AdjustFileName missing")
  }
  if (is_null(config$betaFileName)) {
    stop("betaFileName missing")
  }
  if (is_null(config$EWAScatalogFileName)) {
    stop("EWAScatalogFileName missing")
  }
  if (is_null(config$MultiModProbesFileName)) {
    error(multiModProbesFileNameMissing)
  }
  if (is_null(config$baseDir)) {
    stop("baseDir not set")
  }
  if (is_null(config$winsorTrim)) {
    stop("winsorTrim not set")
  }
  if (is_null(config$mergeAttribut)) {
    stop("mergeAttribut not set")
  }
  if (is_null(config$genderAttribut)) {
    stop("genderAttribut not set")
  }
  if (is_null(config$genderFemaleValue)) {
    stop("genderFemaleValue not set")
  }
  if (is_null(config$genderMaleValue)) {
    stop("genderMaleValue not set")
  }
  if (is_null(config$debugMode)) {
    stop("debugMode not set")
  }
}

getTraitsDFLong <- function(sessionVariables) {
  if (dir.exists(sessionVariables$folder)) {
    tryCatch({
      setwd(sessionVariables$folder)
    }, error=function(err){
      message(paste0("unable to open work folder ", sessionVariables$folder))
      return()
    });
    tryCatch({
    }, error=function(err){
      #      print(paste0("unable open work folder ", dataDir))
      errortext = paste0("unable open read config.yml in ", sessionVariables$folder)
      message(errortext)
      id <- showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    });
    tryCatch({
      PHENOFileName <- config$PHENOFileName
      sessionVariables$dataFileName = PHENOFileName
      PHENO<-fread(file=PHENOFileName,sep="\t",dec=".",header=TRUE)
      PHENO<-as.data.frame(PHENO)
      colnames(PHENO)<-gsub(" ",".",colnames(PHENO)) # replace " " with "." for compatibility to filenames
      colnames(PHENO)<-gsub("-",".",colnames(PHENO)) # replace "-" with "." for compatibility to filenames
      AdjustFileName <- config$AdjustFileName
      PHENO1<-fread(file=AdjustFileName,sep="\t",dec=".")
      P<-subset(PHENO1, select=c(config$mergeAttribut, config$genderAttribut))
      P = as.data.frame(P)
      P1 = base::merge(PHENO, P, by.x = config$mergeAttribut, by.y = config$mergeAttribut, all.x = FALSE, all.y=FALSE)
      PHENO<-P1
      firstPHENOVar <- config$firstPHENOVar + ncol(P) - 1
      lastPHENOVar <- config$lastPHENOVar + ncol(P) - 1
#      PHENO = addXToName(PHENO,firstPHENOVar,lastPHENOVar)
#      PHENOWinsorized = winsorize(PHENO,0.1,firstPHENOVar,lastPHENOVar)
#browser()
      PHENOWinsorized = winsorize(PHENO,config$winsorTrim,firstPHENOVar,lastPHENOVar)
      PHENO = PHENOWinsorized
      return (PHENO)
    }, error=function(err){
#      print(paste0("unable open work folder ", dataDir))
      errortext = paste0("unable open work folder ", sessionVariables$folder)
      message(errortext)
      id <- showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    });
  }
}

loadObjects <- function(){
#browser()
  if (dir.exists(dataDir)) {
  #   PHENO <- getTraitsDFLong(dataDir)
  # 
  #   assign("PHENO",PHENO,envir=globalenv())
  # #  typeof(PHENO[15,15])
  #   
  #   PHENOFull<-PHENO
  #   assign("PHENOFull",PHENOFull,envir=globalenv())
  #   
    #load list of multimodal CpG
    MultiModProbesFileName <- config$MultiModProbesFileName
    #MultiModCpG<-fread(file="y:/Gruppen/immu/Epigenetik/450k/MultiModalCpG.csv",sep="\t",dec=".")
    MultiModProbes<-fread(file=MultiModProbesFileName,sep="\t",dec=".")
    assign("MultiModProbes",MultiModProbes,envir=globalenv())
    
    betaFileName <- config$betaFileName
    #beta <- fread("f:/roeder/methylation_transposed2.csv",stringsAsFactors=FALSE,header=TRUE,sep="\t")
    beta <- fread(betaFileName,stringsAsFactors=FALSE,header=TRUE,sep="\t")
    
    beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
    #beta<-as.data.table(beta)
    
    beta_wo_outliers<-removeOutliers3IQR(as.matrix(beta))
    beta_wo_outliers<-as.data.frame(beta_wo_outliers[[1]])
    beta<-beta_wo_outliers
#browser()
    beta <- removeMultiModelCpGFromBeta(beta,MultiModProbes)

    nonVariableProbesFileName <- config$nonVariableProbesFileName
#    nonVariableProbes <- fread(file = nonVariableProbesFileName,sep = "\t",dec = ".")
#    nonVariableProbes <- readRDS(file="f:/roeder/nonVariableProbes.RDS")
    nonVariableProbes <- readRDS(file = nonVariableProbesFileName)
    assign("nonVariableProbes",nonVariableProbes,envir=globalenv())
    beta <- removeNonVariableProbes(beta,nonVariableProbes)

    assign("beta",beta,envir=globalenv())
    beta.t<-t(beta)
    assign("beta.t",beta.t,envir=globalenv())
    
    annotation <- meffil.get.features("450k")
    annotation$relation.to.island = as.factor(annotation$relation.to.island)
    #remove unmeasured or multimodal probeIDs from annotation 
    annotation = annotation[which(annotation$name %in% rownames(beta)),]
    assign("annotation",annotation,envir=globalenv())
    
    EWAScatalogFileName <- config$EWAScatalogFileName
    EWAScatalog = fread(file=paste0(appDir,EWAScatalogFileName), sep = "\t")
    assign("EWAScatalog",EWAScatalog,envir=globalenv())
    EWAScatalogCount = EWAScatalog %>% group_by(CpG) %>% tally(!is.na(CpG))
    assign("EWAScatalogCount",EWAScatalogCount,envir=globalenv())
#    threshold = 10
  }
}

setBaseSciPen<-function(){
  options(scipen=-10, digits=5)
}

resetSciPen<-function(){
  options(scipen=0, digits=7)
}

# addXToName <- function(traitDF,firstPHENOVar,lastPHENOVar) {
# #browser()
#   foreach(i = firstPHENOVar:lastPHENOVar, .combine = cbind, .verbose=FALSE) %do% {
#     tryCatch({
# #browser()
#       if (!is.na(as.numeric(colnames(traitDF)[i]))) { #check whether name is numeric only
#         colnames(traitDF)[i] = paste0("X",colnames(traitDF)[i])
#       }
#     }, warning=function(warn){
#        #catch unwanted warning messages
# #      print(paste0("Warning: ", warn$message))
#     }, error=function(err){
# #browser()
#       print(paste0("Error: ", err$message))
#     })
#   }
#   return (traitDF)
# }
  
winsorize <- function(traitDF,trim,firstPHENOVar,lastPHENOVar) {
  #winsorize ScenarioDF
  foreach(i = firstPHENOVar:lastPHENOVar, .combine = cbind, .verbose=FALSE) %do% {
    tryCatch({
      traitDF[,i]<-winsor(traitDF[,i],trim)
    }, error=function(err){
      print(paste0("Error: ", err$message))
#      i = i+1
    })
  }
  return (traitDF)
}

getGoTerms <- function(geneSymbols) {
  if (!is_empty(geneSymbols)) {
    # also look for "DOSE" database from slides 21.4. 17:35 to find affected diseases on the general levelBiocManager::install("GOSim")
    BiocManager::install("GOSim")
    library(GOSim)
    library(biomaRt)
    biomaRt::g
    GO_tbl <- getGO(organism = "Homo sapiens", 
                    genes    = c("AT1G06090", "AT1G06100"),
                    filters  = "ensembl_gene_id")
    getGoTerms
    library("GO.db")
    
    sG <- sample(keys(GO.db, "GOID"), 8)
    
    gT <- getGOTerm(sG)
    gP <- getGOParents(sG)
    gC <- getGOChildren(sG)
    gcat <- getGOOntology(sG)
    
  }
}

traitDF <- function(sessionVariables) {
  trait = sessionVariables$trait$trait
  return (sessionVariables$traitsDFLong[,c("ID_Kind", "gender", trait)])
}



geneSymbols <- function(currentProbeID) {
  geneSymbols = str_split (annotation[annotation$name == currentProbeID,]$gene.symbol, " ")
  #  geneSymbols = str_split (subset(annotation, annotation$name == currentProbeID), " ")
  return (geneSymbols)
}

loadResultFile<-function(sessionVariables,numberResults = 100){
  trait = sessionVariables$trait$trait
  if(!is.na(as.numeric(substr(trait,1,1)))) {
    trait = paste0("X",trait)
  }
#  PHENO = addXToName(PHENO,firstPHENOVar,lastPHENOVar)
  
  folder = sessionVariables$folder
  fileName <- paste0(folder,trait,".csv")
  if (debugMode == TRUE) {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", nrows = 10000)
  }
  else {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t")
  }
#  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N","Outlying","Skewed","Clumpy","Sparse","Striated","Convex","Skinny","Stringy","Monotonic","scagnosticsScore2"))
#  all.results<-all.results[,1:16]
  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N"))
  all.results<-all.results[,1:7]
  all.results <- dplyr::left_join(all.results, annotation, by = c("probeID" = "name"))
  all.results<-all.results[all.results$chromosome!="chrY",]
  all.results<-all.results[all.results$chromosome!="chrX",]
  all.results$mLog10FDR<-log10(all.results$FDR)*-1
  all.results$mLog10P_VAL = log10(all.results$P_VAL) * -1
  all.results<-all.results[order(all.results$mLog10P_VAL),]
  #duplicated(all.results$probeID)
  rownames(all.results)<-all.results$probeID
  #in case DeltaMeth does not match BETA
  all.results$DeltaMeth[(all.results$BETA < 0 & all.results$DeltaMeth > 0)] <- all.results$DeltaMeth*-1
  return(all.results)
}

getResultDataSingleTrait <- function(sessionVariables, onlySignificant = FALSE, numberResults = 100){
  id <- showNotification("Reading data...", duration = NULL, closeButton = FALSE)
  on.exit(removeNotification(id), add = TRUE)
  trait = sessionVariables$trait$trait
  if (isTruthy(trait)) {
#    dat <- as.data.table(loadResultFile(trait,numberResults))
    dat <- as.data.table(loadResultFile(sessionVariables,numberResults))
    rownames(dat) <- rownames(dat)
    dat = dat[,1:7]
#wenn gefühlt zu wenig Fälle zurückkommen, dann liegt das am Filtern auf Signifikanz
#    if (onlySignificant == TRUE) {
      dat = dat[dat$P_VAL <= 0.05,]
#      dat = dat[dat$P_VAL <= 0.01,]
#    }
    dat$DeltaMeth = round(dat$DeltaMeth, 5)
    dat <- addLinkToEWASDataHub(dat)
    dat <- addLinkToMRCEWASCatalog(dat)
    return (dat)
  }
}

addLinkToEWASDataHub <- function(df){
  #provide link to EWAS data hub in the form: https://bigd.big.ac.cn/ewas/datahub/probe/cg16867657
  df = mutate(df, probeID = stringr::str_replace_all(probeID, ' ', '%20'),
            EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', config$baseURL_EWASDataHub, probeID, '>', probeID,'</a>' ))
return(df)
}

addLinkToMRCEWASCatalog <- function(df){
  #provide link to EWAS data hub in the form: https://bigd.big.ac.cn/ewas/datahub/probe/cg16867657
  df = mutate(df, probeID = stringr::str_replace_all(probeID, ' ', '%20'),
              MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', config$baseURL_MRCEWASCatalog, probeID, '>', probeID,'</a>' ))
  return(df)
}

addLinkToEWASDataHubToHeader <- function(df) {
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', config$baseURL_EWASDataHub, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = EWASDataHub
    }
  }
  return(df)
}

addLinkToMRCEWASCatalogToHeader <- function(df) {
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', config$baseURL_MRCEWASCatalog, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = MRCEWASCatalog
    }
  }
  return(df)
}

removeNonVariableProbes<-function(df,NonVariableProbesList){
  `%notin%` <- Negate(`%in%`)
  df = df[rownames(df) %notin% NonVariableProbesList,]
  return(df)
}

removeMultiModelCpGFromBeta<-function(df,multiModList){
  #row.name to column
  df$CpGName<-row.names(df)
  #merge
  df <- dplyr::inner_join(df, multiModList, by = c("CpGName" = "CpG"))
  row.names(df)<-df$CpGName
  df$CpGName<-NULL
  #select only CpG with NumModes=1
  df<-df[df$NumModes<2,]
  df$NumModes<-NULL
  df$NormalP<-NULL
  return(df)
}

removeOutliers3IQR<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

getlistOfResultsDF <- function(folder) {
  if (dir.exists(folder)) {
    temp <- list.files(path=folder,pattern="*.csv")
    result = list()
    for (i in 1:length(temp)) {
      exposure <- str_sub(temp[i], 1, str_length(temp[i])-4)
      firstlines <- read.table(file = paste0(as.character(exposure),".csv"),sep = "\t", header = T,nrows = 5)
      if (colnames(firstlines)[1] == "probeID") {
        if (nrow(firstlines) >= 5) {
          if (grepl("adj", temp[i], fixed = TRUE) == TRUE) {
            #read results into DF
            resultDF = loadResultFile(exposure)
            #omit unneccesary variables
            resultDF = resultDF[,c("probeID","P_VAL","DeltaMeth")]
            resultDF = list(resultDF)
            names(resultDF) = exposure
            result = append(result,resultDF)
          }
        }
      }
    }
    saveRDS (result,file="listOfResultsDF.RDS")
    return (result)
  }
}

getExposuresWithSummary <- function(directory){
#  if (isTruthy(directory)) {
  if (dir.exists(directory)) {
    traits <- data.frame(Exposure="ex", MaxN = 1, MinP_VAL = 0, MinFDR = 0, MaxBETA = 1, MaxDeltaMeth = 1, MinDeltaMeth = 1, MaxOutlying = 1, MinOutlying = 0, MaxSkewed = 1, MinSkewed = 0, MaxClumpy = 1, MinClumpy = 0, MaxSparse = 1, MinSparse = 0, MaxStriated = 1, MinStriated = 0, MaxConvex = 1, MinConvex = 0, MaxSkinny = 1, MinSkinny = 0, MaxStringy = 1, MinStringy = 0, MaxMonotonic = 1, MinMonotonic = 0, MaxscagnosticsScore2 = 1, MinscagnosticsScore2 = 0)
    tr = traits
    temp <- list.files(path=directory,pattern="*.csv")
    for (i in 1:length(temp)) {
  #  for (i in 1:2) {
      firstlines <- read.table(file = as.character(temp[i]), sep = "\t", header = T, nrows = 5)
      if (colnames(firstlines)[1] == "probeID") {
        if (nrow(firstlines) >= 5) {
          if (grepl("adj", temp[i], fixed = TRUE) == FALSE) {
            fileName <- str_sub(temp[i], 1, str_length(temp[i])-4)
            tr$Exposure = fileName
            fileName <- paste0(fileName,".csv")
            if (debugMode == TRUE) {
              all.results <- fread(fileName,stringsAsFactors = FALSE,header = TRUE, sep = "\t", nrows = 1000)
            }
            else {
              all.results <- fread(fileName,stringsAsFactors = FALSE,header = TRUE, sep = "\t")
            }
            all.results = all.results[order(all.results$N,decreasing = TRUE),]
            tr$MaxN = all.results$N[1]
            all.results = all.results[order(all.results$P_VAL),]
            tr$MinP_VAL = all.results$P_VAL[1]
            all.results = all.results[order(all.results$FDR),]
            tr$MinFDR = all.results$FDR[1]
            all.results = all.results[order(all.results$BETA),]
            tr$MaxBETA = all.results$BETA[1]
            all.results = all.results[order(all.results$DeltaMeth),]
            tr$MaxDeltaMeth = max(all.results$DeltaMeth)
            tr$MinDeltaMeth = min(all.results$DeltaMeth)
            all.results = all.results[order(all.results$Outlying),]
            tr$MaxOutlying = max(all.results$Outlying)
            tr$MinOutlying = min(all.results$Outlying)
            all.results = all.results[order(all.results$Skewed),]
            tr$MaxSkewed = max(all.results$Skewed)
            tr$MinSkewed = min(all.results$Skewed)
            all.results = all.results[order(all.results$Clumpy),]
            tr$MaxClumpy = max(all.results$Clumpy)
            tr$MinClumpy = min(all.results$Clumpy)
            all.results = all.results[order(all.results$Sparse),]
            tr$MaxSparse = max(all.results$Sparse)
            tr$MinSparse = min(all.results$Sparse)
            all.results = all.results[order(all.results$Striated),]
            tr$MaxStriated = max(all.results$Striated)
            tr$MinStriated = min(all.results$Striated)
            all.results = all.results[order(all.results$Convex),]
            tr$MaxConvex = max(all.results$Convex)
            tr$MinConvex = min(all.results$Convex)
            all.results = all.results[order(all.results$Skinny),]
            tr$MaxSkinny = max(all.results$Skinny)
            tr$MinSkinny = min(all.results$Skinny)
            all.results = all.results[order(all.results$Stringy),]
            tr$MaxStringy = max(all.results$Stringy)
            tr$MinStringy = min(all.results$Stringy)
            all.results = all.results[order(all.results$Monotonic),]
            tr$MaxMonotonic = max(all.results$Monotonic)
            tr$MinMonotonic = min(all.results$Monotonic)
            all.results = all.results[order(all.results$scagnosticsScore2),]
            tr$MaxscagnosticsScore2 = max(all.results$scagnosticsScore2)
            tr$MinscagnosticsScore2 = min(all.results$scagnosticsScore2)
            traits = rbind(traits, tr)
    #        all.results$Exposure = Ex$Exposure
    # #        if (is_empty(g.all.results)) {
    #         if (!exists("g.all.results")) {
    #           g.all.results <- data.frame()
    #         }
    #         g.all.results <- rbind(g.all.results,all.results)
    #         assign("g.all.results",g.all.results,envir=globalenv())
          }
        }
      }
    }
    #remove dummy exposure
    traits <- traits[-c(1), ] 
    return(traits)
  }
}

getDMRs <- function(exposure){
  if (istruthy(exposure)) {
    tryCatch({

      fileNameSelect = paste0(exposure,"_dmrff_dmrs_select",".rda")
      load(fileNameSelect) # into df dmrs_select originally from dmrff.cohort
      DMRs = dmrs_select
      DMRs$Exposure = exposure
    }, error=function(err){
      print(paste0("unable find DMR for ",exposure))
    });
    return(DMRs)
  }
}

getDMRSummary <- function(directory){
#  if (isTruthy(directory)) {
  if (dir.exists(directory)) {
    #template for DMR df
    DMR <- data.frame()#(DMR="DMR")
    DMRs = DMR
    filterSubstring = "_dmrff_dmrs_select_sites"
    temp <- list.files(path=directory,pattern="*.rda")
    for (i in 1:length(temp)) {
      if (grepl(filterSubstring, temp[i], fixed = TRUE) == TRUE) {
        Exposure <- str_sub(temp[i], 1, str_length(temp[i]) - str_length(filterSubstring) - str_length(".rda"))
        DMRs = rbind(DMRs,getDMRs(Exposure))
      }
    }
    #remove dummy exposure
    DMRs <- DMRs[-c(1), ]
    return(DMRs)
  }
}

reducedAnnotation <- function(){
    a = annotation
    a$type = NULL
    a$target = NULL
    a$meth.dye = NULL
    a$chromosome = as.factor(a$chromosome)
    levels(a$chromosome)[levels(a$chromosome)=="chrX"] <-"chr23" #only for sorting
    levels(a$chromosome)[levels(a$chromosome)=="chrY"] <-"chr24" #only for sorting
    a$chromosomeNum = as.factor(as.numeric(gsub("chr","",a$chromosome)))
    a = a[order(a$chromosomeNum,a$position),]
    a$globalPosition <- seq_len(nrow(a))
    
    # minPosChr1 = min(a[which(a$chromosomeNum == 1),]$position)
    # maxPosChr1 = max(a[which(a$chromosomeNum == 1),]$position)
    # lenChr1 = maxPosChr1-minPosChr1

return (a)
}

resultDataSingleScenarioWithAnnotation <- function(df){
  a = reducedAnnotation()
  a$gene.symbolShort = str_sub(a$gene.symbol, 1, 20) #NULL
  a$gene.accession = NULL
  a$gene.region = NULL
  a$cpg.island.name = NULL
  a$relation.to.island = NULL
  a$snp.exclude = NULL
  return (dplyr::left_join(df, a, by = c("probeID" = "name")))
}

resultDataSingleScenarioWithAnnotationEWAScatalogCount <- function(df){
  return (dplyr::left_join(df, EWAScatalogCount, by = c("probeID" = "CpG")))
}

empty_plot <- function(title = NULL){
  plot <- plotly_empty(type = "scatter", mode = "markers") %>%
    config(
      displayModeBar = FALSE
    ) %>%
    layout(
      title = list(
        text = title,
        yref = "paper",
        y = 0.5
      )
    )
  return(plot)
} 
