#' @import data.table
#' @import foreach

utils::globalVariables(c("globalVariables", "debugMode"))

#' very first function
#' @description very first function during package load
#' @importFrom magrittr "%>%"
#' @param libname library name
#' @param pkgname package name
#' @return nothing
#' .onAttach()
.onAttach <- function(libname, pkgname) {
  globalVariables <- list()
  packageStartupMessage("start loading package")
  base::loadNamespace("EpiVisR")
  base::packageStartupMessage("end loading package")
}

# very last function
# @description very last function during package unload
# @param libPath library path
# @return nothing
# .onUnload <- function(libPath) {
#   if (base::isNamespaceLoaded("EpiVisR")) {
#     base::unloadNamespace("EpiVisR")
#   }
#   base::packageStartupMessage("end unloading package")
# }

#' starts the app
#' @description function to start the app
#' @export
EpiVisRApp <- function() {
  shiny::shinyApp(ui, server)
}

# checkconfigVariables <- function(globalVariables) {
#   if (is.null(globalVariables$config$firstPHENOVar)) {
#     stop("firstPHENOVar missing")
#   }
#   if (is.null(globalVariables$config$lastPHENOVar)) {
#     stop("lastPHENOVar missing")
#   }
#   if (is.null(globalVariables$config$AdjustFileName)) {
#     stop("AdjustFileName missing")
#   }
#   if (is.null(globalVariables$config$betaFileName)) {
#     stop("betaFileName missing")
#   }
#   if (is.null(globalVariables$config$EWAScatalogFileName)) {
#     stop("EWAScatalogFileName missing")
#   }
#   if (is.null(globalVariables$config$MultiModProbesFileName)) {
#     stop("multiModProbesFileName missing")
#   }
#   if (is.null(globalVariables$config$baseDir)) {
#     stop("baseDir not set")
#   }
#   if (is.null(globalVariables$config$winsorTrim)) {
#     stop("winsorTrim not set")
#   }
#   if (is.null(globalVariables$config$mergeAttribut)) {
#     stop("mergeAttribut not set")
#   }
#   if (is.null(globalVariables$config$genderAttribut)) {
#     stop("genderAttribut not set")
#   }
#   if (is.null(globalVariables$config$genderFemaleValue)) {
#     stop("genderFemaleValue not set")
#   }
#   if (is.null(globalVariables$config$genderMaleValue)) {
#     stop("genderMaleValue not set")
#   }
#   if (is.null(globalVariables$config$debugMode)) {
#     stop("debugMode not set")
#   }
# }

getTraitsDFLong <- function(globalVariables) {
  if (dir.exists(globalVariables$config$dataDir)) {
    # tryCatch({
    #   setwd(globalVariables$folder)
    # }, error=function(err){
    #   message(paste0("unable to open work folder ", globalVariables$folder))
    #   return()
    # });
    # tryCatch({
    # }, error=function(err){
    #   #      print(paste0("unable open work folder ", dataDir))
    #   errortext = paste0("unable open read config.yml in ", globalVariables$folder)
    #   message(errortext)
    #   id <- shiny::showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    # });
    tryCatch({
#browser()
      traitFileName <- globalVariables$config$traitFileName
      sessionVariables$dataFileName = traitFileName
      traitDFLong<-fread(file=traitFileName,sep="\t",dec=".",header=TRUE, data.table = FALSE)
      traitDFLong<-as.data.frame(traitDFLong)
      rownames(traitDFLong) = traitDFLong[,globalVariables$config$mergeAttribut]
      colnames(traitDFLong)<-gsub(" ",".",colnames(traitDFLong)) # replace " " with "." for compatibility to filenames
      colnames(traitDFLong)<-gsub("-",".",colnames(traitDFLong)) # replace "-" with "." for compatibility to filenames
      genderFileName <- globalVariables$config$genderFileName
      Gender<-fread(file=genderFileName, sep="\t", dec=".", data.table = FALSE)
      Gender<-base::subset(Gender, select=c(globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut))
#      Gender<-Gender[, c(globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut)]
      Gender = as.data.frame(Gender)
      traitDFLong = base::merge(traitDFLong, Gender, by.x = globalVariables$config$mergeAttribut, by.y = globalVariables$config$mergeAttribut, all.x = FALSE, all.y=FALSE)
      # firstPHENOVar <- globalVariables$config$firstPHENOVar + ncol(P) - 1
      # lastPHENOVar <- globalVariables$config$lastPHENOVar + ncol(P) - 1
#      traitDFLongWinsorized = winsorize(traitDFLong,globalVariables$config$winsorTrim,firstPHENOVar,lastPHENOVar)
      # traitDFLongWinsorized = winsorize(traitDFLong,globalVariables$config$winsorTrim,1,ncol(traitDFLong))
      # traitDFLong = traitDFLongWinsorized
      return (traitDFLong)
    }, error=function(err){
#      print(paste0("unable open work folder ", dataDir))
      errortext = paste0("unable open work folder ", globalVariables$folder)
      message(errortext)
      id <- shiny::showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    });
  }
}

#' loadObjects
#' loads globally needed objects (methylation matrix with beta values, annotation, etc.)
#' @param globalVariables contains all global available Objects
#' @return globalVariables
# examples EpiVisR::loadObjects(globalVariables)
loadObjects <- function(globalVariables){
  if (dir.exists(globalVariables$config$dataDir)) {
    print(paste0(Sys.time(), " load beta."))
    betaFileName <- globalVariables$config$betaFileName
#browser()
    if (globalVariables$config$debugMode == FALSE) {
      beta <- data.table::fread(betaFileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", data.table = FALSE)
    }
    else {
      beta <- data.table::fread(betaFileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", nrows = 1000, data.table = FALSE)
    }
    beta <- as.data.frame(beta)
#    beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
    rownames(beta) <- beta$PROBEID
    beta$PROBEID <- NULL
    # message(paste0(Sys.time(), " remove outliers from beta."))
    # beta_wo_outliers<-removeOutliers3IQR(as.matrix(beta))
    # beta_wo_outliers<-as.data.frame(beta_wo_outliers[[1]])
    # beta<-beta_wo_outliers

    #load list of multimodal CpG
    print(paste0(Sys.time(), " load multimodal probes."))
    MultiModProbesFileName <- globalVariables$config$MultiModProbesFileName
    MultiModProbes<-fread(file=MultiModProbesFileName, sep="\t", dec=".", data.table = FALSE)
#    assign("MultiModProbes",MultiModProbes,envir=globalenv())
    beta <- removeMultiModelCpGFromBeta(beta,MultiModProbes)

    # message(paste0(Sys.time(), " loading non-variable probeIDs."))
    # nonVariableProbesFileName <- pkg.env$config$nonVariableProbesFileName
    # nonVariableProbes <- fread(file = nonVariableProbesFileName,sep = "\t",dec = ".")
    # nonVariableProbes <- readRDS(file = nonVariableProbesFileName)
    # assign("nonVariableProbes",nonVariableProbes,envir=globalenv())
    # message(paste0(Sys.time(), " removing non-variable probeIDs."))
    # beta <- removeNonVariableProbes(beta,nonVariableProbes)
    globalVariables$beta = beta
#    assign("beta",beta,envir=globalenv())
    print(paste0(Sys.time(), " transposing beta."))
    beta.t<-t(beta)
    rownames(beta)
    colnames(beta)
    colnames(beta.t) <- rownames(beta)
    colnames(beta.t)
    rownames(beta.t)
    globalVariables$beta.t = beta.t
#    assign("beta.t",beta.t,envir=globalenv())

    print(paste0(Sys.time(), " load annotation."))
    annotation <- meffil::meffil.get.features("450k")
    annotation$relation.to.island = as.factor(annotation$relation.to.island)
    #remove unmeasured or multimodal probeIDs from annotation
    annotation = annotation[which(annotation$name %in% rownames(beta)),]
    globalVariables$annotation = annotation
#    assign("annotation",annotation,envir=globalenv())
    print(paste0(Sys.time(), " load EWAS catalog."))
    EWAScatalogFileName <- globalVariables$config$EWAScatalogFileName
    EWAScatalog = fread(file=EWAScatalogFileName, sep = "\t", data.table = FALSE)
    globalVariables$EWAScatalog = EWAScatalog
#    assign("EWAScatalog",EWAScatalog,envir=globalenv())
#    usethis::use_pipe()
    EWAScatalogCount = EWAScatalog %>% dplyr::group_by(CpG) %>% dplyr::tally(!is.na(CpG))
#    EWAScatalogCount = EWAScatalog %>% dplyr::group_by(EWAScatalog$CpG) %>% dplyr::tally(!is.na(EWAScatalog$CpG))
    globalVariables$EWAScatalogCount = EWAScatalogCount
#    assign("EWAScatalogCount",EWAScatalogCount,envir=globalenv())
  }
  return(globalVariables)
}
#
# setBaseSciPen<-function(){
#   options(scipen=-10, digits=5)
# }
#
# resetSciPen<-function(){
#   options(scipen=0, digits=7)
# }

#' winsorize
#' performs winsorizing
#' @param traitDF data.frame to be used
#' @param trim value to be used for winsorizing
#' @param startVar variable to start with
#' @param endVar variable to end at
#' @return data.frame with winsorized variables
#examples winsorize(df, 0.05, 10, 20)
winsorize <- function(traitDF, trim, startVar, endVar) {
  i <- NULL
  foreach(i = startVar:endVar, .combine = cbind, .verbose=FALSE) %do% {
    tryCatch({
      if (is.numeric(traitDF[,i])){
        traitDF[,i]<-psych::winsor(traitDF[,i],trim)
      }
    }, error=function(err){
      print(paste0("Error: ", err$message))
    })
  }
  return (traitDF)
}

# getGoTerms <- function(geneSymbols) {
#   if (!is_empty(geneSymbols)) {
#     # also look for "DOSE" database from slides 21.4. 17:35 to find affected diseases on the general levelBiocManager::install("GOSim")
#     BiocManager::install("GOSim")
#     library(GOSim)
#     library(biomaRt)
#     biomaRt::g
#     GO_tbl <- getGO(organism = "Homo sapiens",
#                     genes    = c("AT1G06090", "AT1G06100"),
#                     filters  = "ensembl_gene_id")
#     getGoTerms
#     library("GO.db")
#
#     sG <- sample(keys(GO.db, "GOID"), 8)
#
#     gT <- getGOTerm(sG)
#     gP <- getGOParents(sG)
#     gC <- getGOChildren(sG)
#     gcat <- getGOOntology(sG)
#
#   }
# }

traitDF <- function(sessionVariables) {
  trait = sessionVariables$trait$trait
  df = sessionVariables$traitsDFLong[,c(globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut, trait)]
  return (df)
}

# geneSymbols <- function(currentProbeID) {
#   geneSymbols = str_split (annotation[annotation$name == currentProbeID,]$gene.symbol, " ")
#   #  geneSymbols = str_split (subset(annotation, annotation$name == currentProbeID), " ")
#   return (geneSymbols)
# }

#' loadResultFile
#' @param globalVariables globalVariables loads files with globally available information (beta, trait)
#' @param sessionVariables sessionVariables
#' @return data.frame with regression results
# examples loadResultFile(globalVariables, 100)
loadResultFile<-function(globalVariables, sessionVariables){
  trait = sessionVariables$trait$trait
  if(!is.na(as.numeric(substr(trait,1,1)))) {
    trait = paste0("X",trait)
  }
#  PHENO = addXToName(PHENO,firstPHENOVar,lastPHENOVar)
  folder = sessionVariables$folder
  fileName <- paste0(folder,trait,".csv")
  if (globalVariables$config$debugMode == TRUE) {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", nrows = 10000, data.table = FALSE)
  }
  else {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", data.table = FALSE)
  }
#  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N","Outlying","Skewed","Clumpy","Sparse","Striated","Convex","Skinny","Stringy","Monotonic","scagnosticsScore2"))
#  all.results<-all.results[,1:16]
  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N"))
  all.results<-all.results[,1:7]
#  all.results <- dplyr::left_join(all.results, globalVariables$annotation, by = c("probeID" = "name"))
#  all.results <- base::merge(all.results, annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = TRUE)
  all.results <- base::merge(all.results, globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y = FALSE)
  all.results <- stats::na.omit(all.results)
  all.results<-all.results[all.results$chromosome!="chrY",]
  all.results<-all.results[all.results$chromosome!="chrX",]
  all.results$mLog10FDR<-log10(all.results$FDR)*-1
  all.results$mLog10P_VAL = log10(all.results$P_VAL) * -1
  all.results<-all.results[order(all.results$mLog10P_VAL),]
  #duplicated(all.results$probeID)
  rownames(all.results)<-all.results$probeID
  #in case DeltaMeth does not match BETA
#  all.results$DeltaMeth[(all.results$BETA < 0 & all.results$DeltaMeth > 0)] <- all.results$DeltaMeth*-1
  return(all.results)
}

getResultDataSingleTrait <- function(globalVariables, sessionVariables, onlySignificant = FALSE) {
  id <- shiny::showNotification("Reading data...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(id), add = TRUE)
  trait = sessionVariables$trait$trait
  if (shiny::isTruthy(trait)) {
    dat <- loadResultFile(globalVariables, sessionVariables)
#    rownames(dat) <- rownames(dat)
    dat = dat[,1:7]
# if there are too few cases, then filtering for significant values is the reason
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
  #provide link to EWAS data hub
  # df = dplyr::mutate(df, probeID = stringr::str_replace_all(df$probeID, ' ', '%20'),
  #           EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_EWASDataHub, probeID, '>', df$probeID,'</a>' ))
  df$EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_EWASDataHub, df$probeID, '>', df$probeID,'</a>' )
  return(df)
}

addLinkToMRCEWASCatalog <- function(df){
  #provide link to MRC EWAS catalog
  # df = dplyr::mutate(df, probeID = stringr::str_replace_all(df$probeID, ' ', '%20'),
  #             MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_MRCEWASCatalog, probeID, '>', df$probeID,'</a>' ))
  df$MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_MRCEWASCatalog, df$probeID, '>', df$probeID,'</a>' )
  return(df)
}

addLinkToEWASDataHubToHeader <- function(df) {
  i <- NULL
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_EWASDataHub, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = EWASDataHub
    }
  }
  return(df)
}

addLinkToMRCEWASCatalogToHeader <- function(df) {
  i <- NULL
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_MRCEWASCatalog, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = MRCEWASCatalog
    }
  }
  return(df)
}

# removeNonVariableProbes<-function(df,NonVariableProbesList){
#   `%notin%` <- Negate(`%in%`)
#   df = df[rownames(df) %notin% NonVariableProbesList,]
#   return(df)
# }

removeMultiModelCpGFromBeta<-function(df,multiModList){
  #row.name to column
  df$CpGName<-rownames(df)
  #merge
  df <- base::merge(df, multiModList, by.x= "CpGName", by.y = "CpG")
  #select only CpG with NumModes=1
  df<-df[df$NumModes<2,]
  df$NumModes<-NULL
  df$NormalP<-NULL
  rownames(df)<-df$CpGName
  df$CpGName<-NULL
  return(df)
}

# removeOutliers3IQR<-function(probes){
#   require(matrixStats)
#   if(nrow(probes) < ncol(probes)) {
#     warning("expect probes are rows.")
#   }
#   rowIQR <- rowIQRs(probes, na.rm = T)
#   row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
#   maskL <- probes < row2575[,1] - 3 * rowIQR
#   maskU <- probes > row2575[,2] + 3 * rowIQR
#   initial_NAs<-rowSums(is.na(probes))
#   probes[maskL] <- NA
#   removed_lower <- rowSums(is.na(probes))-initial_NAs
#   probes[maskU] <- NA
#   removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
#   N_for_probe<-rowSums(!is.na(probes))
#   Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
#   return(list(probes, Log))
# }

# getlistOfResultsDF <- function(folder) {
#   if (dir.exists(folder)) {
#     temp <- list.files(path=folder,pattern="\\.csv$")
#     result = list()
#     for (i in 1:length(temp)) {
#       exposure <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i])-4)
# #      firstlines <- utils::read.table(file = paste0(as.character(exposure),".csv"),sep = "\t", header = T,nrows = 5)
#       firstlines <- data.table::fread(file = paste0(as.character(exposure),".csv"),sep = "\t", header = T,nrows = 5)
#       if (colnames(firstlines)[1] == "probeID") {
#         if (nrow(firstlines) >= 5) {
#           if (grepl("adj", temp[i], fixed = TRUE) == TRUE) {
#             #read results into DF
# #            resultDF = loadResultFile(exposure)
#             resultDF = loadResultFile(globalVariables = globalVariables)
#             #omit unneccesary variables
#             resultDF = resultDF[,c("probeID","P_VAL","DeltaMeth")]
#             resultDF = list(resultDF)
#             names(resultDF) = exposure
#             result = append(result,resultDF)
#           }
#         }
#       }
#     }
#     saveRDS (result,file="listOfResultsDF.RDS")
#     return (result)
#   }
# }
#
# getExposuresWithSummary <- function(globalVariables, directory){
#   if (dir.exists(directory)) {
#     traits <- data.frame(Exposure="ex", MaxN = 1, MinP_VAL = 0, MinFDR = 0, MaxBETA = 1, MaxDeltaMeth = 1, MinDeltaMeth = 1, MaxOutlying = 1, MinOutlying = 0, MaxSkewed = 1, MinSkewed = 0, MaxClumpy = 1, MinClumpy = 0, MaxSparse = 1, MinSparse = 0, MaxStriated = 1, MinStriated = 0, MaxConvex = 1, MinConvex = 0, MaxSkinny = 1, MinSkinny = 0, MaxStringy = 1, MinStringy = 0, MaxMonotonic = 1, MinMonotonic = 0, MaxscagnosticsScore2 = 1, MinscagnosticsScore2 = 0)
#     tr = traits
#     temp <- list.files(path=directory,pattern="\\.csv$")
#     for (i in 1:length(temp)) {
#   #  for (i in 1:2) {
# #      firstlines <- utils::read.table(file = as.character(temp[i]), sep = "\t", header = T, nrows = 5)
#       firstlines <- data.table::fread(file = as.character(temp[i]), sep = "\t", header = T, nrows = 5)
#       if (colnames(firstlines)[1] == "probeID") {
#         if (nrow(firstlines) >= 5) {
#           if (grepl("adj", temp[i], fixed = TRUE) == FALSE) {
#             fileName <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i])-4)
#             tr$Exposure = fileName
#             fileName <- paste0(fileName,".csv")
#             if (globalVariables$config$debugMode == TRUE) {
#               all.results <- fread(fileName,stringsAsFactors = FALSE,header = TRUE, sep = "\t", nrows = 1000)
#             }
#             else {
#               all.results <- fread(fileName,stringsAsFactors = FALSE,header = TRUE, sep = "\t")
#             }
#             all.results = all.results[order(all.results$N,decreasing = TRUE),]
#             tr$MaxN = all.results$N[1]
#             all.results = all.results[order(all.results$P_VAL),]
#             tr$MinP_VAL = all.results$P_VAL[1]
#             all.results = all.results[order(all.results$FDR),]
#             tr$MinFDR = all.results$FDR[1]
#             all.results = all.results[order(all.results$BETA),]
#             tr$MaxBETA = all.results$BETA[1]
#             all.results = all.results[order(all.results$DeltaMeth),]
#             tr$MaxDeltaMeth = max(all.results$DeltaMeth)
#             tr$MinDeltaMeth = min(all.results$DeltaMeth)
#             all.results = all.results[order(all.results$Outlying),]
#             tr$MaxOutlying = max(all.results$Outlying)
#             tr$MinOutlying = min(all.results$Outlying)
#             all.results = all.results[order(all.results$Skewed),]
#             tr$MaxSkewed = max(all.results$Skewed)
#             tr$MinSkewed = min(all.results$Skewed)
#             all.results = all.results[order(all.results$Clumpy),]
#             tr$MaxClumpy = max(all.results$Clumpy)
#             tr$MinClumpy = min(all.results$Clumpy)
#             all.results = all.results[order(all.results$Sparse),]
#             tr$MaxSparse = max(all.results$Sparse)
#             tr$MinSparse = min(all.results$Sparse)
#             all.results = all.results[order(all.results$Striated),]
#             tr$MaxStriated = max(all.results$Striated)
#             tr$MinStriated = min(all.results$Striated)
#             all.results = all.results[order(all.results$Convex),]
#             tr$MaxConvex = max(all.results$Convex)
#             tr$MinConvex = min(all.results$Convex)
#             all.results = all.results[order(all.results$Skinny),]
#             tr$MaxSkinny = max(all.results$Skinny)
#             tr$MinSkinny = min(all.results$Skinny)
#             all.results = all.results[order(all.results$Stringy),]
#             tr$MaxStringy = max(all.results$Stringy)
#             tr$MinStringy = min(all.results$Stringy)
#             all.results = all.results[order(all.results$Monotonic),]
#             tr$MaxMonotonic = max(all.results$Monotonic)
#             tr$MinMonotonic = min(all.results$Monotonic)
#             all.results = all.results[order(all.results$scagnosticsScore2),]
#             tr$MaxscagnosticsScore2 = max(all.results$scagnosticsScore2)
#             tr$MinscagnosticsScore2 = min(all.results$scagnosticsScore2)
#             traits = rbind(traits, tr)
#           }
#         }
#       }
#     }
#     #remove dummy exposure
#     traits <- traits[-c(1), ]
#     return(traits)
#   }
# }

# getDMRs <- function(exposure){
#   if (shiny::isTruthy(exposure)) {
#     tryCatch({
#       fileNameSelect = paste0(exposure,"_dmrff_dmrs_select",".rda")
#       load(fileNameSelect) # into df dmrs_select originally from dmrff.cohort
# browser() # check reason for dmrs_select
#       DMRs = dmrs_select
#       DMRs$Exposure = exposure
#     }, error=function(err){
#       print(paste0("unable find DMR for ",exposure))
#     });
#     return(DMRs)
#   }
# }

# getDMRSummary <- function(directory){
#   if (dir.exists(directory)) {
#     #template for DMR df
#     DMR <- data.frame()#(DMR="DMR")
#     DMRs = DMR
#     filterSubstring = "_dmrff_dmrs_select_sites"
#     temp <- list.files(path=directory,pattern="*.rda")
#     for (i in 1:length(temp)) {
#       if (grepl(filterSubstring, temp[i], fixed = TRUE) == TRUE) {
#         Exposure <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i]) - stringr::str_length(filterSubstring) - stringr::str_length(".rda"))
#         DMRs = rbind(DMRs,getDMRs(Exposure))
#       }
#     }
#     #remove dummy exposure
#     DMRs <- DMRs[-c(1), ]
#     return(DMRs)
#   }
# }

reducedAnnotation <- function(globalVariables){
#browser()
  a = globalVariables$annotation
  a$type = NULL
  a$target = NULL
  a$meth.dye = NULL
  a$chromosome = as.factor(a$chromosome)
  levels(a$chromosome)[levels(a$chromosome)=="chrX"] <-"chr23" #only for sorting
  levels(a$chromosome)[levels(a$chromosome)=="chrY"] <-"chr24" #only for sorting
  a$chromosomeNum = as.factor(as.numeric(gsub("chr","",a$chromosome)))
  a = a[order(a$chromosomeNum,a$position),]
  a$globalPosition <- seq_len(nrow(a))

  return (a)
}

resultDataSingleScenarioWithAnnotation <- function(globalVariables, df){
  a = reducedAnnotation(globalVariables)
  a$gene.symbolShort = stringr::str_sub(a$gene.symbol, 1, 20) #NULL
  a$gene.accession = NULL
  a$gene.region = NULL
  a$cpg.island.name = NULL
  a$relation.to.island = NULL
  a$snp.exclude = NULL
#browser()
  df = dplyr::left_join(df, a, by = c("probeID" = "name"))
#  df = base::merge(df, a, by.x = "probeID", by.y = "name", all.x = TRUE, all.y=FALSE)
  return (df)
}

resultDataSingleScenarioWithAnnotationEWAScatalogCount <- function(globalVariables, df){
#browser()
  df = dplyr::left_join(df, globalVariables$EWAScatalogCount, by = c("probeID" = "CpG"))
# df = base::merge(df, globalVariables$EWAScatalogCount, by.x = "probeID", by.y = "CpG", all.x = TRUE, all.y=FALSE)
  return (df)
}

empty_plot <- function(title = NULL){
  plot <- plotly::plotly_empty(type = "scatter", mode = "markers") %>%
    plotly::config(
      displayModeBar = FALSE
    ) %>%
    plotly::layout(
      title = list(
        text = title,
        yref = "paper",
        y = 0.5
      )
    )
  return(plot)
}
