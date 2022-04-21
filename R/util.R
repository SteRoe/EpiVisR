#' @import data.table
#' @import foreach

utils::globalVariables(c("globalVariables", "debugMode"))

#' very first function
#' @description very first function during package load
#' @importFrom magrittr "%>%"
#'
#' @param libname library name
#' @param pkgname package name
#'
#' @return nothing
#' @keywords internal
#' .onAttach()
.onAttach <- function(libname, pkgname) {
  globalVariables <- list()
  base::packageStartupMessage("start loading package")
  base::loadNamespace("EpiVisR")
  base::packageStartupMessage("end loading package")
}

# very last function
# description very last function during package unload
# param libPath library path
# return nothing
# keywords internal
# noRd
# .onUnload <- function(libPath) {
#   globalVariables <- NULL
#   if (base::isNamespaceLoaded("EpiVisR")) {
#     base::unloadNamespace("EpiVisR")
#   }
#   base::packageStartupMessage("end unloading package")
# }

#' Starts the App
#' @description Function to start the App. Details on how to work with this interactive package are given in the package vignette.
#' @export
EpiVisRApp <- function() {
  shiny::shinyApp(ui, server)
}

#' gets back the currently selected trait together with gender information and traitnames were replaced with filename compatible characters
#'
#' @param globalVariables contains all global available Objects
#' @return data.frame
#' @keywords internal
#' @noRd
getTraitsDFLong <- function(globalVariables, sessionVariables) {
#  if (dir.exists(globalVariables$config$dataDir)) {
  if (dir.exists(sessionVariables$folder)) {
    tryCatch({
      traitFileName = globalVariables$config$traitFileName #sessionVariables$dataFileName
      traitDFLong<-fread(file=traitFileName,sep="\t",dec=".",header=TRUE, data.table = FALSE)
      traitDFLong<-as.data.frame(traitDFLong)
      rownames(traitDFLong) = traitDFLong[,globalVariables$config$mergeAttribut]
      colnames(traitDFLong)<-gsub(" ",".",colnames(traitDFLong)) # replace " " with "." for compatibility to filenames
      colnames(traitDFLong)<-gsub("-",".",colnames(traitDFLong)) # replace "-" with "." for compatibility to filenames
    }, error=function(err){
      errortext = paste0("unable to open trait file ", globalVariables$config$traitFileName)
      message(errortext)
      id <- shiny::showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    });
    tryCatch({
      genderFileName <- globalVariables$config$genderFileName
      Gender<-fread(file=genderFileName, sep="\t", dec=".", data.table = FALSE)
      Gender<-base::subset(Gender, select=c(globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut))
      Gender = as.data.frame(Gender)
      traitDFLong = base::merge(traitDFLong, Gender, by.x = globalVariables$config$mergeAttribut, by.y = globalVariables$config$mergeAttribut, all.x = FALSE, all.y=FALSE)
      return (traitDFLong)
    }, error=function(err){
      errortext = paste0("unable to open gender file ", globalVariables$config$genderFileName)
      message(errortext)
      id <- shiny::showNotification(errortext, duration = NULL, type = "error", closeButton = TRUE)
    });
  }
}

addPackagePathToConfig <- function(config, packagePath){
  packagePath = paste0(packagePath, "/")
  if (config$debugModeLocalInst == TRUE) {
    packagePath = paste0(packagePath, "inst/")
  }
  if (base::startsWith(config$betaFileName, "./inst/")) {
    config$betaFileName = stringr::str_replace(config$betaFileName, "./inst/", packagePath)
  }
  if (base::startsWith(config$MultiModProbesFileName, "./inst/")) {
    config$MultiModProbesFileName = stringr::str_replace(config$MultiModProbesFileName, "./inst/", packagePath)
  }
  if (base::startsWith(config$traitFileName, "./inst/")) {
    config$traitFileName = stringr::str_replace(config$traitFileName, "./inst/", packagePath)
  }
  if (base::startsWith(config$genderFileName, "./inst/")) {
    config$genderFileName = stringr::str_replace(config$genderFileName, "./inst/", packagePath)
  }
  if (base::startsWith(config$dataDir, "./inst/")) {
    config$dataDir = stringr::str_replace(config$dataDir, "./inst/", packagePath)
  }
  if (base::startsWith(config$EWAScatalogFileName, "./inst/")) {
    config$EWAScatalogFileName = stringr::str_replace(config$EWAScatalogFileName, "./inst/", packagePath)
  }
  return(config)
}

#' loadObjects
#' loads globally needed objects (methylation matrix with beta values, annotation, etc.)
#' @param globalVariables contains all global available Objects
#' @return globalVariables
#' @keywords internal
#' @noRd
loadObjects <- function(globalVariables){
  if (dir.exists(globalVariables$config$dataDir)) {
#  if (dir.exists(sessionVariables$folder)) {
    print(paste0(Sys.time(), " load beta."))
    betaFileName <- globalVariables$config$betaFileName
    if (globalVariables$config$debugMode == FALSE) {
      beta <- data.table::fread(betaFileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", data.table = FALSE)
    }
    else {
      beta <- data.table::fread(betaFileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", nrows = 1000, data.table = FALSE)
    }
    beta <- as.data.frame(beta)
#    beta<-data.frame(column_to_rownames(beta, var = "PROBEID"))
#    rownames(beta) <- beta$probeID
    rownames(beta) <- beta[,globalVariables$config$probeAttribut]
    beta[,globalVariables$config$probeAttribut] <-NULL
    #load list of multimodal CpG
    print(paste0(Sys.time(), " load multimodal probes."))
    MultiModProbesFileName <- globalVariables$config$MultiModProbesFileName
    MultiModProbes<-fread(file=MultiModProbesFileName, sep="\t", dec=".", data.table = FALSE)
    if (!nrow(removeMultiModeCpGFromBeta(beta,MultiModProbes)) == 0) {
      beta <- removeMultiModeCpGFromBeta(beta,MultiModProbes)
    }
    globalVariables$beta = beta
    print(paste0(Sys.time(), " transposing beta."))
    beta.t<-t(beta)
    #rownames(beta) #check
    #colnames(beta) #check
    colnames(beta.t) <- rownames(beta)
    #colnames(beta.t) #check
    #rownames(beta.t) #check
    globalVariables$beta.t = beta.t
    print(paste0(Sys.time(), " load annotation."))
    annotation <- meffil::meffil.get.features("450k")
    annotation$relation.to.island = as.factor(annotation$relation.to.island)
    #remove unmeasured or multimodal probeIDs from annotation
    #annotation = annotation[which(annotation$name %in% rownames(beta)),]
    globalVariables$annotation = annotation
    print(paste0(Sys.time(), " Load EWAS catalog."))
    EWAScatalogFileName <- globalVariables$config$EWAScatalogFileName
    if (base::file.exists(EWAScatalogFileName)) {
      EWAScatalog = data.table::fread(file=EWAScatalogFileName, sep = "\t", data.table = FALSE)
    }
    else {
      EWAScatalog = readTxtGzFromURL(URL=EWAScatalogFileName)
    }
    globalVariables$EWAScatalog = EWAScatalog
    print(paste0(Sys.time(), " Calculating EWASCatalogCount."))
    EWAScatalogCount = data.frame(table(EWAScatalog$CpG))
    colnames(EWAScatalogCount)[1] = "CpG"
    colnames(EWAScatalogCount)[2] = "n"
    globalVariables$EWAScatalogCount = EWAScatalogCount
    print(paste0(Sys.time(), " Done calculating EWASCatalogCount."))
  }
  return(globalVariables)
}

#' readTxtGzFromURL
#' reads structured txt file from URL
#' @param URL URL to be used
#' @return data.frame with contents from URL
#' @keywords internal
#' @noRd
readTxtGzFromURL <- function(URL) {
  con <- base::gzcon(base::url(URL))
  txt <- base::readLines(con)
  return(utils::read.csv(textConnection(txt),sep="\t",dec="."))
}

#' winsorize
#' performs winsorizing
#' @param traitDF data.frame to be used
#' @param trim value to be used for winsorizing
#' @param startVar variable to start with
#' @param endVar variable to end at
#' @return data.frame with winsorized variables
#' @examples
#' EpiVisR::winsorize(df, 0.05, 10, 20)
#' @export
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
#' traitDF
#' @param sessionVariables sessionVariables
#' @param mergeAttribut mergeAttribut
#' @param genderAttribut genderAttribut
#' @return data.frame with trait data, merge attribute and gender attribute
#' @keywords internal
#' @noRd
traitDF <- function(sessionVariables, mergeAttribut, genderAttribut) {
  trait = sessionVariables$trait$trait
  df = sessionVariables$traitsDFLong[,c(mergeAttribut, genderAttribut, trait)]
  return (df)
}

#' loadResultFile
#' @param globalVariables globalVariables loads files with globally available information (beta, trait)
#' @param sessionVariables sessionVariables
#' @return data.frame with regression results
#' @keywords internal
#' @noRd
loadResultFile<-function(globalVariables, sessionVariables){
  trait = sessionVariables$trait$trait
  if(!is.na(as.numeric(substr(trait,1,1)))) {
    trait = paste0("X",trait)
  }
#  PHENO = addXToName(PHENO,firstPHENOVar,lastPHENOVar)
  folder = sessionVariables$folder
#  folder = globalVariables$config$dataDir
  fileName <- paste0(folder,trait,".csv")
  if (globalVariables$config$debugMode == TRUE) {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", nrows = 10000, data.table = FALSE)
  }
  else {
    all.results <- fread(fileName, stringsAsFactors=FALSE, header=TRUE, sep="\t", data.table = FALSE)
  }
  all.results<-setcolorder(all.results, c("probeID","BETA","SE", "P_VAL", "FDR","DeltaMeth","N"))
  all.results<-all.results[,1:7]
  all.results <- base::merge(all.results, globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = FALSE, all.y = FALSE) #was all.x = TRUE, all.y = FALSE)
  all.results <- stats::na.omit(all.results)
  all.results<-all.results[all.results$chromosome!="chrY",]
  all.results<-all.results[all.results$chromosome!="chrX",]
  all.results$mLog10FDR<-log10(all.results$FDR)*-1
  all.results$mLog10P_VAL = log10(all.results$P_VAL) * -1
  all.results<-all.results[order(all.results$mLog10P_VAL),]
  #duplicated(all.results$probeID)
  #rownames(all.results)<-all.results$probeID
  rownames(all.results)<-all.results[,globalVariables$config$probeAttribut]
  #in case DeltaMeth does not match BETA
#  all.results$DeltaMeth[(all.results$BETA < 0 & all.results$DeltaMeth > 0)] <- all.results$DeltaMeth*-1
  return(all.results)
}

#' getResultDataSingleTrait
#' gets back the currently selected trait
#' @param globalVariables contains all global available Objects
#' @param sessionVariables contains all session objects
#' @param significanceBorder border for selecting cases
#' @return data.frame
#' @keywords internal
#' @noRd
# examples getResultDataSingleTrait(globalVariables, sessionVariables, onlySignificant)
getResultDataSingleTrait <- function(globalVariables, sessionVariables, significanceBorder = 0.05) {
  id <- shiny::showNotification("Reading data...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(id), add = TRUE)
  trait = sessionVariables$trait$trait
  if (shiny::isTruthy(trait)) {
#    dat <- NULL
    dat <- loadResultFile(globalVariables, sessionVariables)
#    rownames(dat) <- rownames(dat)
    dat = dat[,1:7]
# if there are too few cases, then filtering for significant values is the reason
#    dat = dat[dat$P_VAL <= 0.05,]
    if (!nrow(dat[dat$P_VAL <= significanceBorder,]) < 1)
    {
      dat = dat[dat$P_VAL <= significanceBorder,]
#      dat = dat[dat$P_VAL <= 0.01,]
    }
    dat$DeltaMeth = round(dat$DeltaMeth, 5)
    dat <- addLinkToEWASDataHub(dat, globalVariables$config$baseURL_EWASDataHub, globalVariables$config$probeAttribut)
    dat <- addLinkToMRCEWASCatalog(dat, globalVariables$config$baseURL_MRCEWASCatalog, globalVariables$config$probeAttribut)
    return (dat)
  }
}

#' addLinkToEWASDataHub
#' adds links to EWASDataHub to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @param probeAttribut string describing the name of the probe variable
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToEWASDataHub(data.frame, baseURL)
addLinkToEWASDataHub <- function(df, baseURL, probeAttribut){
  #provide link to EWAS data hub
  df$EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df[,probeAttribut],'</a>' )
  return(df)
}

#' addLinkToMRCEWASCatalog
#' adds links to MRC EWAS catalog to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @param probeAttribut string describing the name of the probe variable
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToMRCEWASCatalog(data.frame)
addLinkToMRCEWASCatalog <- function(df, baseURL, probeAttribut){
  #provide link to MRC EWAS catalog
  df$MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df[,probeAttribut],'</a>' )
  return(df)
}

#' addLinkToEWASDataHubToHeader
#' adds links to EWAS data hub to a data.frame into first line
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToEWASDataHubToHeader(data.frame)
addLinkToEWASDataHubToHeader <- function(df, baseURL) {
  i <- NULL
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = EWASDataHub
    }
  }
  return(df)
}

#' addLinkToMRCEWASCatalogToHeader
#' adds links to MRC EWAS catalog to a data.frame into first line
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToMRCEWASCatalogToHeader(data.frame)
addLinkToMRCEWASCatalogToHeader <- function(df, baseURL) {
  i <- NULL
  foreach(i=1:ncol(df)) %do% {
    if (grepl("cg", colnames(df)[i], fixed = TRUE) == TRUE) {
      probeID = colnames(df)[i]
      MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, probeID, '>', probeID,'</a>' )
      colnames(df)[i] = MRCEWASCatalog
    }
  }
  return(df)
}

#' removeMultiModeCpGFromBeta
#' removes multimodal CpG from data.frame
#' @param df data.frame to which links should be added
#' @param multiModList list with multimodal CpG
#' @return data.frame
#' @keywords internal
#' @noRd
# examples removeMultiModeCpGFromBeta(data.frame, multiModList)
removeMultiModeCpGFromBeta<-function(df, multiModList){
  #row.name to column
  df$CpGName <- rownames(df)
  #merge
  df <- dplyr::left_join(df, multiModList, by = c("CpGName" = "CpG"))
  #replace NA with 0 in NumModes
  df <- df %>% dplyr::mutate(NumModes = replace(df$NumModes, is.na(df$NumModes), 0))
  #select only CpG with NumModes=1
  df <- df[df$NumModes<2,]
  df$NumModes <- NULL
  df$NormalP <- NULL
  rownames(df) <- df$CpGName
  df$CpGName <- NULL
  return(df)
}

#' reducedAnnotation
#' expects a data.frame with annotation for Illuminas 450k array like provided by meffil package
#' removes chromosome X and Y from annotation data.frame and sorts it by chromosome number and position
#' removes attributes type, target and meth.dye from data.frame
#' @param a data.frame to which links should be added
#' @return data.frame
#' @keywords internal
#' @noRd
# examples reducedAnnotation(data.frame)
#reducedAnnotation <- function(globalVariables){
reducedAnnotation <- function(a){
#  a = globalVariables$annotation
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

#' resultDataSingleScenarioWithAnnotation
#' expects annotation and data.frame
#' removes long stromgs from gene.symbol
#' removes attributes accession, region, cpg.island.name, relation.to.island and snp.exclude from data.frame
#' @param annotation data.frame to which links should be added
#' @param df data.frame to which links should be added
#' @return data.frame
#' @keywords internal
#' @noRd
# examples resultDataSingleScenarioWithAnnotation(data.frame)
resultDataSingleScenarioWithAnnotation <- function(annotation, df){
#  a = reducedAnnotation(globalVariables)
  a = reducedAnnotation(annotation)
  a$gene.symbolShort = stringr::str_sub(a$gene.symbol, 1, 20) #NULL
  a$gene.accession = NULL
  a$gene.region = NULL
  a$cpg.island.name = NULL
  a$relation.to.island = NULL
  a$snp.exclude = NULL
  df = dplyr::left_join(df, a, by = c("probeID" = "name"))
  return (df)
}

#' resultDataSingleScenarioWithAnnotationEWAScatalogCount
#' @param globalVariables globalVariables
#' @param df data.frame to which links should be added
#' @return data.frame
#' @keywords internal
#' @noRd
resultDataSingleScenarioWithAnnotationEWAScatalogCount <- function(globalVariables, df){
  df = dplyr::left_join(df, globalVariables$EWAScatalogCount, by = c("probeID" = "CpG"))
  df$n[is.na(df$n)] = 1
  return (df)
}

#' empty_plot
#' @param title title for empty plot
#' @return plot empty plot
#' @keywords internal
#' @noRd
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

replaceBackslashes <- function(directory) {
  print(paste0(Sys.time(), " replace \\ in folder name ", directory, "."))
  directory = gsub("\\\\","/",directory) # replace "\" with "/"
  print(paste0(Sys.time(), " end with /."))
  if (!grepl("/$",directory)) { #does not end with /
    directory = paste0(directory,"/")
  }
  return(directory)
}
