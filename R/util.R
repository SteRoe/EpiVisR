#' @import data.table
#' @import foreach

utils::globalVariables(c("globalVariables", "debugMode"))

#' very first function
#' @description very first function during package load
#' @importFrom magrittr "%>%"
#' @param libname library name
#' @param pkgname package name
#' @return nothing
#' @keywords internal
#' @noRd
#' .onAttach()
.onAttach <- function(libname, pkgname) {
  globalVariables <- list()
  base::packageStartupMessage("start loading package")
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

#' getTraitsDFLong
#' gets back the currently selected trait together with gender information and traitnames were replaced with filename compatible characters
#' @param globalVariables contains all global available Objects
#' @return data.frame
#' @keywords internal
#' @noRd
# examples EpiVisR::getTraitsDFLong(globalVariables)
getTraitsDFLong <- function(globalVariables) {
  if (dir.exists(globalVariables$config$dataDir)) {
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

#' loadObjects
#' loads globally needed objects (methylation matrix with beta values, annotation, etc.)
#' @param globalVariables contains all global available Objects
#' @return globalVariables
#' @keywords internal
#' @noRd
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

#' winsorize
#' performs winsorizing
#' @param traitDF data.frame to be used
#' @param trim value to be used for winsorizing
#' @param startVar variable to start with
#' @param endVar variable to end at
#' @return data.frame with winsorized variables
#' @keywords internal
#' @noRd
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
#' traitDF
#' @param sessionVariables sessionVariables
#' @param mergeAttribut mergeAttribut
#' @param genderAttribut genderAttribut
#' @return data.frame with trait data, merge attribute and gender attribute
#' @keywords internal
#' @noRd
#examples traitDF(sessionVariables, "ID", "gender")
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
#examples loadResultFile(globalVariables, 100)
loadResultFile<-function(globalVariables, sessionVariables){
  trait = sessionVariables$trait$trait
  if(!is.na(as.numeric(substr(trait,1,1)))) {
    trait = paste0("X",trait)
  }
#  PHENO = addXToName(PHENO,firstPHENOVar,lastPHENOVar)
#  folder = sessionVariables$folder
  folder = globalVariables$config$dataDir
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

#' getResultDataSingleTrait
#' gets back the currently selected trait
#' @param globalVariables contains all global available Objects
#' @param sessionVariables contains all session objects
#' @param significeBorder border for selecting cases
#' @return data.frame
#' @keywords internal
#' @noRd
# examples getResultDataSingleTrait(globalVariables, sessionVariables, onlySignificant)
getResultDataSingleTrait <- function(globalVariables, sessionVariables, significeBorder = 0.05) {
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
    dat = dat[dat$P_VAL <= significeBorder,]
#      dat = dat[dat$P_VAL <= 0.01,]
#    }
    dat$DeltaMeth = round(dat$DeltaMeth, 5)
    dat <- addLinkToEWASDataHub(dat, globalVariables$config$baseURL_EWASDataHub)
    dat <- addLinkToMRCEWASCatalog(dat, globalVariables$config$baseURL_MRCEWASCatalog)
    return (dat)
  }
}

#' addLinkToEWASDataHub
#' adds links to EWASDataHub to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToEWASDataHub(data.frame, baseURL)
addLinkToEWASDataHub <- function(df, baseURL){
  #provide link to EWAS data hub
  # df = dplyr::mutate(df, probeID = stringr::str_replace_all(df$probeID, ' ', '%20'),
  #           EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_EWASDataHub, probeID, '>', df$probeID,'</a>' ))
  df$EWASDataHub = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df$probeID,'</a>' )
  return(df)
}

#' addLinkToMRCEWASCatalog
#' adds links to MRC EWAS catalog to a data.frame as separate column
#' @param df data.frame to which links should be added
#' @param baseURL string describing link to be included
#' @return data.frame
#' @keywords internal
#' @noRd
# examples addLinkToMRCEWASCatalog(data.frame)
addLinkToMRCEWASCatalog <- function(df, baseURL){
  #provide link to MRC EWAS catalog
  # df = dplyr::mutate(df, probeID = stringr::str_replace_all(df$probeID, ' ', '%20'),
  #             MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', globalVariables$config$baseURL_MRCEWASCatalog, probeID, '>', df$probeID,'</a>' ))
  df$MRCEWASCatalog = paste0('<a target=_blank rel="noopener noreferrer" href=', baseURL, df$probeID, '>', df$probeID,'</a>' )
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

# removeNonVariableProbes<-function(df,NonVariableProbesList){
#   `%notin%` <- Negate(`%in%`)
#   df = df[rownames(df) %notin% NonVariableProbesList,]
#   return(df)
# }

#' removeMultiModelCpGFromBeta
#' removes multimodal CpG from data.frame
#' @param df data.frame to which links should be added
#' @param multiModList list with multimodal CpG
#' @return data.frame
#' @keywords internal
#' @noRd
# examples removeMultiModelCpGFromBeta(data.frame, multiModList)
removeMultiModelCpGFromBeta<-function(df, multiModList){
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
#browser()
  df = dplyr::left_join(df, a, by = c("probeID" = "name"))
#  df = base::merge(df, a, by.x = "probeID", by.y = "name", all.x = TRUE, all.y=FALSE)
  return (df)
}

#' resultDataSingleScenarioWithAnnotationEWAScatalogCount
#' @param globalVariables
#' @param df data.frame to which links should be added
#' @return data.frame
#' @keywords internal
#' @noRd
resultDataSingleScenarioWithAnnotationEWAScatalogCount <- function(globalVariables, df){
#browser()
  df = dplyr::left_join(df, globalVariables$EWAScatalogCount, by = c("probeID" = "CpG"))
# df = base::merge(df, globalVariables$EWAScatalogCount, by.x = "probeID", by.y = "CpG", all.x = TRUE, all.y=FALSE)
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
