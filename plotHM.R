plotHM_UI <- function(id){
  ns <- NS(id)
  tagList(
###because <InteractiveComplexHeatmap> does not support modules, we have to do this in server.R
###nevertheless, we refer to the functions defined in plotHM.R

     #shiny::actionButton(ns("btnplotHeatMap"), label = "Plot HeatMap")
     # fluidRow(
     #   column(width = 12,
     #          InteractiveComplexHeatmapOutput(ns("heatmap_1"),
     #                                          height1 = 800, height2 = 800),
     #   )
     # )
  )
}

plotHM_SERVER <- function(id, sessionVariables) {
  moduleServer(id, function(input, output, session) {
    moduleVariables <- reactiveValues(numberResults = 10, df = data.frame())

    observeEvent(input$btnplotHeatMap, ignoreInit = TRUE, {
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, interactiveHeatMapP_Val(), "heatmap_1")
    }, ignoreNULL = FALSE)
  })
}

getResultDF <- function(listOfResultDF,P_D) {
  #    if (missing(P_D)) P_D = "P"
  merged_P_VAL = data.frame()
  merged_DeltaMeth = data.frame()
  foreach(i=1:length(listOfResultDF)) %do% {
    DF = listOfResultDF[i]
    DF = as.data.frame(DF)
    colnames(DF) = c("probeID","P_VAL","DeltaMeth")
    #browser()
    DF_P_VAL = DF
    DF_DeltaMeth = DF
    DF_P_VAL$DeltaMeth = NULL
    DF_DeltaMeth$P_VAL = NULL
    if (i == 1) {
      merged_P_VAL = DF_P_VAL
      merged_DeltaMeth = DF_DeltaMeth
      rownames(merged_P_VAL) = DF_P_VAL$probeID
      rownames(merged_DeltaMeth) = DF_DeltaMeth$probeID
      #        colnames(merged_P_VAL)[i+1] = names(listOfResultDF)[i]
      #        colnames(merged_DeltaMeth)[i+1] = names(listOfResultDF)[i]
    }
    else {
      #merge
      merged_P_VAL = dplyr::inner_join(merged_P_VAL, DF_P_VAL, by = c("probeID" = "probeID"))
      merged_DeltaMeth = dplyr::inner_join(merged_DeltaMeth, DF_DeltaMeth, by = c("probeID" = "probeID"))
    }
    colnames(merged_P_VAL)[i+1] = names(listOfResultDF)[i]
    colnames(merged_DeltaMeth)[i+1] = names(listOfResultDF)[i]
  }
  rownames(merged_P_VAL) = merged_P_VAL$probeID
  merged_P_VAL$probeID = NULL
  rownames(merged_DeltaMeth) = merged_DeltaMeth$probeID
  merged_DeltaMeth$probeID = NULL
  result = list(merged_P_VAL,merged_DeltaMeth)
  return (result)
}

interactiveHeatMapP_Val <- function() {
browser()
    fileNameHM = "ComplexHeatmapP_VAL.RDS"
    if (file_test("-f", paste0(sessionVariables$folder,fileNameHM)) == TRUE && debugMode != TRUE) {
      complexhmP_Val = readRDS (file=paste0(sessionVariables$folder,fileNameHM))
    }
    else {
      #load or generate DF with P_Val from all results merged together
      fileNameLPV = "listResultP_Val_DeltaMeth.RDS"
      if (file_test("-f", paste0(sessionVariables$folder,fileNameLPV)) == TRUE) {
        listResultP_Val_DeltaMeth = readRDS (file=paste0(sessionVariables$folder,fileNameLPV))
      }
      else {
        #load or generate list of DFs with results
        fileNameLR = "listOfResultsDF.RDS"
        if (file_test("-f", paste0(sessionVariables$folder,fileNameLR)) == TRUE) {
          listOfResultDF = readRDS (file=paste0(sessionVariables$folder,fileNameLR))
        }
        else {
          #        listOfResultDF = getlistOfResultsDF(dataDir())
          listOfResultDF = getlistOfResultsDF(sessionVariables$folder)
          saveRDS(listOfResultDF, file = paste0(sessionVariables$folder,fileNameLR))
        }
        listResultP_Val_DeltaMeth = getResultDF (listOfResultDF)
        saveRDS(listResultP_Val_DeltaMeth, file = paste0(sessionVariables$folder,fileNameLPV))
      }

      ResultDFP_Val = as.data.frame(listResultP_Val_DeltaMeth [1])
      #      ResultDFDeltaMeth = as.data.frame(listResultP_Val_DeltaMeth [2])
      ResultDFP_Val = as.matrix(ResultDFP_Val)
      #      ResultDFDeltaMeth = as.matrix(ResultDFDeltaMeth)

      #choose fewer CpG
      if (debugMode == TRUE) {
        ResultDFP_Val = ResultDFP_Val[sample(1:nrow(ResultDFP_Val), 100),]
        #        ResultDFDeltaMeth=ResultDFDeltaMeth[sample(1:nrow(ResultDFDeltaMeth), 10),]
      }
      else {
        ResultDFP_Val = ResultDFP_Val[sample(1:nrow(ResultDFP_Val), 10000),]
        #        ResultDFDeltaMeth=ResultDFDeltaMeth[sample(1:nrow(ResultDFDeltaMeth), 10000),]
      }
      # #omit nas
      # ResultDFP_Val[ResultDFP_Val > 0.05] <- NA
      browser()
    #this works:
      ComplexHeatmap::Heatmap(ResultDFP_Val, name = "heatmap_1", use_raster = TRUE)
    #this not:
    #  complexhmP_Val = ComplexHeatmap::Heatmap(ResultDFP_Val, name = "heatmap_1", use_raster = TRUE)
    #  ComplexHeatmap::Heatmap(as.matrix(mtcars), name = "heatmap_1", use_raster = TRUE)
      # if (debugMode != TRUE) {
      #   saveRDS(complexhmP_Val,file = paste0(dataDir, fileNameHM))
      # }
    }
}

interactiveHeatMapP_Val_BAK <- function() {
  ComplexHeatmap::Heatmap(as.matrix(mtcars), name = "heatmap_1", use_raster = TRUE)
  browser()
#   fileNameHM = "ComplexHeatmapP_VAL.RDS"
#   if (file_test("-f", paste0(dataDir,fileNameHM)) == TRUE && debugMode != TRUE) {
#     complexhmP_Val = readRDS (file=paste0(dataDir,fileNameHM))
#   }
#   else {
#     #load or generate DF with P_Val from all results merged together
#     fileNameLPV = "listResultP_Val_DeltaMeth.RDS"
#     if (file_test("-f", paste0(dataDir,fileNameLPV)) == TRUE) {
#       listResultP_Val_DeltaMeth = readRDS (file=paste0(dataDir,fileNameLPV))  
#     }
#     else {
#       #load or generate list of DFs with results
#       fileNameLR = "listOfResultsDF.RDS"
#       if (file_test("-f", paste0(dataDir,fileNameLR)) == TRUE) {
#         listOfResultDF = readRDS (file=paste0(dataDir,fileNameLR))
#       }
#       else {
#         #        listOfResultDF = getlistOfResultsDF(dataDir())
#         listOfResultDF = getlistOfResultsDF(dataDir)
#         saveRDS(listOfResultDF, file = paste0(dataDir,fileNameLR))
#       }
#       listResultP_Val_DeltaMeth = getResultDF (listOfResultDF)
#       saveRDS(listResultP_Val_DeltaMeth, file = paste0(dataDir,fileNameLPV))
#     }
#     
#     ResultDFP_Val = as.data.frame(listResultP_Val_DeltaMeth [1])
#     #      ResultDFDeltaMeth = as.data.frame(listResultP_Val_DeltaMeth [2])
#     ResultDFP_Val = as.matrix(ResultDFP_Val)
#     #      ResultDFDeltaMeth = as.matrix(ResultDFDeltaMeth)
#     
#     #choose fewer CpG
#     if (debugMode == TRUE) {
#       ResultDFP_Val = ResultDFP_Val[sample(1:nrow(ResultDFP_Val), 100),]
#       #        ResultDFDeltaMeth=ResultDFDeltaMeth[sample(1:nrow(ResultDFDeltaMeth), 10),]
#     }
#     else {
#       ResultDFP_Val = ResultDFP_Val[sample(1:nrow(ResultDFP_Val), 10000),]
#       #        ResultDFDeltaMeth=ResultDFDeltaMeth[sample(1:nrow(ResultDFDeltaMeth), 10000),]
#     }
#     
#     # #omit nas
#     # ResultDFP_Val[ResultDFP_Val > 0.05] <- NA
#     #complexhmP_Val = ComplexHeatmap::Heatmap(ResultDFP_Val, name = "ResultDFP_Val", use_raster = TRUE)
# browser()
# ComplexHeatmap::Heatmap(as.matrix(mtcars), name = "heatmap_1", use_raster = TRUE)
# #    ComplexHeatmap::Heatmap(ResultDFP_Val, name = "heatmap_1", use_raster = TRUE)
# #    complexhmP_Val = draw(complexhmP_Val)
#     if (debugMode != TRUE) {
#       saveRDS(complexhmP_Val,file = paste0(dataDir, fileNameHM))
#     }
#     # complexhmDeltaMeth = ComplexHeatmap::Heatmap(ResultDFDeltaMeth,use_raster = TRUE)
#     # saveRDS(complexhmDeltaMeth,file = paste0(dataDir, "ComplexHeatmapDeltaMeth.RDS"))
#   }
# #  return (complexhmP_Val)
}