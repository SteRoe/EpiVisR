library(shiny)
library(shinyBS)

inputTrait_UI <- function(id){
  ns <- NS(id)
  
#  tagList(
    bsCollapse(id = ns("collapseExample"), open = c(ns("Folder"),ns("Traits")), multiple = TRUE,
      bsCollapsePanel(ns("Folder"),
        fluidRow(
          column(6, #style = "margin-top: 25px;", #align="center", id="txtinput",
                 shiny::textInput(ns("directory"), "Folder for Traits", dataDir, placeholder = TRUE)
          ),
          column(6, style = "margin-top: 25px;", #align="center", id="button",
#                 shiny::actionButton(ns("btnReadFolder"), label = "Select Folder")
                  shiny::verbatimTextOutput(ns("txtDataFile"), placeholder = TRUE)
          )
        ),
        bsCollapsePanel(ns("Traits"),
          tabsetPanel(
            tabPanel("Short Select",
              fluidRow(
                column(6,
                  shiny::selectInput(ns("selSelectedTrait"),"trait","",selectize = FALSE)
                ),
                column(6, style = "margin-top: 25px;",
                  shiny::actionButton(ns("btnSelectTraitShort"), label = "Select Trait")
                )
              )
            ),
            tabPanel("Long Select",
              DT::dataTableOutput(ns("traits"))
#              shiny::actionButton(ns("btnSelectTraitLong"), label = "Select Trait")
            )
          ),
          fluidRow(
            column(6, tags$html("last selected trait"),
                   verbatimTextOutput(ns("txtSelectedTrait"), placeholder = TRUE)
            ),
            column(6, tags$html("selected traits"),
                   verbatimTextOutput(ns("txtSelectedTraits"), placeholder = TRUE)
            )
          )
        ), label = "")
      )
#  )
}

inputTrait_SERVER <- function(id, sessionVariables) {
#browser()
  moduleServer(id, function(input, output, session) {
      selectedTrait <- reactiveValues(trait = "", traits = list())
      moduleVariables <- reactiveValues(traitsWithSummary = data.frame())

      observeEvent(input$btnReadFolder, { #ignoreInit = TRUE, 
#browser()
        sessionVariables$folder = input$directory
        updateTraitsTable(session, output, sessionVariables, moduleVariables)
      }, ignoreNULL = FALSE)
      
      observeEvent(input$directory, { #ignoreInit = TRUE,
#browser()
        sessionVariables$folder = input$directory
        message(paste0(Sys.time(), " reading folder ", sessionVariables$folder))
        updateTraitsTable(session, output, sessionVariables, moduleVariables)
      })
      
      observeEvent(input$traits_cell_clicked, {
        if (isTruthy(input$traits_row_last_clicked)) {
          selectedTrait$trait = moduleVariables$traitsWithSummary[input$traits_row_last_clicked,]$trait
          if (substr(selectedTrait$trait, 1, 1) == "X") {
            selectedTrait$trait<-gsub("X","",selectedTrait$trait) # remove "X" in traitname for compatibility to filenames
          }
          selectedTrait$traits = moduleVariables$traitsWithSummary[input$traits_rows_selected,]$trait
          sessionVariables$trait <- selectedTrait
          shiny::validate(shiny::need(selectedTrait$trait,"Select at least one trait."))
          output$txtSelectedTrait = renderText({sessionVariables$trait$trait})
          output$txtSelectedTraits = renderText({sessionVariables$trait$traits})
        }
      })
      
      observeEvent(input$btnSelectTraitShort, {
        if(isTruthy(input$selSelectedTrait)) {
#browser()
          selectedTrait$trait = input$selSelectedTrait
          sessionVariables$trait <- selectedTrait
          output$txtSelectedTrait = renderText({sessionVariables$trait$trait})
        }
      })
      return(selectedTrait)
    }
  )
}

updateTraitsTable <- function(session, output, sessionVariables, moduleVariables) {
  id <- shiny::showNotification("Filling data table...", duration = NULL, closeButton = FALSE)
  on.exit(removeNotification(id), add = TRUE)
#browser()
  output$txtDataFile = shiny::renderText(sessionVariables$dataFileName)
  sessionVariables$traitsDFLong = getTraitsDFLong(sessionVariables)
  tryCatch({
    traits <- getTraits(sessionVariables$folder)
    }, error=function(err){
      print(paste0("unable to read folder ", sessionVariables$folder))
  });
  updateSelectInput(
    session = session,
    inputId = "selSelectedTrait",
    choices = traits
  )
  
  moduleVariables$traitsWithSummary <- getTraitsWithSummary(sessionVariables$folder)
  output$traits <- DT::renderDataTable({
    id <- showNotification("Printing data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    datatable(moduleVariables$traitsWithSummary, extensions = list("Scroller"), style = "bootstrap", class = "compact", width = "100%",
              options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE))
  }, server = FALSE)
}

getTraits<-function(directory){
  if (dir.exists(directory)) {
    tryCatch({
      #  testthat::expect_is(df, "data.frame")
#      testthat::expect_is(directory, "character")
      #  testthat::expect_equal(length(column), 1)
  
      fileNames=list()
      temp <- list.files(path=directory,pattern="*.csv")
      for (i in 1:length(temp)) {
        firstlines <- read.table(file = paste0(directory,as.character(temp[i])), sep="\t", header = T, nrows = 5)
        if (colnames(firstlines)[1] == "probeID") {
          if (nrow(firstlines) >= 5) {
            fileName <- str_sub(temp[i], 1, str_length(temp[i])-4)
            #        fileNames<-c(fileNames,fileName)
            fileNames<-append(fileNames,fileName)
          }
        }
      }
      return(fileNames)
      }, error=function(err){
#browser()
        print(paste0("unable to load ",directory))
      return(NULL)
    });
  }
}

getTraitsWithSummary <- function(directory){
  #template for exposure df
#  if (isTruthy(directory)) {
  message(paste0(Sys.time(), " start reading traits"))
  if (dir.exists(directory)) {
    traits <- data.frame(trait="ex", MaxN = 1, MinP_VAL = 0, MinFDR = 0, MaxBETA = 1, MaxDeltaMeth = 1, MinDeltaMeth = 1, MaxOutlying = 1, MinOutlying = 0, MaxSkewed = 1, MinSkewed = 0, MaxClumpy = 1, MinClumpy = 0, MaxSparse = 1, MinSparse = 0, MaxStriated = 1, MinStriated = 0, MaxConvex = 1, MinConvex = 0, MaxSkinny = 1, MinSkinny = 0, MaxStringy = 1, MinStringy = 0, MaxMonotonic = 1, MinMonotonic = 0)
    tr = traits
    message(paste0(Sys.time(), " start loading folder"))
    tryCatch({
      temp <- list.files(path=directory,pattern="*.csv")
    }, error=function(err){
      print(paste0("unable to load '",directory,"'"))
    });
    tryCatch({
      if (debugMode == TRUE) {
        max <- min(50:length(temp))  
      }
      else {
        max <- length(temp)
      }
      fileNameTraitsSummary = "traitsSummary.RDS" #check for summary file
      RDSfileName = paste0(directory,fileNameTraitsSummary)
      if (file_test("-f", RDSfileName) == TRUE) {
        
        message(paste0(Sys.time(), " reading RDS", RDSfileName))
        traits = readRDS (file = RDSfileName)
      }
      else {
        for (i in 1:max) {
          firstlines <- read.table(file = paste0(directory,as.character(temp[i])),sep = "\t", header = T,nrows = 5)
          if (colnames(firstlines)[1] == "probeID") {
            if (nrow(firstlines) >= 5) {
              if (grepl("adj", temp[i], fixed = TRUE) == FALSE) {

                fileName <- str_sub(temp[i], 1, str_length(temp[i])-4)
                tr$trait = fileName
                fileName <- paste0(directory,fileName,".csv")
                message(paste0(Sys.time(), " reading trait file ", fileName))
                if (deepDebugMode == TRUE) {
                  all.results <- fread(fileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t", nrows = 1000)
                }
                else {
                  all.results <- fread(fileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
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
                traits = rbind(traits, tr)
              }
            }
          }
        }
        #remove dummy exposure
        traits <- traits[-c(1), ]
        #save DF for faster later access
        if (debugMode == FALSE) {
          message(paste0(Sys.time(), "writing RDS", RDSfileName))
          saveRDS(traits, file = RDSfileName)
        }
      }
      message(paste0(Sys.time(), " finished reading traits"))
      return(traits)
      
    }, error=function(err){
      print(paste0("unable to load ",fileName, " err: ", err))
    });
  }
}
