inputTrait_UI <- function(id){
  ns <- shiny::NS(id)

    shinyBS::bsCollapse(id = ns("collapse"), open = c(ns("DataFolderandFilePanel"), ns("TraitsPanel")), multiple = TRUE,
      shinyBS::bsCollapsePanel(ns("DataFolderandFilePanel"),
      shiny::fluidRow(
          shiny::column(5, style = "margin-top: 25px;",
                shiny::textInput(ns("txtDirectory"), label = "Folder for Traits", value = "", placeholder = "folder with methylation .csv-files")
          ),
          shiny::column(5, style = "margin-top: 25px;",
                shiny::textInput(ns("txtDataFile"), label = "File for Traits", value = "", placeholder = "file with trait values")
          ),
          shiny::column(2, style = "margin-top: 50px;",
                shiny::actionButton(ns("btnSelectFolderFile"), label = "Select Data Folder & File")
          )
        ),
        shinyBS::bsCollapsePanel(id = ns("TraitsPanel"),
          shiny::fluidRow(
            shiny::column(12, htmltools::tags$html(tags$body(h4('Traits for selection'))),
              DT::dataTableOutput(ns("traits"))
            )
          ),
          shiny::fluidRow(
            shiny::column(6, htmltools::tags$html("last selected trait (check only)"),
                   shiny::verbatimTextOutput(ns("txtSelectedTrait"), placeholder = TRUE)
            )
          )
        ), label = "")
      )
}

inputTrait_SERVER <- function(id, globalVariables, sessionVariables) {
  print(paste0(Sys.time(), " starting input trait SERVER."))
  shiny::moduleServer(id, function(input, output, session) {
      print(paste0(Sys.time(), " creating reactive values."))
      selectedTrait <- shiny::reactiveValues(trait = "", traits = list()) #selectedTrait$trait: last single selected trait; selectedTrait$traits: list of all selected traits
      moduleVariables <- shiny::reactiveValues(traitsWithSummary = data.frame())
      print(paste0(Sys.time(), " finished creating reactive values."))
      globalVariables$config$dataDir <- replaceBackslashes(globalVariables$config$dataDir)
      shiny::updateTextInput(session, "txtDirectory", value = globalVariables$config$dataDir)
      shiny::updateTextInput(session, "txtDataFile", value = globalVariables$config$traitFileName)

      shiny::observeEvent(input$btnSelectFolderFile, {
        print(paste0(Sys.time(), " fired directory."))
        directory <- replaceBackslashes(input$txtDirectory)
        shiny::updateTextInput(session, inputId="txtDirectory", label="Folder for Traits", value=directory, placeholder = "folder with methylation .csv-files")
        sessionVariables$folder <- input$txtDirectory
        dataFile <- replaceBackslashes(input$txtDataFile)
        shiny::updateTextInput(session, inputId="txtDataFile", label="File for Traits", value=dataFile)
        globalVariables$config$traitFileName <- input$txtDataFile
        print(paste0(Sys.time(), " reading folder ", sessionVariables$folder)) #globalVariables$config$dataDir))
        updateTraitsTable(session, output, globalVariables, sessionVariables, moduleVariables)
        print(paste0(Sys.time(), " finished: reading folder."))
      }, ignoreInit = TRUE)

      # shiny::observeEvent(input$directory, { #ignoreInit = TRUE,
      #   print(paste0(Sys.time(), " fired directory."))
      #   directory = replaceBackslashes(input$directory)
      #   shiny::updateTextInput(session, "directory", directory)
      #   sessionVariables$folder = input$directory
      #   print(paste0(Sys.time(), " reading folder ", globalVariables$config$dataDir))
      #   updateTraitsTable(session, output, globalVariables, sessionVariables, moduleVariables)
      #   print(paste0(Sys.time(), " finished reading folder."))
      # }, ignoreInit = TRUE)

      shiny::observeEvent(input$traits_cell_clicked, {
        print(paste0(Sys.time(), " fired traits_cell_clicked."))
        if (shiny::isTruthy(input$traits_row_last_clicked)) {
          print(paste0(Sys.time(), " select traits."))
          selectedTrait$trait <- moduleVariables$traitsWithSummary[input$traits_row_last_clicked,]$trait
          if (substr(selectedTrait$trait, 1, 1) == "X") {
            selectedTrait$trait<-gsub("X","",selectedTrait$trait) # remove "X" in traitname for compatibility to filenames
          }

          # if (substr(selectedTrait$trait, nchar(selectedTrait$trait)-2, nchar(selectedTrait$trait)) == "adj") {
          #   selectedTrait$trait<-gsub("adj","",selectedTrait$trait) # remove "adj" in traitname for compatibility to filenames
          # }

          selectedTrait$traits <- moduleVariables$traitsWithSummary[input$traits_rows_selected,]$trait
          sessionVariables$trait <- selectedTrait
          shiny::validate(shiny::need(selectedTrait$trait,"Select at least one trait."))
          print(paste0(Sys.time(), " before renderText."))
          output$txtSelectedTrait <- shiny::renderText({sessionVariables$trait$trait})
#          output$txtSelectedTraits = shiny::renderText({sessionVariables$trait$traits})
          print(paste0(Sys.time(), " finished: select trait."))
        }
      }, ignoreInit = TRUE)

    }
  )
}

updateTraitsTable <- function(session, output, globalVariables, sessionVariables, moduleVariables) {
  print(paste0(Sys.time(), " filling TraitsTable."))
  id <- shiny::showNotification("Filling data table...", duration = NULL, closeButton = FALSE)
  on.exit(shiny::removeNotification(id), add = TRUE)
  print(paste0(Sys.time(), " before renderText."))
  sessionVariables$traitsDFLong <- getTraitsDFLong(globalVariables, sessionVariables)
  output$txtDataFile <- shiny::renderText(globalVariables$config$traitFileName)
  print(paste0(Sys.time(), " assign traitsDFLong."))
  tryCatch({
    print(paste0(Sys.time(), " getTraits."))
    traits <- getTraits(globalVariables, sessionVariables$folder)
    }, error=function(err){
      message(Sys.time(), paste0("unable to read folder ", sessionVariables$folder))
  });
  print(paste0(Sys.time(), " updateSelectInput."))
  # shiny::updateSelectInput(
  #   session = session,
  #   inputId = "selSelectedTrait",
  #   choices = traits
  # )
  print(paste0(Sys.time(), " getTraitsWithSummary."))
  moduleVariables$traitsWithSummary <- getTraitsWithSummary(globalVariables = globalVariables, directory = sessionVariables$folder, exludeMultiModalProbes = TRUE)

  print(paste0(Sys.time(), " before renderDataTable."))
  output$traits <- DT::renderDataTable({
    id <- shiny::showNotification("printing data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
#    tableData <- moduleVariables$traitsWithSummary %>% dplyr::mutate_if(is.numeric, signif, digits = 3)
    tableData <- moduleVariables$traitsWithSummary
#    colnames(tableData) <- stringr::str_to_title(colnames(tableData))
    DT::datatable(tableData, selection = "single", editable = FALSE, extensions = list("Scroller"), style = "bootstrap", class = "compact", width = "100%",
              options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE)) %>%
    DT::formatSignif(2:ncol(tableData), digits = 2)
  }, server = FALSE)
}

#%>% DT::formatSignif(c("MinP_Val", "MinFDR", "MaxDeltaMeth", "MaxOutlying", "MinOutlying", "MaxSkewed", "MinSkewed", "MaxClumpy", "MinClumpy", "MaxSparse", "MinSparse", "axStriated", "MinStriated", "MaxConvex", "MinConvex", "MaxSkinny", "MinSkinny", "MaxStringy", "MinStringy", "MaxMonotonic", "MinMonotonic"), digits = 5)

getTraits<-function(globalVariables, directory){
  if (dir.exists(directory)) {
    tryCatch({
      testthat::expect_type(directory, "character")
      print(paste0(Sys.time(), " start getTraits()."))
      fileNames=list()
      temp <- list.files(path=directory,pattern="\\.csv$")
      for (i in 1:length(temp)) {
        firstlines <- data.table::fread(file = paste0(directory,as.character(temp[i])), sep="\t", header = T, nrows = 5, data.table = FALSE)
        testthat::expect_type(firstlines, "list")
        if (colnames(firstlines)[1] == globalVariables$config$probeAttribut) {
          if (nrow(firstlines) >= 5) {
            fileName <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i])-4)
            fileNames<-append(fileNames,fileName)
          }
        }
      }
      return(fileNames)
      }, error=function(err){
        message(paste0(Sys.time(), " unable to load ",directory))
      return(NULL)
    });
  }
  else {
    message(paste0(Sys.time(), " directory ",directory, " does not exist."))
  }
}

getTraitsWithSummary <- function(globalVariables, directory, exludeMultiModalProbes = TRUE){
  print(paste0(Sys.time(), " start reading traits"))
  if (dir.exists(directory)) {
    traits <- data.frame(Trait="ex", MaxN = 1, MinP_VAL = 0, MinFDR = 0, MaxBETA = 1, MaxDeltaMeth = 1, MinDeltaMeth = 1, MaxOutlying = 1, MinOutlying = 0, MaxSkewed = 1, MinSkewed = 0, MaxClumpy = 1, MinClumpy = 0, MaxSparse = 1, MinSparse = 0, MaxStriated = 1, MinStriated = 0, MaxConvex = 1, MinConvex = 0, MaxSkinny = 1, MinSkinny = 0, MaxStringy = 1, MinStringy = 0, MaxMonotonic = 1, MinMonotonic = 0)
    tr <- traits
    print(paste0(Sys.time(), " start loading folder"))
    tryCatch({
      temp <- list.files(path=directory,pattern="\\.csv$")
    }, error=function(err){
      message(paste0(Sys.time(), "unable to load '",directory,"'"))
    });
    tryCatch({
      if (globalVariables$config$debugMode == TRUE) {
        max <- min(50:length(temp))
      }
      else {
        max <- length(temp)
      }
      fileNameTraitsSummary <- "traitsSummary.RDS" #check for summary file
      RDSfileName <- paste0(directory,fileNameTraitsSummary)
      if (utils::file_test("-f", RDSfileName) == TRUE) {

        print(paste0(Sys.time(), " reading RDS ", RDSfileName))
        traits = readRDS (file = RDSfileName)
      }
      else {
        for (i in 1:max) {
          #firstlines <- utils::read.table(file = paste0(directory,as.character(temp[i])), sep = "\t", header = T, nrows = 5)
          firstlines <- data.table::fread(file = paste0(directory,as.character(temp[i])), sep = "\t", header = T, nrows = 5, data.table = FALSE)
          if (colnames(firstlines)[1] == globalVariables$config$probeAttribut) {
            if (nrow(firstlines) >= 5) {
                fileName <- stringr::str_sub(temp[i], 1, stringr::str_length(temp[i])-4)
                tr$trait <- fileName
                fileName <- paste0(directory,fileName,".csv")
                print(paste0(Sys.time(), " reading trait file ", fileName))

                all.results <- data.table::fread(fileName, stringsAsFactors = FALSE, header = TRUE, sep = "\t", data.table = FALSE)

                if (exludeMultiModalProbes == TRUE) {
                  #merge with annotation and exclude multimodal probes...
                  print(paste0(Sys.time(), " before excluding multimodal probes"))
                  #all.results <- base::merge(all.results, globalVariables$annotation, by.x = "probeID", by.y = "name", all.x = FALSE, all.y = FALSE)
                  all.results <- base::merge(all.results, globalVariables$annotation, by.x = globalVariables$config$probeAttribut, by.y = "name", all.x = FALSE, all.y = FALSE)
                  all.results <- all.results[,1:16]
                  print(paste0(Sys.time(), " finished excluding multimodal probes"))
                }

                testthat::expect_type(all.results, "list")
                all.results <- all.results[order(all.results$N,decreasing = TRUE),]
                tr$MaxN <- all.results$N[1]
                all.results <- all.results[order(all.results$P_VAL),]
                tr$MinP_VAL <- all.results$P_VAL[1]
                all.results <- all.results[order(all.results$FDR),]
                tr$MinFDR <- all.results$FDR[1]
                all.results <- all.results[order(all.results$BETA),]
                tr$MaxBETA <- all.results$BETA[1]
                all.results <- all.results[order(all.results$DeltaMeth),]
                tr$MaxDeltaMeth <- max(all.results$DeltaMeth)
                tr$MinDeltaMeth <- min(all.results$DeltaMeth)
                tryCatch({
                  all.results <- all.results[order(all.results$Outlying),]
                  tr$MaxOutlying <- max(all.results$Outlying)
                  tr$MinOutlying <- min(all.results$Outlying)
                  all.results <- all.results[order(all.results$Skewed),]
                  tr$MaxSkewed <- max(all.results$Skewed)
                  tr$MinSkewed <- min(all.results$Skewed)
                  all.results <- all.results[order(all.results$Clumpy),]
                  tr$MaxClumpy <- max(all.results$Clumpy)
                  tr$MinClumpy <- min(all.results$Clumpy)
                  all.results <- all.results[order(all.results$Sparse),]
                  tr$MaxSparse <- max(all.results$Sparse)
                  tr$MinSparse <- min(all.results$Sparse)
                  all.results <- all.results[order(all.results$Striated),]
                  tr$MaxStriated <- max(all.results$Striated)
                  tr$MinStriated <- min(all.results$Striated)
                  all.results <- all.results[order(all.results$Convex),]
                  tr$MaxConvex <- max(all.results$Convex)
                  tr$MinConvex <- min(all.results$Convex)
                  all.results <- all.results[order(all.results$Skinny),]
                  tr$MaxSkinny <- max(all.results$Skinny)
                  tr$MinSkinny <- min(all.results$Skinny)
                  all.results <- all.results[order(all.results$Stringy),]
                  tr$MaxStringy <- max(all.results$Stringy)
                  tr$MinStringy <- min(all.results$Stringy)
                  all.results <- all.results[order(all.results$Monotonic),]
                  tr$MaxMonotonic <- max(all.results$Monotonic)
                  tr$MinMonotonic <- min(all.results$Monotonic)
                }, error=function(err){
                  message(paste0(Sys.time(), " no scagnostics information found"))
                });
                traits <- rbind(traits, tr)
            }
          }
          else {
            message(paste0(Sys.time(), " attribut", globalVariables$config$probeAttribut, " from <config.yml> not found in data file ", temp[i] ))
          }
        }
        #remove dummy trait
        traits <- traits[-c(1), ]
        #save DF for faster later access
        if (globalVariables$config$debugMode == FALSE) {
          print(paste0(Sys.time(), "writing RDS ", RDSfileName))
          saveRDS(traits, file = RDSfileName)
        }
      }
      print(paste0(Sys.time(), " finished reading traits"))
      return(traits)

    }, error=function(err){
      message(paste0(Sys.time(), "unable to load ",fileName, " err: ", err))
    });
  }
}

