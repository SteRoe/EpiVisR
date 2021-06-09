library(shiny)
library(shinyBS)

inputFolder_UI <- function(id){
  ns <- NS(id)
  
#  tagList(
    bsCollapse(id = ns("collapseExample"), open = c(ns("Folder"),ns("Traits")), multiple = TRUE,
      bsCollapsePanel(ns("Folder"),
        shiny::textInput(ns("directory"), "Folder for Traits", dataDir, placeholder = TRUE),
        shiny::actionButton(ns("btnReadFolder"), label = "Select Folder")
        )
      )
#  )
}

inputFolder_SERVER <- function(id, sessionVariables) {
  moduleServer(
    id,
    function(input, output, session) {
      selectedTrait <- reactiveValues(trait = "", traits = list())
      observeEvent(input$btnReadFolder, { #ignoreInit = TRUE,
        sessionVariables.Folder = input$directory
#         id <- shiny::showNotification("Filling data table...", duration = NULL, closeButton = FALSE)
#         on.exit(removeNotification(id), add = TRUE)
#         traits <- getTraits(input$directory)
#         updateSelectInput(
#           session = session,
#           inputId = "selSelectedTrait",
#           choices = traits
#         )
#                                     
#         traitsWithSummary <- getTraitsWithSummary(input$directory)
#         output$traits <- DT::renderDataTable({
#           datatable(traitsWithSummary, extensions = list("Scroller"), style = "bootstrap", class = "compact", width = "100%",
#                     options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE))
#         }, server = FALSE)
       }, ignoreNULL = FALSE)
#      return(selectedTrait)
    }
  )
}
