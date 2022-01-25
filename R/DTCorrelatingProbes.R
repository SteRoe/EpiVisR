DTCorrelatingProbes_UI <- function(id){
  ns <- shiny::NS(id)
  htmltools::tagList(

    shiny::fluidRow(
      shiny::column(width = 12,
             DT::dataTableOutput(ns("DTCorrelatingProbes"))
      )
    ),
    shiny::fluidRow(
      shiny::column(width = 12,
              shiny::verbatimTextOutput(ns("selectedProbes"), placeholder = TRUE),
              htmltools::tags$head(htmltools::tags$style("#selectedProbes{overflow-y:scroll; height: 50px;}"))
      )
    ),
    shiny::fluidRow(
      shiny::column(width = 12, shiny::wellPanel(
        plotly::plotlyOutput(ns("plotlyCorrelatingProbes"))
      ))
    )
  )
}

DTCorrelatingProbes_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
<<<<<<< HEAD

    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
=======
#    selectedProbe <- reactiveValues(probe = "", probes = list())
#    moduleVariables <- reactiveValues(numberResults = 10)

    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)

>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
    reDFCorrelatingProbes <- shiny::reactive({correlatingProbes(globalVariables, sessionVariables)})

    output$DTCorrelatingProbes <- DT::renderDataTable({
      id <- shiny::showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
<<<<<<< HEAD
      tryCatch({
        print(paste0(Sys.time(), " rendering correlating data table."))
        CP <- reDFCorrelatingProbes()
#        CP <- addLinkToEWASDataHub(CP, globalVariables$config$baseURL_EWASDataHub)
#        CP <- addLinkToMRCEWASCatalog(CP, globalVariables$config$baseURL_MRCEWASCatalog)
        DT::datatable(CP, escape = F, extensions="Scroller", style="bootstrap", class="compact", width="100%",
                  options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      }, error = function(err) {
=======
      #    probeID <- currentProbeID()
      # datatable(shared_df_correlatingProbes, extensions="Scroller", style="bootstrap", class="compact", width="100%",
      #           options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      tryCatch({
        print(paste0(Sys.time(), " rendering correlating data table."))
        #    CP = correlatingProbes(sessionVariables)
        CP <- reDFCorrelatingProbes()
        CP <- addLinkToEWASDataHub(CP, globalVariables$config$baseURL_EWASDataHub)
        CP <- addLinkToMRCEWASCatalog(CP, globalVariables$config$baseURL_MRCEWASCatalog)
        DT::datatable(CP, escape = F, extensions="Scroller", style="bootstrap", class="compact", width="100%",
                  options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      }, errror = function(err) {
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
        shiny::validate(shiny::need(nrow(CP)>0,"No data to show"))
      })
    }, server = TRUE)

<<<<<<< HEAD
    shiny::observeEvent(input$DTCorrelatingProbes_cell_clicked, {
=======
      shiny::observeEvent(input$DTCorrelatingProbes_cell_clicked, {
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
      selected = input$DTCorrelatingProbes_rows_selected
      if (length(selected)) {
        CP <- reDFCorrelatingProbes()
        selectedProbeIDs <- CP[selected,]$probeID
        selectedProbes <- globalVariables$beta.t[,selectedProbeIDs]
        selectedProbes <- as.data.frame(selectedProbes)
        if (ncol(selectedProbes) == 1) {
          colnames(selectedProbes) = selectedProbeIDs
        }
        output$selectedProbes = shiny::renderPrint({
          cat (colnames(selectedProbes))
        })
        if (!is.null(selectedProbes)) {
          result = list()
          i = NULL
          foreach(i=1:ncol(selectedProbes)) %do% {
            traitVar<-traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut)
<<<<<<< HEAD
  #            selectedProbes$ID_Kind <- rownames(selectedProbes)
=======
#            selectedProbes$ID_Kind <- rownames(selectedProbes)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
            selectedProbes$ID <- rownames(selectedProbes)
            selectedProbe <- as.data.frame(selectedProbes[,i])
            colnames(selectedProbe)[1] <- colnames(selectedProbes)[i]
            print(paste0(Sys.time(), " merging correlating probes: ", colnames(selectedProbes)[i]))
            selectedProbe$ID <- rownames(selectedProbes)
            selectedProbeWithExposure = base::merge(traitVar, selectedProbe, by.x = globalVariables$config$mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)
            selectedProbeWithExposure <- stats::na.omit(selectedProbeWithExposure)
            result = c(result,list(selectedProbeWithExposure))
          }

          output$plotlyCorrelatingProbes <- plotly::renderPlotly({

            plotList = list()
            dfList = result
<<<<<<< HEAD
  #            if (!is_null(dfList)) {
=======
#            if (!is_null(dfList)) {
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
            if (!is.null(dfList)) {
              id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
              on.exit(shiny::removeNotification(id), add = TRUE)
              i = NULL
              foreach(i=1:length(dfList)) %do% {
                df = as.data.frame(dfList[i])
                fmla = stats::as.formula(paste0("`", colnames(df)[4], "` ~ `", colnames(df)[3], "`"))
                m <- stats::lm(fmla, data = df)
<<<<<<< HEAD
=======
                #    aug = broom::augment(m,se_fit=TRUE)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
                plot = broom::augment(m,se_fit=TRUE) %>%
                  plotly::plot_ly(x = stats::as.formula(paste0("~ `", colnames(df)[3], "`")), showlegend = FALSE)%>%
                  plotly::add_markers(x = df[,3], y = df[,4], name = colnames(df)[4], showlegend = TRUE)%>%
                  plotly::add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                              ymax = ~.fitted + 1.96 * .se.fit,
                              color = I("gray80"), showlegend = FALSE)%>%
                  plotly::add_lines(y = ~.fitted, color = I("steelblue"), showlegend = FALSE)

                plotList = c(plotList,list(plot))
              }
              plotlyscatter <- plotly::subplot(plotList, shareX = TRUE, nrows = length(dfList))
              return(plotlyscatter)
            }
          })
        }
      }
    })
  })
}

getCorrelatingProbes<-function(globalVariables, probeID, currentData) {
  if (shiny::isTruthy(probeID)) {
    tryCatch({
      probeID <- as.character(probeID)
      correlatingProbes <- array(stats::cor(globalVariables$beta.t,globalVariables$beta.t[,probeID], method = c("pearson")))
      correlatingProbes <- as.data.frame(correlatingProbes)
      colnames(correlatingProbes) <- "corr.coeff"
      correlatingProbes$probeID<-colnames(globalVariables$beta.t)
      rownames(correlatingProbes)<-correlatingProbes$probeID
<<<<<<< HEAD
      #merge with annotation
      correlatingProbes <- dplyr::left_join(correlatingProbes, globalVariables$annotation, by = c("probeID" = "name"))
=======
      #    correlatingProbes <- subset(correlatingProbes, probeID!=probeID)
      #merge with annotation
      correlatingProbes <- dplyr::left_join(correlatingProbes, globalVariables$annotation, by = c("probeID" = "name"))
#      correlatingProbes = base::merge(correlatingProbes, annotation, by.x = "probeID", by.y = "name", all.x = TRUE, all.y=FALSE)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
      correlatingProbes$type <- NULL
      correlatingProbes$target <- NULL
      correlatingProbes$meth.dye <- NULL
      #add DeltaMeth
      DeltaMeth <- currentData[,c("probeID","DeltaMeth","P_VAL","FDR")]
      correlatingProbes <- dplyr::inner_join(DeltaMeth, correlatingProbes, by = c("probeID" = "probeID"))
<<<<<<< HEAD
=======
#      correlatingProbes = base::merge(DeltaMeth, correlatingProbes, by.x = "probeID", by.y = "probeID", all.x = TRUE, all.y=TRUE)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
      #sort by corr.coeff
      correlatingProbes <- correlatingProbes[order(abs(correlatingProbes[,"corr.coeff"]),decreasing=TRUE),]
      return (correlatingProbes)
    }, error=function(err){
      message(paste0(Sys.time(), "unable to correlate ", probeID))
    })
  }
}

correlatingProbes <- function(globalVariables, sessionVariables){
  print(paste0(Sys.time(), " getting correlating probes."))
  probeID = sessionVariables$probe$probe
  if (!is.null(probeID)) {
    df <- sessionVariables$resultDataSingleTrait
    correlatingProbes <- getCorrelatingProbes(globalVariables, probeID, df)
<<<<<<< HEAD
    if (nrow(correlatingProbes) > 100) {
      correlatingProbes <- correlatingProbes[1:100,]
    }
=======
    correlatingProbes <- correlatingProbes[1:100,]
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
    return (correlatingProbes)
  }
}
