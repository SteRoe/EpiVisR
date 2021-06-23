DTCorrelatingProbes_UI <- function(id){
  ns <- NS(id)
  tagList(

    fluidRow(
      column(width = 12,
             DT::dataTableOutput(ns("DTCorrelatingProbes"))
      )
    ),
    fluidRow(
      column(width = 12,
             verbatimTextOutput(ns("selectedProbes"), placeholder = TRUE),
             tags$head(tags$style("#selectedProbes{overflow-y:scroll; height: 50px;}"))
      )
    ),
    fluidRow(
      column(width = 12, wellPanel(
        plotly::plotlyOutput(ns("plotlyCorrelatingProbes"))
      ))
    )
  )
}

correlatingProbes <- function(sessionVariables){
#browser()
  message(paste0(Sys.time(), " getting correlating probes."))
  probeID = sessionVariables$probe$probe
  if (!is_empty(probeID)) {
#browser()
#    if(isTruthy(sessionVariables$resultDataSingleTrait)) {
      df <- sessionVariables$resultDataSingleTrait
    # }
    # else {
    #   df <- isolate(getResultDataSingleTrait(t,moduleVariables$numberResults))
    # }
#    resultDataSingleTrait = resultDataSingleTrait(sessionVariables$trait$trait)
    correlatingProbes <- getCorrelatingProbes(probeID, df)
    correlatingProbes <- correlatingProbes[1:100,]
    return (correlatingProbes)
  }
}

DTCorrelatingProbes_SERVER <- function(id, sessionVariables) {
  moduleServer(id, function(input, output, session) {
#    selectedProbe <- reactiveValues(probe = "", probes = list())
#    moduleVariables <- reactiveValues(numberResults = 10)
    
    id <- showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    reDFCorrelatingProbes <- reactive({correlatingProbes(sessionVariables)})
    
    output$DTCorrelatingProbes <- DT::renderDataTable({
      id <- showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      #    probeID <- currentProbeID()
      # datatable(shared_df_correlatingProbes, extensions="Scroller", style="bootstrap", class="compact", width="100%",
      #           options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      tryCatch({
        message(paste0(Sys.time(), " rendering correlating data table."))
        #    CP = correlatingProbes(sessionVariables)
        CP <- reDFCorrelatingProbes()
        CP <- addLinkToEWASDataHub(CP)
        CP <- addLinkToMRCEWASCatalog(CP)
        datatable(CP, escape = F, extensions="Scroller", style="bootstrap", class="compact", width="100%",
                  options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      }, errror = function(err) {
        validate(need(nrow(CP>0,"No data to show")))
      })
    }, server = TRUE)

      observeEvent(input$DTCorrelatingProbes_cell_clicked, {
      selected = input$DTCorrelatingProbes_rows_selected
      if (length(selected)) {
        CP <- reDFCorrelatingProbes()
        selectedProbeIDs <- CP[selected,]$probeID
        selectedProbes <- beta.t[,selectedProbeIDs]
        selectedProbes <- as.data.frame(selectedProbes)
        if (ncol(selectedProbes) == 1) {
          colnames(selectedProbes) = selectedProbeIDs
        }
        output$selectedProbes = renderPrint({
          cat (colnames(selectedProbes))
        })
        if (!is_empty(selectedProbes)) {
          result = list()
          foreach(i=1:ncol(selectedProbes)) %do% {
            traitVar<-traitDF(sessionVariables)
#            selectedProbes$ID_Kind <- rownames(selectedProbes)
            selectedProbes$ID <- rownames(selectedProbes)
            selectedProbe <- as.data.frame(selectedProbes[,i])
            colnames(selectedProbe)[1] <- colnames(selectedProbes)[i]
            message(paste0(Sys.time(), " merging correlating probes: ", colnames(selectedProbes)[i]))
            selectedProbe$ID <- rownames(selectedProbes)
            selectedProbeWithExposure = base::merge(traitVar, selectedProbe, by.x = config$mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)
            selectedProbeWithExposure <- na.omit(selectedProbeWithExposure)
            result = c(result,list(selectedProbeWithExposure))
          }
          
          output$plotlyCorrelatingProbes <- plotly::renderPlotly({

            plotList = list()
            dfList = result
            if (!is_null(dfList)) {
              id <- showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
              on.exit(removeNotification(id), add = TRUE)
              foreach(i=1:length(dfList)) %do% {
                df = as.data.frame(dfList[i])
                fmla = as.formula(paste0("`", colnames(df)[4], "` ~ `", colnames(df)[3], "`"))
                m <- lm(fmla, data = df)
                #    aug = broom::augment(m,se_fit=TRUE)
                plot = broom::augment(m,se_fit=TRUE) %>%
                  plot_ly(x = as.formula(paste0("~ `", colnames(df)[3], "`")), showlegend = FALSE)%>%
#                  add_markers(x = df[,3], y = df[,4], color = I("black"), name = colnames(df)[4], showlegend = TRUE)%>%
                  add_markers(x = df[,3], y = df[,4], name = colnames(df)[4], showlegend = TRUE)%>%
                  add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                              ymax = ~.fitted + 1.96 * .se.fit,
                              color = I("gray80"), showlegend = FALSE)%>%
                  add_lines(y = ~.fitted, color = I("steelblue"), showlegend = FALSE)                

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

getCorrelatingProbes<-function(probeID, currentData) {
  if (isTruthy(probeID)) {
    tryCatch({
      probeID <- as.character(probeID)
      #    setBaseSciPen()
      correlatingProbes <- array(cor(beta.t,beta.t[,probeID], method = c("pearson")))
      correlatingProbes <- as.data.frame(correlatingProbes)
      colnames(correlatingProbes) <- "corr.coeff"
      correlatingProbes$probeID<-colnames(beta.t)
      rownames(correlatingProbes)<-correlatingProbes$probeID
      #    correlatingProbes <- subset(correlatingProbes, probeID!=probeID)
      #merge with annotation
      correlatingProbes <- dplyr::left_join(correlatingProbes, annotation, by = c("probeID" = "name"))
      correlatingProbes$type <- NULL
      correlatingProbes$target <- NULL
      correlatingProbes$meth.dye <- NULL
      #add DeltaMeth
      DeltaMeth <- currentData[,c("probeID","DeltaMeth","P_VAL","FDR")]
      correlatingProbes <- dplyr::inner_join(DeltaMeth, correlatingProbes, by = c("probeID" = "probeID"))
      #sort by corr.coeff
      correlatingProbes <- correlatingProbes[order(abs(correlatingProbes[,"corr.coeff"]),decreasing=TRUE),]
      #    resetSciPen()
      return (correlatingProbes)
    }, error=function(err){
      print(paste0("unable to correlate ",probeID))
    })
  }
}