library(shiny)
library(plotly)

plotManhattanVolcano_UI <- function(id){
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::fluidRow(
      shiny::verbatimTextOutput(ns("trait"), placeholder = TRUE)
    ),
    shinyBS::bsCollapse(id = "collapsePlot", open = c("plot"), multiple = TRUE,
          shinyBS::bsCollapsePanel("plot", label = "Select Probe",
            shiny::fluidRow(
              shiny::sliderInput(ns("numberResults"),"top n of results",
                          10, 1000, 1000, 1, width = "100%")
            ),
          shiny::tabsetPanel(
            shiny::tabPanel("Visualisation",
#            fluidRow(
                shiny::column(6, plotly::plotlyOutput(ns("plotManhattan"))),
                shiny::column(6, plotly::plotlyOutput(ns("plotVolcano")))
#            )
              ),
              shiny::tabPanel("Table",
                DT::dataTableOutput(ns("dt"))
              ),
              shiny::tabPanel("Link to Pathway Analysis",
                 shiny::fluidRow(
                   htmltools::tags$html("Gene symbols for PathwayAanalysis"),
                   DT::dataTableOutput(ns("DTGeneSymbolOut")),
                 )
              )
            ),
            shiny::fluidRow(
              shiny::column(6, htmltools::tags$html("last selected probe"),
                     shiny::verbatimTextOutput(ns("txtSelectedProbe"), placeholder = TRUE)
              ),
              shiny::column(6, htmltools::tags$html("selected probes"),
                     shiny::verbatimTextOutput(ns("txtSelectedProbes"), placeholder = TRUE)
              )
            )
          )
    )
  )
}

plotManhattanVolcano_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
    #update sliderInput
    shiny::updateSliderInput(
      session = session,
      "numberResults",
      max = nrow(sessionVariables$resultDataSingleTrait),
      value = min(250,nrow(sessionVariables$resultDataSingleTrait))
    )
    reDFManhattanVolcano <- shiny::reactive({getDFforManhattanVolcano(globalVariables, sessionVariables, input$numberResults)})
    output$plotManhattan <- plotly::renderPlotly(plotlyManhattanVolcano(globalVariables,reDFManhattanVolcano(),"M"))
    output$plotVolcano <- plotly::renderPlotly(plotlyManhattanVolcano(globalVariables,reDFManhattanVolcano(),"V"))

    output$dt <- DT::renderDataTable({
      id <- shiny::showNotification("printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      tryCatch({
        print(paste0(Sys.time(), " render data table Manhattan/ volcano."))
        DT::datatable(reDFManhattanVolcano(), escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
                      options = list(searching = TRUE, pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv', 'excel','pdf')))
      }, error = function(err) {
        shiny::validate(shiny::need(nrow(df)>0,"No data to show"))
      })
    }, server = FALSE)

    output$DTGeneSymbolOut <- DT::renderDataTable(server=FALSE,{
      id <- shiny::showNotification("printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      tryCatch({
        # Load data
        data <- getDFforPathwayAnalysis(reDFManhattanVolcano())
        # Show data
        print(paste0(Sys.time(), " render data table pathway."))
        DT::datatable(data, extensions = 'Buttons', height = 400,
                      options = list(scrollY = TRUE, scroller = TRUE, searching = TRUE,
                                     ordering = TRUE, dom = 'ftBS', buttons = c('copy', 'csv', 'excel','pdf')))
      }, error = function(err) {
        shiny::validate(shiny::need(nrow(df)>0,"No data to show"))
      })
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = "plotlyManhattan"), suspended = FALSE, {
      sessionVariables$probe$probe <- as.character(plotly::event_data("plotly_click", source = "plotlyManhattan")$key[1])
      id <- shiny::showNotification(paste0("Selected probe: ", sessionVariables$probe$probe), duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      output$txtSelectedProbe <- shiny::renderText(sessionVariables$probe$probe)
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$dt_cell_clicked, {
      df = reDFManhattanVolcano()
      selected = input$dt_rows_selected
      if (length(selected)) {
        #selectedProbeIDs <- df[selected,]$probeID
        selectedProbeIDs <- df[selected,globalVariables$config$probeAttribut]
        sessionVariables$probe$probe <- selectedProbeIDs[1]
        sessionVariables$probe$probes <- selectedProbeIDs
        output$txtSelectedProbe <- shiny::renderText(sessionVariables$probe$probe)
        output$txtSelectedProbes <- shiny::renderText(sessionVariables$probe$probes)
      }
    })

  })
}

getDFforManhattanVolcano <- function(globalVariables, sessionVariables, n) {
  df <- sessionVariables$resultDataSingleTrait
  df <- df[order(df$P_VAL,decreasing=FALSE),]
  df <- df[1:n,]
  df <- resultDataSingleScenarioWithAnnotation(globalVariables$annotation, df)
  df <- resultDataSingleScenarioWithAnnotationEWAScatalogCount(globalVariables, df)
  df=stats::na.omit(df,globalVariables$config$probeAttribut)
  return(df)
}


plotlyManhattanVolcano <- function(globalVariables, DF, M_V) {
  tryCatch({
    print(paste0(Sys.time(), " plot manhattan / volcano."))
    if(missing(M_V)) M_V = "V"
    colors = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12,'Paired'))
    plot = plotly::plot_ly(data = DF, source = "plotlyManhattan")
    if (M_V == "M") {
      plot = plot %>% plotly::add_trace(x = ~globalPosition, y = ~P_VAL, #x = sharedDF$data()$n, y = sharedDF$data()$P_VAL,
                                color = ~chromosomeNum, colors = colors(24),
                                type = 'scatter', mode = 'markers',
                                #                                marker = list(size = sharedDF$data()$n, opacity = 0.5),
                                marker = list(opacity = 0.5, sizemode = 'diameter'),
                                size = DF$n,
                                fill = ~'',
                                #text = paste0(DF$probeID,"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                text = paste0(DF[,globalVariables$config$probeAttribut],"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                hoverinfo = 'text', key = ~probeID)
      plot = plot %>% plotly::layout(xaxis = list(type = "lin"),
                             yaxis = list(type = "log", autorange = "reversed"))
    }
    else {
      plot = plot %>% plotly::add_trace(x = ~DeltaMeth, y = ~P_VAL,
                                color = ~chromosomeNum, colors = colors(24),
                                type = 'scatter', mode = 'markers',
                                marker = list(opacity = 0.5, sizemode = 'diameter'),
                                size = DF$n,
                                fill = ~'',
                                #text = paste0(DF$probeID,"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                text = paste0(DF[,globalVariables$config$probeAttribut],"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                hoverinfo = 'text', key = ~probeID)
      plot = plot %>% plotly::layout(xaxis = list(type = "lin", zeroline = FALSE, showline = TRUE,
                                          showticklabels = TRUE, showgrid = FALSE),
                             yaxis = list(type = "log", autorange = "reversed"))
    }
    return (plot)
  }, error=function(err){
    message(paste0(Sys.time(), " unable to print manhattan plot; ", err$message))
    return(empty_plot(err$message))
  })

}

getDFforPathwayAnalysis <- function (df) {
  #subset probeID, Gene.symbol      logFC    adj.P.Val
  df <- subset(df, select=c("probeID", "gene.symbol", "DeltaMeth", "P_VAL"))
  df$Gene.symbol <- ""
  i = NULL
  foreach(i=1:nrow(df)) %do% {
    df[i,]$Gene.symbol <- unlist(strsplit(df[i,]$gene.symbol,";"))[1]
  }
  df$gene.symbol <- NULL
  df <- subset(df, !is.na(df$Gene.symbol))
  setcolorder(df, c("probeID","Gene.symbol","DeltaMeth","P_VAL"))
  return(df)
}
