library(shiny)
library(plotly)

plotManhattanVolcano_UI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      verbatimTextOutput(ns("trait"), placeholder = TRUE)
    ),
    bsCollapse(id = "collapsePlot", open = c("plot"), multiple = TRUE,
          bsCollapsePanel("plot", label = "Select Probe",
            fluidRow(
              column(10,
              sliderInput(ns("numberResults"),"top n of results",
                          10, 1000, 1000, 1, width = "100%")
              ),
              column(2, style = "margin-top: 25px;",
                shiny::actionButton(ns("btnPlot"), label = "Update # of Probes")
              )
            ),
            tabsetPanel(
              tabPanel("Visualisation",
#            fluidRow(
                column(6, plotly::plotlyOutput(ns("plotManhattan"))),
                column(6, plotly::plotlyOutput(ns("plotVolcano")))
#            )
              ),
              tabPanel("Table",
                DT::dataTableOutput(ns("dt"))
              ),
              tabPanel("Link to Pathway Analysis",
                 fluidRow(
                   tags$html("Gene symbols for PathwayAanalysis"),
                   DT::dataTableOutput(ns("DTGeneSymbolOut")),
                   
                   # textOutput(ns("txtGeneSymbol")),
                   # tags$html("Entrez Gene ID"),
                   # textOutput(ns("txtEntrezGeneID")),
                   # tags$html("KEGG ID"),
                   # textOutput(ns("txtKEGGID")),
                   # verbatimTextOutput(ns("vtxtPathway"))
                 )
              )
            ),
            fluidRow(
              column(6, tags$html("last selected probe"),
                     verbatimTextOutput(ns("txtSelectedProbe"), placeholder = TRUE)
              ),
              column(6, tags$html("selected probes"),
                     verbatimTextOutput(ns("txtSelectedProbes"), placeholder = TRUE)
              )
            )
          )
    )
  )
}

plotManhattanVolcano_SERVER <- function(id, sessionVariables) {
  moduleServer(id, function(input, output, session) {
    moduleVariables <- reactiveValues(df = data.frame())
    #update sliderInput
    updateSliderInput(
      session = session,
      "numberResults",
      max = nrow(sessionVariables$resultDataSingleTrait),
      value = min(250,nrow(sessionVariables$resultDataSingleTrait))
    )
    plotManhattanVolcano(input, output, sessionVariables, moduleVariables)
    
    observeEvent(event_data("plotly_click", source = "plotlyManhattan"), suspended = FALSE, {
      sessionVariables$probe$probe <- event_data("plotly_click", source = "plotlyManhattan")$key
      id <- showNotification(paste0("Selected probe: ", sessionVariables$probe$probe), duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      output$txtSelectedProbe <- renderText(sessionVariables$probe$probe)
    }, ignoreNULL = FALSE)
    
    observeEvent(input$dt_cell_clicked, {
      df = moduleVariables$df
      selected = input$dt_rows_selected
      if (length(selected)) {
        selectedProbeIDs <- df[selected,]$probeID
        sessionVariables$probe$probe <- selectedProbeIDs[1]
        sessionVariables$probe$probes <- selectedProbeIDs
        output$txtSelectedProbe <- renderText(sessionVariables$probe$probe)
        output$txtSelectedProbes <- renderText(sessionVariables$probe$probes)
      }
    })
    
    observeEvent(input$btnPlot, ignoreInit = FALSE, {
      plotManhattanVolcano(input, output, sessionVariables, moduleVariables)
    }, ignoreNULL = FALSE)
  })
}

plotManhattanVolcano <- function (input, output, sessionVariables, moduleVariables) {
  output$trait <- renderText(sessionVariables$trait$trait)
  df <- sessionVariables$resultDataSingleTrait
  df = df[order(df$P_VAL,),]
  moduleVariables$df = df
  df <- df %>% slice_min(P_VAL, n = input$numberResults)
  df <- resultDataSingleScenarioWithAnnotation(df)
  df <- resultDataSingleScenarioWithAnnotationEWAScatalogCount(df)
  output$plotManhattan <- plotly::renderPlotly(plotlyManhattanVolcano(df,"M"))
  output$plotVolcano <- plotly::renderPlotly(plotlyManhattanVolcano(df,"V"))
  
  output$dt <- DT::renderDataTable({
    id <- showNotification("Printing data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    tryCatch({
    DT::datatable(df, escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
              options = list(searching = TRUE, pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv', 'excel','pdf')))
    }, errror = function(err) {
      validate(need(nrow(df>0,"No data to show")))
    })
  }, server = FALSE)
  output$DTGeneSymbolOut <- DT::renderDataTable(server=FALSE,{
    id <- showNotification("Printing data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    tryCatch({
    # Load data
      data <- getDFforPathwayAnalysis(df)
      # Show data
      datatable(data, extensions = 'Buttons', height = 400,
                options = list(scrollY = TRUE, scroller = TRUE, searching = TRUE, 
                               ordering = TRUE, dom = 'ftBS', buttons = c('copy', 'csv', 'excel','pdf')))
    }, errror = function(err) {
      validate(need(nrow(df>0,"No data to show")))
    })
  })
  # gene.symbol <- extractGeneSymbols(df)
  # entrezID <- extractEntrezID(gene.symbol)
  # KEGGID <- extractKEGGID(entrezID)
  
  #      output$txtPathway <- renderText(t)
  #      output$vtxtPathway <- renderText(t)
  
  # output$txtGeneSymbol <- renderText(gene.symbol)
  # output$txtEntrezGeneID <- renderText(entrezID)
  # output$txtKEGGID <-renderText(KEGGID)
}

plotlyManhattanVolcano <- function(DF, M_V) {
  tryCatch({
    #browser()
    if(missing(M_V)) M_V = "V"
    colors = colorRampPalette(brewer.pal(12,'Paired'))
    plot = plot_ly(data = DF, source = "plotlyManhattan")
    if (M_V == "M") {
      plot = plot %>% add_trace(x = ~globalPosition, y = ~P_VAL, #x = sharedDF$data()$n, y = sharedDF$data()$P_VAL,
                                color = ~chromosomeNum, colors = colors(24),
                                type = 'scatter', mode = 'markers',
                                #                                marker = list(size = sharedDF$data()$n, opacity = 0.5),
                                marker = list(opacity = 0.5, sizemode = 'diameter'),
                                size = DF$n,
                                fill = ~'',
                                text = paste0(DF$probeID,"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                hoverinfo = 'text', key = ~probeID)
      plot = plot %>% layout(xaxis = list(type = "lin"),
                             yaxis = list(type = "log", autorange = "reversed"))
    }
    else {
      plot = plot %>% add_trace(x = ~DeltaMeth, y = ~P_VAL,
                                color = ~chromosomeNum, colors = colors(24), #colors = colorRampPalette('Paired'),
                                type = 'scatter', mode = 'markers',
                                marker = list(opacity = 0.5, sizemode = 'diameter'),
                                size = DF$n,
                                fill = ~'',
                                text = paste0(DF$probeID,"\nDeltaMeth: ", DF$DeltaMeth,"\n Gene symbol: ",DF$gene.symbolShort,"\nn:",DF$n),
                                hoverinfo = 'text', key = ~probeID)
      plot = plot %>% layout(xaxis = list(type = "lin", zeroline = FALSE, showline = TRUE,
                                          showticklabels = TRUE, showgrid = FALSE),
                             yaxis = list(type = "log", autorange = "reversed"))
    }
    return (plot)
    #    event_register("plotlyManhattan", 'plotly_click')
    # d <- event_data("plotly_click", source = "plotlyManhattan")
    # #if (is.null(d)) "Click events appear here (double-click to clear)" else d
    # if (!is.null(d)) {
    #   browser()     
    # }
  }, error=function(err){
    print(paste0("unable to print manhattan plot; ", err$message))
    return(empty_plot(err$message))
  })
  
}

getDFforPathwayAnalysis <- function (df) {
  #subset probeID, Gene.symbol      logFC    adj.P.Val
  df <- subset(df, select=c("probeID", "gene.symbol", "DeltaMeth", "P_VAL"))
  df$Gene.symbol <- ""
  foreach(i=1:nrow(df)) %do% {
    df[i,]$Gene.symbol <- unlist(strsplit(df[i,]$gene.symbol,";"))[1]
  }
  df <- df[,-"gene.symbol"]
  df <- subset(df, !is.na(df$Gene.symbol))
  setcolorder(df, c("probeID","Gene.symbol","DeltaMeth","P_VAL"))
  ##Get the Entrez gene IDs associated with those symbols
  #EG_IDs = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
  
  ##Then get the KEGG IDs associated with those entrez genes.
  #KEGG_IDs = mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA)

  #make DF wit gene.symbol, logFoldChange and P_VAL for pathway analysis like in pathfindR
  #decompose

  #extractEntrezID
  #return (mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA))
  #extractKEGGID
  #return (mget(as.character(EG_IDs), org.Hs.egPATH,ifnotfound=NA))
  return(df)
}
