plotDNAmProfile_UI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 10,
             tags$html("Window size"),
             sliderInput(ns("DMRWindow"),"DMR window size",
                         1, 50, 5, 1, width = "100%")
      ),
      column(width = 2, style = "margin-top: 25px;",
             shiny::actionButton(ns("btnDMRWindow"), label = "Window Size"),
      )
    ),
    tabsetPanel(
      tabPanel("Visualisation",
         fluidRow(
           column(width = 10,
                  plotly::plotlyOutput(ns("PlotlyPcPDMPNearRange"))
           ),
           # fluidRow(
           #   tags$html("color code: range of:"),
           #   textOutput(ns("txtTrait")),
           # ),
           # fluidRow(
           column(width = 2,
                  plotly::plotlyOutput(ns("PlotlyViolinDMPNearRange"))
                  # )
           )
         )
      ),
      tabPanel("Table",
          DT::dataTableOutput(ns("PcPDMPNearRangeData"))
      )
    )
  )
}

getDMPNearRangeprobeID <- function(DMP, range) {
  a = annotation[order(annotation$chromosome, annotation$position),]
  position = which(a$name == DMP)
  positionStart = position - range
  positionEnd = position + range 
  probeIDs = a[positionStart:positionEnd,]$name
}

getDMPNearRange <- function(sessionVariables, DMP, range) {
  tryCatch({
#browser()
    trait = sessionVariables$trait$trait
    exposure<-gsub("adj","",trait)
    #get exposure data
    traitVar = traitDF(sessionVariables)
    #get DMP range data
    probeIDs = getDMPNearRangeprobeID(DMP,range)
    DMPNearRangeData = as.data.frame(beta.t[,probeIDs])
#    DMPNearRangeData$ID_Kind = rownames(DMPNearRangeData)
    DMPNearRangeData$ID = rownames(DMPNearRangeData)
    #merge all
#    traitDMPNearRangeData = base::merge(traitVar, DMPNearRangeData, by.x = config$mergeAttribut, by.y = config$mergeAttribut, all.x = FALSE, all.y=FALSE)
    traitDMPNearRangeData = base::merge(traitVar, DMPNearRangeData, by.x = config$mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)
    
#    rownames(traitDMPNearRangeData) = traitDMPNearRangeData$ID_Kind
    rownames(traitDMPNearRangeData) = traitDMPNearRangeData$ID
  }, error=function(err){
    print(paste0("unable find near range for ", DMP, exposure, ". - ", err$message))
  });
  return(traitDMPNearRangeData)
}

plotDNAmProfile_SERVER <- function(id, sessionVariables) {
  moduleServer(id, function(input, output, session) {
#    selectedProbe <- reactiveValues(probe = "", probes = list())
#    moduleVariables <- reactiveValues(numberResults = 10)
    
    id <- showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    plotPcPForDMP(input, output, sessionVariables)
#    observeEvent(input$btnDMRWindow, ignoreInit = FALSE, {
    observeEvent(input$btnDMRWindow, ignoreInit = TRUE, {
      plotPcPForDMP(input, output, sessionVariables)
    })
    
#    return(selectedProbe)
  })
}

plotPcPForDMP <- function(input, output, sessionVariables) {
#browser()
  trait = sessionVariables$trait$trait
  probe = sessionVariables$probe$probe
  output$txtTrait = renderText(trait)
  DMPNearRange <- getDMPNearRange(sessionVariables, probe, input$DMRWindow)
  output$PlotlyPcPDMPNearRange <- renderPlotly({
    if (!is_empty(DMPNearRange)) {
      id <- showNotification("Plotting near DMP data...", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      plotlyPcPForDMP(DMPNearRange, sessionVariables, input$DMRWindow)
    }
  })
  
  output$PlotlyViolinDMPNearRange <- plotly::renderPlotly({
    if (!is_empty(DMPNearRange)) {
      id <- showNotification("Plotting vertical violin plot...", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(id), add = TRUE)
      plotlyViolinForDMP(DMPNearRange)
    }
  })
#browser()
#  DMPNearRange <- addLinkToEWASDataHubToHeader(DMPNearRange)
  DMPNearRange <- addLinkToMRCEWASCatalogToHeader(DMPNearRange)
  output$PcPDMPNearRangeData <- DT::renderDataTable({
    id <- showNotification("Printing data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    datatable(DMPNearRange, escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
              options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv', 'excel','pdf')))
  }, server = FALSE)
  
  # datatable(DMPNearRange, extensions = 'Buttons', height = 400, options = list(scrollY = TRUE, scroller = TRUE, searching = TRUE, ordering = TRUE, dom = 'tBS', buttons = c('copy', 'csv', 'excel','pdf')))
  # datatable(DMPNearRange, extensions = c('Scroller', 'Buttons'), height = 400, options = list(scrollY = TRUE, scroller = TRUE, searching = TRUE, ordering = TRUE, dom = 'tBS', buttons = c('copy', 'csv', 'excel','pdf')))
  
}

plotlyPcPForDMP <- function(DMPNearRange, sessionVariables, DMRWindow) {
  tryCatch({
    traitName = colnames(DMPNearRange)[3]
    DMP = sessionVariables$probe$probe
    df <- sessionVariables$resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(df)
  gene.symbol = df[which(df$probeID == DMP)]$gene.symbol
#DMPNearRange <- getDMPNearRange(sessionVariables, DMP, input$DMRWindow)
DMPNearRange <- getDMPNearRange(sessionVariables, DMP, DMRWindow)
    DMPNearRange = na.omit(DMPNearRange)
    DMPNearRangeShort = DMPNearRange[,4:ncol(DMPNearRange)]
    dimensionsList=list()
    P_VAL = sessionVariables$resultDataSingleTrait$P_VAL[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]
    DeltaMeth = sessionVariables$resultDataSingleTrait$DeltaMeth[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]
    for (i in 1:ncol(DMPNearRangeShort)) {
      dimension = list(label = colnames(DMPNearRangeShort)[i], values = DMPNearRangeShort[,i],
                       range = c(0, 1))
      dimensionsList = append(dimensionsList,list(dimension))
    }
#browser()
    plot = plot_ly(data = DMPNearRange)
    plot = plot %>% add_trace(type = 'parcoords',
                              line = list(shape = 'spline',
                                          color =  DMPNearRange[,3],
                                          colorscale = 'Jet',
                                          showscale = TRUE,
                                          reversescale = TRUE,
                                          cmin = min(DMPNearRange[,3],na.rm=TRUE),
                                          cmax = max(DMPNearRange[,3],na.rm=TRUE)),
                              dimensions = dimensionsList
    )
    plot = plot %>% layout(
      title = paste0(traitName, " vs. ", DMP, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
      xaxis = list(
        title = "Location",
        showgrid = F,
        tickangle = 45),
      yaxis = list(
        title = "% Methylation"),
      coloraxis = list(
#      coloraxis.colorbar = list(
        title = traitName)
    )
    return (plot)
  }, error=function(err){
    print(paste0("unable to print pc plot; ", err$message))
    return(empty_plot(err$message))
  });
}

plotlyViolinForDMP <- function(DMPNearRange) {
  tryCatch({
#browser()
    DMPNearRange = na.omit(DMPNearRange)
    dens <- density(DMPNearRange[,3], bw = "sj")
#    femaleDMPNearRange = DMPNearRange[(DMPNearRange$gender == 'f' || DMPNearRange$gender == 'w'),]
    femaleDMPNearRange = DMPNearRange[DMPNearRange$gender == config$genderFemaleValue,] #femaleDMPNearRange = DMPNearRange[DMPNearRange$gender == 'w',]
    maleDMPNearRange = DMPNearRange[DMPNearRange$gender == config$genderMaleValue,] #maleDMPNearRange = DMPNearRange[DMPNearRange$gender == 'm',]
    femaleDens <- density(femaleDMPNearRange[,3], bw = "sj")
    maleDens <- density(maleDMPNearRange[,3], bw = "sj")
    plot = plot_ly()
    plot = plot %>% add_trace(x = femaleDens$y, y = femaleDens$x, type = 'scatter', mode = 'spline', color = I('deeppink'), fill = 'tozerox', name = 'female')
    plot = plot %>% add_trace(x = maleDens$y * -1, y = maleDens$x, type = 'scatter', mode = 'spline', color = I('blue'), fill = 'tozerox', name = 'male')
    #omit x axis
    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
    plot <- plot %>% layout(xaxis = ax)
    #    plot
    return (plot)
  }, error=function(err){
#browser()
    print(paste0("unable to plot violin; ", err$message))
    return(empty_plot(err$message))
  });
}
