plotDNAmProfile_UI <- function(id){
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::fluidRow(
      shiny::column(12, htmltools::tags$html(tags$body(h4('DMR window size'))),
       shiny::sliderInput(ns("DMRWindow"), "", #"DMR window size",
                   1, 50, 5, 1, width = "100%")
      )
    ),
    shiny::fluidRow(
      shiny::column(width = 12,
        shiny::tabsetPanel(
          shiny::tabPanel("Visualisation",
            shiny::fluidRow(
              shiny::column(width = 10,
                      plotly::plotlyOutput(ns("PlotlyPcPDMPNearRange"),
                                           width = "100%",
                                           height = "800px")
               ),
              shiny::column(width = 2,
                      plotly::plotlyOutput(ns("PlotlyViolinDMPNearRange"),
                                           width = "100%",
                                           height = "800px")
               )
             )
          ),
          shiny::tabPanel("Table",
              DT::dataTableOutput(ns("PcPDMPNearRangeData"))
          )
        )
      )
    )
  )
}

getDMPNearRangeprobeID <- function(globalVariables, DMP, range) {
  a = globalVariables$annotation[order(globalVariables$annotation$chromosome, globalVariables$annotation$position),]
  position = which(a$name == DMP)
  positionStart = position - range
  positionEnd = position + range
  probeIDs = a[positionStart:positionEnd,]$name
}

getDMPNearRange <- function(globalVariables, sessionVariables, range) {
  tryCatch({
    trait = sessionVariables$trait$trait
    DMP = sessionVariables$probe$probe
    trait<-gsub("adj","",trait)
    #get trait data
    traitVar = traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut)
    #get DMP range data
    probeIDs = getDMPNearRangeprobeID(globalVariables, DMP,range)
    DMPNearRangeData = as.data.frame(globalVariables$beta.t[,probeIDs])
    DMPNearRangeData$ID = rownames(DMPNearRangeData)
    #merge all
    traitDMPNearRangeData = base::merge(traitVar, DMPNearRangeData, by.x = globalVariables$config$mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)

    rownames(traitDMPNearRangeData) = traitDMPNearRangeData$ID
  }, error=function(err){
    print(paste0("unable find near range for ", DMP, " <-> ", trait, ". - ", err$message))
  });
  return(traitDMPNearRangeData)
}

plotDNAmProfile_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    reDMPNearRange <- shiny::reactive({getDMPNearRange(globalVariables, sessionVariables, input$DMRWindow)})

    output$PlotlyPcPDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        id <- shiny::showNotification("Plotting near DMP data...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(id), add = TRUE)
        plotlyPcPForDMP(globalVariables, sessionVariables, DMPNearRange)
      }
    })

    output$PlotlyViolinDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        id <- shiny::showNotification("Plotting vertical violin plot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(id), add = TRUE)
        plotlyViolinForDMP(globalVariables, DMPNearRange)
      }
    })

    output$PcPDMPNearRangeData <- DT::renderDataTable({
      id <- shiny::showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      DMPNearRange = reDMPNearRange()
      DMPNearRange <- addLinkToMRCEWASCatalogToHeader(DMPNearRange, globalVariables$config$baseURL_MRCEWASCatalog)
      DT::datatable(DMPNearRange, escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
                options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv')))
    }, server = FALSE)
  })
}

plotlyPcPForDMP <- function(globalVariables, sessionVariables, DMPNearRange) {
  tryCatch({
    traitName <- colnames(DMPNearRange)[3]
    DMP <- sessionVariables$probe$probe
    df <- sessionVariables$resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(globalVariables$annotation, df)
    gene.symbol <- df[which(df$probeID == DMP),]$gene.symbol
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    dimensionsList=list()
    P_VAL <- sessionVariables$resultDataSingleTrait$P_VAL[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]
    DeltaMeth <- sessionVariables$resultDataSingleTrait$DeltaMeth[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]
    for (i in 1:ncol(DMPNearRangeShort)) {
      lblCpG = colnames(DMPNearRangeShort)[i]
      lblP = signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$P_VAL,3)
      if (shiny::isTruthy(lblP)) {
        lblP = paste0(",\n p: ", lblP)
        lblDM = paste0(",\n d: ", signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$DeltaMeth,3))
      }
      else {
        lblP = ",\n p: n.s."
        lblDM = ""
      }
      lblSym = paste0(",\n sbl:", globalVariables$annotation[which(globalVariables$annotation$name == colnames(DMPNearRangeShort)[i]),]$gene.symbol)
      lblPos = paste0(",\n pos:", globalVariables$annotation[which(globalVariables$annotation$name == colnames(DMPNearRangeShort)[i]),]$position)
      label = paste0(lblCpG, lblP, lblDM, lblSym, lblPos)
      dimension = list(label = label, values = DMPNearRangeShort[,i],
                       range = c(0, 1))
      dimensionsList = append(dimensionsList,list(dimension))
    }
    plot <- plotly::plot_ly(data = DMPNearRange)
    plot <- plot %>% plotly::add_trace(type = 'parcoords',
                              line = list(shape = 'spline',
                                          color =  DMPNearRange[,3],
                                          colorscale = 'Jet',
                                          showscale = TRUE,
                                          reversescale = TRUE,
                                          cmin = min(DMPNearRange[,3],na.rm=TRUE),
                                          cmax = max(DMPNearRange[,3],na.rm=TRUE)),
                              dimensions = dimensionsList
    )
    plot <- plot %>% plotly::layout(
      title = paste0(traitName, " vs. ", DMP, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
      xaxis = list(
        title = "Location",
        showgrid = F,
        tickangle = 45),
      yaxis = list(
        title = "Methylation [%]"),
      coloraxis = list(
        title = traitName)
    )
    return (plot)
  }, error=function(err){
    print(paste0("unable to print pc plot; ", err$message))
    return(empty_plot(err$message))
  });
}

plotlyViolinForDMP <- function(globalVariables, DMPNearRange) {
  tryCatch({
    DMPNearRange = stats::na.omit(DMPNearRange)
    min <- min(DMPNearRange[,3])
    max <- max(DMPNearRange[,3])
    dens <- stats::density(DMPNearRange[,3], bw = "sj")
    femaleDMPNearRange <- DMPNearRange[DMPNearRange$gender == globalVariables$config$genderFemaleValue,] #femaleDMPNearRange = DMPNearRange[DMPNearRange$gender == 'w',]
    maleDMPNearRange <- DMPNearRange[DMPNearRange$gender == globalVariables$config$genderMaleValue,] #maleDMPNearRange = DMPNearRange[DMPNearRange$gender == 'm',]
    femaleDens <- stats::density(femaleDMPNearRange[,3], bw = "sj")
    maleDens <- stats::density(maleDMPNearRange[,3], bw = "sj")
    plot <- plotly::plot_ly()
    plot <- plot %>% plotly::add_trace(x = femaleDens$y, y = femaleDens$x, type = 'scatter', mode = 'spline', color = I('deeppink'), fill = 'tozerox', name = 'female')
    plot <- plot %>% plotly::add_trace(x = maleDens$y * -1, y = maleDens$x, type = 'scatter', mode = 'spline', color = I('blue'), fill = 'tozerox', name = 'male')
    #omit x axis
    ax <- list(
      title = "Density"
#      zeroline = FALSE,
#      showline = FALSE,
#      showticklabels = FALSE,
#      showgrid = FALSE
    )
    ay <- list(
      title = colnames(DMPNearRange)[3],
      range = c(min,max)
      #      zeroline = FALSE,
      #      showline = FALSE,
      #      showticklabels = FALSE,
      #      showgrid = FALSE
    )
    plot <- plot %>% plotly::layout(xaxis = ax, yaxis = ay)
    #    plot
    return (plot)
  }, error=function(err){
    print(paste0("unable to plot violin; ", err$message))
    return(empty_plot(err$message))
  });
}
