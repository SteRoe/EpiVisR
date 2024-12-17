plotDNAmProfile_UI <- function(id){
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::fluidRow(
      shiny::column(12, htmltools::tags$html(tags$body(h4('DMR window size'))),
       shiny::sliderInput(ns("DMRWindow"), "", #"DMR window size",
                   1, 50, 5, 1, width = "100%"),
       shinyBS::bsTooltip(id = ns("DMRWindow"), title = "size of window around selected CpG", placement = "right", trigger = "hover")
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

# getDMPNearRangeprobeIDSessionVariables <- function(globalVariables, DMP, range) {
#   probeID = globalVariables$annotation[order(globalVariables$annotation$chromosome, globalVariables$annotation$position),]
#   result <- getDMPNearRangeprobeID(probeID, DMP, range)
#   return(result)
# }

# getDMPNearRangeprobeID <- function(probeID, DMP, range) {
#   #probeID = globalVariables$annotation[order(globalVariables$annotation$chromosome, globalVariables$annotation$position),]
#   position = which(probeID$name == DMP)
#   positionStart = position - range
#   positionEnd = position + range
#   probeIDs = probeID[positionStart:positionEnd,]$name
#   return(probeIDs)
# }

#' gets a data frame with the near range probeIDs of a DMR
#' @description gets a data structure with the near range of a DMR
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @param DMP CpG to find the range around
#' @param range number of CpG to the left and to the right to include to range (truncated at chromosome end)
#' @return data.frame
#' @export
getDMPNearRangeprobeID <- function(annotation, DMP, range) {
  probeID = annotation[order(annotation$chromosome, annotation$position),]
  position = which(probeID$name == DMP)
  positionStart = position - range
  positionEnd = position + range
  probeIDs = probeID[positionStart:positionEnd,]$name
  return(probeIDs)
}

getDMPNearRangeSessionVariables <- function(globalVariables, sessionVariables, range) {
  DMP = sessionVariables$probe$probe
  traitVar <- traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$sexAttribut)
  probeIDs <- getDMPNearRangeprobeID(globalVariables$annotation, DMP, range)
  DMPNearRangeData <- as.data.frame(globalVariables$beta.t[, probeIDs])
  mergeAttribut <- globalVariables$config$mergeAttribut
  result <- getDMPNearRange(DMP, DMPNearRangeData, traitVar, mergeAttribut, range)
  return(result)
}

#' gets a data structure with the near range of a DMR
#' @description gets a data structure with the near range of a DMR
#' @param DMP #tbc()a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe #tbc  probe (in the center), that the parallel coordinate plot should show; cg number
#' @param range trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @return data.frame
#' @export
#getDMPNearRange <- function(globalVariables, sessionVariables, range) {
getDMPNearRange <- function(DMP, DMPNearRangeData, traitVar, mergeAttribut, range) {
  tryCatch({
    #trait = sessionVariables$trait$trait
    #DMP = sessionVariables$probe$probe
    #trait<-gsub("adj","",trait)
    #get trait data
    #traitVar = traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$sexAttribut)
    #get DMP range data
    #probeIDs = getDMPNearRangeprobeIDSessionVariables(globalVariables, DMP,range)
    #probeIDs = getDMPNearRangeprobeID(globalVariables$annotation, DMP,range)
    #DMPNearRangeData = as.data.frame(globalVariables$beta.t[,probeIDs])
    DMPNearRangeData$ID <- rownames(DMPNearRangeData)
    #merge all
    traitDMPNearRangeData <- base::merge(traitVar, DMPNearRangeData, by.x = mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)

    rownames(traitDMPNearRangeData) = traitDMPNearRangeData$ID
  }, error=function(err){
    #print(paste0("unable find near range for ", DMP, " <-> ", trait, ". - ", err$message))
    print(paste0("unable find near range for ", DMP, " <-> ", traitVar, ". - ", err$message))
  });
  return(traitDMPNearRangeData)
}

#' Plots the near range of a DMR
#' @description plots the near range of a DMR defined inside session variables
plotDNAmProfile_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    reDMPNearRange <- shiny::reactive({getDMPNearRangeSessionVariables(globalVariables, sessionVariables, input$DMRWindow)})

    output$PlotlyPcPDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        id <- shiny::showNotification("Plotting near DMP data...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(id), add = TRUE)
        plotlyPcPForDMPSessionVariables(globalVariables, sessionVariables, DMPNearRange)
      }
    })

    output$PlotlyViolinDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        id <- shiny::showNotification("Plotting vertical violin plot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(id), add = TRUE)
        plotlyViolinForDMPSessionVariables(globalVariables, DMPNearRange)
      }
    })

    output$PcPDMPNearRangeData <- DT::renderDataTable({
      id <- shiny::showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      DMPNearRange = reDMPNearRange()
      DMPNearRange <- addLinkToMRCEWASCatalogToHeader(DMPNearRange, globalVariables$config$baseURL_MRCEWASCatalog)
#      colnames(DMPNearRange) <- stringr::str_to_title(colnames(DMPNearRange))
      DT::datatable(DMPNearRange, escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
                options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv'))) %>%
      DT::formatSignif(3:ncol(DMPNearRange), digits = 2)
    }, server = FALSE, escape = FALSE)
  })
}

#' Plots the near range of a DMR
#' @description plots the near range of a DMR and gives a plotly object back using session variables
plotlyPcPForDMPSessionVariables <- function(globalVariables, sessionVariables, DMPNearRange) {
  probe <- sessionVariables$probe$probe
  resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
  annotation <- globalVariables$annotation
  result <- plotlyPcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation)
  #result <- plotlyPcPForDMP(DMPNearRange, probe, P_VAL, DeltaMeth, annotation)
  return(result)
}

#' Plots the near range of a DMR
#' @description plots the near range of a DMR in a parallel coordinate plot and gives a plotly object back
#' @param DMPNearRange a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe probe (in the center), that the parallel coordinate plot should show; cg number
#' @param resultDataSingleTrait trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @param selection a list of cg numbers to mark within the methylation profile plot
#' @param shortlabel if TRUE only the cg number is plotted into the methylation profile plot
#' @return plotly object
#' @export
plotlyPcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation, selection, shortlabel = TRUE) {
  tryCatch({
    traitName <- colnames(DMPNearRange)[3]
    df <- resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(annotation, df)
    gene.symbol <- df[which(df$probeID == probe),]$gene.symbol
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    dimensionsList=list()
    P_VAL <- resultDataSingleTrait$P_VAL[resultDataSingleTrait$probeID == probe]
    DeltaMeth <- resultDataSingleTrait$DeltaMeth[resultDataSingleTrait$probeID == probe]
    for (i in 1:ncol(DMPNearRangeShort)) {
      lblCpG <- colnames(DMPNearRangeShort)[i]
      if (!missing(shortlabel) && shortlabel != TRUE) {
        lblP <- signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$P_VAL,3)
        if (shiny::isTruthy(lblP)) {
          lblP <- paste0(",\n p: ", lblP)
          lblDM <- paste0(",\n d: ", signif(df[which(df$probeID == colnames(DMPNearRangeShort)[i]),]$DeltaMeth,3))
        }
        else {
          lblP <- ",\n p: n.s."
          lblDM <- ""
        }
        lblSym <- paste0(",\n sbl:", annotation[which(annotation$name == colnames(DMPNearRangeShort)[i]),]$gene.symbol)
        lblPos <- paste0(",\n pos:", annotation[which(annotation$name == colnames(DMPNearRangeShort)[i]),]$position)
        label <- paste0(lblCpG, lblP, lblDM, lblSym, lblPos)
      }
      else {
        label <- lblCpG
      }
      markCpG <- FALSE
      if (!missing(selection)) {
        if (lblCpG %in% selection) {
          markCpG <- TRUE
        }
      }
      if (markCpG) {
        dimension <- list(label = label,
                          values = DMPNearRangeShort[,i],
                          range = c(0, 1),
                          constraintrange = c(0, 1))
      }
      else {
        dimension <- list(label = label,
                          values = DMPNearRangeShort[,i],
                          range = c(0, 1))
      }
      dimensionsList <- append(dimensionsList, list(dimension))
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
                                       dimensions = dimensionsList,
                                       labelangle = 270
    )
    plot <- plot %>% plotly::layout(
      title = paste0(traitName, " vs. ", probe, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
      xaxis = list(
        title = "Location",
        showgrid = F,
        tickangle = 45),
      yaxis = list(
        title = "Methylation [%]"),
      coloraxis = list(
        title = traitName)
    )
    plot <- plot %>%
      plotly::add_annotations(
        text = "Methylation [%]",
        x = 0,
        y = 0.5,
        yref = "paper",
        xref = "paper",
        xanchor = "center",
        yanchor = "center",
        xshift = -50,
        showarrow = FALSE,
        textangle = 270,
        font = list(size = 15)
      )
    plot <- plot %>%
      plotly::add_annotations(
        text = "globalArrayPosition",
        x = 0.5,
        y = 0,
        yref = "paper",
        xref = "paper",
        xanchor = "center",
        yanchor = "center",
        yshift = -35,
        showarrow = FALSE,
        textangle = 0,
        font = list(size = 15)
      )
    return (plot)
  }, error=function(err){
    print(paste0("unable to print pc plot: ", err$message))
    return(empty_plot(err$message))
  });
}

plotlyViolinForDMPSessionVariables <- function(globalVariables, DMPNearRange) {
  sexFemaleValue <- globalVariables$config$sexFemaleValue
  sexMaleValue <- globalVariables$config$sexMaleValue
  result <- plotlyViolinForDMP(sexFemaleValue, sexMaleValue, DMPNearRange)
  return(result)
}

#' Plots the violin plot for both sexes of a near range of a DMR
#' @description plots the violin plot for both sexes of a near range of a DMR and gives a plotly object back
#' @export
plotlyViolinForDMP <- function(sexFemaleValue, sexMaleValue, DMPNearRange) {
#plotlyViolinForDMP <- function(globalVariables, DMPNearRange) {
  tryCatch({
    DMPNearRange = stats::na.omit(DMPNearRange)
    min <- min(DMPNearRange[,3])
    max <- max(DMPNearRange[,3])
    dens <- stats::density(DMPNearRange[,3], bw = "sj")
    femaleDMPNearRange <- DMPNearRange[DMPNearRange$sex == sexFemaleValue,]
    maleDMPNearRange <- DMPNearRange[DMPNearRange$sex == sexMaleValue,]
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
