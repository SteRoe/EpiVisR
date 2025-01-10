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
              shiny::fluidRow( #insert numerical overview and table headers here
                shiny::column(width = 12,
                  shiny::verbatimTextOutput(ns("PcPTitle"), placeholder = TRUE) #text with title here
                )
              ),
              shiny::fluidRow(
              shiny::column(width = 10,
                shiny::fluidRow(
                  shiny::column(width = 2, #distance to start of chr
                                shiny::verbatimTextOutput(ns("PcPTitleLeft"), placeholder = TRUE)
                  ),
                  shiny::column(width = 6, #table headers
                                shiny::verbatimTextOutput(ns("PcPTitleMiddle"), placeholder = TRUE) #text with title here
                  ),
                  shiny::column(width = 2, #distance to end of chr
                                shiny::verbatimTextOutput(ns("PcPTitleRight"), placeholder = TRUE)
                  )
                ),
                shiny::fluidRow(
                      plotly::plotlyOutput(ns("PlotlyPcPDMPNearRange"),
                                           width = "100%",
                                           height = "800px")
                )
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

#' gets length in base pairs from a certain chromosome from hg38
#' @description gets length in base pairs from a certain chromosome from hg38
#' @param chr chromosome name as string, e.g. "chr1" or "chrX"
#' @return integer
getChromosomeLength <- function(chr) {
  #library(GenomicFeatures)
  chromosomes <- getChromInfoFromUCSC("hg38")
  chromosomes <- chromosomes[which(chromosomes$assembled == TRUE), ]
  result <- chromosomes[which(chromosomes$chrom == chr), ]$size
}

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
#' @param DMP a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param DMPNearRangeData a data.frame containing the CpGs around a selected CpG
#' @param traitVar a data.frame containing the real world measurements for a certain trait
#' @param mergeAttribut the attribut from traitVar, that is used for merging with DMPNearRangeData
#' @param range trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @return data.frame
#' @export
getDMPNearRange <- function(DMP, DMPNearRangeData, traitVar, mergeAttribut, range) {
  tryCatch({
    DMPNearRangeData$ID <- rownames(DMPNearRangeData)
    #merge all
    traitDMPNearRangeData <- base::merge(traitVar, DMPNearRangeData, by.x = mergeAttribut, by.y = "ID", all.x = FALSE, all.y=FALSE)

    rownames(traitDMPNearRangeData) = traitDMPNearRangeData$ID
  }, error=function(e){
    print(paste0("unable find near range for ", DMP, " <-> ", traitVar, ". - ", e$message))
  });
  return(traitDMPNearRangeData)
}

#' Plots the near range of a DMR
#' @description plots the near range of a DMR and gives a plotly object back using session variables
#' @param globalVariables package global variables
#' @param sessionVariables package session variables
#' @param DMPNearRange range to plot
plotlyPcPForDMPSessionVariables <- function(globalVariables, sessionVariables, DMPNearRange) {
  probe <- sessionVariables$probe$probe
  resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
  annotation <- globalVariables$annotation
  result <- plotlyPcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation)
  return(result)
}

#' generates the left title for a PC plot
#' @param DMPNearRange a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe probe (in the center), that the parallel coordinate plot should show; cg number
#' @param resultDataSingleTrait trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @return text
#' @export
leftTitlePcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation) {
  tryCatch({
    df <- resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(annotation, df)
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    firstProbe <- colnames(DMPNearRangeShort)[1]
    pos <- annotation[which(annotation$name == firstProbe),]$position
    title = paste0(pos, " bp from chr start")
    return(title)
  }, warning = function(w) {
    print(paste0("warning while printing left title for pc plot: ", w$message))
    return(empty_plot(w$message))
  }, error = function(e) {
    print(paste0("unable to print left title for pc plot: ", e$message))
    return(empty_plot(e$message))
  });
}

#' generates the middle title for a PC plot
#' @param DMPNearRange a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe probe (in the center), that the parallel coordinate plot should show; cg number
#' @param resultDataSingleTrait trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @return text
#' @export
middleTitlePcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation) {
  tryCatch({
    df <- resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(annotation, df)
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    firstProbe <- colnames(DMPNearRangeShort)[1]
    lastProbe <- colnames(DMPNearRangeShort)[ncol(DMPNearRangeShort)]
    posStart <- annotation[which(annotation$name == firstProbe),]$position
    posEnd <- annotation[which(annotation$name == lastProbe),]$position
    chrF <- factor(annotation[which(annotation$name == probe),]$chromosome)
    chrLength <- getChromosomeLength(chrF)
    chrRange <- base::round((posEnd - posStart)*100/chrLength, digits = 2) #selected area comprises x% of chromosome...
    title = base::paste0(chrRange, " % bp of chr")
    return(title)
  }, warning = function(w) {
    print(paste0("warning while printing middle title for pc plot: ", w$message))
    return(empty_plot(w$message))
  }, error = function(e) {
    print(paste0("unable to print middle title for pc plot: ", e$message))
    return(empty_plot(e$message))
  });
}

#' generates the right title for a PC plot
#' @param DMPNearRange a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe probe (in the center), that the parallel coordinate plot should show; cg number
#' @param resultDataSingleTrait trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @return text
#' @export
rightTitlePcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation) {
  tryCatch({
    df <- resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(annotation, df)
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    lastProbe <- colnames(DMPNearRangeShort)[ncol(DMPNearRangeShort)]
    pos <- annotation[which(annotation$name == lastProbe),]$position
    chrF <- factor(annotation[which(annotation$name == probe),]$chromosome)
    chrLength <- getChromosomeLength(chrF)
    pos <- chrLength - pos
    title = paste0(pos, " bp from chr end")
    return(title)
  }, warning = function(w) {
    print(paste0("warning while printing right title for pc plot: ", w$message))
    return(empty_plot(w$message))
  }, error = function(e) {
    print(paste0("unable to print right title for pc plot: ", e$message))
    return(empty_plot(e$message))
  });
}

#' generates the title for a PC plot
#' @param DMPNearRange a data.frame containing row.names(ID), IDs (ID_Kind), sex, trait and cgs around probe
#' @param probe probe (in the center), that the parallel coordinate plot should show; cg number
#' @param resultDataSingleTrait trait original methylation data as data.frame: probeID, BETA, SE, P_VAL, FDR, DeltaMeth, N
#' @param annotation annotation data (pre loaded); the annotation from meffil
#' @return text
#' @export
titlePcPForDMP <- function(DMPNearRange, probe, resultDataSingleTrait, annotation) {
  tryCatch({
    traitName <- colnames(DMPNearRange)[3]
    df <- resultDataSingleTrait
    df <- resultDataSingleScenarioWithAnnotation(annotation, df)
    gene.symbol <- df[which(df$probeID == probe),]$gene.symbol
    chr <- df[which(df$probeID == probe),]$chromosome
    pos <- df[which(df$probeID == probe),]$position
    DMPNearRange <- stats::na.omit(DMPNearRange)
    DMPNearRangeShort <- DMPNearRange[,4:ncol(DMPNearRange)]
    P_VAL <- resultDataSingleTrait$P_VAL[resultDataSingleTrait$probeID == probe]
    DeltaMeth <- resultDataSingleTrait$DeltaMeth[resultDataSingleTrait$probeID == probe]
    title = paste0(traitName, " vs. ", probe, " gene.symbol: ", gene.symbol, " chromosome: ", chr, " position: ", pos, " P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth)
    return(title)
  }, warning = function(w) {
    print(paste0("warning while printing title for pc plot: ", w$message))
    return(empty_plot(w$message))
  }, error = function(e) {
    print(paste0("unable to print title for pc plot: ", e$message))
    return(empty_plot(e$message))
  });
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
      CpGidentical <- FALSE
      if (probe == lblCpG) {
        CpGidentical <- TRUE
      }
      #this does not work as expected - with fractions as constraintrange the colorscale in plotly::add_trace is inactive; maybe error in plotly library?
      # if (CpGidentical) {
      #   dimension <- list(label = label,
      #                     values = DMPNearRangeShort[,i],
      #                     range = c(0, 1),
      #                     constraintrange = list(c(0.2,0.4),c(0.6,0.8)))
      # }
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
      #title = paste0(traitName, " vs. ", probe, " gene.symbol: ", gene.symbol ," P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth),
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
#' @param sexFemaleValue value to mark for females
#' @param sexMaleValue value to mark for males
#' @param DMPNearRange plots a violin plot with the distribution of a certain trait contained in DMPNearRange
#' @export
plotlyViolinForDMP <- function(sexFemaleValue, sexMaleValue, DMPNearRange) {
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
    )
    ay <- list(
      title = colnames(DMPNearRange)[3],
      range = c(min,max)
    )
    plot <- plot %>% plotly::layout(xaxis = ax, yaxis = ay)
    return (plot)
  }, error=function(e){
    print(paste0("unable to plot violin; ", e$message))
    return(empty_plot(e$message))
  });
}

#' Plots the near range of a DMR
#' @description plots the near range of a DMR defined inside session variables
#' @param id id to identify against UI
#' @param globalVariables package global variables
#' @param sessionVariables package session variables
plotDNAmProfile_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
    shinyId <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(shinyId), add = TRUE)
    reDMPNearRange <- shiny::reactive({getDMPNearRangeSessionVariables(globalVariables, sessionVariables, input$DMRWindow)})

    output$PcPTitleLeft <- shiny::renderText({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Printing left title for PCPlot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        probe <- sessionVariables$probe$probe
        resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
        annotation <- globalVariables$annotation
        return(leftTitlePcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation))
      }
    })

    output$PcPTitleMiddle <- shiny::renderText({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Printing middle title for PCPlot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        probe <- sessionVariables$probe$probe
        resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
        annotation <- globalVariables$annotation
        return(middleTitlePcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation))
      }
    })

    output$PcPTitleRight <- shiny::renderText({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Printing right title for PCPlot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        probe <- sessionVariables$probe$probe
        resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
        annotation <- globalVariables$annotation
        return(rightTitlePcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation))
      }
    })

    output$PcPTitle <- shiny::renderText({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Printing title for PCPlot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        probe <- sessionVariables$probe$probe
        resultDataSingleTrait <- sessionVariables$resultDataSingleTrait
        annotation <- globalVariables$annotation
        return(titlePcPForDMP(DMPNearRange, probe, resultDataSingleTrait, annotation))
      }
    })

    output$PlotlyPcPDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Plotting near DMP data...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        plotlyPcPForDMPSessionVariables(globalVariables, sessionVariables, DMPNearRange)
      }
    })

    output$PlotlyViolinDMPNearRange <- plotly::renderPlotly({
      DMPNearRange = reDMPNearRange()
      if (!is.null(DMPNearRange)) {
        shinyId <- shiny::showNotification("Plotting vertical violin plot...", duration = NULL, closeButton = FALSE)
        on.exit(shiny::removeNotification(shinyId), add = TRUE)
        plotlyViolinForDMPSessionVariables(globalVariables, DMPNearRange)
      }
    })

    output$PcPDMPNearRangeData <- DT::renderDataTable({
      sinyId <- shiny::showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(shinyId), add = TRUE)
      DMPNearRange = reDMPNearRange()
      DMPNearRange <- addLinkToMRCEWASCatalogToHeader(DMPNearRange, globalVariables$config$baseURL_MRCEWASCatalog)
#      colnames(DMPNearRange) <- stringr::str_to_title(colnames(DMPNearRange))
      DT::datatable(DMPNearRange, escape = F, extensions = c('Scroller', 'Buttons'), style = "bootstrap", class = "compact", width = "100%",
                options = list(pageLength = 10, deferRender = TRUE, scrollY = 300, scrollX = TRUE, scroller = TRUE, dom = 'ftBS', buttons = c('copy', 'csv'))) %>%
      DT::formatSignif(3:ncol(DMPNearRange), digits = 2)
    }, server = FALSE, escape = FALSE)
  })
}
