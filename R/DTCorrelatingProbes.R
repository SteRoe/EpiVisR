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

    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    reDFCorrelatingProbes <- shiny::reactive({correlatingProbes(globalVariables, sessionVariables)})

    output$DTCorrelatingProbes <- DT::renderDataTable({
      id <- shiny::showNotification("Printing data...", duration = NULL, closeButton = FALSE)
      on.exit(shiny::removeNotification(id), add = TRUE)
      tryCatch({
        print(paste0(Sys.time(), " rendering correlating data table."))
        CP <- reDFCorrelatingProbes()
        DT::datatable(CP, escape = F, extensions="Scroller", style="bootstrap", class="compact", width="100%",
                  options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE)) %>%
          DT::formatSignif(3:5, digits = 2)
        #if JavaScript error occurs, this has something to do with temp space on server: restart R session
        # DT::datatable(CP, escape = F, extensions="Scroller", style="bootstrap", class="compact", width="100%",
        #               callback = JS("$.fn.dataTable.ext.errMode = 'none';"),
        #               options=list(pageLength = 5, deferRender=TRUE, scrollY=300, scroller=TRUE))
      }, error = function(err) {
        shiny::validate(shiny::need(nrow(CP)>0,"No data to show"))
      })
    }, server = TRUE)

    shiny::observeEvent(input$DTCorrelatingProbes_cell_clicked, {
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
  #            if (!is_null(dfList)) {
            if (!is.null(dfList)) {
              id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
              on.exit(shiny::removeNotification(id), add = TRUE)
              i = NULL
              foreach(i=1:length(dfList)) %do% {
                df = as.data.frame(dfList[i])
                fmla = stats::as.formula(paste0("`", colnames(df)[4], "` ~ `", colnames(df)[3], "`"))
                m <- stats::lm(fmla, data = df)
                plot <- broom::augment(m,se_fit=TRUE) %>%
                  plotly::plot_ly(x = stats::as.formula(paste0("~ `", colnames(df)[3], "`")), showlegend = FALSE) %>%
                  plotly::add_markers(x = df[,3], y = df[,4], name = colnames(df)[4], showlegend = TRUE) %>%
                  plotly::add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                              ymax = ~.fitted + 1.96 * .se.fit,
                              color = I("gray80"), showlegend = FALSE) %>%
                  plotly::add_lines(y = ~.fitted, color = I("steelblue"), showlegend = FALSE) %>%
                  plotly::layout(
                    yaxis = list(title = '%', range = c(0,1))
#                    yaxis = list(range = c(0,1))
                  )
                plotList = c(plotList,list(plot))
              }
              # m <- list(
              #   l = 100,
              #   r = 0,
              #   b = 0,
              #   t = 0
              # )
              plotlyscatter <- plotly::subplot(plotList, shareX = TRUE, shareY = TRUE, nrows = length(dfList)) %>%
              # plotly::layout(autosize = T, margin = m) %>%
              # plotly::layout(annotations = list(
              #     list(x = -0.15, y = 0.5, text = "Methylation [%]",
              #          xshift = -50,
              #          textangle = 270,
              #          showarrow = F, xref='paper', yref='paper', size=48)
              #   ))
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
      #merge with annotation
      correlatingProbes <- dplyr::left_join(correlatingProbes, globalVariables$annotation, by = c("probeID" = "name"))
      correlatingProbes$type <- NULL
      correlatingProbes$target <- NULL
      correlatingProbes$meth.dye <- NULL
      #add DeltaMeth
      DeltaMeth <- currentData[,c("probeID","DeltaMeth","P_VAL","FDR")]
      correlatingProbes <- dplyr::inner_join(DeltaMeth, correlatingProbes, by = c("probeID" = "probeID"))
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
    if (nrow(correlatingProbes) > 100) {
      correlatingProbes <- correlatingProbes[1:100,]
    }
    return (correlatingProbes)
  }
}
