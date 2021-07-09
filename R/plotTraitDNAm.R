plotTraitDNAm_UI <- function(id){
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::fluidRow(
      plotly::plotlyOutput(ns("plotlyOneProbe"), height = 400),
      plotly::plotlyOutput(ns("plotlyHorizontalViolin"), height = 100)
    )
  )
}

plotTraitDNAm_SERVER <- function(id, globalVariables, sessionVariables) {
  shiny::moduleServer(id, function(input, output, session) {
    id <- shiny::showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)

    trait = sessionVariables$trait$trait
    probe = sessionVariables$probe$probe
    output$plotlyOneProbe <- plotly::renderPlotly(printScatterPlotlyForOneProbeID(globalVariables, sessionVariables))
    output$plotlyHorizontalViolin <- plotly::renderPlotly(plotlyHorizontalViolin(traitDF(sessionVariables)))
  })
}

printScatterPlotlyForOneProbeID<-function(globalVariables, sessionVariables){
  tryCatch({
    probeID = sessionVariables$probe$probe
    trait = sessionVariables$trait$trait
    Name<-as.character(trait)
    Name<-gsub("adj","",Name)
    probeID<-as.character(probeID)
    #if we have a data.table here, rownames are missing and beta[probeID,] does not work
    #solution: use data.frame instead
    beta_single<-globalVariables$beta[probeID,]
    beta_single<-t(beta_single)
    beta_single<-(as.data.frame(beta_single))
    beta_single$ID = rownames(beta_single)
    #pick variable from traits dataframe
    traitVar<-traitDF(sessionVariables)
    traitName = sessionVariables$trait$trait
    DMP = sessionVariables$probe$probe
    P_VAL = sessionVariables$resultDataSingleTrait$P_VAL[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]
    DeltaMeth = sessionVariables$resultDataSingleTrait$DeltaMeth[sessionVariables$resultDataSingleTrait$probeID == sessionVariables$probe$probe]

    #merge with pheno
#    plotData = base::merge(beta_single, traitVar, by.x = config$mergeAttribut, by.y = config$mergeAttribut, all.x = FALSE, all.y=FALSE)
    plotData = base::merge(beta_single, traitVar, by.x = "ID", by.y = globalVariables$config$mergeAttribut, all.x = FALSE, all.y=FALSE)
    plotData <- stats::na.omit(plotData)
    gender <- "factor(gender)"
    plotDataFemale = subset(plotData, plotData$gender == globalVariables$config$genderFemaleValue)
    plotDataMale = subset(plotData, plotData$gender == globalVariables$config$genderMaleValue)

    fmla = stats::as.formula(paste0("`", probeID, "` ~ `", Name,"`"))
    m <- stats::lm(fmla, data = plotData)
    plot = broom::augment(m,se_fit=TRUE) %>%
      plotly::plot_ly(x = stats::as.formula(paste0("~ `",Name,"`")), showlegend = FALSE) %>%
      plotly::add_markers(y = stats::as.formula(paste0("~ `",probeID,"`")), color = I("black"), name = "all", showlegend = TRUE)%>%
      plotly::add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                  ymax = ~.fitted + 1.96 * .se.fit,
                  color = I("gray80"), showlegend = FALSE)%>%
      plotly::add_lines(y = ~.fitted, color = I("steelblue"), showlegend = FALSE) %>%
      plotly::layout(
        title = paste0(traitName, " vs. ", DMP, " P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth)
      )
    mFemale <- stats::lm(fmla, data = plotDataFemale)
    plotFemale = broom::augment(mFemale,se_fit=TRUE) %>%
      plotly::plot_ly(x = stats::as.formula(paste0("~ `", Name, "`")), showlegend = FALSE) %>%
      plotly::add_markers(y = stats::as.formula(paste0("~ `", probeID, "`")), color = I("deeppink"), name = "female", showlegend = TRUE)%>%
      plotly::add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                  ymax = ~.fitted + 1.96 * .se.fit,
                  color = I("gray80"), showlegend = FALSE)%>%
      plotly::add_lines(y = ~.fitted, color = I("deeppink"), showlegend = FALSE)
    mMale <- stats::lm(fmla, data = plotDataMale)
    plotMale = broom::augment(mMale,se_fit=TRUE) %>%
      plotly::plot_ly(x = stats::as.formula(paste0("~ `", Name, "`")), showlegend = FALSE) %>%
      plotly::add_markers(y = stats::as.formula(paste0("~ `", probeID, "`")), color = I("blue"), name = "male", showlegend = TRUE)%>%
      plotly::add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                  ymax = ~.fitted + 1.96 * .se.fit,
                  color = I("gray80"), showlegend = FALSE)%>%
      plotly::add_lines(y = ~.fitted, color = I("blue"), showlegend = FALSE)
    plot = plotly::subplot(plot,plotFemale,plotMale,nrows=3)
    return(plot)
  }, error=function(err){
    message(paste0(Sys.time(), "unable to plot ", probeID, " vs. ", Name, "; ", err$message))
    return(empty_plot(err$message))
  });
}

plotlyHorizontalViolin <- function(traitDF) {
  tryCatch({
    traitDF = stats::na.omit(traitDF)
    dens <- stats::density(traitDF[,3], bw = "sj")
    femaletraitDF = traitDF[(traitDF$gender == globalVariables$config$genderFemaleValue),]
    maletraitDF = traitDF[traitDF$gender == globalVariables$config$genderMaleValue,]
    femaleDens <- stats::density(femaletraitDF[,3], bw = "sj")
    maleDens <- stats::density(maletraitDF[,3], bw = "sj")
    plot = plotly::plot_ly()
    plot = plot %>% plotly::add_trace(x = femaleDens$x, y = femaleDens$y, type = 'scatter', mode = 'spline', color = I('deeppink'), fill = 'tozeroy', name = 'female')
    plot = plot %>% plotly::add_trace(x = maleDens$x, y = maleDens$y * -1, type = 'scatter', mode = 'spline', color = I('blue'), fill = 'tozeroy', name = 'male')
    #omit x axis
    ay <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
    plot <- plot %>% plotly::layout(yaxis = ay)
    #    plot
    return (plot)
  }, error=function(err){
    print(paste0("unable to print pc plot; ",err$message))
    return(empty_plot(err$message))
  });
}
