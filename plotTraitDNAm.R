plotTraitDNAm_UI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      plotly::plotlyOutput(ns("plotlyOneProbe"), height = 400),
      plotly::plotlyOutput(ns("plotlyHorizontalViolin"), height = 100)
    )
  )
}

plotTraitDNAm_SERVER <- function(id, sessionVariables) {
  moduleServer(id, function(input, output, session) {

    id <- showNotification("Plotting data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)

    trait = sessionVariables$trait$trait
    probe = sessionVariables$probe$probe
    output$plotlyOneProbe <- plotly::renderPlotly(printScatterPlotlyForOneProbeID(sessionVariables))
#    output$plotlyHorizontalViolin <- plotly::renderPlotly(plotlyHorizontalViolin(traitDF(trait)))
    output$plotlyHorizontalViolin <- plotly::renderPlotly(plotlyHorizontalViolin(traitDF(sessionVariables)))
  })
}

printScatterPlotlyForOneProbeID<-function(sessionVariables){
  tryCatch({
    probeID = sessionVariables$probe$probe
    trait = sessionVariables$trait$trait
    Name<-as.character(trait)
    Name<-gsub("adj","",Name)
    probeID<-as.character(probeID)
    #das funktioniert nicht, weil aus irgendeinem Grund schon wieder die rownames fehlen!!!
    #Lösung: es war data.table, dann geht beta[probeID,] nicht
    beta_single<-beta[probeID,]
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
    plotData = base::merge(beta_single, traitVar, by.x = "ID", by.y = config$mergeAttribut, all.x = FALSE, all.y=FALSE)
    plotData <- na.omit(plotData)
    gender <- "factor(gender)"
    plotDataFemale = subset(plotData, plotData$gender == config$genderFemaleValue)
    plotDataMale = subset(plotData, plotData$gender == config$genderMaleValue)
    
    fmla = as.formula(paste0("`", probeID, "` ~ `", Name,"`"))
    m <- lm(fmla, data = plotData)
    plot = broom::augment(m,se_fit=TRUE) %>%
      plot_ly(x = as.formula(paste0("~ `",Name,"`")), showlegend = FALSE) %>%
      add_markers(y = as.formula(paste0("~ `",probeID,"`")), color = I("black"), name = "all", showlegend = TRUE)%>%
      add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                  ymax = ~.fitted + 1.96 * .se.fit,
                  color = I("gray80"), showlegend = FALSE)%>%
      add_lines(y = ~.fitted, color = I("steelblue"), showlegend = FALSE) %>%
      layout(
        title = paste0(traitName, " vs. ", DMP, " P_VAL: ", P_VAL, " DeltaMeth: ", DeltaMeth)
      )
    mFemale <- lm(fmla, data = plotDataFemale)
    plotFemale = broom::augment(mFemale,se_fit=TRUE) %>%
      plot_ly(x = as.formula(paste0("~ `", Name, "`")), showlegend = FALSE) %>%
      add_markers(y = as.formula(paste0("~ `", probeID, "`")), color = I("deeppink"), name = "female", showlegend = TRUE)%>%
      add_ribbons(ymin = ~.fitted - 1.96 * .se.fit, 
                  ymax = ~.fitted + 1.96 * .se.fit, 
                  color = I("gray80"), showlegend = FALSE)%>%
      add_lines(y = ~.fitted, color = I("deeppink"), showlegend = FALSE)
    mMale <- lm(fmla, data = plotDataMale)
    plotMale = broom::augment(mMale,se_fit=TRUE) %>%
      plot_ly(x = as.formula(paste0("~ `", Name, "`")), showlegend = FALSE) %>%
      add_markers(y = as.formula(paste0("~ `", probeID, "`")), color = I("blue"), name = "male", showlegend = TRUE)%>%
      add_ribbons(ymin = ~.fitted - 1.96 * .se.fit,
                  ymax = ~.fitted + 1.96 * .se.fit,
                  color = I("gray80"), showlegend = FALSE)%>%
      add_lines(y = ~.fitted, color = I("blue"), showlegend = FALSE)
    plot = subplot(plot,plotFemale,plotMale,nrows=3)
    return(plot)
  }, error=function(err){
    print(paste0("unable to plot ",probeID," vs. ", Name, "; ",err$message))
    return(empty_plot(err$message))
  });
}

plotlyHorizontalViolin <- function(traitDF) {
  tryCatch({
    traitDF = na.omit(traitDF)
    dens <- density(traitDF[,3], bw = "sj")
    femaletraitDF = traitDF[(traitDF$gender == config$genderFemaleValue),]
    maletraitDF = traitDF[traitDF$gender == config$genderMaleValue,]
    femaleDens <- density(femaletraitDF[,3], bw = "sj")
    maleDens <- density(maletraitDF[,3], bw = "sj")
    plot = plot_ly()
    plot = plot %>% add_trace(x = femaleDens$x, y = femaleDens$y, type = 'scatter', mode = 'spline', color = I('deeppink'), fill = 'tozeroy', name = 'female')
    plot = plot %>% add_trace(x = maleDens$x, y = maleDens$y * -1, type = 'scatter', mode = 'spline', color = I('blue'), fill = 'tozeroy', name = 'male')
    #omit x axis
    ay <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      showticklabels = FALSE,
      showgrid = FALSE
    )
    plot <- plot %>% layout(yaxis = ay)
    #    plot
    return (plot)
  }, error=function(err){
    print(paste0("unable to print pc plot; ",err$message))
    return(empty_plot(err$message))
  });
}