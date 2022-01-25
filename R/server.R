server <- function(input, output, session) {
  globalVariables <- list()
  print(paste0(Sys.time(), " loading configuration."))
  globalVariables$config <- config::get(file = "config.yml")
<<<<<<< HEAD
#  output$txtdebug <- shiny::renderText(globalVariables$config$debugMode)
  #  checkconfigVariables()
  sessionVariables <- shiny::reactiveValues(folder = "", trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), traitDF = data.frame())
=======
  output$txtdebug <- shiny::renderText(globalVariables$config$debugMode)
  #  checkconfigVariables()
  #  sessionVariables <- reactiveValues(folder = "", dataFileName = "", trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), traitDF = data.frame(), traitsDFLong = data.frame())
  sessionVariables <- shiny::reactiveValues(folder = "", trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), traitDF = data.frame())
#  sessionVariables$dataFileName = globalVariables$config$traitFileName
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80

  print(paste0(Sys.time(), " loading objects."))
  globalVariables = loadObjects(globalVariables)
#  assign("globalVariables", globalVariables, envir = .GlobalEnv)
#  checkConfigVariables()
  print(paste0(Sys.time(), " making session variables."))
  print(paste0(Sys.time(), " reading inputTrait_SERVER."))
  inputTrait_SERVER("trait", globalVariables, sessionVariables)

  shiny::observeEvent(input$btnSelectTrait, { # ignoreInit = TRUE,
    id <- shiny::showNotification(paste0("Selected trait: ",sessionVariables$trait$trait), duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    #load DF from selected trait
    if (shiny::isTruthy(sessionVariables$trait$trait)) {
      print(paste0(Sys.time(), " plotting Manhattan/Volcano plots."))
      df = getResultDataSingleTrait(globalVariables, sessionVariables, 0.05)
      sessionVariables$resultDataSingleTrait = df
      sessionVariables$traitDF <- traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$genderAttribut) #if error occurs here, then the wrong Traits file is referenced
      plotManhattanVolcano_SERVER("M", globalVariables, sessionVariables)
    }
  })

  shiny::observeEvent(input$btnSelectProbe, ignoreInit = TRUE, {
    if (shiny::isTruthy(sessionVariables$probe$probe)) {
      shiny::updateTextInput(session, "txtSelectedProbes", value = sessionVariables$probe$probe)
<<<<<<< HEAD

=======
      # print(paste0(Sys.time(), " plotting trait DNAm."))
      # plotTraitDNAm_SERVER("TraitDNAm", globalVariables, sessionVariables)
      # print(paste0(Sys.time(), " plotting DNAm profile."))
      # plotDNAmProfile_SERVER("DNAmProfile", globalVariables, sessionVariables)
      # print(paste0(Sys.time(), " plotting correlating probes."))
      # DTCorrelatingProbes_SERVER("CorrProbes", globalVariables, sessionVariables)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
    }
  })
  shiny::observeEvent(input$txtSelectedProbes, ignoreInit = FALSE, {
    if (shiny::isTruthy(input$txtSelectedProbes)) {
      sessionVariables$probe$probe = input$txtSelectedProbes
<<<<<<< HEAD
=======
#      output$txtSelectedProbesOld <- shiny::renderText(sessionVariables$probe$probe)
#      shiny::updateTextInput(session, "txtSelectedProbes", value = sessionVariables$probe$probe)
>>>>>>> 293e544a2d61b0de48d77a832b4da81bba457b80
      print(paste0(Sys.time(), " plotting trait DNAm."))
      plotTraitDNAm_SERVER("TraitDNAm", globalVariables, sessionVariables)
      print(paste0(Sys.time(), " plotting DNAm profile."))
      plotDNAmProfile_SERVER("DNAmProfile", globalVariables, sessionVariables)
      print(paste0(Sys.time(), " plotting correlating probes."))
      DTCorrelatingProbes_SERVER("CorrProbes", globalVariables, sessionVariables)
    }
  })
}
