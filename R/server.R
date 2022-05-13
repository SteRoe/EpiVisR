server <- function(input, output, session) {
  globalVariables <- list()
  print(paste0(Sys.time(), " loading configuration."))
  packagePath <- find.package("EpiVisR", lib.loc=NULL, quiet = TRUE)
  configFileName <- paste0(packagePath,"/","config.yml")
  #globalVariables$config <- config::get(file = "config.yml")
  globalVariables$config <- config::get(file = configFileName)
  globalVariables$config <- addPackagePathToConfig(globalVariables$config, packagePath)
  sessionVariables <- shiny::reactiveValues(folder = "", trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), traitDF = data.frame())
  print(paste0(Sys.time(), " loading objects."))
  globalVariables = loadObjects(globalVariables)
  print(paste0(Sys.time(), " making session variables."))
  print(paste0(Sys.time(), " reading inputTrait_SERVER."))
  inputTrait_SERVER("trait", globalVariables, sessionVariables)

  shiny::observeEvent(input$btnSelectTrait, ignoreInit = TRUE, {
    id <- shiny::showNotification(paste0("Selected trait: ",sessionVariables$trait$trait), duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    #load DF from selected trait
    if (shiny::isTruthy(sessionVariables$trait$trait)) {
      print(paste0(Sys.time(), " plotting Manhattan/Volcano plots."))
      df = getResultDataSingleTrait(globalVariables, sessionVariables, 0.05)
      sessionVariables$resultDataSingleTrait = df
      sessionVariables$traitDF <- traitDF(sessionVariables, globalVariables$config$mergeAttribut, globalVariables$config$sexAttribut) #if error occurs here, then the wrong Traits file is referenced
      plotManhattanVolcano_SERVER("M", globalVariables, sessionVariables)
    }
  })

  shiny::observeEvent(input$btnSelectProbe, ignoreInit = TRUE, {
    id <- shiny::showNotification(paste0("Selected probe: ",sessionVariables$probe$probe), duration = NULL, closeButton = FALSE)
    on.exit(shiny::removeNotification(id), add = TRUE)
    if (shiny::isTruthy(sessionVariables$probe$probe)) {
#      shiny::updateTextInput(session, "txtSelectedProbe", value = sessionVariables$probe$probe)
#      sessionVariables$probe$probe = input$txtSelectedProbes
      print(paste0(Sys.time(), " plotting trait DNAm."))
      plotTraitDNAm_SERVER("TraitDNAm", globalVariables, sessionVariables)
      print(paste0(Sys.time(), " plotting DNAm profile."))
      plotDNAmProfile_SERVER("DNAmProfile", globalVariables, sessionVariables)
      print(paste0(Sys.time(), " plotting correlating probes."))
      DTCorrelatingProbes_SERVER("CorrProbes", globalVariables, sessionVariables)

    }
  })
#  shiny::observeEvent(input$txtSelectedProbe, ignoreInit = FALSE, {
#    if (shiny::isTruthy(input$txtSelectedProbes)) {
    #   sessionVariables$probe$probe = input$txtSelectedProbes
    #   print(paste0(Sys.time(), " plotting trait DNAm."))
    #   plotTraitDNAm_SERVER("TraitDNAm", globalVariables, sessionVariables)
    #   print(paste0(Sys.time(), " plotting DNAm profile."))
    #   plotDNAmProfile_SERVER("DNAmProfile", globalVariables, sessionVariables)
    #   print(paste0(Sys.time(), " plotting correlating probes."))
    #   DTCorrelatingProbes_SERVER("CorrProbes", globalVariables, sessionVariables)
#    }
#  })
}
