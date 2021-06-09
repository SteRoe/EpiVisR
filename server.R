# reactive Create a variable that can be changed over time by user inputs, evaluates "lazy" meaning only when called.
# 
# observe Continually monitor reactive events and variables, whenever ANY reactive variable is changed in the environment (the observed environment), the code is evaluated. Can change values of previously defined reactive variables, cannot create/return variables.
# 
# observeEvent (Domino Effect) Continually monitor ONE defined reactive variable/event (the trigger) and run the code when the the trigger is activated by change/input of that trigger. Can change values of previously defined reactive variables, cannot create/return variables.
# 
# eventReactive Create a variable, with a defined trigger similar to observeEvent. Use this when you want a reactive variable that evaluates due to a trigger instead of when it is called.
# 
# I hope this helps, and if I am mistaken in my understanding or there could be more clarification, feel free to edit this answer.

#library(crosstalk)

#create dummy object for shared data frame
# fc <- c("a","b")
# sc <- c("c", "d")
# shared_df_correlatingProbes <- data.frame(fc, sc)
# assign("shared_df_correlatingProbes",shared_df_correlatingProbes,envir=globalenv())

#source(paste0(appDir,"inputFolder.R"))
source(paste0(appDir,"inputTrait.R"))
source(paste0(appDir,"plotManhattanVolcano.R"))
source(paste0(appDir,"plotTraitDNAm.R"))
source(paste0(appDir,"plotDNAmProfile.R"))
source(paste0(appDir,"DTCorrelatingProbes.R"))
source(paste0(appDir,"plotHM.R"))

server <- function(input, output, session) {
###because <InteractiveComplexHeatmap> does not support modules, we have to do this here outside the heatmap module plotHM.R
###nevertheless, we refer to the functions defined in plotHM.R

  checkConfigVariables()

  sessionVariables <- reactiveValues(folder = "", dataFileName = "", trait = list(), probe = list(), df_plot = data.frame(), resultDataSingleTrait = data.frame(), traitDF = data.frame(), traitsDFLong = data.frame())
  inputTrait_SERVER("trait", sessionVariables)

  observeEvent(input$btnSelectTrait, { # ignoreInit = TRUE,
    id <- showNotification(paste0("Selected trait: ",sessionVariables$trait$trait), duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    #load DF from selected trait
    if (isTruthy(sessionVariables$trait$trait)) {
      df = getResultDataSingleTrait(sessionVariables,100)
      sessionVariables$resultDataSingleTrait = df
#browser()
      sessionVariables$traitDF <- traitDF(sessionVariables)
      
      plotManhattanVolcano_SERVER("M", sessionVariables)
    }
#    return(sessionVariables$probe)
  })
  
  observeEvent(input$btnSelectProbe, ignoreInit = TRUE, {
    if (isTruthy(sessionVariables$probe$probe)) {
      output$txtSelectedProbes <- shiny::renderText(sessionVariables$probe$probe)
      plotTraitDNAm_SERVER("TraitDNAm", sessionVariables)
      plotDNAmProfile_SERVER("DNAmProfile", sessionVariables)
      DTCorrelatingProbes_SERVER("CorrProbes", sessionVariables)
    }
  })
}
