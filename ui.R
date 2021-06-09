library(shiny)
library(shinyBS)

#source(paste0(appDir,"inputFolder.R"))
source(paste0(appDir,"inputTrait.R"))
source(paste0(appDir,"plotManhattanVolcano.R"))
source(paste0(appDir,"plotTraitDNAm.R"))
source(paste0(appDir,"plotDNAmProfile.R"))
source(paste0(appDir,"DTCorrelatingProbes.R"))
source(paste0(appDir,"plotHM.R"))

ui <- shinyUI(
  fluidPage(
    # Some custom CSS for a smaller font for preformatted text
    tags$head(
      tags$style(HTML("      pre, table.table {        font-size: smaller;      }    "))
    ),

    tabsetPanel(
      tabPanel("Trait",
#         inputFolder_UI("folder"),
         inputTrait_UI("trait"),
         shiny::actionButton("btnSelectTrait", label = "Select Trait"),
         plotManhattanVolcano_UI("M"),
         shiny::actionButton("btnSelectProbe", label = "Select Probe"),
         verbatimTextOutput("txtSelectedProbes", placeholder = TRUE),
         tabsetPanel(
           tabPanel("Trait vs DNAm",
              fluidRow(
                plotTraitDNAm_UI("TraitDNAm")
              )
           ),
           tabPanel("Range DNAm",
                    fluidRow(
                      plotDNAmProfile_UI("DNAmProfile")
                    )
           ),
           tabPanel("Correlating Probes",
                    fluidRow(
                      DTCorrelatingProbes_UI("CorrProbes")
                    )
           )
         )
      ),
      # tabPanel("HM",
      #   shiny::actionButton("btnplotHM", label = "Plot HM"),
      #   plotHM_UI("HeatMap"),
      #   fluidRow(
      #     column(width = 12,
      #            InteractiveComplexHeatmapOutput("heatmap_1",
      #                                            height1 = 800, height2 = 800),
      #     )
      #   ),
      # ),
      #    fluidRow(
      #      plotHM_UI("HeatMap")
      #    )
      # ),
      tabPanel("Settings",
         fluidRow(
           tags$html("Settings are now in <config.yml>"),
           # shiny::textInput("inpTraitsFile", "File with Traits", "FileName", placeholder = TRUE),
           # shiny::textInput("inpDNAmFile", "File with DNAm", "FileName", placeholder = TRUE),
           # shiny::textInput("inpAnnotation1File", "File with Annotation1", "FileName", placeholder = TRUE),
           # shiny::textInput("inpAnnotation2File", "File with Annotation2", "FileName", placeholder = TRUE),
           # shiny::textInput("inpWinsor", "Winsorize Traits %", "%", placeholder = TRUE)
         )
      )
    )
  ) 
)
