ui <- shinyUI(
  fluidPage(
    # Some custom CSS for a smaller font for preformatted text
    htmltools::tags$head(
      htmltools::tags$style(HTML("      pre, table.table {        font-size: smaller;      }    "))
    ),

    inputTrait_UI("trait"),
    shiny::actionButton("btnSelectTrait", label = "Select Trait"),
    plotManhattanVolcano_UI("M"),
#    shiny::verbatimTextOutput("txtSelectedProbesOld", placeholder = TRUE),
    shiny::textInput("txtSelectedProbes", "probeID", "", placeholder = TRUE),
    shiny::actionButton("btnSelectProbe", label = "Select Probe"),
    shiny::tabsetPanel(
      shiny::tabPanel("Trait vs DNAm",
                      shiny::fluidRow(
          plotTraitDNAm_UI("TraitDNAm")
        )
     ),
     shiny::tabPanel("Range DNAm",
                     shiny::fluidRow(
          plotDNAmProfile_UI("DNAmProfile")
        )
     ),
     shiny::tabPanel("Correlating Probes",
                     shiny::fluidRow(
          DTCorrelatingProbes_UI("CorrProbes")
        )
      )
    )
#,
#  shiny::verbatimTextOutput("txtdebug", placeholder = TRUE),
  )
)
