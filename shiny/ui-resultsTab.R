riskTab <- tabPanel("Risk",
  uiOutput("riskPlots")
)

survivalTab <- tabPanel("Survival",
  uiOutput("survPlots")
)

forestTab <- tabPanel("Forest",
  uiOutput("forestPlots")
)

resultsTab <- tabPanel("Results",
  hr(),
  uiOutput("geneHeader"),
  fluidRow(column(12,
    tabsetPanel(id = "resultsPage",
      riskTab,
      survivalTab,
      forestTab
      )
    )
    
  )
)
