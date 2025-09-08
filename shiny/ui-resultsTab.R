riskTab <- tabPanel("Risk",
  uiOutput("riskPlots")
)

survivalTab <- tabPanel("Survival",
  uiOutput("survPlots")
)

forestTab <- tabPanel("Forest",
                   
)

resultsTab <- tabPanel("Results",
  hr(),
  fluidRow(column(12,
    tabsetPanel(id = "resultsPage",
      riskTab,
      survivalTab,
      forestTab
      )
    )
    
  )
)
