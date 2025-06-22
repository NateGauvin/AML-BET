riskTab <- tabPanel("Risk",
  uiOutput("riskPlots")
)

survivalTab <- tabPanel("Survival",
  
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
