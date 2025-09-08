library(ggplot2)

source("functions/queryMongoGeneratePlots.R")

# Make "query_mongo_km" to expand plot-generating to include kaplan-meier

plotMultiple <- function() {
  plots <- query_mongo(input$inputGene, "risk", datasets = "all", input$inputTest)
  plots_km <- query_mongo_km(input$inputGene, datasets = "all")
  
  
  output$riskPlots <- renderUI({
    lapply(1:length(plots), function(i) {
      if (is.na(plots[[i]][1])) (return())
      id <- paste0("risk_", i)
      col <- column(width = 3, plotOutput(id))
      output[[id]] <- renderPlot({
        plots[[i]]
      })
      col
      
    }
    )
  })
  
  output$survPlots <- renderUI({
    lapply(1:length(plots_km), function(i) {
      if (is.na(plots_km[[i]][1])) (return())
      id <- paste0("surv_", i)
      col <- column(width = 3, plotOutput(id))
      output[[id]] <- renderPlot({
        plots_km[[i]]
      })
      col
      
    }
    )
  })
}