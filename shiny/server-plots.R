library(ggplot2)

source("functions/queryMongoGeneratePlots.R")

plotMultiple <- function() {
  plots <- query_mongo(input$inputGene, "risk", datasets = "all", input$inputTest)
  
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
}