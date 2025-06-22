library(shiny)

source("shiny/ui-resultsTab.R")

function(input, output, session) {
  
  insertTab("mainPage", resultsTab, "Home", position = "after")
  hideTab("mainPage", "Results")
  
  source("shiny/server-plots.R", local = TRUE)
  
  m <- connect_mongo("TCGA_expr")
  genes <- m$find()$gene

  geneButton <- observeEvent(input$geneSearchButton, {
    print("Event Activated")
    
    plotMultiple()
    showTab("mainPage", "Results", select = TRUE)
  })
  
  updateSelectizeInput(session, "inputGene", choices = genes,
                       selected = character(0), server = TRUE,
                       options = list(maxOptions = 20, placeholder = "Enter gene..."))
}