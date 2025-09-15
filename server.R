# TODO
# Upload all datasets
# Make collection that contains all available genes
# Add KM plotting capabilities (get HR using coxph, median is cutpoint)

library(shiny)

source("./shiny/ui-resultsTab.R")

function(input, output, session) {
  
  insertTab("mainPage", resultsTab, "Home", position = "after")
  hideTab("mainPage", "Results")
  
  source("./shiny/server-plots.R", local = TRUE)
  
  # Expand this to include genes in all datasets
  # Make a collection that contains all genes, pull from there
  m <- connect_mongo("genes")
  genes <- m$find()$gene
  genes <- genes[order(genes)][15:length(genes)]

  geneButton <- observeEvent(input$geneSearchButton, {
    print("Event Activated")
    
    plotMultiple()
    showTab("mainPage", "Results", select = TRUE)
  })
  
  updateSelectizeInput(session, "inputGene", choices = genes,
                       selected = character(0), server = TRUE,
                       options = list(maxOptions = 20, placeholder = "Enter gene..."))
}