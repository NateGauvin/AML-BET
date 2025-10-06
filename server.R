# TODO
# Upload all datasets

library(shiny)

selected_gene <<- "No Gene Selected"

source("./shiny/ui-resultsTab.R")

function(input, output, session) {
  
  insertTab("mainPage", resultsTab, "Home", position = "after")
  hideTab("mainPage", "Results")
  
  source("./shiny/server-plots.R", local = TRUE)
  
  m <- connect_mongo("genes")
  genes <- m$find()$gene
  genes <- genes[order(genes)][15:length(genes)]

  geneButton <- observeEvent(input$geneSearchButton, {
    print("Calculating...")
    selected_gene <<- input$inputGene
    output$geneHeader <- renderUI({
      h1(selected_gene)
    })
    print(selected_gene)
    plotMultiple()
    showTab("mainPage", "Results", select = TRUE)
  })
  
  updateSelectizeInput(session, "inputGene", choices = genes,
                       selected = character(0), server = TRUE,
                       options = list(maxOptions = 20, placeholder = "Enter gene..."))
}