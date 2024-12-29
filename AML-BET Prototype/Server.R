function(input, output) {
  
  output$main_plot <- renderPlot({
    
    riskByExpression(beat, input$geneid)
  })
}
