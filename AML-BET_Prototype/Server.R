function(input, output) {
  
  output$main_plot <- renderPlot({
    
    #riskByExpression(beat, input$geneid)
    #ggplot(data = iris, aes(Sepal.Length, Sepal.Width)) + geom_point()
    plot(1:10)
    
  })
}
