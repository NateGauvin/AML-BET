bootstrapPage(
  
  selectInput(inputId = "geneid",
              label = "Which gene to graph?",
              choices = c(rownames(beat$X)[1:20]),
              selected = "TSPAN6"),
  
  plotOutput(outputId = "main_plot", height = "600px"),
  
)
