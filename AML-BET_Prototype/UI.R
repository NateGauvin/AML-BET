bootstrapPage(
  
  h1("Test application"),
  
  selectInput(inputId = "geneid",
              label = "Which gene to graph?",
              #choices = c(rownames(beat$X)[1:20]),
              choices = c("TSPAN6", "CD24"),
              selected = "TSPAN6"),
  
  plotOutput(outputId = "main_plot", height = "600px"),
  
)
