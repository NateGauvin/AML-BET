selectingSingleGene <- list(
  
  h4("Select a Gene:"),
  
  fluidRow(
    column(3,
      selectizeInput("inputGene", label = "Begin typing or select gene",
                     choices = sort(c("A1BG")))
    ),
    column(3,
      selectInput("inputTest", label = "Select a testing format",
                     choices = c("t-test", "wilcox", "anova"))
      ),
    column(2,
      div(
        actionButton("geneSearchButton", "Evaluate Gene")
      )
    )
  )
  
)

homePage <- list(
  
  h1("Acute Myeloid Leukemia Biomarker Evaluation Tool (AML-BET)"),
  hr(color = "blue", style = "height: 2px"),
  
  selectingSingleGene
)