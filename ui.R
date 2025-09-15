library(shiny)
library(bslib)

source("./shiny/ui-home.R")
source("./shiny/ui-about.R")

navbarPage(
  
  title = HTML("AML-BET"),
  windowTitle = "Acute Myeloid Leukemia Biomarker Evaluation Tool",
  id = "mainPage",
  
  tabPanel("Home",
           div(id = 'homepage',
               homePage
               )
           ),
  tabPanel("About",
           div(id = "aboutPage",
               aboutPage
               )
           )
  
)

