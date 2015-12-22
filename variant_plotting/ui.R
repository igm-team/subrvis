library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
#  navbarPage("subRVIS",
#             tabPanel("Plot",

  # Application title
  titlePanel("Plot Variants Across subRVIS Regions"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      textInput("gene", label = h3("Gene Name"), 
                value = "MAPT"),
      
      tags$style(type="text/css", "textarea {width:100%}"),
      tags$textarea(id = 'variants', 
                    placeholder = 'Variants (build 37), e.g. chr1 10003', 
                    rows = 8, ""),
      
      actionButton("plot", label = "Plot")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("genePlot")
    )),
#  ),
  h3("Terms of use"),
  p("The content of the subRVIS site is intended strictly for educational 
              and research purposes. The data derived from this website may not be used 
              for any commercial purpose. The data from this website may not be 
              replicated on any other website without written consent."),
  h3("Citation"),
  p("Gussow AB, Petrovski S, Wang Q, Allen AS, Goldstein DB. 
     The intolerance to functional genetic variation of protein domains predicts the localization of pathogenic mutations within genes. 
     Under review."),
  h3("Contact"),
  p("For bug reports, please contact Ayal Gussow (abg2164@columbia.edu).")
))
#  tabPanel("About", ""),
#  tabPanel("Terms", 
#           )
#  )
