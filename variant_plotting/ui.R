library(shiny)

# Define UI
shinyUI(fluidPage(

  # Application title
  titlePanel("Plot Variants Across subRVIS Regions"),

  # Sidebar
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
  p("This webtool can be used to plot the domain subRVIS scores of a given gene. 
     To the right of the plot is a legend denoting each aligned conserved domain, with the CDD
     ID in parentheses next to the domain name."),
  p("Underneath the plot are two subRVIS measures, per gene, intended to aid in the interpretation
     of the results. The first, labelled \"subRVIS P-value\", is the gene's FDR corrected p-value from the subRVIS paper, 
     denoting whether the gene's domain subRVIS scores corresponded to the distribution of 
     previously reported pathogenic variants within the gene. Please note that an insignificant p-value
     may reflect a lack of power, and should not be taken to indicate homogeneity of intolerance 
     throughout the gene."),
  p("The second measure, labelled \"subRVIS SDP\", is a score intended to reflect the degree to which genes 
     vary in intolerance
     across their regions. This score is the standard deviation of the gene's domain subRVIS scores, 
     converted
     to percentiles across all genes. A higher SDP indicates a higher degree of variation in the 
     intolerance scores across the gene's regions. Though these scores can be useful in predicting 
     whether we expect pathogenic mutations to cluster in specific sub regions within a given gene, 
     currently there is not a relationship between these scores and whether known pathogenic 
     mutations actually do so."),
  h3("About"),
  p("SubRVIS is a gene sub region based score intended to help in the intepretation of human
    sequence data. The sub regions are defined as the regions within the gene that aligned 
    to the Conserved Domain Database (http://www.ncbi.nlm.nih.gov/cdd/) and the unaligned regions 
    between each conserved domain alignment."),
  p("In its current form, subRVIS is based upon allele frequency as represented 
    in whole exome sequence data from the NHLBI-ESP6500 data set. The score is designed to 
    rank the genic sub regions based on their intolerance to functional variation. A lower score
    indicates a more intolerant sub region, while a higher score indicates a more tolerant sub region."),
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

