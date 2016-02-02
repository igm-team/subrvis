library(shiny)

# Define UI
shinyUI(fluidPage(
  
  # Google Analytics
  tags$head(includeScript("google-analytics.js")),

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
      radioButtons("rgn", NULL,
                   choices=c("Domains", "Exons"),
                   selected = NULL, 
                   inline = T, 
                   width = NULL),
      radioButtons("score.type", NULL,
                   choices=c("Raw", "Percentiles"),
                   selected = NULL, 
                   inline = T, 
                   width = NULL),
      actionButton("plot", label = "Plot")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel('Plot',
                 plotOutput("genePlot")),
        tabPanel('Table',
                 dataTableOutput("var.table"))
      )
    )),
  p("This online tool can be used to plot the domain or exon subRVIS scores of a given gene."),
  p("To the right of the plot is a legend. For domains, the legend denotes each aligned 
     conserved domain, with the CDD ID in parentheses next to the domain name. For exons,
     the legend denotes the exon number. Note that genes that are on the minus strand
     are plotted in reverse, i.e. regardless of its strand the first base of the gene
     is at position zero at the leftmost point of the plot."),
  p("There is the option to plot either the raw subRVIS scores or the subRVIS percentiles. The raw subRVIS
     scores are the scores obtained from directly applying the subRVIS methodology (see ", a("Gussow et al", href="http://www.genomebiology.com/2016/17/1/9", target="_blank"), ").
     The subRVIS percentiles are a scaled version of the score, where the raw scores were converted to 
     percentiles across the genome. Though the raw scores offer superior resolution, the percentiles
     are easier to interpret as they are put into the context of the rest of the sub regions
     in the genome. Lower subRVIS scores and percentiles correspond to more intolerant regions. 
     We have found that 35% and lower is a reasonable threshold for the subRVIS percentiles - in other words, we are more inclined to believe
     that a mutation is likely to cause disease when it falls in a sub region in the lower 35th percentile of
     subRVIS scores. However, by no means does that mean that mutations falling in sub regions below this percentile
	 are necessarily pathogenic, nor does it mean that mutations falling in sub regions above this percentile
	 are necessarily benign."),
  p("Underneath the plot are two subRVIS related measures, per gene, intended to aid in the 
     interpretation of the results. The first, labelled 
    \"subRVIS P-value for known pathogenic mutations\", 
    assesses the relationship between subRVIS scores and the distribution of known 
    pathogenic mutations (see ", a("Gussow et al", href="http://www.genomebiology.com/2016/17/1/9", target="_blank"), "). The subRVIS p-value for known pathogenic 
    mutations obviously is only reported for known disease causing genes."),
  p("The second measure, labelled \"subRVIS SDP\", is intended as a general score 
     to help assess whether the location of a mutation in a gene is an important factor 
     in determining whether the mutation is pathogenic. Unlike the p-value for 
     pathogenic mutations, this score can be constructed for both genes already 
     implicated in diseases and for genes not currently implicated. 
     This score is the standard deviation of the gene's domain or exon raw subRVIS scores, converted to percentiles.
     A higher SDP indicates a higher degree of variation in the intolerance scores across the 
     gene's regions. Though these scores can be useful in predicting whether we expect 
     pathogenic mutations to cluster in specific sub regions within a given gene, 
     currently there is not a relationship between these scores and whether known 
     pathogenic mutations actually do so."),
  h3("About"),
  p("SubRVIS is a gene sub region based score intended to help in the intepretation of human
    sequence data. The sub regions are defined as the regions within the gene that aligned 
    to the Conserved Domain Database (http://www.ncbi.nlm.nih.gov/cdd/) and the unaligned regions 
    between each conserved domain alignment."),
  p("In its current form, subRVIS is based upon allele frequency as represented 
    in whole exome sequence data from the NHLBI-ESP6500 data set. The score is designed to 
    rank the genic sub regions based on their intolerance to functional variation. A lower score
    indicates a more intolerant sub region, while a higher score indicates a more tolerant sub region."),
  p("The paper can be found ", a("here.", href="http://www.genomebiology.com/2016/17/1/9", target="_blank")),
  h3("Terms of use"),
  p("The content of the subRVIS site is intended strictly for educational 
     and research purposes. The data derived from this website may not be used 
     for any commercial purpose. The data from this website may not be 
     replicated on any other website without written consent."),
  h3("Citation"),
  p("Gussow AB, Petrovski S, Wang Q, Allen AS, Goldstein DB. 
     The intolerance to functional genetic variation of protein domains predicts the localization of pathogenic mutations within genes. 
     Genome Biology 2016, 17:9. doi:10.1186/s13059-016-0869-4."),
  h3("Contact"),
  p("For bug reports, please contact Ayal Gussow (abg2164@columbia.edu).")
))

