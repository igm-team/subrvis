library(shiny)
source("src/geneplots.R")

subrgn.data <- ReadSubrgnData("data/domains_table.txt")
subrgn.bed <- ReadSubrgnBed("data/domains.bed")
domain.map <- ReadDomainMap("data/cdd_table.txt")
strand.map <- ReadStrandData("data/strand_data.txt")
gene.pvals <- ReadGenePvals("data/domains_gene_pvals.txt")
gene.sds <- ReadGeneSds("data/domains_gene_sds.txt")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot

  GetGene <- eventReactive(input$plot, {
    input$gene
  })
    
  GetVariants <- eventReactive(input$plot, {
    
    if (is.null(input$variants) || input$variants == "") {
      return(NULL)
    }
    
    
    variants <- tryCatch({
      read.table(textConnection(input$variants), header=F)
      }, error=function(e) {
        stop("Unrecognized format for variants.")
      })
    
    if (ncol(variants) != 2) {
      stop("Variants should have two columns: chromosome and position.")
    }
    
    
    colnames(variants) <- c("chrom", "pos")
    if (!is.numeric(variants$pos)){
      stop("Variant position should be numeric.")
    }
    
    return(variants)
  })
  
  
  output$genePlot <- renderPlot({
    gene <- GetGene()
    if (!(toupper(gene) %in% toupper(subrgn.data$gene))) {
      err.msg <- paste("Gene", gene, "cannot be found.")
      stop(err.msg)
    }
    
    strand <- strand.map[toupper(rownames(strand.map)) == toupper(gene), ]
    if((length(strand) != 1) || (!(strand %in% c("+", "-")))) {
      err.msg <- paste("Gene", gene, "has no strand data.")
      stop(err.msg)      
    }
    
    
    variants <- GetVariants()
    
    gene.subrgns <- GetGeneSubrgn(subrgn.data, subrgn.bed, domain.map, gene)
    gene.name <- unique(gene.subrgns$gene)
    stopifnot(length(gene.name)==1)
    gene.plot <- PlotGeneSkel(gene.subrgns, gene.name)
        
    gene.plot <- PlotConnectors(gene.plot, gene.subrgns, strand)
    if (strand == "-") {
      gene.plot <- gene.plot + scale_x_reverse(label=function(x) -(x-max(gene.subrgns$end)))
    }
    
    if (!is.null(variants)) {
      gene.plot <- AddVariants(gene.plot, gene.subrgns, subrgn.bed, variants, strand)
    }

    p.val <- gene.pvals[toupper(gene.name), ]
    gene.sd <- gene.sds[toupper(gene.name), "percentile"]
    plot(gene.plot + labs(x=paste0("Gene Position\nsubRVIS P-value: ", 
                                   round(p.val, digits=3))))
  })
})