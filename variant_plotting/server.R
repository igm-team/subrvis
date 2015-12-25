library(shiny)
source("src/geneplots.R")

## functions
Scores2Percentiles <- function(scores) {
  stopifnot(is.numeric(scores))
  return((rank(scores) / length(scores)) * 100)
}

# for domains
dmn.data <- ReadSubrgnData("data/domains_table.txt")
dmn.data$perc <- Scores2Percentiles(dmn.data$subRVIS)
dmn.bed <- ReadSubrgnBed("data/domains.bed")
dmn.map <- ReadDomainMap("data/cdd_table.txt")
dmn.pvals <- ReadGenePvals("data/domains_gene_pvals.txt")
dmn.sds <- ReadGeneSds("data/domains_gene_sds.txt")

# for exons
exn.data <- ReadSubrgnData("data/exons_table.txt")
exn.data$perc <- Scores2Percentiles(exn.data$subRVIS)
exn.bed <- ReadSubrgnBed("data/exons.bed")
exn.map <- NULL
exn.pvals <- ReadGenePvals("data/exons_gene_pvals.txt")
exn.sds <- ReadGeneSds("data/exons_gene_sds.txt")

# for all
strand.map <- ReadStrandData("data/strand_data.txt")


# Define server logic
shinyServer(function(input, output) {

  # Expression that generates the plot. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot

  GetGene <- eventReactive(input$plot, {
    input$gene
  })
  
  GetRgn <- eventReactive(input$plot, {
    input$rgn
  })
    
  GetScoreType <- eventReactive(input$plot, {
    input$score.type
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
    rgn <- GetRgn()
    score.type <- GetScoreType()
    
    if (rgn == "Domains") {
      subrgn.data <- dmn.data
      subrgn.bed <- dmn.bed
      subrgn.map <- dmn.map
      gene.pvals <- dmn.pvals
      gene.sds <- dmn.sds
    } else {
      subrgn.data <- exn.data
      subrgn.bed <- exn.bed
      subrgn.map <- exn.map
      gene.pvals <- exn.pvals
      gene.sds <- exn.sds
    }
    
    gene.subrgns <- GetGeneSubrgn(subrgn.data, subrgn.bed, subrgn.map, gene, rgn, score.type)
    gene.name <- unique(gene.subrgns$gene)
    stopifnot(length(gene.name)==1)
    gene.plot <- PlotGeneSkel(gene.subrgns, gene.name, rgn, score.type)
        
    gene.plot <- PlotConnectors(gene.plot, gene.subrgns, strand)
    if (strand == "-") {
      gene.plot <- gene.plot + scale_x_reverse(label=function(x) -(x-max(gene.subrgns$end)))
    }
    
    if (!is.null(variants)) {
      gene.plot <- AddVariants(gene.plot, gene.subrgns, subrgn.bed, variants, strand)
    }

    p.val <- gene.pvals[toupper(gene.name), ]
    gene.sd <- gene.sds[toupper(gene.name), "percentile"]
    
    ylab.txt <- ifelse(score.type=="Percentiles", "subRVIS Percentile", "subRVIS Raw Score")
    plot(gene.plot + labs(x=paste0("Gene Position\n",
                                   "subRVIS P-value for known pathogenic mutations: ", 
                                   round(p.val, digits=3),
                                   "\nsubRVIS SDP: ",
                                   round(gene.sd, digits=1),
                                   "%"),
                          y=ylab.txt))
  })
})