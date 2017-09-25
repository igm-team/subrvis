library(shiny)
source("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/src/geneplots.R")

## functions
Scores2Percentiles <- function(scores) {
        stopifnot(is.numeric(scores))
        return((rank(scores) / length(scores)) * 100)
}

# for domains
dmn.data <- ReadSubrgnData("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/domains_table.txt")
dmn.data$perc <- Scores2Percentiles(dmn.data$subRVIS)
dmn.bed <- ReadSubrgnBed("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/domains.bed")
dmn.map <- ReadDomainMap("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/cdd_table.txt")
dmn.pvals <- ReadGenePvals("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/domains_gene_pvals.txt")
dmn.sds <- ReadGeneSds("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/domains_gene_sds.txt")

# for exons
exn.data <- ReadSubrgnData("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/exons_table.txt")
exn.data$perc <- Scores2Percentiles(exn.data$subRVIS)
exn.bed <- ReadSubrgnBed("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/exons.bed")
exn.map <- NULL
exn.pvals <- ReadGenePvals("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/exons_gene_pvals.txt")
exn.sds <- ReadGeneSds("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/exons_gene_sds.txt")

# for all
strand.map <- ReadStrandData("C:/Users/m119979/Desktop/subrvis-master/variant_plotting/data/strand_data.txt")

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
        
        GetPlotData <- eventReactive(input$plot, {
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
                
                gene.subrgns <- GetGeneSubrgn(subrgn.data, subrgn.bed, subrgn.map, gene, rgn, score.type)
                gene.name <- unique(gene.subrgns$gene)
                stopifnot(length(gene.name)==1)
                
                p.val <- gene.pvals[toupper(gene.name), ]
                gene.sd <- gene.sds[toupper(gene.name), "percentile"]
                gene.bed <- subrgn.bed[toupper(subrgn.bed$gene) == toupper(gene), ]
                
                return(list(gene.subrgns = gene.subrgns, 
                            subrgn.bed = subrgn.bed,
                            gene.bed = gene.bed,
                            gene.name = gene.name, 
                            p.val = p.val, 
                            gene.sd = gene.sd, 
                            variants = variants, 
                            rgn = rgn, 
                            score.type = score.type, 
                            strand = strand))
        })
        
        
        output$genePlot <- renderPlot({
                
                plot.data <- GetPlotData()
                
                gene.plot <- PlotGeneSkel(plot.data$gene.subrgns, plot.data$gene.name, plot.data$rgn, plot.data$score.type)
                
                gene.plot <- PlotConnectors(gene.plot, plot.data$gene.subrgns, plot.data$strand)
                if (plot.data$strand == "-") {
                        gene.plot <- gene.plot + scale_x_reverse(label=function(x) -(x-max(plot.data$gene.subrgns$end)))
                }
                
                if (!is.null(plot.data$variants)) {
                        gene.plot <- AddVariants(gene.plot, plot.data$gene.subrgns, plot.data$subrgn.bed, plot.data$variants, plot.data$strand)
                }
                
                ylab.txt <- ifelse(plot.data$score.type=="Percentiles", "subRVIS Percentile", "subRVIS Raw Score")
                plot(gene.plot + labs(x=paste0("Gene Position\n",
                                               "subRVIS P-value for known pathogenic mutations: ", 
                                               round(plot.data$p.val, digits=3),
                                               "\nsubRVIS SDP: ",
                                               round(plot.data$gene.sd, digits=1),
                                               "%"),
                                      y=ylab.txt))
        })
        
        output$var.table = renderDataTable({
                plot.data <- GetPlotData()
                
                plot.table <- plot.data$gene.subrgns
                plot.table$var.count <- rep(0, nrow(plot.table))
                
                gene.chrom <- unique(plot.data$gene.bed$chrom)
                stopifnot(length(gene.chrom) == 1)
                
                variants <- plot.data$variants
                variants <- variants[variants$chrom == gene.chrom, ]
                
                positions <- variants$pos - 1 # convert to 0-based BED coords
                
                for (pos in positions) {
                        subrgn.name <- plot.data$gene.bed[pos >= plot.data$gene.bed$start & pos < plot.data$gene.bed$end, "subrgn"]
                        stopifnot(length(subrgn.name) <= 1)
                        
                        if(length(subrgn.name) == 1) {
                                plot.table[subrgn.name, "var.count"] <- plot.table[subrgn.name, "var.count"] + 1
                        }
                        
                }
                
                plot.table$name <- rownames(plot.table)
                
                cols.to.plot <- c("name", "start", "end", "size", "coverage_percent", "subRVIS", plot.data$rgn, "var.count")
                
                plot.table <- plot.table[, cols.to.plot]
                plot.table$start <- plot.table$start + 1
                
                if((plot.data$strand == "-") && (nrow(plot.table) != 1)) {
                        # Need to reverse coordinates
                        plot.table <- plot.table[order(-plot.table$start), ]
                        plot.table$start[1] <- 1
                        plot.table$end[1] <- plot.table$size[1]
                        for (i in 2:nrow(plot.table)) {
                                plot.table[i, "start"] <- plot.table$end[i-1] + 1
                                plot.table$end[i] <- plot.table$size[i] + plot.table$start[i] - 1
                        }
                }
                
                colnames(plot.table) <- c("name", "start", "end", "size",
                                          "percent_covered", "score",
                                          "sub-region", "variant_count")
                
                return(plot.table)
        })
        
        output$downloadData <- downloadHandler(
                filename = function() {
                        paste('subRVIS table','.csv', sep='')
                },
                content = function(file){
                        plot.data <- GetPlotData()
                        
                        plot.table <- plot.data$gene.subrgns
                        plot.table$var.count <- rep(0, nrow(plot.table))
                        
                        gene.chrom <- unique(plot.data$gene.bed$chrom)
                        stopifnot(length(gene.chrom) == 1)
                        
                        variants <- plot.data$variants
                        variants <- variants[variants$chrom == gene.chrom, ]
                        
                        positions <- variants$pos - 1 # convert to 0-based BED coords
                        
                        for (pos in positions) {
                                subrgn.name <- plot.data$gene.bed[pos >= plot.data$gene.bed$start & pos < plot.data$gene.bed$end, "subrgn"]
                                stopifnot(length(subrgn.name) <= 1)
                                
                                if(length(subrgn.name) == 1) {
                                        plot.table[subrgn.name, "var.count"] <- plot.table[subrgn.name, "var.count"] + 1
                                }
                                
                        }
                        
                        plot.table$name <- rownames(plot.table)
                        
                        cols.to.plot <- c("name", "start", "end", "size", "coverage_percent", "subRVIS", plot.data$rgn, "var.count")
                        
                        plot.table <- plot.table[, cols.to.plot]
                        plot.table$start <- plot.table$start + 1
                        
                        if((plot.data$strand == "-") && (nrow(plot.table) != 1)) {
                                # Need to reverse coordinates
                                plot.table <- plot.table[order(-plot.table$start), ]
                                plot.table$start[1] <- 1
                                plot.table$end[1] <- plot.table$size[1]
                                for (i in 2:nrow(plot.table)) {
                                        plot.table[i, "start"] <- plot.table$end[i-1] + 1
                                        plot.table$end[i] <- plot.table$size[i] + plot.table$start[i] - 1
                                }
                        }
                        
                        colnames(plot.table) <- c("name", "start", "end", "size",
                                                  "percent_covered", "score",
                                                  "sub-region", "variant_count")
                        
                        table_data <- plot.table
                        
                        write.csv(table_data, file)
                })  
        
        
        output$downloadPlot <- downloadHandler(
                filename = function() {
                        paste("subRVIS plot", ".svg", sep='')
                },
                content = function(file) {
                        plot.data <- GetPlotData()
                        
                        gene.plot <- PlotGeneSkel(plot.data$gene.subrgns, plot.data$gene.name, plot.data$rgn, plot.data$score.type)
                        
                        gene.plot <- PlotConnectors(gene.plot, plot.data$gene.subrgns, plot.data$strand)
                        if (plot.data$strand == "-") {
                                gene.plot <- gene.plot + scale_x_reverse(label=function(x) -(x-max(plot.data$gene.subrgns$end)))
                        }
                        
                        if (!is.null(plot.data$variants)) {
                                gene.plot <- AddVariants(gene.plot, plot.data$gene.subrgns, plot.data$subrgn.bed, plot.data$variants, plot.data$strand)
                        }
                        
                        ylab.txt <- ifelse(plot.data$score.type=="Percentiles", "subRVIS Percentile", "subRVIS Raw Score")
                        
                        plot(gene.plot + labs(x=paste0("Gene Position\n",
                                                               "subRVIS P-value for known pathogenic mutations: ", 
                                                               round(plot.data$p.val, digits=3),
                                                               "\nsubRVIS SDP: ",
                                                               round(plot.data$gene.sd, digits=1),
                                                               "%"),
                                                      y=ylab.txt))
                        
                        ggsave(file, plot = last_plot (), dpi = 300, width = 12, height = 6, units = "in")
                })
})