library(ggplot2)

ReadSubrgnData <- function(data.path) {
  # read sub region data
  return(read.table(data.path, header=T, row.names=1, stringsAsFactors=F, comment.char=""))
}

ReadStrandData <- function(strand.path) {
  # read the strand data per gene.
  return(read.table(strand.path, header=T, row.names=1, comment.char=""))
}

ReadSubrgnBed <- function(bed.path) {
  # read sub region BED file
  ret <- read.table(bed.path, header=F, stringsAsFactors=F, comment.char="")
  colnames(ret) <- c("chrom", "start", "end", "subrgn")
  ret$gene <- sapply(strsplit(ret$subrgn, ":"), function (x) x[[1]])
  
  return(ret)
}

ReadGenePvals <- function(pvals.path) {
  # read each gene's pval
  ret <- read.table(pvals.path, header=T, row.names=1, stringsAsFactors=F, comment.char="")
  rownames(ret) <- toupper(rownames(ret))
  return(ret)
}

ReadGeneSds <- function(sds.path) {
  # read each gene's SD
  ret <- read.table(sds.path, header=T, row.names=1, stringsAsFactors=F, comment.char="")
  rownames(ret) <- toupper(rownames(ret))
  return(ret)
}

ReadDomainMap <- function(domain.map.path) {
  # read domain map file
  ret <- read.table(domain.map.path, 
                    header=F, 
                    stringsAsFactors=F, 
                    sep="\t", 
                    comment.char="",
                    quote="",
                    row.names=1)
  ret["-", ] <- c("-", "No domain aligned")
  colnames(ret) <- c("origid", "name")

  return(ret)
}

GetGeneSubrgn <- function(subrgn.data, subrgn.bed, domains.map, gene) {
  # extract sub region gene data
  gene.data <- subrgn.data[toupper(subrgn.data$gene)==toupper(gene), ]
  gene.bed <- subrgn.bed[toupper(subrgn.bed$gene)==toupper(gene), ]
  gene.size <- sum(gene.bed$end - gene.bed$start)
  
  for (subrgn in rownames(gene.data)) {
    subrgn.gene.bed <- gene.bed[gene.bed$subrgn==subrgn, ]
    
    gene.data[subrgn, "size"] <- sum(subrgn.gene.bed$end - subrgn.gene.bed$start)
    
    # add 1 to start as BED is 0-based, non-inclusive
    gene.data[subrgn, "start"] <- min(subrgn.gene.bed$start) + 1
  }
  
  gene.data <- gene.data[order(gene.data$start), ]
  
  if (length(gene.data$start) != 1) {
    gene.data$start <- c(0, gene.data$size[1:(length(gene.data$size)-1)])
    
    for (i in 2:length(gene.data$start)) {
      gene.data$start[i] <- gene.data$start[i] + gene.data$start[i-1]
    }
  } else {
    gene.data$start <- 0
  }
  
  gene.data$end <- gene.data$start + gene.data$size
  
  gene.data$Domain <- sapply(strsplit(rownames(gene.data), ":"), function (x) x[[2]])
  gene.data$Domain <- paste0(domains.map[gene.data$Domain, "name"], 
                             " (", gene.data$Domain, ")")
  
  return(gene.data)
}

PlotGeneSkel <- function(gene.subrgns, title) {
  # plot gene skeleton without variants
  gene.skel <- ggplot(data=gene.subrgns, aes(x=0, ymin=subRVIS, ymax=subRVIS, color=Domain)) + 
    geom_segment(size=5, aes(y=min(subRVIS) - 1, 
                                              x=start, 
                                              xend=end, 
                                              yend=min(subRVIS) - 1))
  
  # plot subRVIS lines
  gene.skel <- gene.skel + 
    labs(x="Gene Position", y="subRVIS Score", title=title) +
    geom_segment(aes(y=subRVIS, 
                     x=start, 
                     xend=end, 
                     yend=subRVIS), color="purple")
  
  return(gene.skel)
  
}

PlotConnectors <- function(gene.skel, gene.subrgns, strand) {
  # plot connectors for subRVIS lines
  if(nrow(gene.subrgns) == 1) {
    return(gene.skel)
  }
  
  # more than one region
  for(i in 1:(nrow(gene.subrgns) - 1)) {
    gene.skel <- gene.skel + geom_segment(y=gene.subrgns$subRVIS[i],
                                          x=ifelse(strand=="+",
                                                   gene.subrgns$end[i], 
                                                   -gene.subrgns$end[i]),
                                          xend=ifelse(strand=="+",
                                                      gene.subrgns$end[i], 
                                                      -gene.subrgns$end[i]),
                                          yend=gene.subrgns$subRVIS[i+1],
                                          color="purple")
  } 
  
  return(gene.skel)
}


AddVariants <- function(gene.plot, gene.subrgns, subrgn.bed, variants, strand) {
  # add variants to gene plot
  gene.name <- unique(gene.subrgns$gene)
  stopifnot(length(gene.name)==1)
  gene.bed <- subrgn.bed[toupper(subrgn.bed$gene)==toupper(gene.name), ]
  gene.chrom <- unique(gene.bed$chrom)
  stopifnot(length(gene.chrom)==1)
  
  for (i in 1:nrow(variants)) {
    variant <- variants[i, ]
    if (variant$chrom == gene.chrom) { 
      gene.plot <- PlotVariant(variant, gene.subrgns, gene.bed, strand, gene.plot)
    }
  }
  
  return(gene.plot)
}

PlotVariant <- function(variant, gene.subrgns, gene.bed, strand, gene.plot) {
  # plot the variant on the exonic coding sequence based on genomic coordinates.
  gene.chrom <- unique(gene.bed$chrom)
  pos <- variant$pos - 1 # convert to 0-based BED coords
  
  if (!(gene.chrom==variant$chrom)) {
    return(gene.plot)
  }
  
  if (!(pos >= min(gene.bed$start) && pos < max(gene.bed$end))) {
    return(gene.plot)
  }
  
  # variant is within gene boundaries
  gene.subrgn.names <- rownames(gene.subrgns)
  
  plot.pos <- 0
  for (gene.subrgn.name in gene.subrgn.names) {
    gene.subrgn.bed <- gene.bed[gene.bed$subrgn == gene.subrgn.name, ]
    
    for (i in 1:nrow(gene.subrgn.bed)) {
      exn <- gene.subrgn.bed[i, ]
      
      if (pos >= exn$start && pos < exn$end) { #in range
        plot.pos <- plot.pos + (pos - exn$start)
        variant.height <- (((max(gene.subrgns$subRVIS) - min(gene.subrgns$subRVIS)) / 2) + 
                           min(gene.subrgns$subRVIS))
        
        if (strand == "+") {
          return(gene.plot + geom_point(x=plot.pos,
                                        y=variant.height, 
                                        color="dodgerblue3", 
                                        size=4))
        } else {
          return(gene.plot + geom_point(x=-plot.pos,
                                        y=variant.height, 
                                        color="dodgerblue3", 
                                        size=4))
          }
        
      } else {
        plot.pos <- plot.pos + (exn$end - exn$start)
      }
    }
  }
  
  return(gene.plot)
}

