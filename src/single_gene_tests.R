#!/usr/bin/Rscript
###############################################################################
# Purpose:  Run single-gene tests across subRVIS genes.                       #
# Author:   Ayal Gussow                                                       #
# Input:    counts file, scores file.                                         #
# Output:   P-value per gene.                                                 #
###############################################################################

### Libraries -----------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(aod))

### Functions -----------------------------------------------------------------
GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("-s", "--scores"), action="store", 
                default="scores.txt", dest="scores.path", 
                help="File with scores, [default %default]"),

    make_option(c("-c", "--counts"), action="store", 
                default="counts.txt", dest="counts.path", 
                help="File with counts / mrt, [default %default]"),

    make_option(c("--perms"), action="store", 
                dest="nperms", type="integer",
                help="Number of permutations"),

     make_option(c("-m", "--min-regions"), action="store", 
                dest="min.regions", type="integer",
                help="Minimum regions"),

    make_option(c("-n", "--min-patho-count"), action="store", 
                dest="min.patho.count", type="integer",
                help="Minimum pathogenic variants"),

    make_option(c("-o", "--out"), action="store", 
                default="out.txt", dest="out.path", 
                help="Output path, [default %default]")

  )

  return(option.list)
}


FormatRegionNames <- function(top.names, sub.names) {
  # Formats top.names, sub.names into region format.
  #
  # Args:
  #   top.names: Vector of upper level names, e.g. c(A1BG)
  #   name: Vector of sub names, in format XX_[0-9], e.g. c(-_9)
  #
  # Returns:
  #   Region names in correct format, e.g. c(A1BG:-:-_9)
  middle.names <- sapply(strsplit(sub.names, "_"), function(x) x[1])
  return(paste(top.names, middle.names, sub.names, sep=":"))
}

ExtractTopNames <- function(region.names) {
  # Extracts top level gene name from region name.
  #
  # Args:
  #   region.names: Vector of region names, e.g. c(A1BG:143228:143228_0)
  # 
  # Returns:
  #   Vector of top level names, e.g. c(A1BG)
  return(sapply(strsplit(region.names, ":"), function (x) x[1]))
}

ReadScores <- function(scores.path) {
  # Reads in scores.
  #
  # Args:
  #   scores.path: path to scores file.
  #
  # Returns:
  #   Scores in data.frame.
  ret <- read.table(scores.path, stringsAsFactors=F, header=F, row.names=1)
  colnames(ret) <- c("score")
  ret$gene <- ExtractTopNames(rownames(ret))

  return(ret)
}

ReadCounts <- function(counts.path) {
  # Reads in counts.
  #
  # Args:
  #   counts.path: path to counts file.
  #
  # Returns:
  #   Counts in data.frame.
  ret <- read.table(counts.path, stringsAsFactors=F, header=F)
  colnames(ret) <- c("gene", "sub.name", "count", "mrt")
  rownames(ret) <- FormatRegionNames(ret$gene, ret$sub.name)
  ret$sub.name <- NULL

  return(ret)
}

CreateModelDataframe <- function(scores, counts) {
  # Creates a dataframe of input to the single gene test for all genes.
  #
  # Args:
  #   scores: Scores data frame. 
  #   counts: Counts data frame.
  #
  # Returns:
  #   Data frame for input to the single gene test.
  ret <- merge(scores, counts, by.x=0, by.y=0)
  rownames(ret) <- ret$Row.names
  ret$Row.names <- NULL
  if (!identical(ret$gene.x, ret$gene.y)) {
    stop("Gene vectors unequal!")
  }

  ret$gene <- ret$gene.x
  ret$gene.x <- NULL

  return(ret)
}

SingleGeneTest <- function(model.data, nperm) {
  # Runs single gene test on input data.
  #
  # Args:
  #   model.data: Data frame with count, score and mrt.
  #   nperm: Number of permutations.
  #
  # Returns:
  #   Permutation p-value for score.

  exp.counts <- sum(model.data$count) * (model.data$mrt / sum(model.data$mrt))

  dep.counts <- model.data$count - exp.counts

  dep.score.cov <- cov(dep.counts, model.data$score)

  cov.perms <- rep(NA, nperm)
  for(i in 1:nperm){
    count.perm <- rmultinom(1, sum(model.data$count), model.data$mrt / sum(model.data$mrt))
    dep.perm <- count.perm - exp.counts
    cov.perms[i] <- cov(dep.perm, model.data$score)
  }
  
  p.value <- (sum(cov.perms <= dep.score.cov) + 1) / (nperm + 1)

  return(p.value)
}

ExtractGeneEntries <- function(region.df, gene) {
  # Extracts gene entries from a region data frame.
  #
  # Args:
  #   region.df: Region data frame with gene column.
  #
  # Returns:
  #   Entries from region.df corresponding to gene.
  return(region.df[region.df$gene == gene, ])
}

### Main ----------------------------------------------------------------------
args <- commandArgs(trailing=T)

option.list <- GetOptList()
opt <- parse_args(OptionParser(option_list=option.list),
                    args=args)
## read data
counts <- ReadCounts(opt$counts.path)
scores <- ReadScores(opt$scores.path)
all.model.data <- CreateModelDataframe(scores, counts)
genes <- sort(unique(all.model.data$gene))

## iterate over genes
ret <- data.frame()

i <- 1
for (gene in genes) {
  i <- i + 1
  model.data <- ExtractGeneEntries(all.model.data, gene)
  regions.count <- nrow(model.data)
  if ((regions.count >= opt$min.regions) && 
      (sum(model.data$count) >= opt$min.patho.count)) {
    pval <- SingleGeneTest(model.data, opt$nperms)
    ret[gene, "score"] <- pval
  }
}

ret$score <- p.adjust(ret$score, method="fdr")
colnames(ret) <- c("score_fdr_pval")
write.table(ret, file=opt$out.path, quote=F, sep="\t")