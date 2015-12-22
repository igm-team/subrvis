#!/usr/bin/Rscript
###############################################################################
# Purpose:  Run genome-wide subRVIS test.                                     #
# Author:   Ayal Gussow                                                       #
# Input:    scores table, gene list.                                          #
# Output:   Genome-Wide model results.                                        #
###############################################################################

### Libraries -----------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(MASS))

### Functions ------------------------------------------------------------------
ReadData <- function(path) {
  # Read the data from <path>.
  return(read.table(path, header=T, row.names=1))
}

RunModel <- function(model.data) {
  # Return model.
  return(glm(as.numeric(count > 0) ~ . - gene - mrt + log(mrt), 
             data=model.data, family="binomial"))
}

FormatData <- function(model.data, genes, predictors) {
  # Format data for RunModel.
  ret <- data.frame(row.names=rownames(model.data))
  for (predictor in predictors) {
    ret[[predictor]] <- model.data[[predictor]]
  }
  
  ret$count <- model.data$count
  ret$gene <- model.data$gene
  
  ret <- ret[model.data$gene %in% genes, ]
  ret <- na.omit(ret)

  for (predictor in predictors) {
    if (predictor != "mrt") {
      ret[[predictor]] <- ret[[predictor]]  / sd(ret[[predictor]])
    }
  }

  return(ret)
}

### Args -----------------------------------------------------------------------
GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--table"), action="store", 
                default="", dest="table.path", 
                help="sub regional table"),

    make_option(c("--covars"), action="store",
                default="", dest="covars",
                help="covars to include, comma separated."),

    make_option(c("--genes"), action="store",
                default="", dest="genes.path",
                help="gene list to assess."),

    make_option(c("--out"), action="store",
                default="", dest="out.txt",
                help="Output report.")
  )
}

### Main ----------------------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

option.list <- GetOptList()
opt <- parse_args(OptionParser(option_list=option.list), args=args)

genes = scan(opt$genes.path, what=character(), quiet=T)
covars <- strsplit(opt$covars, ",")[[1]]

all.model.data <- ReadData(opt$table.path)
model.data <- FormatData(all.model.data, genes, c("mrt", covars))
lr.model <- RunModel(model.data)

## write results
write.table(coef(summary(lr.model)),
            file=opt$out.txt)