# SubRVIS
This repository contains the software used in the subRVIS publication (http://bit.ly/2vY1KP5).

Citation:

Gussow AB, Petrovski S, Wang Q, Allen AS, Goldstein DB. The intolerance to functional genetic variation of protein domains predicts the localization of pathogenic mutations within genes. Genome Biology. 2016; 17(1):9.

## Regional Score Gene-Specific Prediction Test Model
To run the single gene assessment used in the manuscript, the usage is:

    Rscript src/single_gene_tests.R -s <scores> -c <counts> --perms <number of permutations> -n <minimum pathogenic variants> -m <minimum regions> -o <out file path>

For example:

    Rscript src/single_gene_tests.R -s data/domains.rvis -c data/domains_counts.txt --perms 20000 -n 1 -m 2 -o subrvis_single_gene_test_results.txt

## Regional Score Genome-Wide Prediction Test Model
To run the genome-wide assessment used in the manuscript, the usage is:

    Rscript src/genome_wide_test.R --table <scores table> --genes <gene list> --covars <predictor names> --out <out file path>

For example:

    Rscript src/genome_wide_test.R --table data/domains_table.txt --genes data/genes.txt --covars subRVIS --out subrvis_genome_wide_test_results.txt

## subRVIS Variant Plotting
The variant_plotting directory contains the subRVIS website tool (www.subrvis.org), which is based on the Shiny framework (http://shiny.rstudio.com).
