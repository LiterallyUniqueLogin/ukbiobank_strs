The contents of this directory that are part of my analysis and not raw downloaded data
(each directory has its own README):

TODO check order
Main analysis:
workflow/ - contains a file describing the various steps in running the analyses in this project that
    snakemake can execute
sample_qc/ - contains code and results from filtering samples for various quality metrics
snpstr/ - Shubham's SNP-STR reference panel and some analyses done on that reference panel
pre_imputation_qc/ - analyses done before imputing STRs into the dataset
str_imputed/ - contains the imputed STRs and the code that generated them 
traits/ - contains the phenotype data and transformations of it
association/ - contains code and results from association testing the phenotypes in traits/ against the 
    the STRs in str_imputed/ and the SNPs in array_imputed/
signals/ - contains code and results from analyzing the results in associations/ and grouping them
    into various regions
finemapping/ - contains code and results from performing fine-mapping on the results in association/
post_finemapping/ - contains code and results from interpreting the results in finemapping/

Supporting directories:
side_analyses/ - various side analyses
logs/ - outputs of SLURM jobs from any part of the analysis, sorted by name
scratch/ - a directory for writing scratch files
utilities/ - contains various installed software and scripts I wrote used in parts of the analysis

Deprecated:
exome/ - code for a side analysis that was never completed, calling STRs in UKB exome sequencing CRAMs
report - TODO

TODO
export
export_scripts
soft.ver
ukb.env: the command used to generate the default environment with the software 
    tools used in this project
