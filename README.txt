This is the code repository for most of the analyses and figures described in the paper:
Polymorphic short tandem repeats make widespread contributions to blood and serum traits
Jonathan Margoliash, Shai Fuchs, Yang Li, Xuan Zhang, Arya Massarat, Alon Goren, Melissa Gymrek
bioRxiv 2022.08.01.502370; doi: https://doi.org/10.1101/2022.08.01.502370

workflow/expanse_wdl/figures.wdl is the entry point for generating most of this.

Code directories (some scripts are old have not been migrated to the WDL pipeline): 
association - running the STR GWAS. Written before splitting associaTR into its own package.
export_scripts - miscellaneous old scripts
finemapping - for fine-mapping STRs with SNPs
post_finemapping - for analyzing results for the study after fine-mapping and all other analyses are done
pre_imputation_qc - old scripts contributing to information used to motivate imputation choices
sample_qc - for filtering samples and dividing ethnicities
side_analyses - variaous side analyses, not organized
signals - for calculating peaks and fine-mapping trait regions
str_imputed - old scripts used to run STR imputation with Beagle
traits - for extracting and transforming phenotypes and covariates from the main dataset
utilities - various utilities
wgs - scripts for dealing with WGS data on the UKB RAP or dealing with summaries generated from that data
workflow - workflow scripts for orchestrating the full analysis

Dataset directories (actual contents not all stored in the repository):
array_imputed - imputed SNP + indel genotypes 
main_dataset - UKB main dataset available on the UKB showcase
microarray - phased hard-called array genotypes
misc_data - miscillaneous additional datasets
snpstr - SNP-STR reference panel related to Saini et al.
