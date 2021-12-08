table.tab contains one row for each variant that was not filtered prior to association testing and either was reported as an association in literature, or I manually marked it as interesting or both had association p-value >= 5e-8 and FINEMAP posterior probability of causaility >= 0.05

table.tab contains the following columns:
chrom
signal_region - {start bp}_{end bp} of region of association the STR resides in, 1-based, inclusive
start_pos - 1-based, inclusive
end_pos - 1-based, inclusive
SNPSTR_start_pos - the start position of the STR in the SNP-STR reference panel (may be smaller that start_pos if HipSTR called this STR with a physically phased SNP upstream of the STR)
SNPSTR_ID - ID of the STR in the SNPSTR reference panel
period - as specified in the SNPSTR reference panel
repeat_unit - as inferred by the algorithm described in the paper methods
alleles - possible # copies of the repeat unit
reference_allele
total_per_allele_dosages - sum of imputed dosage of each allele across both chromsomes in all samples
total_hardcall_alleles - number of each allele in the population using Beagle's called most probable phased genotypes
total_hardcall_genotypes - number of each unphased genotype in the population using Beagle's called most probable phased genotypes
subset_total_per_allele_dosages - as total_per_allele_dosages, but calculated on the subset of samples used in the association
subset_total_hardcall_alleles - as total_hardcall_alleles , but calculated on the subset of samples used in the association
subset_total_hardcall_genotypes - as total_hardcall_genotypes, but calculated on the subset of samples used in the association
subset_multiallelicness - % of dosage attributed to alleles that are not the first or second most frequent alleles, by dosage, within the subset of samples used in the association
subset_heterozygosity - the expected percent of heterozygous samples in the association sample subset as calculated from allele frequencies
subset_entropy - the entropy in the association sample subset calculated from allele frequencies
subset_HWEP - the Hardy-Weinberg p-value in the association sample subset comparing hardcall genotypes and allele frequencies
subset_allele_dosage_r2 - a metric of imputation accuracy, the per-allele R^2 between dosage-weighted allele length and allele hardcalls across samples
association_p_value - linear regression of rank inverse normalized phenotype values vs repeat length on QC'ed sample subset, with covariates
direction_of_association - + if an increase in STR length causes an increase in phenotype, - if the reverse, NaN if p-value > 0.05
Δphenotype_per_additional_repeat_unit - linear regression on raw phenotypes vs repeat length on QC'ed sample subset, no covariates included, phenotype measured in mmol/mol
Δphenotype_per_s.d._increase_in_repeat_size - linear regression on raw phenotypes vs repeat length on QC'ed sample subset, no covariates included, phenotype measured in mmol/mol
pcausal - FINEMAP posterior probability of causality
mentioned_in_literature - Whether or not this we know of a citation saying this STR is likely causal for this trait. False here means we have not checked whether or not this is reported in the literature, not that it has not.
literature_inclusion_url - If mentioned_in_literature, then the corresponding ULR. Otherwise NA. 
included_only_due_to_literature - NA if not mentioned in literature. True if the locus would have been filtered, but is included due to literature. False otherwise. 
curated - Whether or not manual examination was used to decide this was an exciting locus worthy of extra attention. 
included_only_due_to_curation - NA if not curated. True if the locus would have been filtered, but is included due to manual curation. False if the locus was curated but would have been included regardless (not filtered). 
nearby_exons - A comma separated list of distance:'upstream'/'downstream'/'inside':exon-start:exon-end:gene-name:gene-type, where upstream means upstream of the exon in that gene's direction. All exons within a 1000bp radius, or 'none'. (See https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.) 
closest_gene - distance:'upstream'/'downstream'/'inside':gene-name:gene-type. Possibly a comma separated list if multiple genes are tied by distance. (See https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.)
nearby_genes - A comma separated list of entries like that in closest_gene. All genes within a 100kbp radius, or 'none'.
relation_to_gene - If this STR intersects a single gene with a single category, then 'CDS'/'5_prime_UTR'/'3_prime_UTR'/'unannotated_UTR'/'exon_other'/'intron':gene-name:gene-type. If it intersects multiple categories of the same gene, then those categories are comma separated before the gene and gene type. If it intersects multiple genes, then 'multigene;' followed by the categorization for each of the genes as above, separated by semicolons. If it intersects no genes, then 'intergenic'. If it partially intersects at least one gene, then 'intergenic;' followed by the information above. Note that STRs which have only exonic annotations might still be introns in relation to other transcripts of the gene they overlap. (See https://www.gencodegenes.org/pages/biotypes.html for gene type meanings.) 
transcribed - Does this intersect a transcript? If so, 'transcribed;' followed by a comma separated list of transcript-name:transcript-type:transcript-support-level. Otherwise 'untranscribed'.(See https://www.gencodegenes.org/pages/biotypes.html for transcript-type meanings and https://www.gencodegenes.org/pages/data_format.html for transcript-support-level meanings.)

