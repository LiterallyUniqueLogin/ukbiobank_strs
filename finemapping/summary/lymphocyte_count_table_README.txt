table.tab contains one row for each variant that was not filtered prior to association testing and either was previously reported as an association or both hadassociation p-value >= 5e-8 andFINEMAP posterior probability of causaility >= 0.05

table.tab contains the following columns:
chrom
signal_region - {start bp}_{end bp} of region of association the STR resides in, 1-based, inclusive
start_pos - 1-based, inclusive
end_pos - 1-based, inclusive
SNPSTR_start_pos - the start position of the STR in the SNP-STR reference panel (may be smaller that start_pos if HipSTR called this STR with a physically phased SNP upstream of the STR)
SNPSTR_ID - ID of the STR in the SNPSTR reference panel
period - as specified in the SNPSTR reference panel
repeat_unit - as inferred by TRTools from the STR sequence and the period
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
Δphenotype_per_additional_repeat_unit - linear regression on raw phenotypes vs repeat length on QC'ed sample subset, no covariates included, phenotype measured in 10^9 cells/Litre
Δphenotype_per_s.d._increase_in_repeat_size - linear regression on raw phenotypes vs repeat length on QC'ed sample subset, no covariates included, phenotype measured in 10^9 cells/Litre
pcausal - FINEMAP posterior probability causality
previously_reported_association - Whether or not this row was included in the table because it was reported previously in the literature

