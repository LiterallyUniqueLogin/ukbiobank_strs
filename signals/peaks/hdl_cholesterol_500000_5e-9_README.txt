Considering all variants from input files signals/peak_inputs/hdl_cholesterol_str_5e-9.tab and signals/peak_inputs/hdl_cholesterol_imputed_snp_5e-9.tab
Choosing association peak variants in the following order (per chromosome): Round all variants with p < 1e-300 to p=0. Take the variant with the lowest p-value that isn't within 500000 bp from any variant already selected, with the following tiebreakers (tiebreakers should only frequently occur for varaints with p rounded down to 0 and it is unlikely the second or later tiebreakers will ever be used.):
* choose the variant with the smallest starting bp
* choose SNPs over STRs
* choose SNPs with shorter reference alleles
* choose SNPs with lexicographically earlier alternate alleles
This is continued until all variants are examined. STR/SNP peak varaints within 500000 of variants of the other type that pass the p-value threshold are marked as being tagged by the other type of variants.