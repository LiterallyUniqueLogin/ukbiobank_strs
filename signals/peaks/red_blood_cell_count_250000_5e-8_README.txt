Considering all variants from input files association/results/red_blood_cell_count/my_str/results.tab and association/results/red_blood_cell_count/plink_snp/results.tab
Choosing association peak variants in the following order (per chromosome): Round all variants with p < 1e-300 to p=0. Take the variant with the lowest p-value that isn't within 250000 bp from any variant already selected, with the following tiebreakers (tiebreakers should only frequently occur for varaints with p rounded down to 0 and it is unlikely the second or later tiebreakers will ever be used.):
* choose the variant with the smallest starting bp
* choose SNPs over STRs
* choose SNPs with shorter reference alleles
* choose SNPs with lexicographically earlier alternate alleles
This is continued until all variants with p_value < 5e-08 are examined.
STR/SNP peak varaints within 250000 of variants of the other type that pass the p-value threshold are marked as being tagged by the other type of variants.