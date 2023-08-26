Columns:
chromosome
base_pair_location: beginning of the repeat (hg19, 1-indexed, inclusive)
alleles:
  lengths alleles in the population measured in number of repeat units. For example, the allele 5 for an
  AC repeat implies the bases "ACACACACAC" (possibly with some impurity).
  Occasionally the repeat unit will be listed as none, but the period of the repeat will still be given.
  In that case, the length of an allele in base pairs can still be calculated by multiplying the allele by the period.
beta:
  measured effect size of the linear association of the rank-inverse-normalized phenotype against the 
  length-dosages of unnormalized STR genotypes, measured in number of repeat units.
  As phenotypes are rank-inverse-normalized, this is only comparable to other betas from this study.
  p-values are more appropriate for comparison across studies.
standard_error: See caveats for beta
allele_frequencies
p_value: p-values less than 1e-300 exceeded our software's numeric precision and are listed as 0
ref_allele: measured in number of repeat units
repeat_unit: the standardized repeat unit of this STR, or none if there was no one clear repeat unit
period: the length of the repeat unit
end_pos (hg19): end of the repeat (1-indexed, inclusive)
start_pos (hg38)
end_pos (hg38)
n:
  this study worked with imputed calls and no call-level filters, as such n will be equivalent for each
  variant associated with the same phenotype
number_of_common_alleles:
  the number of alleles in the population with frequency >= 1%
mean_{phenotype}_per_summed_gt:
  the mean phenotype value for each sum of allele lengths, where each participant's contribution to the
  phenotype mean for each length-sum is weighted by the imputed probability of their true genotype sum being
  equal to that length-sum. Can be used for plotting graphs of mean phenotype value vs summed-length.
  Summed gts are measured in number of summed repeat units.
summed_0.05_significance_CI:
  The 95% symmetric confidence interval for each of the means above
summed_5e-8_significance_CI:
  The (1 - 5e-8) symmetric confidence interval for each of the means above
mean_{phenotype}_per_paired_gt:
  the mean phenotype value for each unordered pair of allele lengths, where each participant's contribution to the
  phenotype mean for each pair  is weighted by the imputed probability of their true genotype pair being
  equal to that pair. Can be used for plotting graphs of mean phenotype value vs length pair.
  Each gt in each pair is measured in number of repeat units.
paired_0.05_significance_CI:
  The 95% symmetric confidence interval for each of the means above
paired_5e-8_significance_CI:
  The (1 - 5e-8) symmetric confidence interval for each of the means above
