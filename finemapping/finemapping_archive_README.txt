* One tarball finemapping_${trait}.tgz for each trait

* The tarball contains two summary files
finemapping_first_pass.tab
finemapping_followup.tab
* In these files
 - positions are hg19 and 1-indexed
 - Regions are denoted by ${chrom}_${start_pos_inclusive}_${end_pos_inclusive}
 - variant names for SNPS are SNP_{position}_{ref}_{alt}
 - variant names for STRs are STR_{start_position}
 - p_val, coeff and se are from the association of that variant with the trait being fine-mapped
 - susie_cs values denote which pure cs a variant is present in, or -1 if its in no pure cs
 - finemap_pip and susie_alpha values were used as CP values in the paper
 - susie_cs values of -1 imply that the susie_alpha value should be ignored (treated as zero)
 - only variants included in fine-mapping are present in these files
 - some variants lack susie values if they were not included in the susie fine-mapping (probably because they had too high a p-value)
* The followup files only contain variants from the regions subjectted to follow-up fine-mapping conditions.
 - traits with no regions that were followed-up on are empty.
 - In addition to the columns in the first pass file (except for coeff or se), there are additional FINEMAP or SuSiE columns for each extra condition.
 - the best_guess column corresponds to the use of best guess genotypes from imputation for fine-mapping
 - the ratio columns correspond to the prior of favoring SNP over STR causality by a 4-to-1 ratio
 - the repeat column corresponds to the repeat FINEMAP run with no changed settings
 - the total_prob column corresponds to the FINEMAP run with a prior of there being 4 total causal variants
 - the prior_std_derived and prior_std_low columns correspond to the priors for the effect sizes of causal variants of 0.05% and 0.0025%, respectively
 - the conv_tol column corresponds to the flag --prob-conv-sss-tol 0.0001
 - the mac column corresponds to the non-major allele dosage threshold of 100
 - the p_thresh column corresponds to the p_value threshold of 1e-4
 - values for the additional columns will be missing for variants which were not fine-mapped in specifically those conditions (say, a variant with p-value of 1e-3 in the 1e-4 threshold column)

* In addition to the two summary files, the tarballs contain raw fine-mapping output files for each region.
* The files for FINEMAP are ( described at http://christianbenner.com/ ):
FINEMAP_first_pass_${region}_finemap_output.log
FINEMAP_first_pass_${region}_finemap_output.snp
FINEMAP_first_pass_${region}_finemap_output.config
FINEMAP_first_pass_${region}_finemap_output.credX
* and the following files for SuSiE :
SuSiE_first_pass_${region}_alpha.tab
SuSiE_first_pass_${region}_colnames.txt
SuSiE_first_pass_${region}_csX.txt
SuSiE_first_pass_${region}_lbf.tab
SuSiE_first_pass_${region}_lbf_variable.tab
SuSiE_first_pass_${region}_lfsr.tab
SuSiE_first_pass_${region}_sigma2.txt
SuSiE_first_pass_${region}_V.tab
* Except for the colnames and cs files, these are arrays written from the output fields of the susie() function described at https://stephenslab.github.io/susieR/reference/susie.html
* colnames contains one variant name per line, each line corresponding to one column of the alpha array
* csX contains three rows:
 - the first contains the 1-indexed numbers of the variables included in the credible set, in ascending order
 - the second contains the coverage of the credible set (Always greater than 0.9 which was the requested coverage)
 - the third contains three numbers, the min, mean and median absolute correlations between each variable in the credible set
   ( see susie_get_cs() at https://stephenslab.github.io/susieR/reference/susie_get_methods.html )
 - the 1-indexed number in the filename corresponds to the corresponding row of the alpha.tab array
* Additionally, for each of the regions we followed up on, this tarball contains the same files as above, but with
  the following prefixes, corresponding to the follow-up fine-mapping condition being tested:
FINEMAP_derived_effect_size_prior - (effect size prior of 0.05%)
FINEMAP_low_effect_size_prior - (effect size prior of 0.0025%)
FINEMAP_mac_threshold_100 - (non-major allele dosage threshold of 100)
FINEMAP_prior_4_signals - (prior of 4 causal variants per region)
FINEMAP_prior_snps_over_strs - (prior of favoring SNP over STR causality by a 4-to-1 ratio)
FINEMAP_pval_threshold_1e4 - (p_value threshold of 1e-4)
FINEMAP_stricter_stopping_threshold - (using the flag --prob-conv-sss-tol 0.0001)
SuSiE_prior_snps_over_strs - (prior of favoring SNP over STR causality by a 4-to-1 ratio)
SuSiE_best_guess_genotypes - (the use of best guess genotypes from imputation for fine-mapping)
* See the paper and Supplementary Note 3 for more details on each follow-up condition.
* We do not have raw output files for FINEMAP_repeat runs.

Notes:
* For FINEMAP, I have focused on the posterior probabilities in the finemap_output.snp files
* For SuSiE I focused on the cs files with min absolute correlation > 0.8,
  and then the values in alpha.tab with rows identified by the CS number and columns identified by the variants in the first row of the CS file
* the log files for FINEMAP saying v1.4.1 seem to be buggy, the actual FINEMAP output to the command line indicates that v1.4.2 was run.
