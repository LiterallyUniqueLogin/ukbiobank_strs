version 1.0

import "finemapping_tasks.wdl"
import "retryable_susie_load_gts.wdl"
import "escalating_susie_run.wdl"

workflow susie_one_region {

  input {
    String script_dir

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    File? shared_covars
    File phenotype_samples
    File transformed_phenotype_data

    File? my_str_gwas
    File plink_snp_gwas

    String phenotype_name
    region bounds

    File all_samples_list

    Float? tol
    Float? snp_p_over_str_p
    Float? res_var
    Float? prior_var
    Boolean best_guess = false
  }

  # call SuSiE 

  Int chrom_minus_one = bounds.chrom-1
 
  call finemapping_tasks.susie_choose_vars { input :
    script_dir = script_dir,
    str_assoc_results = my_str_gwas,
    snp_assoc_results = plink_snp_gwas,
    variants_to_filter = snp_vars_to_filter_from_finemapping[chrom_minus_one],
    phenotype_name = phenotype_name,
    bounds = bounds,
  }

  call retryable_susie_load_gts.retryable_susie_load_gts { input :
    script_dir = script_dir,
    strs = str_vcfs[chrom_minus_one],
    snps = imputed_snp_bgens[chrom_minus_one],
    all_samples = all_samples_list,
    phenotype_samples = phenotype_samples,
    pheno_data = transformed_phenotype_data,
    shared_covars = shared_covars,
    colnames = susie_choose_vars.colnames,
    phenotype_name = phenotype_name,
    bounds = bounds,
    best_guess = best_guess
  }

  call escalating_susie_run.escalating_susie_run { input :
    script_dir = script_dir,
    gts_h5 = retryable_susie_load_gts.gts_h5,
    pheno_residuals_h5 = retryable_susie_load_gts.pheno_residuals_h5,
    varnames_file = susie_choose_vars.colnames,
    tol = tol,
    snp_p_over_str_p = snp_p_over_str_p,
    res_var = res_var,
    prior_var = prior_var,
    colnames = susie_choose_vars.colnames
  }

  output {
    SuSiE_output susie_output = escalating_susie_run.susie_output
  }
}
