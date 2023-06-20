version 1.0

import "expanse_files.wdl"
import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/finemapping_workflow.wdl"

workflow finemapping {

  input {
    String script_dir  = "."

    String phenotype_name
  }

  call expanse_files.files

  call gwas_tasks.phenotype_names

  Int phenotype_idx = phenotype_names.idxs[phenotype_name]

  call finemapping_workflow.finemapping { input :
    script_dir = script_dir,
    finemap_command = "finemap",

    chr_lens = files.chr_lens,

    str_vcfs = files.str_vcfs,
    imputed_snp_bgens = files.imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = files.snps_to_filter,

    shared_covars = files.shared_covars,
    phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
    transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][phenotype_idx],
    
    my_str_gwas = files.str_gwas_results[phenotype_idx],
    ethnic_my_str_gwass =  files.pheno_to_ethnic_to_str_gwas_results[phenotype_idx],
    plink_snp_gwas = files.snp_gwas_results[phenotype_idx],

    phenotype_name = phenotype_name,
    is_binary = false,

    all_samples_list = files.all_samples_list
  }

  output {
    File first_pass_regions_tsv = finemapping.first_pass_regions_tsv
    File first_pass_regions_readme = finemapping.first_pass_regions_readme

    Array[serializable_SuSiE_output] original_susie = finemapping.original_susie
    Array[Array[File]] original_susie_CSs = finemapping.original_susie_CSs
    Array[serializable_FINEMAP_output] original_finemap = finemapping.original_finemap
    Array[Array[File]] original_finemap_creds = finemapping.original_finemap_creds

    File first_pass_df = finemapping.first_pass_df
    File susie_min_abs_corrs = finemapping.susie_min_abs_corrs

    File followup_regions_tsv = finemapping.followup_regions_tsv

    Array[serializable_SuSiE_output] best_guess_susie = finemapping.best_guess_susie
    Array[Array[File]] best_guess_susie_CSs = finemapping.best_guess_susie_CSs
    Array[serializable_FINEMAP_output] total_prob_finemap = finemapping.total_prob_finemap
    Array[Array[File]] total_prob_finemap_creds = finemapping.total_prob_finemap_creds
    Array[serializable_FINEMAP_output] derived_prior_std_finemap = finemapping.derived_prior_std_finemap
    Array[Array[File]] derived_prior_std_finemap_creds = finemapping.derived_prior_std_finemap_creds
    Array[serializable_FINEMAP_output] prob_conv_sss_tol_finemap = finemapping.prob_conv_sss_tol_finemap
    Array[Array[File]] prob_conv_sss_tol_finemap_creds = finemapping.prob_conv_sss_tol_finemap_creds
    Array[serializable_FINEMAP_output] mac_finemap = finemapping.mac_finemap
    Array[Array[File]] mac_finemap_creds = finemapping.mac_finemap_creds
    Array[serializable_FINEMAP_output] threshold_finemap = finemapping.threshold_finemap
    Array[Array[File]] threshold_finemap_creds = finemapping.threshold_finemap_creds

    Array[serializable_SuSiE_output] ratio_susie = finemapping.ratio_susie
    Array[Array[File]] ratio_susie_CSs = finemapping.ratio_susie_CSs
    Array[serializable_FINEMAP_output] ratio_finemap = finemapping.ratio_finemap
    Array[Array[File]] ratio_finemap_creds = finemapping.ratio_finemap_creds
    Array[serializable_FINEMAP_output] low_prior_std_finemap = finemapping.low_prior_std_finemap
    Array[Array[File]] low_prior_std_finemap_creds = finemapping.low_prior_std_finemap_creds

    File? followup_df = finemapping.followup_df
  }
}
