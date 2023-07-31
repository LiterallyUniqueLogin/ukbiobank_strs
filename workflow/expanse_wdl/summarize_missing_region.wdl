version 1.0
import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"
import "../finemapping_wdl/finemap_one_region_workflow.wdl"
import "../finemapping_wdl/susie_one_region_workflow.wdl"
import "expanse_files.wdl"

workflow summarize_missing_region {

  input {
    String phenotype_name = "mean_platelet_volume"
  }

  call gwas_tasks.phenotype_names

  Int mpv_idx = phenotype_names.idxs[phenotype_name]

  call expanse_files.files 

  String script_dir = "."
  String finemap_command = "finemap"

  # one per chrom
  Array[VCF]+ str_vcfs = files.str_vcfs
  Array[bgen]+ imputed_snp_bgens = files.imputed_snp_bgens
  Array[File] snp_vars_to_filter_from_finemapping = files.snps_to_filter

  File shared_covars = files.shared_covars
  File phenotype_samples = files.unrelated_samples_for_ethnicity_for_phenotype[0][mpv_idx]
  File transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][mpv_idx]

  File my_str_gwas = files.str_gwas_results[mpv_idx]
  Array[File] ethnic_my_str_gwass =  files.pheno_to_ethnic_to_str_gwas_results[mpv_idx]
  File plink_snp_gwas = files.snp_gwas_results[mpv_idx]

  File all_samples_list = files.all_samples_list

  region followup_bounds = {
    "chrom": 17,
    "start": 2341352,
    "end": 2710113
  }
  String followup_regions = "~{followup_bounds.chrom}_~{followup_bounds.start}_~{followup_bounds.end}"
  Int followup_chroms = followup_bounds.chrom

  call finemapping_tasks.followup_finemapping_conditions_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    original_susie_outputs = [{
      "V": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/V.tab",
      "alpha": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/alpha.tab",
      "colnames": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/colnames.txt",
      "converged": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/converged.txt",
      "lbf": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/lbf.tab",
      "lbf_variable": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/lbf_variable.tab",
      "lfsr": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/lfsr.tab",
      "requested_coverage": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/requested_coverage.txt",
      "sigma2": "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/sigma2.txt",
    }],
    original_susie_CSs = [[
      "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/cs1.txt",
      "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/cs2.txt",
      "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/cs4.txt",
      "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/cs5.txt",
      "finemapping/susie_results/mean_platelet_volume/17_2341352_2710113/cs6.txt",
    ]],
    derived_prior_std_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_derived_effect_size_prior_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_derived_effect_size_prior_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_derived_effect_size_prior_17_2341352_2710113_finemap_output.snp",
    }],
    derived_prior_std_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_derived_effect_size_prior_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_derived_effect_size_prior_17_2341352_2710113_finemap_output.cred3",
    ]],
    original_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_first_pass_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_first_pass_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_first_pass_17_2341352_2710113_finemap_output.snp",
    }],
    original_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_first_pass_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_first_pass_17_2341352_2710113_finemap_output.cred3",
    ]],
    low_prior_std_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.snp",
    }],
    low_prior_std_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.cred3",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_low_effect_size_prior_17_2341352_2710113_finemap_output.cred4",
    ]],
    mac_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_mac_threshold_100_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_mac_threshold_100_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_mac_threshold_100_17_2341352_2710113_finemap_output.snp",
    }],
    mac_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_mac_threshold_100_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_mac_threshold_100_17_2341352_2710113_finemap_output.cred3",
    ]],
    total_prob_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_4_signals_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_4_signals_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_4_signals_17_2341352_2710113_finemap_output.snp",
    }],
    total_prob_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_4_signals_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_4_signals_17_2341352_2710113_finemap_output.cred3",
    ]],
    ratio_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_snps_over_strs_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_snps_over_strs_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_snps_over_strs_17_2341352_2710113_finemap_output.snp",
    }],
    ratio_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_snps_over_strs_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_prior_snps_over_strs_17_2341352_2710113_finemap_output.cred3",
    ]],
    threshold_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_pval_threshold_1e4_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_pval_threshold_1e4_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_pval_threshold_1e4_17_2341352_2710113_finemap_output.snp",
    }],
    threshold_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_pval_threshold_1e4_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_pval_threshold_1e4_17_2341352_2710113_finemap_output.cred3",
    ]],
    repeat_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_repeat_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_repeat_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_repeat_17_2341352_2710113_finemap_output.snp",
    }],
    repeat_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_repeat_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_repeat_17_2341352_2710113_finemap_output.cred3",
    ]],
    conv_tol_finemap_outputs = [{
      "config": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_stricter_stopping_threshold_17_2341352_2710113_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_stricter_stopping_threshold_17_2341352_2710113_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_stricter_stopping_threshold_17_2341352_2710113_finemap_output.snp",
    }],
    conv_tol_finemap_creds = [[
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_stricter_stopping_threshold_17_2341352_2710113_finemap_output.cred2",
      "wdl_cache/finemapping/mean_platelet_volume_FINEMAP_stricter_stopping_threshold_17_2341352_2710113_finemap_output.cred3",
    ]],
    best_guess_susie_outputs = [{
      "V": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_V.tab",
      "alpha": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_alpha.tab",
      "colnames": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_colnames.txt",
      "converged": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_converged.txt",
      "lbf": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_lbf.tab",
      "lbf_variable": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_lbf_variable.tab",
      "lfsr": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_lfsr.tab",
      "requested_coverage": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_requested_coverage.txt",
      "sigma2": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_sigma2.txt",
    }],
    best_guess_susie_CSs = [[
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs1.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs2.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs3.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs4.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs5.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_best_guess_genotypes_17_2341352_2710113_cs6.txt",
    ]],
    ratio_susie_outputs = [{
      "V": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_V.tab",
      "alpha": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_alpha.tab",
      "colnames": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_colnames.txt",
      "converged": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_converged.txt",
      "lbf": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_lbf.tab",
      "lbf_variable": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_lbf_variable.tab",
      "lfsr": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_lfsr.tab",
      "requested_coverage": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_requested_coverage.txt",
      "sigma2": "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_sigma2.txt",
    }],
    ratio_susie_CSs = [[
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_cs1.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_cs2.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_cs4.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_cs5.txt",
      "wdl_cache/finemapping/mean_platelet_volume_SuSiE_prior_snps_over_strs_17_2341352_2710113_cs6.txt",
    ]],
    original_regions = [followup_regions],
    original_chroms = [followup_chroms],
    followup_regions = [followup_regions],
    followup_chroms = [followup_chroms]
  }

  output {
#    Array[serializable_SuSiE_output] best_guess_susie = best_guess_susie_output_
#    Array[Array[File]] best_guess_susie_CSs = best_guess_susie_CSs_
#    Array[serializable_FINEMAP_output] repeat_finemap = repeat_finemap_output_
#    Array[Array[File]] repeat_finemap_creds = repeat_finemap_creds_
#    Array[serializable_FINEMAP_output] total_prob_finemap = total_prob_finemap_output_
#    Array[Array[File]] total_prob_finemap_creds = total_prob_finemap_creds_
#    Array[serializable_FINEMAP_output] derived_prior_std_finemap = prior_std_derived_finemap_output_
#    Array[Array[File]] derived_prior_std_finemap_creds = prior_std_derived_finemap_creds_
#    Array[serializable_FINEMAP_output] prob_conv_sss_tol_finemap = prob_conv_sss_tol_finemap_output_
#    Array[Array[File]] prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_creds_
#    Array[serializable_FINEMAP_output] mac_finemap = mac_finemap_output_
#    Array[Array[File]] mac_finemap_creds = mac_finemap_creds_
#    Array[serializable_FINEMAP_output] threshold_finemap = threshold_finemap_output_
#    Array[Array[File]] threshold_finemap_creds = threshold_finemap_creds_
#
#    Array[serializable_SuSiE_output] ratio_susie = ratio_susie_output_
#    Array[Array[File]] ratio_susie_CSs = ratio_susie_CSs_
#    Array[serializable_FINEMAP_output] ratio_finemap = ratio_finemap_output_
#    Array[Array[File]] ratio_finemap_creds = ratio_finemap_creds_
#    Array[serializable_FINEMAP_output] low_prior_std_finemap = prior_std_low_finemap_output_
#    Array[Array[File]] low_prior_std_finemap_creds = prior_std_low_finemap_creds_

    File followup_df = followup_finemapping_conditions_df.df
  }
}
