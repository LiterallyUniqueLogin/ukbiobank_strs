version 1.0
import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"
import "../finemapping_wdl/finemap_one_region_workflow.wdl"
import "../finemapping_wdl/susie_one_region_workflow.wdl"
import "expanse_files.wdl"

workflow ldl_apob_region {

  input {
    String phenotype_name = "ldl_cholesterol_direct"
  }

  call gwas_tasks.phenotype_names

  Int ldl_idx = phenotype_names.idxs[phenotype_name]

  call expanse_files.files 

  String script_dir = "."
  String finemap_command = "finemap"

  # one per chrom
  Array[VCF]+ str_vcfs = files.str_vcfs
  Array[bgen]+ imputed_snp_bgens = files.imputed_snp_bgens
  Array[File] snp_vars_to_filter_from_finemapping = files.snps_to_filter

  File shared_covars = files.shared_covars
  File phenotype_samples = files.unrelated_samples_for_ethnicity_for_phenotype[0][ldl_idx]
  File transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][ldl_idx]

  File my_str_gwas = files.str_gwas_results[ldl_idx]
  Array[File] ethnic_my_str_gwass =  files.pheno_to_ethnic_to_str_gwas_results[ldl_idx]
  File plink_snp_gwas = files.snp_gwas_results[ldl_idx]

  File all_samples_list = files.all_samples_list

  region followup_bounds = {
    "chrom": 2,
    "start": 20225783,
    "end": 21819451
  }
  String followup_regions = "~{followup_bounds.chrom}_~{followup_bounds.start}_~{followup_bounds.end}"
  Int followup_chroms = followup_bounds.chrom

  call finemap_one_region_workflow.finemap_one_region as repeat_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    prefix = "~{phenotype_name}_FINEMAP_repeat_~{followup_regions}_",
    cache_breaker = 1
  }
  serializable_FINEMAP_output repeat_finemap_output_ = repeat_finemap_.finemap_output.subset
  Array[File] repeat_finemap_creds_ = repeat_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as ratio_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    snp_str_ratio = 4,
    prefix = "~{phenotype_name}_FINEMAP_prior_snps_over_strs_~{followup_regions}_"
  }
  serializable_FINEMAP_output ratio_finemap_output_ = ratio_finemap_.finemap_output.subset
  Array[File] ratio_finemap_creds_ = ratio_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as total_prob_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    total_prob = 4,
    prefix = "~{phenotype_name}_FINEMAP_prior_4_signals_~{followup_regions}_"
  }
  serializable_FINEMAP_output total_prob_finemap_output_ = total_prob_finemap_.finemap_output.subset
  Array[File] total_prob_finemap_creds_ = total_prob_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as prior_std_derived_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    prior_std = 0.0224,
    prefix = "~{phenotype_name}_FINEMAP_derived_effect_size_prior_~{followup_regions}_"
  }
  serializable_FINEMAP_output prior_std_derived_finemap_output_ = prior_std_derived_finemap_.finemap_output.subset
  Array[File] prior_std_derived_finemap_creds_ = prior_std_derived_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as prior_std_low_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    prior_std = 0.005,
    prefix = "~{phenotype_name}_FINEMAP_low_effect_size_prior_~{followup_regions}_"
  }
  serializable_FINEMAP_output prior_std_low_finemap_output_ = prior_std_low_finemap_.finemap_output.subset
  Array[File] prior_std_low_finemap_creds_ = prior_std_low_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as prob_conv_sss_tol_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    prob_conv_sss_tol = 0.0001,
    prefix = "~{phenotype_name}_FINEMAP_stricter_stopping_threshold_~{followup_regions}_"
  }
  serializable_FINEMAP_output prob_conv_sss_tol_finemap_output_ = prob_conv_sss_tol_finemap_.finemap_output.subset
  Array[File] prob_conv_sss_tol_finemap_creds_ = prob_conv_sss_tol_finemap_.finemap_output.creds

  Int followup_chroms_minus_one = followup_chroms - 1
  call gwas_tasks.imputed_snp_frequencies { input :
    script_dir = script_dir,
    imputed_snp_bgen = imputed_snp_bgens[followup_chroms_minus_one],
    samples = phenotype_samples,
    all_samples = all_samples_list,
    bounds = followup_bounds,
  }

  call finemap_one_region_workflow.finemap_one_region as mac_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    mac = 100,
    snp_macs = imputed_snp_frequencies.counts,
    prefix = "~{phenotype_name}_FINEMAP_mac_threshold_100_~{followup_regions}_"
  }
  serializable_FINEMAP_output mac_finemap_output_ = mac_finemap_.finemap_output.subset
  Array[File] mac_finemap_creds_ = mac_finemap_.finemap_output.creds

  call finemap_one_region_workflow.finemap_one_region as threshold_finemap_ { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    inclusion_threshold = 0.0005,
    prefix = "~{phenotype_name}_FINEMAP_pval_threshold_1e4_~{followup_regions}_"
  }
  serializable_FINEMAP_output threshold_finemap_output_ = threshold_finemap_.finemap_output.subset
  Array[File] threshold_finemap_creds_ = threshold_finemap_.finemap_output.creds

  call susie_one_region_workflow.susie_one_region as best_guess_susie_ { input :
    script_dir = script_dir,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    shared_covars = shared_covars,
    phenotype_samples = phenotype_samples,
    transformed_phenotype_data = transformed_phenotype_data,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    best_guess = true,
    prefix = "~{phenotype_name}_SuSiE_best_guess_genotypes_~{followup_regions}_"
  }
  serializable_SuSiE_output best_guess_susie_output_ = best_guess_susie_.susie_output.subset
  Array[File] best_guess_susie_CSs_ = best_guess_susie_.susie_output.CSs

  call susie_one_region_workflow.susie_one_region as ratio_susie_ { input :
    script_dir = script_dir,
    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    shared_covars = shared_covars,
    phenotype_samples = phenotype_samples,
    transformed_phenotype_data = transformed_phenotype_data,
    my_str_gwas = my_str_gwas,
    plink_snp_gwas = plink_snp_gwas,
    phenotype_name = phenotype_name,
    bounds = followup_bounds,
    all_samples_list = all_samples_list,
    snp_p_over_str_p = 4,
    prefix = "~{phenotype_name}_SuSiE_prior_snps_over_strs_~{followup_regions}_"
  }
  serializable_SuSiE_output ratio_susie_output_ = ratio_susie_.susie_output.subset
  Array[File] ratio_susie_CSs_ = ratio_susie_.susie_output.CSs

   call finemapping_tasks.followup_finemapping_conditions_df { input :
    script_dir = ".",
    phenotype_name = phenotype_name,
    snp_assoc_results = files.snp_gwas_results[ldl_idx],
    str_assoc_results = files.str_gwas_results[ldl_idx],
    ethnic_str_assoc_results = files.pheno_to_ethnic_to_str_gwas_results[ldl_idx],
    original_susie_outputs = [{
      "V": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/V.tab",
      "alpha": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/alpha.tab",
      "colnames": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/colnames.txt",
      "converged": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/converged.txt",
      "lbf": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/lbf.tab",
      "lbf_variable": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/lbf_variable.tab",
      "lfsr": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/lfsr.tab",
      "requested_coverage": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/requested_coverage.txt",
      "sigma2": "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/sigma2.txt",
    }],
    original_susie_CSs = [[
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs1.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs2.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs3.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs4.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs5.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs6.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs7.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs8.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs9.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs10.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs11.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs12.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs13.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs14.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs15.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs16.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs17.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs18.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs19.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs20.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs21.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs22.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs23.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs24.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs25.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs26.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs27.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs28.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs29.txt",
      "finemapping/susie_results/ldl_cholesterol_direct/2_20225783_21819451/cs30.txt",
    ]],
    original_finemap_outputs = [{
      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.config",
      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.log_sss",
      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.snp",
    }],
    original_finemap_creds = [[
      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.cred15",
      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.cred16",
      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_first_pass_2_20225783_21819451_finemap_output.cred17",
    ]],
    best_guess_susie_outputs = [best_guess_susie_output_],
    best_guess_susie_CSs = [best_guess_susie_CSs_],
    repeat_finemap_outputs = [repeat_finemap_output_],
    repeat_finemap_creds = [repeat_finemap_creds_],
    total_prob_finemap_outputs = [total_prob_finemap_output_],
    total_prob_finemap_creds = [total_prob_finemap_creds_],
    derived_prior_std_finemap_outputs = [prior_std_derived_finemap_output_],
    derived_prior_std_finemap_creds = [prior_std_derived_finemap_creds_],
    conv_tol_finemap_outputs = [prob_conv_sss_tol_finemap_output_],
    conv_tol_finemap_creds = [prob_conv_sss_tol_finemap_creds_],
    mac_finemap_outputs = [mac_finemap_output_],
    mac_finemap_creds = [mac_finemap_creds_],
    threshold_finemap_outputs = [threshold_finemap_output_],
    threshold_finemap_creds = [threshold_finemap_creds_],
    ratio_susie_outputs = [ratio_susie_output_],
    ratio_susie_CSs = [ratio_susie_CSs_],
    ratio_finemap_outputs = [ratio_finemap_output_],
    ratio_finemap_creds = [ratio_finemap_creds_],
    low_prior_std_finemap_outputs = [prior_std_low_finemap_output_],
    low_prior_std_finemap_creds = [prior_std_low_finemap_creds_],
    original_regions = [followup_regions],
    original_chroms = [followup_chroms],
    followup_regions = [followup_regions],
    followup_chroms = [followup_chroms]
  }

#  call finemapping_tasks.followup_finemapping_conditions_df { input :
#    script_dir = script_dir,
#    phenotype_name = phenotype_name,
#    snp_assoc_results = plink_snp_gwas,
#    str_assoc_results = my_str_gwas,
#    ethnic_str_assoc_results = ethnic_my_str_gwass,
#
#    derived_prior_std_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_derived_effect_size_prior_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_derived_effect_size_prior_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_derived_effect_size_prior_2_20225783_21819451_finemap_output.snp",
#    }],
#    derived_prior_std_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_derived_effect_size_prior_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_derived_effect_size_prior_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    
#    low_prior_std_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.snp",
#    }],
#    low_prior_std_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.cred3",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_low_effect_size_prior_2_20225783_21819451_finemap_output.cred4",
#    ]],
#    mac_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_mac_threshold_100_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_mac_threshold_100_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_mac_threshold_100_2_20225783_21819451_finemap_output.snp",
#    }],
#    mac_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_mac_threshold_100_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_mac_threshold_100_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    total_prob_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_4_signals_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_4_signals_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_4_signals_2_20225783_21819451_finemap_output.snp",
#    }],
#    total_prob_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_4_signals_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_4_signals_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    ratio_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_snps_over_strs_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_snps_over_strs_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_snps_over_strs_2_20225783_21819451_finemap_output.snp",
#    }],
#    ratio_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_snps_over_strs_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_prior_snps_over_strs_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    threshold_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_pval_threshold_1e4_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_pval_threshold_1e4_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_pval_threshold_1e4_2_20225783_21819451_finemap_output.snp",
#    }],
#    threshold_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_pval_threshold_1e4_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_pval_threshold_1e4_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    repeat_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_repeat_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_repeat_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_repeat_2_20225783_21819451_finemap_output.snp",
#    }],
#    repeat_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_repeat_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_repeat_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    conv_tol_finemap_outputs = [{
#      "config": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_stricter_stopping_threshold_2_20225783_21819451_finemap_output.config",
#      "log_sss": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_stricter_stopping_threshold_2_20225783_21819451_finemap_output.log_sss",
#      "snp_file": "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_stricter_stopping_threshold_2_20225783_21819451_finemap_output.snp",
#    }],
#    conv_tol_finemap_creds = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_stricter_stopping_threshold_2_20225783_21819451_finemap_output.cred2",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_FINEMAP_stricter_stopping_threshold_2_20225783_21819451_finemap_output.cred3",
#    ]],
#    best_guess_susie_outputs = [{
#      "V": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_V.tab",
#      "alpha": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_alpha.tab",
#      "colnames": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_colnames.txt",
#      "converged": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_converged.txt",
#      "lbf": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_lbf.tab",
#      "lbf_variable": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_lbf_variable.tab",
#      "lfsr": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_lfsr.tab",
#      "requested_coverage": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_requested_coverage.txt",
#      "sigma2": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_sigma2.txt",
#    }],
#    best_guess_susie_CSs = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs1.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs2.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs3.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs4.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs5.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_best_guess_genotypes_2_20225783_21819451_cs6.txt",
#    ]],
#    ratio_susie_outputs = [{
#      "V": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_V.tab",
#      "alpha": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_alpha.tab",
#      "colnames": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_colnames.txt",
#      "converged": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_converged.txt",
#      "lbf": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_lbf.tab",
#      "lbf_variable": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_lbf_variable.tab",
#      "lfsr": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_lfsr.tab",
#      "requested_coverage": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_requested_coverage.txt",
#      "sigma2": "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_sigma2.txt",
#    }],
#    ratio_susie_CSs = [[
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_cs1.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_cs2.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_cs4.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_cs5.txt",
#      "wdl_cache/finemapping/ldl_cholesterol_direct_SuSiE_prior_snps_over_strs_2_20225783_21819451_cs6.txt",
#    ]],
#    original_regions = [followup_regions],
#    original_chroms = [followup_chroms],
#    followup_regions = [followup_regions],
#    followup_chroms = [followup_chroms]
#  }

  output {
    serializable_SuSiE_output best_guess_susie = best_guess_susie_output_
    Array[File] best_guess_susie_CSs = best_guess_susie_CSs_
    serializable_FINEMAP_output repeat_finemap = repeat_finemap_output_
    Array[File] repeat_finemap_creds = repeat_finemap_creds_
    serializable_FINEMAP_output total_prob_finemap = total_prob_finemap_output_
    Array[File] total_prob_finemap_creds = total_prob_finemap_creds_
    serializable_FINEMAP_output derived_prior_std_finemap = prior_std_derived_finemap_output_
    Array[File] derived_prior_std_finemap_creds = prior_std_derived_finemap_creds_
    serializable_FINEMAP_output prob_conv_sss_tol_finemap = prob_conv_sss_tol_finemap_output_
    Array[File] prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_creds_
    serializable_FINEMAP_output mac_finemap = mac_finemap_output_
    Array[File] mac_finemap_creds = mac_finemap_creds_
    serializable_FINEMAP_output threshold_finemap = threshold_finemap_output_
    Array[File] threshold_finemap_creds = threshold_finemap_creds_

    serializable_SuSiE_output ratio_susie = ratio_susie_output_
    Array[File] ratio_susie_CSs = ratio_susie_CSs_
    serializable_FINEMAP_output ratio_finemap = ratio_finemap_output_
    Array[File] ratio_finemap_creds = ratio_finemap_creds_
    serializable_FINEMAP_output low_prior_std_finemap = prior_std_low_finemap_output_
    Array[File] low_prior_std_finemap_creds = prior_std_low_finemap_creds_

    File followup_df = followup_finemapping_conditions_df.df
  }
}
