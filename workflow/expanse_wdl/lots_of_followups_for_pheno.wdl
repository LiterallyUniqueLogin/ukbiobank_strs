version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "expanse_files.wdl"
import "expanse_finemapping_files.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"
import "../finemapping_wdl/finemap_one_region_workflow.wdl"
import "../finemapping_wdl/susie_one_region_workflow.wdl"

workflow lots_of_followups_for_pheno {
  input {
    String phenotype = "mean_platelet_volume"
  }

  call gwas_tasks.phenotype_names
  call expanse_files.files
  call expanse_finemapping_files.finemapping_files

  Int phenotype_idx = phenotype_names.idxs[phenotype]
  Int ldl_idx = phenotype_names.idxs["ldl_cholesterol_direct"]
  Int bilirubin_idx = phenotype_names.idxs["total_bilirubin"]
  if (phenotype_idx > bilirubin_idx) {
    Int phenotype_idx_minus = phenotype_idx - 1
  }
  Int phenotype_idx_temp = select_first([phenotype_idx_minus, phenotype_idx])
  if (phenotype_idx_temp > ldl_idx) {
    Int phenotype_idx_temp_minus = phenotype_idx_temp - 1
  }
  Int phenotype_idx_sans = select_first([phenotype_idx_temp_minus, phenotype_idx_temp])

  File finemapping_region = files.finemapping_regions[phenotype_idx]

  Array[Array[String]] finemapping_regions_tsv = read_tsv(finemapping_region)

  # first passs finemap results 
  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
    region first_pass_bounds = {
      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
    }
    if (!(
      phenotype == "urate" &&
      first_pass_bounds.chrom == 4 &&
      first_pass_bounds.start == 8165642 &&
      first_pass_bounds.end == 11717761
    ) && !(
      phenotype == "total_bilirubin" &&
      first_pass_bounds.chrom == 12 &&
      first_pass_bounds.start == 19976272 &&
      first_pass_bounds.end == 22524428
    ) && !(
      phenotype == "alkaline_phosphatase" &&
      first_pass_bounds.chrom == 1 &&
      first_pass_bounds.start == 19430673 &&
      first_pass_bounds.end == 24309348
    )) {
      String first_pass_regions = "~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
      Int first_pass_chroms = first_pass_bounds.chrom

      call finemap_one_region_workflow.finemap_one_region as original_finemap_ { input :
        script_dir = ".",
        finemap_command = "finemap",
        str_vcfs = files.str_vcfs,
        imputed_snp_bgens = files.imputed_snp_bgens,
        snp_vars_to_filter_from_finemapping = files.snps_to_filter,
        phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
        my_str_gwas = files.str_gwas_results[phenotype_idx],
        plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
        phenotype_name = phenotype,
        bounds = first_pass_bounds,
        all_samples_list = files.all_samples_list,
        prefix = "~{phenotype}_FINEMAP_first_pass_~{first_pass_regions}_"
      }
      serializable_FINEMAP_output original_finemap_output_ = original_finemap_.finemap_output.subset
      Array[File] original_finemap_creds_ = original_finemap_.finemap_output.creds
    }
  }

  call finemapping_tasks.first_pass_finemapping_df { input :
    script_dir = ".",
    phenotype_name = phenotype,
    snp_assoc_results = files.snp_gwas_results[phenotype_idx],
    str_assoc_results = files.str_gwas_results[phenotype_idx],
    ethnic_str_assoc_results = files.pheno_to_ethnic_to_str_gwas_results[phenotype_idx],
    original_finemap_outputs = select_all(original_finemap_output_),
    original_finemap_creds = select_all(original_finemap_creds_),
    original_susie_outputs = finemapping_files.original_susie_snakemake[phenotype_idx_sans],
    original_susie_CSs = finemapping_files.original_susie_CSs_snakemake[phenotype_idx_sans],
    regions = select_all(first_pass_regions),
    chroms = select_all(first_pass_chroms)
  }

  call finemapping_tasks.generate_followup_regions_tsv { input :
    script_dir = ".",
    phenotype = phenotype,
    first_pass_df = first_pass_finemapping_df.all_regions_concordance
  }

  # followup finemapping on all loci
  scatter (followup_bounds in first_pass_bounds) {
    String followup_regions = "~{followup_bounds.chrom}_~{followup_bounds.start}_~{followup_bounds.end}"
    Int followup_chroms = followup_bounds.chrom

    call finemap_one_region_workflow.finemap_one_region as repeat_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      prefix = "~{phenotype}_FINEMAP_repeat_~{followup_regions}_",
      cache_breaker = 1
    }
    serializable_FINEMAP_output repeat_finemap_output_ = repeat_finemap_.finemap_output.subset
    Array[File] repeat_finemap_creds_ = repeat_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as ratio_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      snp_str_ratio = 4,
      prefix = "~{phenotype}_FINEMAP_prior_snps_over_strs_~{followup_regions}_"
    }
    serializable_FINEMAP_output ratio_finemap_output_ = ratio_finemap_.finemap_output.subset
    Array[File] ratio_finemap_creds_ = ratio_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as total_prob_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      total_prob = 4,
      prefix = "~{phenotype}_FINEMAP_prior_4_signals_~{followup_regions}_"
    }
    serializable_FINEMAP_output total_prob_finemap_output_ = total_prob_finemap_.finemap_output.subset
    Array[File] total_prob_finemap_creds_ = total_prob_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as prior_std_derived_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      prior_std = 0.0224,
      prefix = "~{phenotype}_FINEMAP_derived_effect_size_prior_~{followup_regions}_"
    }
    serializable_FINEMAP_output prior_std_derived_finemap_output_ = prior_std_derived_finemap_.finemap_output.subset
    Array[File] prior_std_derived_finemap_creds_ = prior_std_derived_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as prior_std_low_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      prior_std = 0.005,
      prefix = "~{phenotype}_FINEMAP_low_effect_size_prior_~{followup_regions}_"
    }
    serializable_FINEMAP_output prior_std_low_finemap_output_ = prior_std_low_finemap_.finemap_output.subset
    Array[File] prior_std_low_finemap_creds_ = prior_std_low_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as prob_conv_sss_tol_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      prob_conv_sss_tol = 0.0001,
      prefix = "~{phenotype}_FINEMAP_stricter_stopping_threshold_~{followup_regions}_"
    }
    serializable_FINEMAP_output prob_conv_sss_tol_finemap_output_ = prob_conv_sss_tol_finemap_.finemap_output.subset
    Array[File] prob_conv_sss_tol_finemap_creds_ = prob_conv_sss_tol_finemap_.finemap_output.creds

    Int followup_chroms_minus_one = followup_chroms - 1
    call gwas_tasks.imputed_snp_frequencies { input :
      script_dir = ".",
      imputed_snp_bgen = files.imputed_snp_bgens[followup_chroms_minus_one],
      samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      all_samples = files.all_samples_list,
      bounds = followup_bounds,
    }

    call finemap_one_region_workflow.finemap_one_region as mac_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      mac = 100,
      snp_macs = imputed_snp_frequencies.counts,
      prefix = "~{phenotype}_FINEMAP_mac_threshold_100_~{followup_regions}_"
    }
    serializable_FINEMAP_output mac_finemap_output_ = mac_finemap_.finemap_output.subset
    Array[File] mac_finemap_creds_ = mac_finemap_.finemap_output.creds

    call finemap_one_region_workflow.finemap_one_region as threshold_finemap_ { input :
      script_dir = ".",
      finemap_command = "finemap",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      inclusion_threshold = 0.0005,
      prefix = "~{phenotype}_FINEMAP_pval_threshold_1e4_~{followup_regions}_"
    }
    serializable_FINEMAP_output threshold_finemap_output_ = threshold_finemap_.finemap_output.subset
    Array[File] threshold_finemap_creds_ = threshold_finemap_.finemap_output.creds

    call susie_one_region_workflow.susie_one_region as best_guess_susie_ { input :
      script_dir = ".",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      shared_covars = files.shared_covars,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][phenotype_idx],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      best_guess = true,
      prefix = "~{phenotype}_SuSiE_best_guess_genotypes_~{followup_regions}_"
    }
    serializable_SuSiE_output best_guess_susie_output_ = best_guess_susie_.susie_output.subset
    Array[File] best_guess_susie_CSs_ = best_guess_susie_.susie_output.CSs

    call susie_one_region_workflow.susie_one_region as ratio_susie_ { input :
      script_dir = ".",
      str_vcfs = files.str_vcfs,
      imputed_snp_bgens = files.imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = files.snps_to_filter,
      shared_covars = files.shared_covars,
      phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[phenotype_idx][0],
      transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][phenotype_idx],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      plink_snp_gwas = files.snp_gwas_results[phenotype_idx],
      phenotype_name = phenotype,
      bounds = followup_bounds,
      all_samples_list = files.all_samples_list,
      snp_p_over_str_p = 4,
      prefix = "~{phenotype}_SuSiE_prior_snps_over_strs_~{followup_regions}_"
    }
    serializable_SuSiE_output ratio_susie_output_ = ratio_susie_.susie_output.subset
    Array[File] ratio_susie_CSs_ = ratio_susie_.susie_output.CSs
  }

  call finemapping_tasks.followup_finemapping_conditions_df as followup_on_all_regions_df { input :
    script_dir = ".",
    phenotype_name = phenotype,
    snp_assoc_results = files.snp_gwas_results[phenotype_idx],
    str_assoc_results = files.str_gwas_results[phenotype_idx],
    ethnic_str_assoc_results = files.pheno_to_ethnic_to_str_gwas_results[phenotype_idx],
    original_susie_outputs = finemapping_files.original_susie_snakemake[phenotype_idx_sans],
    original_susie_CSs = finemapping_files.original_susie_CSs_snakemake[phenotype_idx_sans],
    original_finemap_outputs = select_all(original_finemap_output_),
    original_finemap_creds = select_all(original_finemap_creds_),
    best_guess_susie_outputs = best_guess_susie_output_,
    best_guess_susie_CSs = best_guess_susie_CSs_,
    repeat_finemap_outputs = repeat_finemap_output_,
    repeat_finemap_creds = repeat_finemap_creds_,
    total_prob_finemap_outputs = total_prob_finemap_output_,
    total_prob_finemap_creds = total_prob_finemap_creds_,
    derived_prior_std_finemap_outputs = prior_std_derived_finemap_output_,
    derived_prior_std_finemap_creds = prior_std_derived_finemap_creds_,
    conv_tol_finemap_outputs = prob_conv_sss_tol_finemap_output_,
    conv_tol_finemap_creds = prob_conv_sss_tol_finemap_creds_,
    mac_finemap_outputs = mac_finemap_output_,
    mac_finemap_creds = mac_finemap_creds_,
    threshold_finemap_outputs = threshold_finemap_output_,
    threshold_finemap_creds = threshold_finemap_creds_,
    ratio_susie_outputs = ratio_susie_output_,
    ratio_susie_CSs = ratio_susie_CSs_,
    ratio_finemap_outputs = ratio_finemap_output_,
    ratio_finemap_creds = ratio_finemap_creds_,
    low_prior_std_finemap_outputs = prior_std_low_finemap_output_,
    low_prior_std_finemap_creds = prior_std_low_finemap_creds_,
    original_regions = select_all(first_pass_regions),
    original_chroms = select_all(first_pass_chroms),
    followup_regions = followup_regions,
    followup_chroms = followup_chroms
  }

  call finemapping_tasks.followup_finemapping_conditions_comparison { input :
    script_dir = ".",
    followup_conditions_tsvs = [followup_on_all_regions_df.df]
  }

  output {
    Array[serializable_FINEMAP_output] original_finemap = select_all(original_finemap_output_)
    Array[Array[File]] original_finemap_creds = select_all(original_finemap_creds_)

    File first_pass_df = first_pass_finemapping_df.all_regions_concordance
    File susie_min_abs_corrs = first_pass_finemapping_df.susie_all_regions_min_abs_corrs

    File followup_regions_tsv = generate_followup_regions_tsv.tsv

    Array[serializable_SuSiE_output] best_guess_susie = best_guess_susie_output_
    Array[Array[File]] best_guess_susie_CSs = best_guess_susie_CSs_
    Array[serializable_FINEMAP_output] repeat_finemap = repeat_finemap_output_
    Array[Array[File]] repeat_finemap_creds = repeat_finemap_creds_
    Array[serializable_FINEMAP_output] total_prob_finemap = total_prob_finemap_output_
    Array[Array[File]] total_prob_finemap_creds = total_prob_finemap_creds_
    Array[serializable_FINEMAP_output] derived_prior_std_finemap = prior_std_derived_finemap_output_
    Array[Array[File]] derived_prior_std_finemap_creds = prior_std_derived_finemap_creds_
    Array[serializable_FINEMAP_output] prob_conv_sss_tol_finemap = prob_conv_sss_tol_finemap_output_
    Array[Array[File]] prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_creds_
    Array[serializable_FINEMAP_output] mac_finemap = mac_finemap_output_
    Array[Array[File]] mac_finemap_creds = mac_finemap_creds_
    Array[serializable_FINEMAP_output] threshold_finemap = threshold_finemap_output_
    Array[Array[File]] threshold_finemap_creds = threshold_finemap_creds_

    Array[serializable_SuSiE_output] ratio_susie = ratio_susie_output_
    Array[Array[File]] ratio_susie_CSs = ratio_susie_CSs_
    Array[serializable_FINEMAP_output] ratio_finemap = ratio_finemap_output_
    Array[Array[File]] ratio_finemap_creds = ratio_finemap_creds_
    Array[serializable_FINEMAP_output] low_prior_std_finemap = prior_std_low_finemap_output_
    Array[Array[File]] low_prior_std_finemap_creds = prior_std_low_finemap_creds_

    File followup_df = followup_on_all_regions_df.df
    
    File susie_best_guess_png = followup_finemapping_conditions_comparison.susie_best_guess_png
    File susie_best_guess_svg = followup_finemapping_conditions_comparison.susie_best_guess_svg
    File finemap_conv_tol_png = followup_finemapping_conditions_comparison.finemap_conv_tol_png
    File finemap_conv_tol_svg = followup_finemapping_conditions_comparison.finemap_conv_tol_svg
    File finemap_total_prob_png = followup_finemapping_conditions_comparison.finemap_total_prob_png
    File finemap_total_prob_svg = followup_finemapping_conditions_comparison.finemap_total_prob_svg
    File finemap_prior_std_derived_png = followup_finemapping_conditions_comparison.finemap_prior_std_derived_png
    File finemap_prior_std_derived_svg = followup_finemapping_conditions_comparison.finemap_prior_std_derived_svg
    File finemap_mac_png = followup_finemapping_conditions_comparison.finemap_mac_png
    File finemap_mac_svg = followup_finemapping_conditions_comparison.finemap_mac_svg
    File finemap_p_thresh_png = followup_finemapping_conditions_comparison.finemap_p_thresh_png
    File finemap_p_thresh_svg = followup_finemapping_conditions_comparison.finemap_p_thresh_svg

    # too conservative
    File finemap_ratio_png = followup_finemapping_conditions_comparison.finemap_ratio_png
    File finemap_ratio_svg = followup_finemapping_conditions_comparison.finemap_ratio_svg
    File susie_ratio_png = followup_finemapping_conditions_comparison.susie_ratio_png
    File susie_ratio_svg = followup_finemapping_conditions_comparison.susie_ratio_svg
    File finemap_prior_std_low_png = followup_finemapping_conditions_comparison.finemap_prior_std_low_png
    File finemap_prior_std_low_svg = followup_finemapping_conditions_comparison.finemap_prior_std_low_svg
  }

}
