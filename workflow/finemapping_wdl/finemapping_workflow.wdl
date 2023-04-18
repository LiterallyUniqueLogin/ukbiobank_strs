version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"
import "susie_one_region_workflow.wdl"

workflow finemapping {

  input {
    String script_dir
    String finemap_command

    File chr_lens

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    File shared_covars
    File phenotype_samples
    File transformed_phenotype_data

    File my_str_gwas
    Array[File] ethnic_my_str_gwass
    File plink_snp_gwas

    String phenotype_name
    Boolean is_binary

    File all_samples_list
  }

  call gwas_tasks.generate_finemapping_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_finemapping_regions.data)

  # finemap each region
  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
    region first_pass_bounds = {
      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
    }

    call finemap_one_region_workflow.finemap_one_region as original_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = first_pass_bounds,
      all_samples_list = all_samples_list
    }
    serializable_FINEMAP_output original_finemap_output = original_finemap_.output.subset
    Array[File] original_finemap_creds = original_finemap_.output.creds

    call susie_one_region_workflow.susie_one_region as original_susie_ { input :
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
      bounds = first_pass_bounds,
      all_samples_list = all_samples_list
    }
    serializable_SuSiE_output original_susie_output = original_susie_.output.subset
    Array[File] original_susie_CSs = original_susie_.output.CSs
  }

  call finemapping_tasks.first_pass_finemapping_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    original_finemap_outputs = original_finemap_.finemap_output.subset,
    original_finemap_creds = original_finemap_.finemap_output.creds,
    original_susie_outputs = original_susie_.susie_output.subset,
    original_susie_CSs = original_susie_.susie_output.CSs
  }

  call finemapping_tasks.generate_followup_regions_tsv { input :
    script_dir = script_dir,
    first_pass_df = first_pass_finemapping_df.all_regions_concordance
  }

  Array[Array[String]] followup_finemapping_regions_tsv = read_tsv(generate_followup_regions_tsv.tsv)

  # finemap each region
  scatter (followup_region_idx in range(length(followup_finemapping_regions_tsv) - 1)) {
    Int followup_region_idx_plus_one = followup_region_idx + 1
    region followup_bounds = {
      "chrom": followup_finemapping_regions_tsv[followup_region_idx_plus_one][1],
      "start": sub(sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one][2], "^[^_]*_", ""), "_[^_]*$", ""),
      "end": sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one][2], "^[^_]*_[^_]*_", ""),
    }

    call finemap_one_region_workflow.finemap_one_region as ratio_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      snp_str_ratio = 4
    }
    serializable_FINEMAP_output ratio_finemap_output = ratio_finemap_.output.subset
    Array[File] ratio_finemap_creds = ratio_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as total_prob_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      total_prob = 4,
    }
    serializable_FINEMAP_output total_prob_finemap_output = total_prob_finemap_.output.subset
    Array[File] total_prob_finemap_creds = total_prob_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as prior_std_derived_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prior_std = 0.0224
    }
    serializable_FINEMAP_output prior_std_derived_finemap_output = prior_std_derived_finemap_.output.subset
    Array[File] prior_std_derived_finemap_creds = prior_std_derived_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as prior_std_low_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prior_std = 0.005
    }
    serializable_FINEMAP_output prior_std_low_finemap_output = prior_std_low_finemap_.output.subset
    Array[File] prior_std_low_finemap_creds = prior_std_low_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as prob_conv_sss_tol_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prob_conv_sss_tol = 0.0001
    }
    serializable_FINEMAP_output prob_conv_sss_tol_finemap_output = prob_conv_sss_tol_finemap_.output.subset
    Array[File] prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as mac_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      mac = 100
    }
    serializable_FINEMAP_output mac_finemap_output = mac_finemap_.output.subset
    Array[File] mac_finemap_creds = mac_finemap_.output.creds

    call finemap_one_region_workflow.finemap_one_region as threshold_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      threshold = 0.0005
    }
    serializable_FINEMAP_output threshold_finemap_output = threshold_finemap_.output.subset
    Array[File] threshold_finemap_creds = threshold_finemap_.output.creds

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
    }
    serializable_SuSiE_output best_guess_susie_output = best_guess_susie_.output.subset
    Array[File] best_guess_susie_CSs = best_guess_susie_.output.CSs

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
      ratio = 4
    }
    serializable_SuSiE_output ratio_susie_output = ratio_susie_.output.subset
    Array[File] ratio_susie_CSs = ratio_susie_.output.CSs
  }

  call finemapping_tasks.followup_finemapping_conditions_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    best_guess_susie_outputs = best_guess_susie_.susie_output.subset,
    best_guess_susie_CSs = best_guess_susie_.susie_output.CSs,
    total_prob_finemap_outputs = total_prob_finemap_.finemap_output.subset,
    total_prob_finemap_creds = total_prob_finemap_.finemap_output.creds,
    derived_prior_std_finemap_outputs = prior_std_derived_finemap_.finemap_output.subset,
    derived_prior_std_finemap_creds = prior_std_derived_finemap_.finemap_output.creds,
    prob_conv_sss_tol_finemap_outputs = prob_conv_sss_tol_finemap_.finemap_output.subset,
    prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_.finemap_output.creds,
    mac_finemap_outputs = mac_finemap_.finemap_output.subset,
    mac_finemap_creds = mac_finemap_.finemap_output.creds,
    threshold_finemap_outputs = threshold_finemap_.finemap_output.subset,
    threshold_finemap_creds = threshold_finemap_.finemap_output.creds,
    ratio_susie_outputs = ratio_susie_.susie_output.subset,
    ratio_susie_CSs = ratio_susie_.susie_output.CSs,
    ratio_finemap_outputs = ratio_finemap_.finemap_output.subset,
    ratio_finemap_creds = ratio_finemap_.finemap_output.creds,
    low_prior_std_finemap_outputs = prior_std_low_finemap_.finemap_output.subset,
    low_prior_std_finemap_creds = prior_std_low_finemap_.finemap_output.creds,
  }

  output {
    File first_pass_regions_tsv = generate_finemapping_regions.data
    File first_pass_regions_readme = generate_finemapping_regions.readme

    Array[serializable_SuSiE_output] original_susie = original_susie_.susie_output.subset
    Array[Array[File]] original_susie_CSs = original_susie_.susie_output.CSs
    Array[serializable_FINEMAP_output] original_finemap = original_finemap_.finemap_output.subset
    Array[Array[File]] original_finemap_creds = original_finemap_.finemap_output.creds

    File first_pass_df = first_pass_finemapping_df.all_regions_concordance
    File susie_min_abs_corrs = first_pass_finemapping_df.susie_all_regions_min_abs_corrs

    File followup_regions_tsv = generate_followup_regions_tsv.tsv

    Array[serializable_SuSiE_output] best_guess_susie = best_guess_susie_.susie_output.subset
    Array[Array[File]] best_guess_susie_CSs = best_guess_susie_.susie_output.CSs
    Array[serializable_FINEMAP_output] total_prob_finemap = total_prob_finemap_.finemap_output.subset
    Array[Array[File]] total_prob_finemap_creds = total_prob_finemap_.finemap_output.creds
    Array[serializable_FINEMAP_output] derived_prior_std_finemap = prior_std_derived_finemap_.finemap_output.subset
    Array[Array[File]] derived_prior_std_finemap_creds = prior_std_derived_finemap_.finemap_output.creds
    Array[serializable_FINEMAP_output] prob_conv_sss_tol_finemap = prob_conv_sss_tol_finemap_.finemap_output.subset
    Array[Array[File]] prob_conv_sss_tol_finemap_creds = prob_conv_sss_tol_finemap_.finemap_output.creds
    Array[serializable_FINEMAP_output] mac_finemap = mac_finemap_.finemap_output.subset
    Array[Array[File]] mac_finemap_creds = mac_finemap_.finemap_output.creds
    Array[serializable_FINEMAP_output] threshold_finemap = threshold_finemap_.finemap_output.subset
    Array[Array[File]] threshold_finemap_creds = threshold_finemap_.finemap_output.creds

    Array[serializable_SuSiE_output] ratio_susie = ratio_susie_.susie_output.subset
    Array[Array[File]] ratio_susie_CSs = ratio_susie_.susie_output.CSs
    Array[serializable_FINEMAP_output] ratio_finemap = ratio_finemap_.finemap_output.subset
    Array[Array[File]] ratio_finemap_creds = ratio_finemap_.finemap_output.creds
    Array[serializable_FINEMAP_output] low_prior_std_finemap = prior_std_low_finemap_.finemap_output.subset
    Array[Array[File]] low_prior_std_finemap_creds = prior_std_low_finemap_.finemap_output.creds

    File followup_df = followup_finemapping_conditions_df.df
  }
}
