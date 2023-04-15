version 1.0

import "../gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"

workflow finemapping {

  input {
    String script_dir
    String finemap_command

    File chr_lens
    File str_loci

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    File shared_covars
    File phenotype_samples
    File transformed_phenotype_data

    File my_str_gwas
    Array[File] ethnic_my_str_gwass,
    File plink_snp_gwas

    String phenotype_name
    Boolean is_binary

    File all_samples_list
  }

  call gwas_tasks.generate_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_regions.data)

  # finemap each region
  scatter (region_idx in range(length(finemapping_regions_tsv) - 1)) {
    region first_pass_bounds = {
      "chrom": finemapping_regions_tsv[region_idx+1][0],
      "start": finemapping_regions_tsv[region_idx+1][1],
      "end": finemapping_regions_tsv[region_idx+1][2],
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
  }

  call finemapping_tasks.first_pass_finemapping_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    original_finemap_outputs = original_finemap_.finemap_output,
    original_susie_outputs = original_susie_.susie_output
  }

  call finemapping_tasks.generate_followup_regions_tsv { input :
    script_dir = script_dir,
    first_pass_df = first_pass_finemapping_df.all_regions_concordance
  }

  Array[Array[String]] followup_finemapping_regions_tsv = read_tsv(generate_followup_regions_tsv.tsv)

  # finemap each region
  scatter (region_idx in range(length(followup_finemapping_regions_tsv) - 1)) {
    region followup_bounds = {
      "chrom": finemapping_regions_tsv[region_idx+1][0],
      "start": finemapping_regions_tsv[region_idx+1][1],
      "end": finemapping_regions_tsv[region_idx+1][2],
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_ratio_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_total_prob_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_prior_std_derived_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_prior_std_low_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_prob_conv_sss_tol_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_mac_ { input :
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

    call finemap_one_region_workflow.finemap_one_region as finemap_threshold_ { input :
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

    call susie_one_region_workflow.susie_one_region as susie_best_guess_ { input :
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

    call susie_one_region_workflow.susie_one_region as susie_ratio_ { input :
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
  }

  call finemapping_tasks.followup_finemapping_conditions_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    original_finemap_outputs = original_finemap.finemap_output,
    original_susie_outputs = original_susie.susie_output,
    total_prob_finemap_outputs = finemap_total_prob_.finemap_output,
    derived_prior_std_finemap_outputs = finemap_prior_std_derived_.finemap_output,
    conv_tol_finemap_outputs = finemap_conv_tol_.finemap_output,
    mac_finemap_outputs = finemap_mac_.finemap_output,
    threshold_finemap_outputs = finemap_threshold_.finemap_output,
    best_guess_susie_outputs = susie_best_guess_.susie_output,
    low_prior_std_finemap_outputs = finemap_prior_std_low_.finemap_output,
    ratio_finemap_outputs = finemap_ratio.finemap_output,
    ratio_susie_outputs = susie_ratio.susie_output,
  }

  output {
    File first_pass_regions_tsv = generate_regions.data
    File first_pass_regions_readme = generate_regions.readme
    Array[SuSiE_output] original_susie = original_susie_.susie_output
    Array[FINEMAP_output] original_finemap = original_finemap_.finemap_output
    File first_pass_df = first_pass_finemapping_df.all_regions_concordance
    File susie_min_abs_corrs = first_pass_finemapping_df.susie_all_regions_min_abs_corrs
    File followup_regions_tsv = generate_followup_regions_tsv.tsv
    Array[SuSiE_output] best_guess_susie = susie_best_guess_.susie_output
    Array[FINEMAP_output] total_prob_finemap = finemap_total_prob_.finemap_output
    Array[FINEMAP_output] derived_prior_std_finemap = finemap_prior_std_dervied_.finemap_output
    Array[FINEMAP_output] prob_conv_sss_tol_finemap = finemap_prob_conv_sss_tol_.finemap_output
    Array[FINEMAP_output] mac_finemap = finemap_mac_.finemap_output
    Array[FINEMAP_output] threshold_finemap = finemap_threshold_.finemap_output
    Array[SuSiE_output] ratio_susie = susie_ratio_.susie_output
    Array[FINEMAP_output] ratio_finemap = finemap_ratio_.finemap_output
    Array[FINEMAP_output] low_prior_std_finemap = finemap_prior_std_low_.finemap_output
    File followup_df = followup_finemapping_conditions_df.df
  }
}
