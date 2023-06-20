version 1.0

import "finemapping.wdl"
import "../gwas_wdl/gwas_tasks.wdl"
import "expanse_files.wdl"
import "expanse_tasks.wdl"

workflow bilirubin_and_ldl_finemapping {

  call gwas_tasks.phenotype_names

  call expanse_files.files

  call finemapping.finemapping as bilirubin_finemapping { input :
    phenotype_name = "total_bilirubin"
  }

  call finemapping.finemapping as ldl_finemapping { input :
    phenotype_name = "ldl_cholesterol_direct"
  }

  call expanse_tasks.rename_file as bilirubin_finemapping_regions { input : f = bilirubin_finemapping.first_pass_regions_tsv, name = "total_bilirubin_finemapping_regions.tsv" }
  call expanse_tasks.rename_file as ldl_finemapping_regions { input : f = ldl_finemapping.first_pass_regions_tsv, name = "ldl_cholesterol_direct_finemapping_regions.tsv" }
  call expanse_tasks.rename_file as bilirubin_finemapping_regions_readme { input : f = bilirubin_finemapping.first_pass_regions_readme, name = "total_bilirubin_finemapping_regions_README.txt" }
  call expanse_tasks.rename_file as ldl_finemapping_regions_readme { input : f = ldl_finemapping.first_pass_regions_readme, name = "ldl_cholesterol_direct_finemapping_regions_README.txt" }

  output {
    File bilirubin_first_pass_regions_tsv = bilirubin_finemapping.first_pass_regions_tsv
    File bilirubin_first_pass_regions_readme = bilirubin_finemapping.first_pass_regions_readme

    Array[serializable_SuSiE_output] bilirubin_original_susie = bilirubin_finemapping.original_susie
    Array[Array[File]] bilirubin_original_susie_CSs = bilirubin_finemapping.original_susie_CSs
    Array[serializable_FINEMAP_output] bilirubin_original_finemap = bilirubin_finemapping.original_finemap
    Array[Array[File]] bilirubin_original_finemap_creds = bilirubin_finemapping.original_finemap_creds

    File bilirubin_first_pass_df = bilirubin_finemapping.first_pass_df
    File bilirubin_susie_min_abs_corrs = bilirubin_finemapping.susie_min_abs_corrs

    File bilirubin_followup_regions_tsv = bilirubin_finemapping.followup_regions_tsv

    Array[serializable_SuSiE_output] bilirubin_best_guess_susie = bilirubin_finemapping.best_guess_susie
    Array[Array[File]] bilirubin_best_guess_susie_CSs = bilirubin_finemapping.best_guess_susie_CSs
    Array[serializable_FINEMAP_output] bilirubin_total_prob_finemap = bilirubin_finemapping.total_prob_finemap
    Array[Array[File]] bilirubin_total_prob_finemap_creds = bilirubin_finemapping.total_prob_finemap_creds
    Array[serializable_FINEMAP_output] bilirubin_derived_prior_std_finemap = bilirubin_finemapping.derived_prior_std_finemap
    Array[Array[File]] bilirubin_derived_prior_std_finemap_creds = bilirubin_finemapping.derived_prior_std_finemap_creds
    Array[serializable_FINEMAP_output] bilirubin_prob_conv_sss_tol_finemap = bilirubin_finemapping.prob_conv_sss_tol_finemap
    Array[Array[File]] bilirubin_prob_conv_sss_tol_finemap_creds = bilirubin_finemapping.prob_conv_sss_tol_finemap_creds
    Array[serializable_FINEMAP_output] bilirubin_mac_finemap = bilirubin_finemapping.mac_finemap
    Array[Array[File]] bilirubin_mac_finemap_creds = bilirubin_finemapping.mac_finemap_creds
    Array[serializable_FINEMAP_output] bilirubin_threshold_finemap = bilirubin_finemapping.threshold_finemap
    Array[Array[File]] bilirubin_threshold_finemap_creds = bilirubin_finemapping.threshold_finemap_creds

    Array[serializable_SuSiE_output] bilirubin_ratio_susie = bilirubin_finemapping.ratio_susie
    Array[Array[File]] bilirubin_ratio_susie_CSs = bilirubin_finemapping.ratio_susie_CSs
    Array[serializable_FINEMAP_output] bilirubin_ratio_finemap = bilirubin_finemapping.ratio_finemap
    Array[Array[File]] bilirubin_ratio_finemap_creds = bilirubin_finemapping.ratio_finemap_creds
    Array[serializable_FINEMAP_output] bilirubin_low_prior_std_finemap = bilirubin_finemapping.low_prior_std_finemap
    Array[Array[File]] bilirubin_low_prior_std_finemap_creds = bilirubin_finemapping.low_prior_std_finemap_creds

    File? bilirubin_followup_df = bilirubin_finemapping.followup_df

    File ldl_first_pass_regions_tsv = ldl_finemapping.first_pass_regions_tsv
    File ldl_first_pass_regions_readme = ldl_finemapping.first_pass_regions_readme

    Array[serializable_SuSiE_output] ldl_original_susie = ldl_finemapping.original_susie
    Array[Array[File]] ldl_original_susie_CSs = ldl_finemapping.original_susie_CSs
    Array[serializable_FINEMAP_output] ldl_original_finemap = ldl_finemapping.original_finemap
    Array[Array[File]] ldl_original_finemap_creds = ldl_finemapping.original_finemap_creds

    File ldl_first_pass_df = ldl_finemapping.first_pass_df
    File ldl_susie_min_abs_corrs = ldl_finemapping.susie_min_abs_corrs

    File ldl_followup_regions_tsv = ldl_finemapping.followup_regions_tsv

    Array[serializable_SuSiE_output] ldl_best_guess_susie = ldl_finemapping.best_guess_susie
    Array[Array[File]] ldl_best_guess_susie_CSs = ldl_finemapping.best_guess_susie_CSs
    Array[serializable_FINEMAP_output] ldl_total_prob_finemap = ldl_finemapping.total_prob_finemap
    Array[Array[File]] ldl_total_prob_finemap_creds = ldl_finemapping.total_prob_finemap_creds
    Array[serializable_FINEMAP_output] ldl_derived_prior_std_finemap = ldl_finemapping.derived_prior_std_finemap
    Array[Array[File]] ldl_derived_prior_std_finemap_creds = ldl_finemapping.derived_prior_std_finemap_creds
    Array[serializable_FINEMAP_output] ldl_prob_conv_sss_tol_finemap = ldl_finemapping.prob_conv_sss_tol_finemap
    Array[Array[File]] ldl_prob_conv_sss_tol_finemap_creds = ldl_finemapping.prob_conv_sss_tol_finemap_creds
    Array[serializable_FINEMAP_output] ldl_mac_finemap = ldl_finemapping.mac_finemap
    Array[Array[File]] ldl_mac_finemap_creds = ldl_finemapping.mac_finemap_creds
    Array[serializable_FINEMAP_output] ldl_threshold_finemap = ldl_finemapping.threshold_finemap
    Array[Array[File]] ldl_threshold_finemap_creds = ldl_finemapping.threshold_finemap_creds

    Array[serializable_SuSiE_output] ldl_ratio_susie = ldl_finemapping.ratio_susie
    Array[Array[File]] ldl_ratio_susie_CSs = ldl_finemapping.ratio_susie_CSs
    Array[serializable_FINEMAP_output] ldl_ratio_finemap = ldl_finemapping.ratio_finemap
    Array[Array[File]] ldl_ratio_finemap_creds = ldl_finemapping.ratio_finemap_creds
    Array[serializable_FINEMAP_output] ldl_low_prior_std_finemap = ldl_finemapping.low_prior_std_finemap
    Array[Array[File]] ldl_low_prior_std_finemap_creds = ldl_finemapping.low_prior_std_finemap_creds

    File? ldl_followup_df = ldl_finemapping.followup_df
  }
}
