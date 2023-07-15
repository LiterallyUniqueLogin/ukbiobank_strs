version 1.0

import "finemapping_tasks.wdl"

workflow post_finemapping {
  input {
    String script_dir

    Array[File] first_pass_dfs
    Array[File] susie_all_min_abs_corrs
    Array[File] followup_dfs
  }

  call finemapping_tasks.first_pass_comparison { input :
    script_dir = script_dir,
    first_pass_dfs = first_pass_dfs,
    susie_all_min_abs_corrs = susie_all_min_abs_corrs
  }
 
  call finemapping_tasks.followup_finemapping_conditions_comparison { input :
    script_dir = script_dir,
    followup_conditions_tsvs = followup_dfs
  }

  output {
    File stat_statements = first_pass_comparison.stat_statements
    File cs_min_abs_corrs_png = first_pass_comparison.cs_min_abs_corrs_png
    File cs_min_abs_corrs_svg = first_pass_comparison.cs_min_abs_corrs_svg
    File susie_alpha_v_pip_png = first_pass_comparison.susie_alpha_v_pip_png
    File susie_alpha_v_pip_svg = first_pass_comparison.susie_alpha_v_pip_svg
    File susie_alpha_histogram_png = first_pass_comparison.susie_alpha_histogram_png
    File susie_alpha_histogram_svg = first_pass_comparison.susie_alpha_histogram_svg
    File finemap_pip_histogram_png = first_pass_comparison.finemap_pip_histogram_png
    File finemap_pip_histogram_svg = first_pass_comparison.finemap_pip_histogram_svg
    File susie_cs_finemap_total_pips_png = first_pass_comparison.susie_cs_finemap_total_pips_png
    File susie_cs_finemap_total_pips_svg = first_pass_comparison.susie_cs_finemap_total_pips_svg
    File finemap_v_susie_consistency_STR_png = first_pass_comparison.finemap_v_susie_consistency_STR_png
    File finemap_v_susie_consistency_STR_svg = first_pass_comparison.finemap_v_susie_consistency_STR_svg
    File finemap_v_susie_consistency_SNP_png = first_pass_comparison.finemap_v_susie_consistency_SNP_png
    File finemap_v_susie_consistency_SNP_svg = first_pass_comparison.finemap_v_susie_consistency_SNP_svg

    File doubly_finemapped_STRs = followup_finemapping_conditions_comparison.doubly_finemapped_STRs
    File confidently_finemapped_STRs = followup_finemapping_conditions_comparison.confidently_finemapped_STRs
    File overconfidently_finemapped_STRs = followup_finemapping_conditions_comparison.overconfidently_finemapped_STRs

    # used for deciding confidently finemap
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
