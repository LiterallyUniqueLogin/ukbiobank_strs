version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/post_finemapping_workflow.wdl"

workflow expanse_post_finemapping {

  call gwas_tasks.phenotype_names

  scatter (phenotype in phenotype_names.n) {
    File first_pass_dfs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/finemapping_all_concordance_~{phenotype}.tab"
    File susie_all_min_abs_corrs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/susie_all_min_abs_corrs_~{phenotype}.npy"
    File followup_dfs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_~{phenotype}.tab"
  }

  call post_finemapping_workflow.post_finemapping { input :
    script_dir = ".",
    first_pass_dfs = first_pass_dfs,
    susie_all_min_abs_corrs = susie_all_min_abs_corrs,
    followup_dfs = followup_dfs
  }

  output {
    File stat_statements = post_finemapping.stat_statements
    File cs_min_abs_corrs_png = post_finemapping.cs_min_abs_corrs_png
    File cs_min_abs_corrs_svg = post_finemapping.cs_min_abs_corrs_svg
    File susie_alpha_v_pip_png = post_finemapping.susie_alpha_v_pip_png
    File susie_alpha_v_pip_svg = post_finemapping.susie_alpha_v_pip_svg
    File susie_alpha_histogram_png = post_finemapping.susie_alpha_histogram_png
    File susie_alpha_histogram_svg = post_finemapping.susie_alpha_histogram_svg
    File finemap_pip_histogram_png = post_finemapping.finemap_pip_histogram_png
    File finemap_pip_histogram_svg = post_finemapping.finemap_pip_histogram_svg
    File susie_cs_finemap_total_pips_png = post_finemapping.susie_cs_finemap_total_pips_png
    File susie_cs_finemap_total_pips_svg = post_finemapping.susie_cs_finemap_total_pips_svg
    File finemap_v_susie_consistency_STR_png = post_finemapping.finemap_v_susie_consistency_STR_png
    File finemap_v_susie_consistency_STR_svg = post_finemapping.finemap_v_susie_consistency_STR_svg
    File finemap_v_susie_consistency_SNP_png = post_finemapping.finemap_v_susie_consistency_SNP_png
    File finemap_v_susie_consistency_SNP_svg = post_finemapping.finemap_v_susie_consistency_SNP_svg

    File doubly_finemapped_STRs = post_finemapping.doubly_finemapped_STRs
    File confidently_finemapped_STRs = post_finemapping.confidently_finemapped_STRs
    File overconfidently_finemapped_STRs = post_finemapping.overconfidently_finemapped_STRs

    # used for deciding confidently finemap
    File susie_best_guess_png = post_finemapping.susie_best_guess_png
    File susie_best_guess_svg = post_finemapping.susie_best_guess_svg
    File finemap_conv_tol_png = post_finemapping.finemap_conv_tol_png
    File finemap_conv_tol_svg = post_finemapping.finemap_conv_tol_svg
    File finemap_total_prob_png = post_finemapping.finemap_total_prob_png
    File finemap_total_prob_svg = post_finemapping.finemap_total_prob_svg
    File finemap_prior_std_derived_png = post_finemapping.finemap_prior_std_derived_png
    File finemap_prior_std_derived_svg = post_finemapping.finemap_prior_std_derived_svg
    File finemap_mac_png = post_finemapping.finemap_mac_png
    File finemap_mac_svg = post_finemapping.finemap_mac_svg
    File finemap_p_thresh_png = post_finemapping.finemap_p_thresh_png
    File finemap_p_thresh_svg = post_finemapping.finemap_p_thresh_svg

    # too conservative
    File finemap_ratio_png = post_finemapping.finemap_ratio_png
    File finemap_ratio_svg = post_finemapping.finemap_ratio_svg
    File susie_ratio_png = post_finemapping.susie_ratio_png
    File susie_ratio_svg = post_finemapping.susie_ratio_svg
    File finemap_prior_std_low_png = post_finemapping.finemap_prior_std_low_png
    File finemap_prior_std_low_svg = post_finemapping.finemap_prior_std_low_svg
  }
}
