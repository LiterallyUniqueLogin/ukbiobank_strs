version 1.0

import "expanse_files.wdl"
import "../finemapping_wdl/simulations_workflow.wdl"
import "../gwas_wdl/gwas_tasks.wdl"

workflow platelet_count_simulations {
  
  call expanse_files.files

  call gwas_tasks.phenotype_names

  Int platelet_count_idx = phenotype_names.idxs["platelet_count"]

  call simulations_workflow.simulations { input :
    script_dir = ".",
    finemap_command = "finemap",
    finemapping_regions = files.finemapping_regions[platelet_count_idx],
    all_regions_finemapping_df = files.finemapping_first_pass_dfs[platelet_count_idx],
    str_vcfs = files.str_vcfs,
    imputed_snp_bgens = files.imputed_snp_bgens,
    imputed_snp_pfiles = files.imputed_snp_pfiles,
    snp_vars_to_filter_from_finemapping = files.snps_to_filter,
    phenotype_name = "platelet_count",
    phenotype_samples = files.unrelated_samples_for_pheno_for_ethnicity[platelet_count_idx][0],
    transformed_phenotype_data = files.ethnic_to_pheno_to_transformed_phenotype_data[0][platelet_count_idx],
    shared_covars = files.shared_covars,
    true_my_str_gwas = files.str_gwas_results[platelet_count_idx],
    true_plink_snp_gwas = files.snp_gwas_results[platelet_count_idx],
    all_samples_list = files.all_samples_list,
    snp_macs = files.platelet_count_snp_macs,
    flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos
  }

  output {
		File bins = simulations.bins
    Array[File] bin_effects = simulations.bin_effects

#    File causal_vars_and_betas_from_finemap_no_strs = simulations.causal_vars_and_betas_from_finemap_no_strs
#    FINEMAP_output finemap_no_strs_simulation_finemap_results = simulations.finemap_no_strs_simulation_finemap_results
#    SuSiE_output finemap_no_strs_simulation_susie_results = simulations.finemap_no_strs_simulation_susie_results
#
#    File causal_vars_and_betas_from_finemap_snps_strs = simulations.causal_vars_and_betas_from_finemap_snps_strs
#    FINEMAP_output finemap_snps_strs_simulation_finemap_results = simulations.finemap_snps_strs_simulation_finemap_results
#    SuSiE_output finemap_snps_strs_simulation_susie_results = simulations.finemap_snps_strs_simulation_susie_results

#    File causal_vars_and_betas_from_susie_no_strs = simulations.causal_vars_and_betas_from_susie_no_strs
#    FINEMAP_output susie_no_strs_simulation_finemap_results = simulations.susie_no_strs_simulation_finemap_results
#    SuSiE_output susie_no_strs_simulation_susie_results = simulations.susie_no_strs_simulation_susie_results
#
#    File causal_vars_and_betas_from_susie_snps_strs = simulations.causal_vars_and_betas_from_susie_snps_strs
#    FINEMAP_output susie_snps_strs_simulation_finemap_results = simulations.susie_snps_strs_simulation_finemap_results
#    SuSiE_output susie_snps_strs_simulation_susie_results = simulations.susie_snps_strs_simulation_susie_results
#
#    File causal_vars_and_betas_from_random_one_var = simulations.causal_vars_and_betas_from_random_one_var
#    FINEMAP_output one_random_causal_simulation_finemap_results = simulations.one_random_causal_simulation_finemap_results
#    SuSiE_output one_random_causal_simulation_susie_results = simulations.one_random_causal_simulation_susie_results
#
#    File causal_vars_and_betas_from_random_two_var = simulations.causal_vars_and_betas_from_random_two_var
#    FINEMAP_output two_random_causal_simulation_finemap_results = simulations.two_random_causal_simulation_finemap_results
#    SuSiE_output two_random_causal_simulation_susie_results = simulations.two_random_causal_simulation_susie_results
#
#    File causal_vars_and_betas_from_random_three_var = simulations.causal_vars_and_betas_from_random_three_var
#    FINEMAP_output three_random_causal_simulation_finemap_results = simulations.three_random_causal_simulation_finemap_results
#    SuSiE_output three_random_causal_simulation_susie_results = simulations.three_random_causal_simulation_susie_results

    Array[File] causal_vars_and_betas_from_susie_no_strs = simulations.causal_vars_and_betas_from_susie_no_strs
    Array[Array[FINEMAP_output]] susie_no_strs_simulation_finemap_results = simulations.susie_no_strs_simulation_finemap_results
    Array[Array[SuSiE_output]] susie_no_strs_simulation_susie_results = simulations.susie_no_strs_simulation_susie_results

    Array[File] causal_vars_and_betas_from_susie_snps_strs = simulations.causal_vars_and_betas_from_susie_snps_strs
    Array[Array[FINEMAP_output]] susie_snps_strs_simulation_finemap_results = simulations.susie_snps_strs_simulation_finemap_results
    Array[Array[SuSiE_output]] susie_snps_strs_simulation_susie_results = simulations.susie_snps_strs_simulation_susie_results

    Array[Array[File]] causal_vars_and_betas_from_random_one_var = simulations.causal_vars_and_betas_from_random_one_var
    Array[Array[FINEMAP_output]] one_random_causal_simulation_finemap_results = simulations.one_random_causal_simulation_finemap_results
    Array[Array[SuSiE_output]] one_random_causal_simulation_susie_results = simulations.one_random_causal_simulation_susie_results

    Array[Array[File]] causal_vars_and_betas_from_random_two_var = simulations.causal_vars_and_betas_from_random_two_var
    Array[Array[FINEMAP_output]] two_random_causal_simulation_finemap_results = simulations.two_random_causal_simulation_finemap_results
    Array[Array[SuSiE_output]] two_random_causal_simulation_susie_results = simulations.two_random_causal_simulation_susie_results

    Array[Array[File]] causal_vars_and_betas_from_random_three_var = simulations.causal_vars_and_betas_from_random_three_var
    Array[Array[FINEMAP_output]] three_random_causal_simulation_finemap_results = simulations.three_random_causal_simulation_finemap_results
    Array[Array[SuSiE_output]] three_random_causal_simulation_susie_results = simulations.three_random_causal_simulation_susie_results

    File simulations_df = simulations.simulations_df
  }
}
