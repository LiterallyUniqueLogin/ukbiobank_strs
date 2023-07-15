version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "rerun_finemap_and_followup_for_pheno.wdl"

workflow rerun_first_half {
  call gwas_tasks.phenotype_names

  Array[String] first_half_phenos = [
    "alanine_aminotransferase",
    "albumin",
    #"alkaline_phosphatase",
    "apolipoprotein_a",
    "apolipoprotein_b",
    "aspartate_aminotransferase",
    "c_reactive_protein",
    "calcium",
    "cholesterol",
    "creatinine",
    #"cystatin_c",
    "eosinophil_count",
    "eosinophil_percent",
    #"gamma_glutamyltransferase",
    "glucose",
    "glycated_haemoglobin",
    "haematocrit",
    "haemoglobin_concentration",
    "hdl_cholesterol",
    "igf_1",
    #"ldl_cholesterol_direct",
    #"lymphocyte_count"
  ]

  scatter (phenotype in first_half_phenos) {
    call rerun_finemap_and_followup_for_pheno.rerun_finemap_and_followup_for_pheno { input :
      phenotype = phenotype
    }
  }

  output {
    Array[Array[serializable_FINEMAP_output]] original_finemap = rerun_finemap_and_followup_for_pheno.original_finemap
    Array[Array[Array[File]]] original_finemap_creds = rerun_finemap_and_followup_for_pheno.original_finemap_creds

    Array[File] first_pass_df = rerun_finemap_and_followup_for_pheno.first_pass_df
    Array[File] susie_min_abs_corrs = rerun_finemap_and_followup_for_pheno.susie_min_abs_corrs

    Array[File] followup_regions_tsv = rerun_finemap_and_followup_for_pheno.followup_regions_tsv

    Array[Array[serializable_SuSiE_output]] best_guess_susie = rerun_finemap_and_followup_for_pheno.best_guess_susie
    Array[Array[Array[File]]] best_guess_susie_CSs = rerun_finemap_and_followup_for_pheno.best_guess_susie_CSs
    Array[Array[serializable_FINEMAP_output]] total_prob_finemap = rerun_finemap_and_followup_for_pheno.total_prob_finemap
    Array[Array[Array[File]]] total_prob_finemap_creds = rerun_finemap_and_followup_for_pheno.total_prob_finemap_creds
    Array[Array[serializable_FINEMAP_output]] derived_prior_std_finemap = rerun_finemap_and_followup_for_pheno.derived_prior_std_finemap
    Array[Array[Array[File]]] derived_prior_std_finemap_creds = rerun_finemap_and_followup_for_pheno.derived_prior_std_finemap_creds
    Array[Array[serializable_FINEMAP_output]] prob_conv_sss_tol_finemap = rerun_finemap_and_followup_for_pheno.prob_conv_sss_tol_finemap
    Array[Array[Array[File]]] prob_conv_sss_tol_finemap_creds = rerun_finemap_and_followup_for_pheno.prob_conv_sss_tol_finemap_creds
    Array[Array[serializable_FINEMAP_output]] mac_finemap = rerun_finemap_and_followup_for_pheno.mac_finemap
    Array[Array[Array[File]]] mac_finemap_creds = rerun_finemap_and_followup_for_pheno.mac_finemap_creds
    Array[Array[serializable_FINEMAP_output]] threshold_finemap = rerun_finemap_and_followup_for_pheno.threshold_finemap
    Array[Array[Array[File]]] threshold_finemap_creds = rerun_finemap_and_followup_for_pheno.threshold_finemap_creds

    Array[Array[serializable_SuSiE_output]] ratio_susie = rerun_finemap_and_followup_for_pheno.ratio_susie
    Array[Array[Array[File]]] ratio_susie_CSs = rerun_finemap_and_followup_for_pheno.ratio_susie_CSs
    Array[Array[serializable_FINEMAP_output]] ratio_finemap = rerun_finemap_and_followup_for_pheno.ratio_finemap
    Array[Array[Array[File]]] ratio_finemap_creds = rerun_finemap_and_followup_for_pheno.ratio_finemap_creds
    Array[Array[serializable_FINEMAP_output]] low_prior_std_finemap = rerun_finemap_and_followup_for_pheno.low_prior_std_finemap
    Array[Array[Array[File]]] low_prior_std_finemap_creds = rerun_finemap_and_followup_for_pheno.low_prior_std_finemap_creds

    Array[File] followup_df = select_all(rerun_finemap_and_followup_for_pheno.followup_df)
  }
}
