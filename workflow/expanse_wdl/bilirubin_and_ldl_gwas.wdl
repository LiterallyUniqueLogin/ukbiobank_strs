version 1.0

import "gwas.wdl"
import "../gwas_wdl/gwas_tasks.wdl"
import "expanse_files.wdl"

workflow bilirubin_and_ldl_gwas {

  call gwas_tasks.phenotype_names

  call expanse_files.files

  Int bili_idx = phenotype_names.idxs["total_bilirubin"]

  call gwas.gwas as bilirubin_gwas { input :
    phenotype_id = 30840,
    categorical_covariate_names = ["total_bilirubin_aliquot"],
    categorical_covariate_ids = [30842],
    phenotype_name = "total_bilirubin",

    cached_unrelated_samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[bili_idx],
    cached_shared_covars = files.shared_covars
  }

  Int ldl_idx = phenotype_names.idxs["ldl_cholesterol_direct"]

  call gwas.gwas as ldl_gwas { input :
    phenotype_id = 30780,
    categorical_covariate_names = ["ldl_cholesterol_direct_aliquot"],
    categorical_covariate_ids = [30782],
    phenotype_name = "ldl_cholesterol_direct",

    cached_unrelated_samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[ldl_idx],
    cached_shared_covars = files.shared_covars
  }

  output {
    Array[File] bilirubin_sample_lists = bilirubin_gwas.sample_lists

    File bilirubin_shared_covars = bilirubin_gwas.shared_covars
    File bilirubin_shared_covar_names = bilirubin_gwas.shared_covar_names

    Array[File] bilirubin_all_samples_for_phenotype = bilirubin_gwas.all_samples_for_phenotype
    Array[File] bilirubin_samples_for_phenotype = bilirubin_gwas.samples_for_phenotype

    Array[File] bilirubin_pheno_data = bilirubin_gwas.pheno_data
    Array[File] bilirubin_transformed_trait_values = bilirubin_gwas.transformed_trait_values
    Array[File] bilirubin_pheno_covar_names = bilirubin_gwas.pheno_covar_names
    Array[File] bilirubin_pheno_readme = bilirubin_gwas.pheno_readme

    File bilirubin_my_str_gwas = bilirubin_gwas.my_str_gwas
    File bilirubin_plink_snp_gwas = bilirubin_gwas.plink_snp_gwas
    File bilirubin_peaks = bilirubin_gwas.peaks
    File bilirubin_peaks_readme = bilirubin_gwas.peaks_readme
    File bilirubin_finemapping_regions = bilirubin_gwas.finemapping_regions
    File bilirubin_finemapping_regions_readme = bilirubin_gwas.finemapping_regions_readme
    Array[File] bilirubin_ethnic_my_str_gwas = bilirubin_gwas.ethnic_my_str_gwas

    Array[File] ldl_sample_lists = ldl_gwas.sample_lists

    File ldl_shared_covars = ldl_gwas.shared_covars
    File ldl_shared_covar_names = ldl_gwas.shared_covar_names

    Array[File] ldl_all_samples_for_phenotype = ldl_gwas.all_samples_for_phenotype
    Array[File] ldl_samples_for_phenotype = ldl_gwas.samples_for_phenotype

    Array[File] ldl_pheno_data = ldl_gwas.pheno_data
    Array[File] ldl_transformed_trait_values = ldl_gwas.transformed_trait_values
    Array[File] ldl_pheno_covar_names = ldl_gwas.pheno_covar_names
    Array[File] ldl_pheno_readme = ldl_gwas.pheno_readme

    File ldl_my_str_gwas = ldl_gwas.my_str_gwas
    File ldl_plink_snp_gwas = ldl_gwas.plink_snp_gwas
    File ldl_peaks = ldl_gwas.peaks
    File ldl_peaks_readme = ldl_gwas.peaks_readme
    File ldl_finemapping_regions = ldl_gwas.finemapping_regions
    File ldl_finemapping_regions_readme = ldl_gwas.finemapping_regions_readme
    Array[File] ldl_ethnic_my_str_gwas = ldl_gwas.ethnic_my_str_gwas
  }
}
