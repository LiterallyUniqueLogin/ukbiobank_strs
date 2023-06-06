version 1.0

import "../finemapping_wdl/first_pass_finemapping_workflow.wdl"

workflow platelet_count_first_pass_finemapping {

	String script_dir  = "."

	File chr_lens = "misc_data/genome/chr_lens.txt"

	File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample" # could instead create a task for 

  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
    bgen imputed_snp_bgens = {
      "bgen": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen",
      "index": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.bgi",
      "bgen_reader_metadata" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata",
      "bgen_reader_metadata2" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata2.mmm",
      "bgen_reader_complex_metadata2" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.complex.metadata2.mmm"
    }
    File snp_vars_to_filter_from_finemapping = "finemapping/str_imp_snp_overlaps/chr~{chrom+1}_to_filter.tab"
  }

  File cached_samples_for_phenotype = "sample_qc/runs/white_brits/platelet_count/combined_unrelated.sample"
  File cached_my_str_gwas = "association/results/platelet_count/my_str/results.tab"
  File cached_plink_snp_gwas = "association/results/platelet_count/plink_snp/results.tab"
  File cached_shared_covars = "traits/shared_covars/shared_covars.npy"
  File cached_transformed_phenotype_data = "traits/subset_transformed_phenotypes/white_brits/platelet_count.npy"

  scatter (ethnicity in ["black", "south_asian", "chinese", "irish", "white_other"]) {
    File cached_ethnic_my_str_gwass = "association/results_finemapped_only/~{ethnicity}/platelet_count/my_str/results.tab"
  }

  call first_pass_finemapping_workflow.first_pass_finemapping { input :
    script_dir = script_dir,
    finemap_command = "finemap",

    chr_lens = chr_lens,

    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,

    phenotype_name = "platelet_count",

    phenotype_samples = cached_samples_for_phenotype,
    all_samples_list = all_samples_list,

    shared_covars = cached_shared_covars,
    transformed_phenotype_data = cached_transformed_phenotype_data,

    my_str_gwas = cached_my_str_gwas,
    plink_snp_gwas = cached_plink_snp_gwas,
    ethnic_my_str_gwass = cached_ethnic_my_str_gwass


  }

  output {
    File regions_tsv = first_pass_finemapping.regions_tsv
    File regions_readme = first_pass_finemapping.regions_readme

    Array[serializable_FINEMAP_output] finemap = first_pass_finemapping.finemap
    Array[Array[File]] finemap_creds = first_pass_finemapping.finemap_creds
    Array[serializable_SuSiE_output] susie = first_pass_finemapping.susie
    Array[Array[File]] susie_CSs = first_pass_finemapping.susie_CSs

    File first_pass_df = first_pass_finemapping.first_pass_df
  }
}
