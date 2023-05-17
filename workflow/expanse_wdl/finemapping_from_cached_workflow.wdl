version 1.0

import "expanse_tasks.wdl"
import "../gwas_wdl/gwas_workflow.wdl"
import "../finemapping_wdl/finemapping_workflow.wdl"

workflow finemapping_from_cached {

  input {
    String script_dir  = "."

    File chr_lens = "misc_data/genome/chr_lens.txt"

    Boolean is_binary = false

    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample" # could instead create a task for downloading this with ukbgene

    String phenotype_name = "platelet_count"
    File cached_shared_covars = "traits/shared_covars/shared_covars.npy"
    File cached_samples_for_phenotype = "sample_qc/runs/white_brits/platelet_count/combined_unrelated.sample"
    File cached_transformed_phenotype = "traits/subset_transformed_phenotypes/white_brits/platelet_count.npy"
    File cached_my_str_gwas = "association/results/platelet_count/my_str/results.tab"
    Array[File] cached_ethnic_str_gwass = [
      "association/results_finemapped_only/black/platelet_count/my_str/results.tab",
      "association/results_finemapped_only/south_asian/platelet_count/my_str/results.tab",
      "association/results_finemapped_only/chinese/platelet_count/my_str/results.tab",
      "association/results_finemapped_only/irish/platelet_count/my_str/results.tab",
      "association/results_finemapped_only/white_other/platelet_count/my_str/results.tab"
    ]
    File cached_plink_snp_gwas = "association/results/platelet_count/plink_snp/results.tab"
  }

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

  call finemapping_workflow.finemapping { input :
    script_dir = script_dir,
    finemap_command = "finemap",

    chr_lens = chr_lens,

    str_vcfs = str_vcfs,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,

    shared_covars = cached_shared_covars,
    phenotype_samples = cached_samples_for_phenotype,
    transformed_phenotype_data = cached_transformed_phenotype,
    
    my_str_gwas = cached_my_str_gwas,
    ethnic_my_str_gwass =  cached_ethnic_str_gwass,
    plink_snp_gwas = cached_plink_snp_gwas,

    phenotype_name = phenotype_name,
    is_binary = is_binary,

    all_samples_list = all_samples_list
  }
}
