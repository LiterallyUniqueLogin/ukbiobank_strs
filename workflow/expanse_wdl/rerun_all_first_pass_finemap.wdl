version 1.0

import "../finemapping_wdl/first_pass_finemap_workflow.wdl"
import "../gwas_wdl/gwas_tasks.wdl"

workflow rerun_all_first_pass_finemap {

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

  call gwas_tasks.phenotype_names

  scatter (name in phenotype_names.n) {
    File cached_samples_for_phenotype = "sample_qc/runs/white_brits/~{name}/combined_unrelated.sample"
    File cached_my_str_gwas = "association/results/~{name}/my_str/results.tab"
    File cached_plink_snp_gwas = "association/results/~{name}/plink_snp/results.tab"

    call first_pass_finemap_workflow.first_pass_finemap { input :
      script_dir = script_dir,
      finemap_command = "finemap",

      chr_lens = chr_lens,

      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,

      phenotype_samples = cached_samples_for_phenotype,
      my_str_gwas = cached_my_str_gwas,
      plink_snp_gwas = cached_plink_snp_gwas,

      phenotype_name = name,

      all_samples_list = all_samples_list
    }
  }

  output {
    # per phenotype
    Array[File] regions_tsv = first_pass_finemap.regions_tsv
    Array[File] regions_readme = first_pass_finemap.regions_readme

    # per phenotype per region
    Array[Array[serializable_FINEMAP_output]] finemap = first_pass_finemap.finemap
    # per phenotype per region many creds
    Array[Array[Array[File]]] finemap_creds = first_pass_finemap.finemap_creds
  }
}
