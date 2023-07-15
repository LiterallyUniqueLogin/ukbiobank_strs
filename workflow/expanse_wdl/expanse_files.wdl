version 1.0

import "../gwas_wdl/gwas_tasks.wdl"

workflow files {

  Array[String] ethnicities_ = [
    "black",
    "south_asian",
    "chinese",
    "irish",
    "white_other",
  ]

  Array[String] all_ethnicities_ = flatten([["white_brits"], ethnicities])

  scatter (ethnicity in all_ethnicities) {
    File qced_samples_for_ethnicity_ = "sample_qc/runs/~{ethnicity}/no_phenotype/combined.sample"
    File unrelated_samples_for_ethnicity_ = "sample_qc/runs/~{ethnicity}/no_phenotype/combined_unrelated.sample"
    scatter (phenotype in phenotype_names.n) {
      File unrelated_samples_for_ethnicity_for_phenotype_  = "sample_qc/runs/~{ethnicity}/~{phenotype}/combined_unrelated.sample"
    }
  }
  scatter (phenotype in phenotype_names.n) {
    scatter (ethnicity in all_ethnicities) {
      File unrelated_samples_for_pheno_for_ethnicity_ = "sample_qc/runs/~{ethnicity}/~{phenotype}/combined_unrelated.sample"
    }
  }

  scatter (chrom in range(22)) {
    VCF str_vcfs_ = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
     PFiles imputed_snp_pfiles_ = {
       "pgen": "array_imputed/pfile_converted/chr~{chrom+1}.pgen",
       "pvar": "array_imputed/pfile_converted/chr~{chrom+1}.pvar",
       "psam": "array_imputed/pfile_converted/chr~{chrom+1}.psam",
     }
    bgen imputed_snp_bgens_ = {
      "bgen": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen",
      "index": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.bgi",
      "bgen_reader_metadata": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata",
      "bgen_reader_metadata2": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata2.mmm",
      "bgen_reader_complex_metadata2": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.complex.metadata2.mmm",
    }
    File snps_to_filter_ = "finemapping/str_imp_snp_overlaps/chr~{chrom+1}_to_filter.tab"

    File intersects_gene_annotation_ = "side_analyses/str_annotations/intersects_gene/chr~{chrom+1}.tab"
    File intersects_exon_annotation_ = "side_analyses/str_annotations/intersects_exon/chr~{chrom+1}.tab"
    File intersects_CDS_annotation_ = "side_analyses/str_annotations/intersects_CDS/chr~{chrom+1}.tab"
    File intersects_five_prime_UTR_annotation_ = "side_analyses/str_annotations/intersects_five_prime_UTR/chr~{chrom+1}.tab"
    File intersects_three_prime_UTR_annotation_ = "side_analyses/str_annotations/intersects_three_prime_UTR/chr~{chrom+1}.tab"
    File intersects_UTR_annotation_ = "side_analyses/str_annotations/intersects_UTR/chr~{chrom+1}.tab"
    File closest_gene_annotation_ = "side_analyses/str_annotations/closest_gene/chr~{chrom+1}.tab"
    File intersects_transcript_support_2_annotation_ = "side_analyses/str_annotations/intersects_transcript_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_CDS_support_2_annotation_ = "side_analyses/str_annotations/intersects_protein_coding_CDS_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_five_prime_UTR_support_2_annotation_ = "side_analyses/str_annotations/intersects_protein_coding_five_prime_UTR_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_three_prime_UTR_support_2_annotation_ = "side_analyses/str_annotations/intersects_protein_coding_three_prime_UTR_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_UTR_support_2_annotation_ = "side_analyses/str_annotations/intersects_protein_coding_UTR_support_2/chr~{chrom+1}.tab"
    File closest_downstream_protein_coding_exon_support_2_annotation_ = "side_analyses/str_annotations/closest_downstream_protein_coding_exon_support_2/chr~{chrom+1}.tab"
    File closest_upstream_protein_coding_exon_support_2_annotation_ = "side_analyses/str_annotations/closest_upstream_protein_coding_exon_support_2/chr~{chrom+1}.tab"
    File closest_downstream_protein_coding_gene_annotation_ = "side_analyses/str_annotations/closest_downstream_protein_coding_gene/chr~{chrom+1}.tab"
    File closest_upstream_protein_coding_gene_annotation_ = "side_analyses/str_annotations/closest_upstream_protein_coding_gene/chr~{chrom+1}.tab"
  }

  call gwas_tasks.phenotype_names

  scatter (ethnicity in all_ethnicities_) {
    scatter (phenotype in phenotype_names.n) {
      if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
        File ethnic_to_pheno_to_covar_names_wdl = "wdl_cache/~{phenotype}_~{ethnicity}_covar_names.txt"
        File ethnic_to_pheno_to_transformed_phenotype_data_wdl = "wdl_cache/~{phenotype}_~{ethnicity}_rank_inverse_normalized.npy"
      }
      if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
        File ethnic_to_pheno_to_covar_names_snakemake = "traits/phenotypes/~{ethnicity}/~{phenotype}_covar_names.txt"
        File ethnic_to_pheno_to_transformed_phenotype_data_snakemake = "traits/subset_transformed_phenotypes/~{ethnicity}/~{phenotype}.npy"
      }
      File ethnic_to_pheno_to_covar_names_ = select_first([ethnic_to_pheno_to_covar_names_snakemake, ethnic_to_pheno_to_covar_names_wdl])
      File ethnic_to_pheno_to_transformed_phenotype_data_ = select_first([ethnic_to_pheno_to_transformed_phenotype_data_snakemake, ethnic_to_pheno_to_transformed_phenotype_data_wdl])
    }
  }

  scatter (ethnicity in ethnicities_) {
    scatter (phenotype in phenotype_names.n) {
      if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
        File ethnic_to_pheno_to_str_gwas_results_wdl = "wdl_cache/~{phenotype}_~{ethnicity}_str_gwas_results.tab"
      }
      if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
        File ethnic_to_pheno_to_str_gwas_results_snakemake = "association/results_finemapped_only/~{ethnicity}/~{phenotype}/my_str/results.tab"
      }
      File ethnic_to_pheno_to_str_gwas_results_ = select_first([ethnic_to_pheno_to_str_gwas_results_snakemake , ethnic_to_pheno_to_str_gwas_results_wdl])
    }
  }


  scatter (phenotype in phenotype_names.n) {
    if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
      File str_gwas_results_wdl = "wdl_cache/~{phenotype}_str_gwas_results.tab"
      File snp_gwas_results_wdl  = "wdl_cache/~{phenotype}_snp_gwas_results.tab"
      File peaks_for_1ef_wdl  = "wdl_cache/~{phenotype}_peaks_250000_5e-8.tab"
      File finemapping_regions_wdl = "wdl_cache/~{phenotype}_finemapping_regions.tab"
      scatter (ethnicity in ethnicities_) {
        File pheno_to_ethnic_to_str_gwas_results_wdl = "wdl_cache/~{phenotype}_~{ethnicity}_str_gwas_results.tab"
      }
    }
    if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
      File str_gwas_results_snakemake = "association/results/~{phenotype}/my_str/results.tab"
      File snp_gwas_results_snakemake = "association/results/~{phenotype}/plink_snp/results.tab"
      File peaks_for_1ef_snakemake = "signals/peaks/~{phenotype}_250000_5e-8.tab"
      File finemapping_regions_snakemake = "signals/regions/~{phenotype}.tab"
      scatter (ethnicity in ethnicities_) {
        File pheno_to_ethnic_to_str_gwas_results_snakemake = "association/results_finemapped_only/~{ethnicity}//~{phenotype}/my_str/results.tab"
      }
    }
    File str_gwas_results_ = select_first([str_gwas_results_wdl, str_gwas_results_snakemake])
    File snp_gwas_results_ = select_first([snp_gwas_results_wdl, snp_gwas_results_snakemake])
    File peaks_for_1ef_ = select_first([peaks_for_1ef_wdl, peaks_for_1ef_snakemake])
    File finemapping_regions_ = select_first([finemapping_regions_wdl, finemapping_regions_snakemake])
    Array[File] pheno_to_ethnic_to_str_gwas_results_ = select_first([pheno_to_ethnic_to_str_gwas_results_wdl, pheno_to_ethnic_to_str_gwas_results_snakemake])
    File finemapping_first_pass_dfs_ = "wdl_cache/finemapping/finemapping_all_regions_concordance_~{phenotype}.tab"
    File susie_min_abs_corrs_ = "wdl_cache/finemapping/susie_all_regions_min_abs_corrs_~{phenotype}.tab"
    File finemapping_followup_dfs_ = "wdl_cache/finemapping/finemapping_followup_concordance_~{phenotype}.tab" # bigger for mean platelet volume, so needs consistent filtering
  }

  output {
    # ------------ cached files directly from UKB
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample"
    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam"
    File withdrawn_sample_list = "sample_qc/common_filters/remove/withdrawn.sample"
    File kinship = "misc_data/ukbgene/ukb46122_rel_s488282.dat" # could create a task for downloading this with ukbgene
    # indexed by chrom
    Array[bgen] imputed_snp_bgens = imputed_snp_bgens_
    String date_of_most_recent_first_occurrence_update = "2021-04-01" # only needed for binary phenotypes

    # -------------- cached other data source files or information I hardcoded
    # -------------- for some of these, could replace caching with WDL download scripts
    File chr_lens = "misc_data/genome/chr_lens.txt"
    File gencode = "misc_data/gencode/gencode.v38lift37.annotation.without_chr.sorted.gene.gff3"
    File pan_ukbb = "misc_data/snp_summary_stats/bilirubin/neale/biomarkers-30840-both_sexes-irnt.tsv"
    File specific_alleles = "association/specific_alleles.tab"
    File eSTR_table = "misc_data/eSTR/eSTRs.csv"
    File CBL_gtex_expression = "misc_data/gtex_yang/CBL_chr11_119206290_GTEX_TPM.tsv"
    File CBL_geuvadis_expression = "misc_data/geuvadis_melissa/CBL_Geuvadis_Exprdata.csv"
    File Liver_SLC2A2_exon4_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon4_psi.tsv"
    File Liver_SLC2A2_exon6_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon6_psi.tsv"
    File RHOT1_geuvadis_expression = "misc_data/geuvadis_melissa/RHOT1_Geuvadis_Exprdata.csv"
    Array[String] ethnicities = ethnicities_
    Array[String] all_ethnicities = all_ethnicities_

    # --------------------- cached files that I generated that aren't yet WDL enabled
    File flank_start_to_start_and_end_pos = "snpstr/flank_trimmed_vcf/vars.tab"
    File repeat_units_table = "snpstr/repeat_units.tab"
    File str_loci = "snpstr/str_loci.txt"
    File str_hg19_pos_bed = "snpstr/str_loci.bed"
    File str_hg38_pos_bed = "snpstr/str_loci.hg38.bed"
    File str_t2t_pos_bed = "snpstr/str_loci.t2tv2.bed"
    # indexed by chrom
    Array[VCF] str_vcfs = str_vcfs_
    Array[PFiles] imputed_snp_pfiles = imputed_snp_pfiles_
    Array[File] snps_to_filter = snps_to_filter_
    Array[File] intersects_gene_annotation = intersects_gene_annotation_
    Array[File] intersects_exon_annotation = intersects_exon_annotation_
    Array[File] intersects_CDS_annotation = intersects_CDS_annotation_
    Array[File] intersects_five_prime_UTR_annotation = intersects_five_prime_UTR_annotation_
    Array[File] intersects_three_prime_UTR_annotation = intersects_three_prime_UTR_annotation_
    Array[File] intersects_UTR_annotation = intersects_UTR_annotation_
    Array[File] closest_gene_annotation = closest_gene_annotation_
    Array[File] intersects_transcript_support_2_annotation = intersects_transcript_support_2_annotation_
    Array[File] intersects_protein_coding_CDS_support_2_annotation = intersects_protein_coding_CDS_support_2_annotation_
    Array[File] intersects_protein_coding_five_prime_UTR_support_2_annotation = intersects_protein_coding_five_prime_UTR_support_2_annotation_
    Array[File] intersects_protein_coding_three_prime_UTR_support_2_annotation = intersects_protein_coding_three_prime_UTR_support_2_annotation_
    Array[File] intersects_protein_coding_UTR_support_2_annotation = intersects_protein_coding_UTR_support_2_annotation_
    Array[File] closest_downstream_protein_coding_exon_support_2_annotation = closest_downstream_protein_coding_exon_support_2_annotation_
    Array[File] closest_upstream_protein_coding_exon_support_2_annotation = closest_upstream_protein_coding_exon_support_2_annotation_
    Array[File] closest_downstream_protein_coding_gene_annotation = closest_downstream_protein_coding_gene_annotation_
    Array[File] closest_upstream_protein_coding_gene_annotation = closest_upstream_protein_coding_gene_annotation_
    # TODO files that are currently cached that MUST be regenerated with WDL
    File eQTL_table = "post_finemapping/gtex_STR2/blessed_qtl_STRs.tab"

    # ------------ cached files that can be regenerated with WDL workflows, though
    # ------------ may change due to code updates, withdrawing samples,
    # ------------ unrelatedness randomness and FINEMAP randomness
    # indexed by all_ethnicities
    # unrelated also implies qced
    Array[File] qced_samples_for_ethnicity = qced_samples_for_ethnicity_
    Array[File] unrelated_samples_for_ethnicity = unrelated_samples_for_ethnicity_
    File unrelated_samples_CBL_hom_not_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    # indexed first by all_ethnicities, then phenotype
    Array[Array[File]] unrelated_samples_for_ethnicity_for_phenotype = unrelated_samples_for_ethnicity_for_phenotype_
    # indexed first by phenotype, then all_ethnicities
    Array[Array[File]] unrelated_samples_for_pheno_for_ethnicity = unrelated_samples_for_pheno_for_ethnicity_

    File shared_covars = "traits/shared_covars/shared_covars.npy"
    # indexed by all ethnicities, then by phenotype
    Array[Array[File]] ethnic_to_pheno_to_covar_names = ethnic_to_pheno_to_covar_names_
    Array[Array[File]] ethnic_to_pheno_to_transformed_phenotype_data = ethnic_to_pheno_to_transformed_phenotype_data_

    # indexed by phenotype names
    Array[File] str_gwas_results = str_gwas_results_
    Array[File] snp_gwas_results = snp_gwas_results_
    File CBL_conditioned_SNP_119080037_A_G_str_results = "association/results/platelet_crit/my_str_conditional/chr11_118447267_119339135_STR__ISNP_119080037_A_G__ASNP.tab"
    File CBL_conditioned_SNP_119080037_A_G_snp_results = "association/results/platelet_crit/plink_snp_conditional/chr11_118447267_119339135_STR__ISNP_119080037_A_G__ASNP/plink2.rin_platelet_crit.glm.linear.done"
    File CBL_conditioned_STR_19077000_SNP_119080037_A_G_str_results = "association/results/platelet_crit/my_str_conditional/chr11_118447267_119339135_STR_119077000__ISNP_119080037_A_G__ASNP.tab"
    File CBL_conditioned_STR_19077000_SNP_119080037_A_G_snp_results = "association/results/platelet_crit/plink_snp_conditional/chr11_118447267_119339135_STR_119077000__ISNP_119080037_A_G__ASNP/plink2.rin_platelet_crit.glm.linear.done"
    # TODO must regenerate
    #File SLC2A2_conditioned_STR_17100913_str_results = "association/results/total_bilirubin/my_str_conditional/chr11_118447267_119339135_STR_17100913__ISNP__ASNP.tab"
    #File SLC2A2_conditioned_STR_17100913_snp_results = "association/results/total_bilirubin/plink_snp_conditional/chr11_118447267_119339135_STR_17100913__ISNP__ASNP/plink2.rin_platelet_crit.glm.linear.done"
    File TAOK1_conditioned_STR_27842010_snp_results = "association/results/mean_platelet_volume/plink_snp_conditional/chr17_25885931_29975306_STR_27842010__ISNP__ASNP/plink2.rin_mean_platelet_volume.glm.linear.done"
    Array[File] peaks_for_1ef = peaks_for_1ef_
    Array[File] finemapping_regions = finemapping_regions_
    # Indexed first by phenotype then ethnicity
    Array[Array[File]] pheno_to_ethnic_to_str_gwas_results = pheno_to_ethnic_to_str_gwas_results_
    # indexed first by ethnicity then phenotype
    Array[Array[File]] ethnic_to_pheno_to_str_gwas_results = ethnic_to_pheno_to_str_gwas_results_ 

    Array[File] finemapping_first_pass_dfs = finemapping_first_pass_dfs_ 
    Array[File] susie_min_abs_corrs = susie_min_abs_corrs_
    Array[File] finemapping_followup_dfs = finemapping_first_pass_dfs_ # bigger for mean platelet volume, so needs consistent filtering
  }
}
