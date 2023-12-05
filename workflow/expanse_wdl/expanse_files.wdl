version 1.0

import "../gwas_wdl/gwas_tasks.wdl"

workflow files {

  inputs {
    String base_location = "/expanse/projects/gymreklab/jmargoli/ukbiobank"
  }

  Array[String] ethnicities_ = [
    "black",
    "south_asian",
    "chinese",
    "irish",
    "white_other",
  ]

  Array[String] all_ethnicities_ = flatten([["white_brits"], ethnicities])

  scatter (ethnicity in all_ethnicities) {
    File qced_samples_for_ethnicity_ = "~{base_location}/sample_qc/runs/~{ethnicity}/no_phenotype/combined.sample"
    File unrelated_samples_for_ethnicity_ = "~{base_location}/sample_qc/runs/~{ethnicity}/no_phenotype/combined_unrelated.sample"
    scatter (phenotype in phenotype_names.n) {
      File unrelated_samples_for_ethnicity_for_phenotype_  = "~{base_location}/sample_qc/runs/~{ethnicity}/~{phenotype}/combined_unrelated.sample"
    }
  }
  scatter (phenotype in phenotype_names.n) {
    scatter (ethnicity in all_ethnicities) {
      File unrelated_samples_for_pheno_for_ethnicity_ = "~{base_location}/sample_qc/runs/~{ethnicity}/~{phenotype}/combined_unrelated.sample"
    }
  }

  scatter (chrom in range(22)) {
    VCF str_vcfs_ = {
      "vcf": "~{base_location}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "~{base_location}/str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
     PFiles imputed_snp_pfiles_ = {
       "pgen": "~{base_location}/array_imputed/pfile_converted/chr~{chrom+1}.pgen",
       "pvar": "~{base_location}/array_imputed/pfile_converted/chr~{chrom+1}.pvar",
       "psam": "~{base_location}/array_imputed/pfile_converted/chr~{chrom+1}.psam",
     }
    bgen imputed_snp_bgens_ = {
      "bgen": "~{base_location}/array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen",
      "index": "~{base_location}/array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.bgi",
      "bgen_reader_metadata": "~{base_location}/array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata",
      "bgen_reader_metadata2": "~{base_location}/array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata2.mmm",
      "bgen_reader_complex_metadata2": "~{base_location}/array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.complex.metadata2.mmm",
    }
    File snps_to_filter_ = "~{base_location}/finemapping/str_imp_snp_overlaps/chr~{chrom+1}_to_filter.tab"

    File intersects_gene_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_gene/chr~{chrom+1}.tab"
    File intersects_exon_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_exon/chr~{chrom+1}.tab"
    File intersects_CDS_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_CDS/chr~{chrom+1}.tab"
    File intersects_five_prime_UTR_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_five_prime_UTR/chr~{chrom+1}.tab"
    File intersects_three_prime_UTR_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_three_prime_UTR/chr~{chrom+1}.tab"
    File intersects_UTR_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_UTR/chr~{chrom+1}.tab"
    File closest_gene_annotation_ = "~{base_location}/side_analyses/str_annotations/closest_gene/chr~{chrom+1}.tab"
    File intersects_transcript_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_transcript_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_CDS_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_protein_coding_CDS_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_five_prime_UTR_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_protein_coding_five_prime_UTR_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_three_prime_UTR_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_protein_coding_three_prime_UTR_support_2/chr~{chrom+1}.tab"
    File intersects_protein_coding_UTR_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/intersects_protein_coding_UTR_support_2/chr~{chrom+1}.tab"
    File closest_downstream_protein_coding_exon_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/closest_downstream_protein_coding_exon_support_2/chr~{chrom+1}.tab"
    File closest_upstream_protein_coding_exon_support_2_annotation_ = "~{base_location}/side_analyses/str_annotations/closest_upstream_protein_coding_exon_support_2/chr~{chrom+1}.tab"
    File closest_downstream_protein_coding_gene_annotation_ = "~{base_location}/side_analyses/str_annotations/closest_downstream_protein_coding_gene/chr~{chrom+1}.tab"
    File closest_upstream_protein_coding_gene_annotation_ = "~{base_location}/side_analyses/str_annotations/closest_upstream_protein_coding_gene/chr~{chrom+1}.tab"
  }

  call gwas_tasks.phenotype_names

  scatter (ethnicity in all_ethnicities_) {
    scatter (phenotype in phenotype_names.n) {
      if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
        File ethnic_to_pheno_to_covar_names_wdl = "~{base_location}/wdl_cache/~{phenotype}_~{ethnicity}_covar_names.txt"
        File ethnic_to_pheno_to_transformed_phenotype_data_wdl = "~{base_location}/wdl_cache/~{phenotype}_~{ethnicity}_rank_inverse_normalized.npy"
      }
      if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
        File ethnic_to_pheno_to_covar_names_snakemake = "~{base_location}/traits/phenotypes/~{ethnicity}/~{phenotype}_covar_names.txt"
        File ethnic_to_pheno_to_transformed_phenotype_data_snakemake = "~{base_location}/traits/subset_transformed_phenotypes/~{ethnicity}/~{phenotype}.npy"
      }
      File ethnic_to_pheno_to_covar_names_ = select_first([ethnic_to_pheno_to_covar_names_snakemake, ethnic_to_pheno_to_covar_names_wdl])
      File ethnic_to_pheno_to_transformed_phenotype_data_ = select_first([ethnic_to_pheno_to_transformed_phenotype_data_snakemake, ethnic_to_pheno_to_transformed_phenotype_data_wdl])
    }

    File wgs_comparison_stats_ = "~{base_location}/wgs/20230907_new_old_~{ethnicity}-locuscompare.tab"
    File wgs_allele_freqs_ = "~{base_location}/wgs/20230908_~{ethnicity}_new_old_WGS_allele_freqs.tab"
  }

  scatter (ethnicity in ethnicities_) {
    scatter (phenotype in phenotype_names.n) {
      if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
        File ethnic_to_pheno_to_str_gwas_results_wdl = "~{base_location}/wdl_cache/~{phenotype}_~{ethnicity}_str_gwas_results.tab"
      }
      if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
        File ethnic_to_pheno_to_str_gwas_results_snakemake = "~{base_location}/association/results_finemapped_only/~{ethnicity}/~{phenotype}/my_str/results.tab"
      }
      File ethnic_to_pheno_to_str_gwas_results_ = select_first([ethnic_to_pheno_to_str_gwas_results_snakemake , ethnic_to_pheno_to_str_gwas_results_wdl])
    }
  }


  scatter (phenotype in phenotype_names.n) {
    if (phenotype == "ldl_cholesterol_direct" || phenotype == "total_bilirubin") {
      File str_gwas_results_wdl = "~{base_location}/wdl_cache/~{phenotype}_str_gwas_results.tab"
      File snp_gwas_results_wdl  = "~{base_location}/wdl_cache/~{phenotype}_snp_gwas_results.tab"
      File peaks_for_1ef_wdl  = "~{base_location}/wdl_cache/~{phenotype}_peaks_250000_5e-8.tab"
      File finemapping_regions_wdl = "~{base_location}/wdl_cache/~{phenotype}_finemapping_regions.tab"
      scatter (ethnicity in ethnicities_) {
        File pheno_to_ethnic_to_str_gwas_results_wdl = "~{base_location}/wdl_cache/~{phenotype}_~{ethnicity}_str_gwas_results.tab"
      }
    }
    if (phenotype != "ldl_cholesterol_direct" && phenotype != "total_bilirubin") {
      File str_gwas_results_snakemake = "~{base_location}/association/results/~{phenotype}/my_str/results.tab"
      File snp_gwas_results_snakemake = "~{base_location}/association/results/~{phenotype}/plink_snp/results.tab"
      File peaks_for_1ef_snakemake = "~{base_location}/signals/peaks/~{phenotype}_250000_5e-8.tab"
      File finemapping_regions_snakemake = "~{base_location}/signals/regions/~{phenotype}.tab"
      scatter (ethnicity in ethnicities_) {
        File pheno_to_ethnic_to_str_gwas_results_snakemake = "~{base_location}/association/results_finemapped_only/~{ethnicity}/~{phenotype}/my_str/results.tab"
      }
    }
    File str_gwas_results_ = select_first([str_gwas_results_wdl, str_gwas_results_snakemake])
    File snp_gwas_results_ = select_first([snp_gwas_results_wdl, snp_gwas_results_snakemake])
    File peaks_for_1ef_ = select_first([peaks_for_1ef_wdl, peaks_for_1ef_snakemake])
    File finemapping_regions_ = select_first([finemapping_regions_wdl, finemapping_regions_snakemake])
    Array[File] pheno_to_ethnic_to_str_gwas_results_ = select_first([pheno_to_ethnic_to_str_gwas_results_wdl, pheno_to_ethnic_to_str_gwas_results_snakemake])
    File finemapping_first_pass_dfs_ = "~{base_location}/wdl_cache/finemapping/finemapping_all_regions_concordance_~{phenotype}.tab"
    File susie_min_abs_corrs_ = "~{base_location}/wdl_cache/finemapping/susie_all_regions_min_abs_corrs_~{phenotype}.npy"
    File finemapping_followup_dfs_ = "~{base_location}/wdl_cache/finemapping/finemapping_followup_concordance_~{phenotype}.tab" # bigger for mean platelet volume, so needs consistent filtering
  }

  scatter (i in range(548)) {
    File platelet_count_snp_macs_ = "~{base_location}/wdl_cache/platelet_count_snp_macs/region_~{i}_snp_macs.tab"
  }

  output {
    # ------------ cached files directly from UKB
    File all_samples_list = "~{base_location}/microarray/ukb46122_hap_chr1_v2_s487314.sample"
    File fam_file = "~{base_location}/microarray/ukb46122_cal_chr1_v2_s488176.fam"
    File withdrawn_sample_list = "~{base_location}/sample_qc/common_filters/remove/withdrawn.sample"
    File kinship = "~{base_location}/misc_data/ukbgene/ukb46122_rel_s488282.dat" # could create a task for downloading this with ukbgene
    # indexed by chrom
    Array[bgen] imputed_snp_bgens = imputed_snp_bgens_
    String date_of_most_recent_first_occurrence_update = "2021-04-01" # only needed for binary phenotypes

    # -------------- cached other data source files or information I hardcoded
    # -------------- for some of these, could replace caching with WDL download scripts
    File chr_lens = "~{base_location}/misc_data/genome/chr_lens.txt"
    File gencode = "~{base_location}/misc_data/gencode/gencode.v38lift37.annotation.without_chr.sorted.gene.gff3"
    File pan_ukbb = "~{base_location}/misc_data/snp_summary_stats/bilirubin/neale/biomarkers-30840-both_sexes-irnt.tsv"
    File specific_alleles = "~{base_location}/association/specific_alleles.tab"
    File eSTR_table = "~{base_location}/misc_data/eSTR/eSTRs.csv"
    File all_eQTL_results = "/expanse/projects/gymreklab/yal084_storage/share_with_Jonathan/eSTR_GB_650pc_combined_fdr10p.csv"
    File all_meQTL_results = "/expanse/projects/gymreklab/yal084_storage/share_with_Jonathan/Combined_DNAm_regression.tsv"
    #File confidently_finemapped_methylation_results = "misc_data/gtex_yang/methylation_confidently_finemapped_strs.tsv"
    File CBL_gtex_expression = "~{base_location}/misc_data/gtex_yang/CBL_chr11_119206290_GTEX_TPM.tsv"
    File CBL_geuvadis_expression = "~{base_location}/misc_data/geuvadis_melissa/CBL_Geuvadis_Exprdata.csv"
    File Liver_SLC2A2_exon4_psi= "~{base_location}/misc_data/gtex_yang/Liver_SLC2A2_exon4_psi.tsv"
    File Liver_SLC2A2_exon6_psi= "~{base_location}/misc_data/gtex_yang/Liver_SLC2A2_exon6_psi.tsv"
    File RHOT1_geuvadis_expression = "~{base_location}/misc_data/geuvadis_melissa/RHOT1_Geuvadis_Exprdata.csv"
    File TAOK1_expression = "~{base_location}/misc_data/gtex_yang/chr17_29514992-ENSG00000160551.11_Thyroid.tsv"
    File TAOK1_residual_expression  = "~{base_location}/misc_data/gtex_yang/chr17_29514992-ENSG00000160551.11_Thyroid_regressing_out_covariates.tsv"
    File cg17724175_methylation = "~{base_location}/misc_data/gtex_yang/chr1_150579759-cg17724175.tsv"
    File cg22674798_methylation = "~{base_location}/misc_data/gtex_yang/chr1_3170058-cg22674798.tsv"
    Array[String] ethnicities = ethnicities_
    Array[String] all_ethnicities = all_ethnicities_

    # from WGS
    Array[File] wgs_comparison_stats = wgs_comparison_stats_
    Array[File] wgs_allele_freqs = wgs_allele_freqs_

    # --------------------- cached files that I generated that aren't yet WDL enabled
    File flank_start_to_start_and_end_pos = "~{base_location}/snpstr/flank_trimmed_vcf/vars.tab"
    File repeat_units_table = "~{base_location}/snpstr/repeat_units.tab"
    File str_loci = "~{base_location}/snpstr/str_loci.txt"
    File str_hg19_pos_bed = "~{base_location}/snpstr/str_loci.bed"
    File str_hg38_pos_bed = "~{base_location}/snpstr/str_loci.hg38.bed"
    File str_t2t_pos_bed = "~{base_location}/snpstr/str_loci.t2tv2.bed"
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

    # ------------ cached files that can be regenerated with WDL workflows, though
    # ------------ may change due to code updates, withdrawing samples,
    # ------------ unrelatedness randomness and FINEMAP randomness
    # indexed by all_ethnicities
    # unrelated also implies qced
    Array[File] qced_samples_for_ethnicity = qced_samples_for_ethnicity_
    Array[File] unrelated_samples_for_ethnicity = unrelated_samples_for_ethnicity_
    File unrelated_samples_CBL_hom_not_begin_C_T_snp = "~{base_location}/sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_begin_C_T_snp = "~{base_location}/sample_qc/subpop_runs/CBL_hom_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    # indexed first by all_ethnicities, then phenotype
    Array[Array[File]] unrelated_samples_for_ethnicity_for_phenotype = unrelated_samples_for_ethnicity_for_phenotype_
    # indexed first by phenotype, then all_ethnicities
    Array[Array[File]] unrelated_samples_for_pheno_for_ethnicity = unrelated_samples_for_pheno_for_ethnicity_

    File shared_covars = "~{base_location}/traits/shared_covars/shared_covars.npy"
    # indexed by all ethnicities, then by phenotype
    Array[Array[File]] ethnic_to_pheno_to_covar_names = ethnic_to_pheno_to_covar_names_
    Array[Array[File]] ethnic_to_pheno_to_transformed_phenotype_data = ethnic_to_pheno_to_transformed_phenotype_data_

    # indexed by phenotype names
    Array[File] str_gwas_results = str_gwas_results_
    Array[File] snp_gwas_results = snp_gwas_results_
    File CBL_conditioned_SNP_119080037_A_G_str_results = "~{base_location}/association/results/platelet_crit/my_str_conditional/chr11_118447267_119339135_STR__ISNP_119080037_A_G__ASNP.tab"
    File CBL_conditioned_SNP_119080037_A_G_snp_results = "~{base_location}/association/results/platelet_crit/plink_snp_conditional/chr11_118447267_119339135_STR__ISNP_119080037_A_G__ASNP/plink2.rin_platelet_crit.glm.linear.done"
    File CBL_conditioned_STR_19077000_SNP_119080037_A_G_str_results = "~{base_location}/association/results/platelet_crit/my_str_conditional/chr11_118447267_119339135_STR_119077000__ISNP_119080037_A_G__ASNP.tab"
    File CBL_conditioned_STR_19077000_SNP_119080037_A_G_snp_results = "~{base_location}/association/results/platelet_crit/plink_snp_conditional/chr11_118447267_119339135_STR_119077000__ISNP_119080037_A_G__ASNP/plink2.rin_platelet_crit.glm.linear.done"
    File RHOT1_conditioned_str_results = "~{base_location}/association/results/red_blood_cell_distribution_width/my_str_conditional/chr17_30287357_30595028_STR_30469467__ISNP__ASNP.tab"
    File RHOT1_conditioned_snp_results = "~{base_location}/association/results/red_blood_cell_distribution_width/plink_snp_conditional/chr17_30287357_30595028_STR_30469467__ISNP__ASNP/plink2.rin_red_blood_cell_distribution_width.glm.linear.done"
    File TAOK1_conditioned_STR_27842010_snp_results = "~{base_location}/association/results/mean_platelet_volume/plink_snp_conditional/chr17_25885931_29975306_STR_27842010__ISNP__ASNP/plink2.rin_mean_platelet_volume.glm.linear.done"
    Array[File] peaks_for_1ef = peaks_for_1ef_
    Array[File] finemapping_regions = finemapping_regions_
    # Indexed first by phenotype then ethnicity
    Array[Array[File]] pheno_to_ethnic_to_str_gwas_results = pheno_to_ethnic_to_str_gwas_results_
    # indexed first by ethnicity then phenotype
    Array[Array[File]] ethnic_to_pheno_to_str_gwas_results = ethnic_to_pheno_to_str_gwas_results_ 

    Array[File] finemapping_first_pass_dfs = finemapping_first_pass_dfs_ 
    Array[File] susie_min_abs_corrs = susie_min_abs_corrs_
    Array[File] finemapping_followup_dfs = finemapping_followup_dfs_ # bigger for mean platelet volume, so needs consistent filtering

    Array[File] platelet_count_snp_macs = platelet_count_snp_macs_

  }
}
