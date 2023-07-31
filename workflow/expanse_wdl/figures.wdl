version 1.0

import "expanse_tasks.wdl"
import "expanse_files.wdl"
import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"
import "../finemapping_wdl/post_finemapping_workflow.wdl"
import "../gwas_wdl/prep_samples_and_phenotype_workflow.wdl"

# also includes tables
workflow figures {

  String script_dir = "."

  # inputs, but written as a task for reusability
  call expanse_files.files
  call gwas_tasks.phenotype_names

  ##### some common task inputs
  call expanse_tasks.extract_field as sc_white_brits { input :
    script_dir = script_dir,
    id = 22006
  }

  call gwas_tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = script_dir,
    sc = sc_white_brits.data
  }

  Int total_bilirubin_idx = phenotype_names.idxs["total_bilirubin"]

  ##### generate figure 1 

  # 1a provided by Rahel

  # 1b WDL not working
#  this crashes because of a missing header in the VCF ##command=Hipstr
#  scatter (chrom in range(22)) {
#    call gwas_tasks.imputed_str_locus_summary as imputed_str_locus_summaries { input :
#      script_dir = script_dir,
#      vcf = str_vcfs[chrom],
#      qced_white_brits = unrelated_white_brits_sample_list
#    }
#  }
#
#  call gwas_tasks.str_multiallelicness_distro as fig_1b { input:
#    script_dir = script_dir,
#    thresh = 0.01,
#    chrom_locus_summaries = imputed_str_locus_summaries.out
#  }

  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    call gwas_tasks.generate_peaks { input :
      script_dir = script_dir,
      snp_assoc_results = files.snp_gwas_results[phenotype_idx],
      str_assoc_results = files.str_gwas_results[phenotype_idx],
      phenotype = phenotype_names.n[phenotype_idx],
      spacing = "250000",
      thresh = "5e-8"
    }
  }

  call gwas_tasks.generate_peaks as bilirubin_overview_manhattan_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = files.snp_gwas_results[total_bilirubin_idx],
    str_assoc_results = files.str_gwas_results[total_bilirubin_idx],
    phenotype = "total_bilirubin", 
    spacing = "10000000",
    thresh = "5e-8"
  }

  call gwas_tasks.overview_manhattan as bilirubin_overview_manhattan { input :
    script_dir = script_dir,
    phenotype_name = "total_bilirubin",
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[total_bilirubin_idx],
    snp_gwas_results = files.snp_gwas_results[total_bilirubin_idx],
    peaks = bilirubin_overview_manhattan_peaks.peaks,
    ext = "png"
  }

  Int platelet_count_idx = phenotype_names.idxs["platelet_count"]
  call gwas_tasks.generate_peaks as platelet_count_overview_manhattan_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = files.snp_gwas_results[platelet_count_idx],
    str_assoc_results = files.str_gwas_results[platelet_count_idx],
    phenotype = "platelet_count", 
    spacing = "10000000",
    thresh = "5e-8"
  }

  call gwas_tasks.overview_manhattan as platelet_count_overview_manhattan { input :
    script_dir = script_dir,
    phenotype_name = "platelet_count",
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[platelet_count_idx],
    snp_gwas_results = files.snp_gwas_results[platelet_count_idx],
    peaks = platelet_count_overview_manhattan_peaks.peaks,
    ext = "png"
  }

  call gwas_tasks.summarize_peaks as fig_1ef { input :
    script_dir = script_dir,
    phenotype_names = phenotype_names.n,
    peak_files = generate_peaks.peaks
  }

  ### haven't WDLed supp fig 1

  ### supp table 1 is manual

  ##### supp fig 2
  call gwas_tasks.compare_bili_to_UKBB as supp_fig_2 { input :
    script_dir = script_dir,
    pan_ukbb = files.pan_ukbb,
    plink_snp = files.snp_gwas_results[total_bilirubin_idx]
  }

  ##### TODO supp table 2

  ##### Supp table 3 not in WDL

  # supp figs 3,4,5,7,8,11,12,13  
  call post_finemapping_workflow.post_finemapping { input :
    script_dir = ".",
    first_pass_dfs = files.finemapping_first_pass_dfs,
    susie_all_min_abs_corrs = files.susie_min_abs_corrs,
    followup_dfs = files.finemapping_followup_dfs
  }

  # TODO supp fig 6, 10

  # supp fig 9
  call finemapping_tasks.susie_finemap_venn_diagram { input :
    script_dir = script_dir,
    first_pass_dfs = files.finemapping_first_pass_dfs
  }

  ######### generate supplementary tables 4 and 5
  call finemapping_tasks.str_tables_for_paper { input :
    script_dir = script_dir,
    flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos,
    str_hg19_pos_bed = files.str_hg19_pos_bed,
    str_hg38_pos_bed = files.str_hg38_pos_bed,
    repeat_units_table = files.repeat_units_table,
    intersects_gene = files.intersects_gene_annotation,
    intersects_exon = files.intersects_exon_annotation,
    intersects_CDS = files.intersects_CDS_annotation,
    intersects_five_prime_UTR = files.intersects_five_prime_UTR_annotation,
    intersects_three_prime_UTR = files.intersects_three_prime_UTR_annotation,
    intersects_UTR = files.intersects_UTR_annotation,
    phenotype_names = phenotype_names.n,
    assocs = files.str_gwas_results,
    black_assocs = files.ethnic_to_pheno_to_str_gwas_results[0],
    south_asian_assocs = files.ethnic_to_pheno_to_str_gwas_results[1],
    chinese_assocs = files.ethnic_to_pheno_to_str_gwas_results[2],
    irish_assocs = files.ethnic_to_pheno_to_str_gwas_results[3],
    white_other_assocs = files.ethnic_to_pheno_to_str_gwas_results[4],
    first_pass_finemapping_dfs = files.finemapping_first_pass_dfs,
    followup_finemapping_dfs = files.finemapping_followup_dfs,
  }

  ######### generate figure 2
  call finemapping_tasks.generate_QTL_table as eQTL_table { input :
    script_dir = script_dir,
    confidently_finemapped_table = str_tables_for_paper.confidently_finemapped_strs_for_paper,
    all_QTL_results = files.all_eQTL_results
  }

  call finemapping_tasks.generate_QTL_table as meQTL_table { input :
    script_dir = script_dir,
    confidently_finemapped_table = str_tables_for_paper.confidently_finemapped_strs_for_paper,
    all_QTL_results = files.all_meQTL_results,
    methylation = true
  }

  call finemapping_tasks.graph_main_hits { input :
    script_dir = script_dir,
    hits_table = str_tables_for_paper.singly_finemapped_strs_for_paper,
    eQTL_table = eQTL_table.QTL_table,
    meQTL_table = meQTL_table.QTL_table, #files.confidently_finemapped_methylation_results,
    closest_gene_annotations = files.closest_gene_annotation
  }

  ######## generate figure 3, supp figure 14 and supp table 6
  call finemapping_tasks.concordance_in_other_ethnicities { input :
    script_dir = script_dir,
    confidently_finemapped_STRs_df = post_finemapping.confidently_finemapped_STRs,
    first_pass_dfs = files.finemapping_first_pass_dfs
  }

  ##### Supp fig 15 provided by Xuan

  ######## generate supp fig 16
  call finemapping_tasks.generate_enrichments_table { input :
    script_dir = script_dir,
    flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos,
    str_loci = files.str_loci,
    repeat_units_table = files.repeat_units_table,
    eSTR_table = files.eSTR_table,
    gencode = files.gencode,
    phenotypes = phenotype_names.n,
    str_assocs = files.str_gwas_results,
    confidently_finemapped_STRs_df = post_finemapping.confidently_finemapped_STRs,

    intersects_transcript_support_2 = files.intersects_transcript_support_2_annotation,
    intersects_protein_coding_CDS_support_2 = files.intersects_protein_coding_CDS_support_2_annotation,
    intersects_protein_coding_five_prime_UTR_support_2 = files.intersects_protein_coding_five_prime_UTR_support_2_annotation,
    intersects_protein_coding_three_prime_UTR_support_2 = files.intersects_protein_coding_three_prime_UTR_support_2_annotation,
    intersects_protein_coding_UTR_support_2 = files.intersects_protein_coding_UTR_support_2_annotation,
    closest_downstream_protein_coding_exon_support_2 = files.closest_downstream_protein_coding_exon_support_2_annotation,
    closest_upstream_protein_coding_exon_support_2 = files.closest_upstream_protein_coding_exon_support_2_annotation,
    closest_downstream_protein_coding_gene = files.closest_downstream_protein_coding_gene_annotation,
    closest_upstream_protein_coding_gene = files.closest_upstream_protein_coding_gene_annotation,
  }

  call finemapping_tasks.calc_enrichments { input :
    script_dir = script_dir,
    enrichment_df = generate_enrichments_table.table
  }

  call finemapping_tasks.graph_enrichments { input :
    script_dir = script_dir,
    enrichment_stats = calc_enrichments.enrichment_stats
  }

  ###### supp fig 17 for Yan

  ####### Supp fig 18 SLC2A2 

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon4 { input :
    script_dir = script_dir,
    individual_tsv = files.Liver_SLC2A2_exon4_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_18a { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2_exon_4",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon4.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon6 { input :
    script_dir = script_dir,
    individual_tsv = files.Liver_SLC2A2_exon6_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_18b { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2_exon_6",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon6.out],
    dosage_threshold = 5
  }

  call expanse_tasks.extract_field as bilirubin_sc { input :
    script_dir = script_dir,
    id = 30840,
  }
  
  call expanse_tasks.extract_field as bilirubin_covariate_sc { input :
    script_dir = script_dir,
    id = 30842,
  }

  call expanse_tasks.extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
  }

  call expanse_tasks.extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }

  call gwas_tasks.load_shared_covars { input:
    script_dir = script_dir,
    fam_file = files.fam_file,
    sc_pcs = pcs.data,
    sc_assessment_ages = assessment_ages.data
  }

  call gwas_tasks.load_continuous_phenotype as bilirubin_all_white_brits { input:
    script_dir = script_dir,
    sc = bilirubin_sc.data,
    qced_sample_list = white_brits_sample_list.data,
    assessment_ages_npy = load_shared_covars.assessment_ages,
    categorical_covariate_names = ["total_bilirubin_aliquot"],
    categorical_covariate_scs = [bilirubin_covariate_sc.data],
  }

  call gwas_tasks.transform_trait_values as transformed_bilirubin { input:
    script_dir = script_dir,
    pheno_data = bilirubin_all_white_brits.data,
    samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[total_bilirubin_idx][0],
    is_binary = false
  }

  call gwas_tasks.str_spot_test as SLC2A2_assoc { input:
    script_dir = script_dir,
    str_vcf = files.str_vcfs[2],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = bilirubin_all_white_brits.data,
    transformed_phenotype = transformed_bilirubin.data,
    all_samples_list = files.all_samples_list,
    is_binary = false,
    chrom = 3,
    pos = 170727702,
    phenotype_name = "total_bilirubin",
  }

  call gwas_tasks.locus_plot as supp_fig_18c { input :
    script_dir = script_dir,
    chrom = 3,
    pos = 170727702,
    phenotype_name = "total_bilirubin",
    unit = phenotype_names.unit[total_bilirubin_idx],
    assoc_results = [SLC2A2_assoc.data],
    dosage_fraction_threshold = 0.001
  }

  call prep_samples_and_phenotype_workflow.prep_samples_and_phenotype as prep_total_bilirubin { input : 
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",
    phenotype_name = "total_bilirubin",
    categorical_covariate_names = ["total_bilirubin_aliquot"],
    categorical_covariate_scs = [bilirubin_covariate_sc.data],
    is_binary = false,
    is_zero_one_neg_nan = false,
    date_of_most_recent_first_occurrence_update = files.date_of_most_recent_first_occurrence_update ,
    fam_file = files.fam_file,
    withdrawn_sample_list = files.withdrawn_sample_list,
    kinship = files.kinship,
    sc_white_brits = sc_white_brits.data,
    sc_ethnicity_self_report = ethnicity_self_report.data,
    sc_sex_aneuploidy = sex_aneuploidy.data,
    sc_genetic_sex = genetic_sex.data,
    sc_reported_sex = reported_sex.data,
    sc_kinship_count = kinship_count.data,
    sc_assessment_ages = assessment_ages.data,
    sc_pcs = pcs.data,
    sc_year_of_birth = year_of_birth.data,
    sc_month_of_birth = month_of_birth.data,
    sc_date_of_death = date_of_death.data,
    sc_phenotype = bilirubin_sc.data,
    cached_unrelated_samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[total_bilirubin_idx],
    cached_shared_covars = files.shared_covars,
  }

  region slc2a2_bounds = {
    "chrom": 3,
    "start": 170424098,
    "end": 170886547
  }

  call gwas_tasks.prep_conditional_input as prep_slc2a2_conditional_input { input :
		script_dir = script_dir,
    all_samples = files.all_samples_list,
    chrom = slc2a2_bounds.chrom,
    str_vcf = files.str_vcfs[2],
    strs = [170727702],
    snps = []
  }

  call gwas_tasks.regional_my_str_gwas as conditional_slc2a2_str_gwas { input :
		script_dir = script_dir,
		str_vcf = files.str_vcfs[2],
		shared_covars = prep_total_bilirubin.shared_covars,
		untransformed_phenotype = prep_total_bilirubin.pheno_data[0],
		transformed_phenotype = prep_total_bilirubin.transformed_trait_values[0],
    conditional_covars = prep_slc2a2_conditional_input.data,
		all_samples_list = files.all_samples_list,
		is_binary = false,
		binary_type = "linear", # won't be used if not binary
		bounds = slc2a2_bounds,
		phenotype_name = "total_bilirubin"
  }

  call gwas_tasks.prep_plink_input as prep_slc2a2_plink_input { input : 
    script_dir = script_dir,
    shared_covars = prep_total_bilirubin.shared_covars,
    shared_covar_names = prep_total_bilirubin.shared_covar_names,
		transformed_phenotype = prep_total_bilirubin.transformed_trait_values[0],
    pheno_covar_names = prep_total_bilirubin.pheno_covar_names[0],
    conditional_genotypes = prep_slc2a2_conditional_input.data,
    conditional_covar_names = prep_slc2a2_conditional_input.varnames,
    is_binary = false,
    binary_type = "linear",
    phenotype_name = "total_bilirubin"
  }

  call gwas_tasks.plink_snp_association as conditional_slc2a2_snp_gwas { input :
    script_dir = script_dir,
    plink_command = "plink2",
    imputed_snp_p_file = files.imputed_snp_pfiles[2],
    pheno_data = prep_slc2a2_plink_input.data,
    chrom = slc2a2_bounds.chrom,
    start = slc2a2_bounds.start,
    end = slc2a2_bounds.end,
    phenotype_name = "total_bilirubin",
    binary_type = "linear"
  }

  call gwas_tasks.manhattan as supp_fig_18d { input :
    script_dir = script_dir,
    phenotype_name = "total_bilirubin",
    unit = phenotype_names.unit[total_bilirubin_idx],
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[total_bilirubin_idx],
    snp_gwas_results = files.snp_gwas_results[total_bilirubin_idx],
    bounds = slc2a2_bounds,
    conditioned_STRs = [170727702],
    conditional_str_results = conditional_slc2a2_str_gwas.data,
    conditional_snp_results = conditional_slc2a2_snp_gwas.data,
    ext = "png",
  }

  ####### generate figure 4 and supp figs 19 - 21 

  call gwas_tasks.cbl_imperfection_ld { input :
    script_dir = script_dir,
    bgen_chr11 = files.imputed_snp_bgens[10],
    all_samples_file = files.all_samples_list,
    phenotype_samples_file = files.unrelated_samples_for_pheno_for_ethnicity[platelet_count_idx][0]
  }

  # fig 4a,b from WGS
#  call gwas_tasks.fig_4a { input :
#    script_dir = script_dir,
#    all_samples_list = files.all_samples_list,
#    white_brits_sample_list = white_brits_sample_list.data,
#    black_sample_list = files.qced_samples_for_ethnicity[1],
#    south_asian_sample_list = files.qced_samples_for_ethnicity[2],
#    chinese_sample_list =  files.qced_samples_for_ethnicity[3],
#    str_vcf_chr_11 = files.str_vcfs[10],
#    specific_alleles = files.specific_alleles,
#  }

  # supp fig 19c-d
  call expanse_tasks.extract_field as platelet_count_sc { input :
    script_dir = script_dir,
    id = 30080,
  }
  
  call expanse_tasks.extract_field as platelet_count_covariate_sc { input :
    script_dir = script_dir,
    id = 30083,
  }

  # all, not qced or subset to unrelated, the sample list for this won't be used, only the data
  call gwas_tasks.load_continuous_phenotype as platelet_count_all_white_brits { input:
    script_dir = script_dir,
    sc = platelet_count_sc.data,
    qced_sample_list = white_brits_sample_list.data,
    assessment_ages_npy = load_shared_covars.assessment_ages,
    categorical_covariate_names = ["platelet_count_device_id"],
    categorical_covariate_scs = [platelet_count_covariate_sc.data],
  }

  call gwas_tasks.transform_trait_values as transformed_platelet_count { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[platelet_count_idx][0],
    is_binary = false
  }

  call gwas_tasks.str_spot_test as imputed_CBL_platelet_count_assoc { input:
    script_dir = script_dir,
    str_vcf = files.str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data,
    transformed_phenotype = transformed_platelet_count.data,
    all_samples_list = files.all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.locus_plot as supp_fig_19c { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [imputed_CBL_platelet_count_assoc.data],
    dosage_fraction_threshold = 0.001
  }

  Int platelet_crit_idx = phenotype_names.idxs["platelet_crit"]
  call expanse_tasks.extract_field as platelet_crit_sc { input :
    script_dir = script_dir,
    id = 30080,
  }
  
  call expanse_tasks.extract_field as platelet_crit_covariate_sc { input :
    script_dir = script_dir,
    id = 30083,
  }

  # all, not qced or subset to unrelated, the sample list for this won't be used, only the data
  call gwas_tasks.load_continuous_phenotype as platelet_crit_all_white_brits { input:
    script_dir = script_dir,
    sc = platelet_crit_sc.data,
    qced_sample_list = white_brits_sample_list.data,
    assessment_ages_npy = load_shared_covars.assessment_ages,
    categorical_covariate_names = ["platelet_crit_device_id"],
    categorical_covariate_scs = [platelet_crit_covariate_sc.data],
  }

  call gwas_tasks.transform_trait_values as transformed_platelet_crit { input:
    script_dir = script_dir,
    pheno_data = platelet_crit_all_white_brits.data,
    samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[platelet_crit_idx][0],
    is_binary = false
  }

  call gwas_tasks.str_spot_test as imputed_CBL_platelet_crit_assoc { input:
    script_dir = script_dir,
    str_vcf = files.str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_crit_all_white_brits.data,
    transformed_phenotype = transformed_platelet_crit.data,
    all_samples_list = files.all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_crit",
  }

  call gwas_tasks.locus_plot as supp_fig_19d { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_crit",
    unit = "10^9 cells/L",
    assoc_results = [imputed_CBL_platelet_crit_assoc.data],
    dosage_fraction_threshold = 0.001
  }

  # fig 4 c-e
  call gwas_tasks.manhattan as fig_4ce { input:
    script_dir = script_dir,
    phenotype_name = "platelet_crit",
    unit = phenotype_names.unit[platelet_crit_idx],
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[platelet_crit_idx],
    snp_gwas_results = files.snp_gwas_results[platelet_crit_idx],
    bounds = {
      "chrom": 11,
      "start": 118447267,
      "end": 119339135
    },
    conditioned_STRs = [119077000],
    conditioned_imputed_SNPs = ["119080037_A_G"],
    conditional_str_results = files.CBL_conditioned_STR_19077000_SNP_119080037_A_G_str_results,
    conditional_snp_results = files.CBL_conditioned_STR_19077000_SNP_119080037_A_G_snp_results,
    ext = "png"
  }

  call gwas_tasks.manhattan as fig_4d { input:
    script_dir = script_dir,
    phenotype_name = "platelet_crit",
    unit = phenotype_names.unit[platelet_crit_idx],
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[platelet_crit_idx],
    snp_gwas_results = files.snp_gwas_results[platelet_crit_idx],
    bounds = {
      "chrom": 11,
      "start": 118447267,
      "end": 119339135
    },
    conditioned_STRs = [],
    conditioned_imputed_SNPs = ["119080037_A_G"],
    conditional_str_results = files.CBL_conditioned_SNP_119080037_A_G_str_results,
    conditional_snp_results = files.CBL_conditioned_SNP_119080037_A_G_snp_results,
    ext = "png"
  }

  # fig 4f from WGS
  # use precomputed sample lists
#  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_not_SNP_samples { input:
#    script_dir = script_dir,
#    pheno_data = platelet_count_all_white_brits.data,
#    samples_for_phenotype = files.unrelated_samples_CBL_hom_not_begin_C_T_snp,
#    is_binary = false
#  }
#
#  call gwas_tasks.str_spot_test as CBL_hom_not_SNP_assoc { input:
#    script_dir = script_dir,
#    str_vcf = files.str_vcfs[10],
#    shared_covars = load_shared_covars.shared_covars, 
#    untransformed_phenotype = platelet_count_all_white_brits.data,
#    transformed_phenotype = platelet_count_transformed_CBL_hom_not_SNP_samples.data, 
#    all_samples_list = files.all_samples_list,
#    is_binary = false,
#    chrom = 11,
#    pos = 119077000,
#    phenotype_name = "platelet_count",
#  }
#
#  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_SNP_samples { input:
#    script_dir = script_dir,
#    pheno_data = platelet_count_all_white_brits.data,
#    samples_for_phenotype = files.unrelated_samples_CBL_hom_begin_C_T_snp,
#    is_binary = false
#  }
#   
#  call gwas_tasks.str_spot_test as CBL_hom_SNP_assoc { input:
#    script_dir = script_dir,
#    str_vcf = files.str_vcfs[10],
#    shared_covars = load_shared_covars.shared_covars, 
#    untransformed_phenotype = platelet_count_all_white_brits.data, 
#    transformed_phenotype = platelet_count_transformed_CBL_hom_SNP_samples.data,
#    all_samples_list = files.all_samples_list,
#    is_binary = false,
#    chrom = 11,
#    pos = 119077000,
#    phenotype_name = "platelet_count",
#  }
#
#  call gwas_tasks.locus_plot as fig_4f { input:
#    script_dir = script_dir,
#    chrom = 11,
#    pos = 119077000,
#    phenotype_name = "platelet_count",
#    unit = "10^9 cells/L",
#    assoc_results = [CBL_hom_not_SNP_assoc.data, CBL_hom_SNP_assoc.data],
#    group_names = ["homozygous (CGG)n", "homozygous CGGTGG(CGG)m"],
#    dosage_threshold = 200,
#  }

  ### plot with data from GTEx Yang

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_CBL_gtex_expression { input :
    script_dir = script_dir,
    individual_tsv = files.CBL_gtex_expression,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "TPM(expression)"
  }

  call gwas_tasks.locus_plot as fig_4g { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "CBL_expression",
    unit = "TPM",
    data_tsvs = [summarized_CBL_gtex_expression.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_CBL_geuvadis_european_expression { input :
    script_dir = script_dir,
    individual_tsv = files.CBL_geuvadis_expression,
    length_sum_column_name = "strgt_sum",
    trait_column_name = "Expr",
    sep = ",",
    subset_field = "superpop",
    subset_value = "EUR"
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_CBL_geuvadis_african_expression { input :
    script_dir = script_dir,
    individual_tsv = files.CBL_geuvadis_expression,
    length_sum_column_name = "strgt_sum",
    trait_column_name = "Expr",
    sep = ",",
    subset_field = "superpop",
    subset_value = "AFR"
  }

  call gwas_tasks.locus_plot as fig_4h { input: 
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "CBL_expression",
    unit = "RPKM",
    data_tsvs = [summarized_CBL_geuvadis_european_expression.out, summarized_CBL_geuvadis_african_expression.out],
    group_names = ["European", "African"],
    dosage_threshold = 5
  }

  ### Supp fig 19a-b from WGS
  ### Supp fig 20a from WGS
  ### Supp fig 20b not in WDL
  ### Supp fig 21 from Melissa

  ### TODO supp fig 22a
  #### Supp fig 22b from WGS

  #### Supp fig 23a from WGS
  ### GWAS for tstat for TAOK1 for Fig 23b

  call expanse_tasks.extract_field as ethnicity_self_report { input :
    script_dir = script_dir,
    id = 21000
  }

  call expanse_tasks.extract_field as sex_aneuploidy { input:
    script_dir = script_dir,
    id = 22019
  }

  call expanse_tasks.extract_field as genetic_sex { input:
    script_dir = script_dir,
    id = 22001
  }

  call expanse_tasks.extract_field as reported_sex { input:
    script_dir = script_dir,
    id = 31
  }

  call expanse_tasks.extract_field as kinship_count { input:
    script_dir = script_dir,
    id = 22021
  }

  call expanse_tasks.extract_field as year_of_birth { input :
    script_dir = script_dir,
    id = 34
  }

  call expanse_tasks.extract_field as month_of_birth { input :
    script_dir = script_dir,
    id = 52
  }

  call expanse_tasks.extract_field as date_of_death { input :
    script_dir = script_dir,
    id = 40000
  }

  Int mean_platelet_volume_idx = phenotype_names.idxs["mean_platelet_volume"]

  call expanse_tasks.extract_field as mean_platelet_volume_sc { input :
    script_dir = script_dir,
    id = phenotype_names.ID[mean_platelet_volume_idx]
  }

	call expanse_tasks.extract_field as mean_platelet_volume_device_id { input :
		script_dir = script_dir,
		id = phenotype_names.ID[mean_platelet_volume_idx] + 3
  }

  call prep_samples_and_phenotype_workflow.prep_samples_and_phenotype as prep_mean_platelet_volume { input : 
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",
    phenotype_name = "mean_platelet_volume",
    categorical_covariate_names = ["mean_platelet_volume_device_id"],
    categorical_covariate_scs = [mean_platelet_volume_device_id.data],
    is_binary = false,
    is_zero_one_neg_nan = false,
    date_of_most_recent_first_occurrence_update = files.date_of_most_recent_first_occurrence_update ,
    fam_file = files.fam_file,
    withdrawn_sample_list = files.withdrawn_sample_list,
    kinship = files.kinship,
    sc_white_brits = sc_white_brits.data,
    sc_ethnicity_self_report = ethnicity_self_report.data,
    sc_sex_aneuploidy = sex_aneuploidy.data,
    sc_genetic_sex = genetic_sex.data,
    sc_reported_sex = reported_sex.data,
    sc_kinship_count = kinship_count.data,
    sc_assessment_ages = assessment_ages.data,
    sc_pcs = pcs.data,
    sc_year_of_birth = year_of_birth.data,
    sc_month_of_birth = month_of_birth.data,
    sc_date_of_death = date_of_death.data,
    sc_phenotype = mean_platelet_volume_sc.data,
    cached_unrelated_samples_for_phenotype = files.unrelated_samples_for_pheno_for_ethnicity[mean_platelet_volume_idx],
    cached_shared_covars = files.shared_covars,
  }

  region taok1_bounds = {
    "chrom": 17,
    "start": 25885931,
    "end": 29975306
  }

  call gwas_tasks.regional_my_str_gwas as taok1_str_gwas { input :
		script_dir = script_dir,
		str_vcf = files.str_vcfs[16],
		shared_covars = prep_mean_platelet_volume.shared_covars,
		untransformed_phenotype = prep_mean_platelet_volume.pheno_data[0],
		transformed_phenotype = prep_mean_platelet_volume.transformed_trait_values[0],
		all_samples_list = files.all_samples_list,
		is_binary = false,
		binary_type = "linear", # won't be used if not binary
		bounds = taok1_bounds,
		phenotype_name = "mean_platelet_volume"
  }

  call gwas_tasks.prep_conditional_input as prep_taok1_conditional_input { input :
		script_dir = script_dir,
    all_samples = files.all_samples_list,
    chrom = taok1_bounds.chrom,
    str_vcf = files.str_vcfs[16],
    strs = [27842010],
    snps = []
  }

  call gwas_tasks.regional_my_str_gwas as conditional_taok1_str_gwas { input :
		script_dir = script_dir,
		str_vcf = files.str_vcfs[16],
		shared_covars = prep_mean_platelet_volume.shared_covars,
		untransformed_phenotype = prep_mean_platelet_volume.pheno_data[0],
		transformed_phenotype = prep_mean_platelet_volume.transformed_trait_values[0],
    conditional_covars = prep_taok1_conditional_input.data,
		all_samples_list = files.all_samples_list,
		is_binary = false,
		binary_type = "linear", # won't be used if not binary
		bounds = taok1_bounds,
		phenotype_name = "mean_platelet_volume"
  }

  call gwas_tasks.manhattan as supp_fig_23b { input :
    script_dir = script_dir,
    phenotype_name = "mean_platelet_volume",
    unit = phenotype_names.unit[mean_platelet_volume_idx],
    chr_lens = files.chr_lens,
    str_gwas_results = taok1_str_gwas.data,
    snp_gwas_results = files.snp_gwas_results[mean_platelet_volume_idx],
    bounds = taok1_bounds,
    conditioned_STRs = [27842016],
    conditional_str_results = conditional_taok1_str_gwas.data,
    conditional_snp_results = files.TAOK1_conditioned_STR_27842010_snp_results,
    ext = "png",
    use_tstat = true
  }

  ### Supp fig 24a,c from WGS

  # from snakemake cache, don't need the WDL
#  call gwas_tasks.manhattan as supp_fig_24b { input:
#    script_dir = script_dir,
#    phenotype_name = "haematocrit",
#    unit = phenotype_names.unit[haematorcrit_idx], #
#    chr_lens = files.chr_lens,
#    str_gwas_results = files.str_gwas_results[haematorcrit_idx],
#    snp_gwas_results = files.snp_gwas_results[haematorcrit_idx],
#    bounds = {
#      "chrom": 11, #
#      "start": 118447267, #
#      "end": 119339135 #
#    },
#    conditioned_STRs = [119077000], #
#    conditioned_imputed_SNPs = ["119080037_A_G"], #
#    conditional_str_results = files.CBL_conditioned_SNP_119080037_A_G_str_results, #
#    conditional_snp_results = files.CBL_conditioned_SNP_119080037_A_G_snp_results, #
#    ext = "png"
#  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_RHOT1_geuvadis_european_expression { input :
    script_dir = script_dir,
    individual_tsv = files.RHOT1_geuvadis_expression,
    length_sum_column_name = "strgt_sum",
    trait_column_name = "Expr",
    sep = ",",
    subset_field = "superpop",
    subset_value = "EUR"
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_RHOT1_geuvadis_african_expression { input :
    script_dir = script_dir,
    individual_tsv = files.RHOT1_geuvadis_expression,
    length_sum_column_name = "strgt_sum",
    trait_column_name = "Expr",
    sep = ",",
    subset_field = "superpop",
    subset_value = "AFR"
  }

  call gwas_tasks.locus_plot as supp_fig_25b { input: 
    script_dir = script_dir,
    chrom = 17,
    pos = 30469471,
    phenotype_name = "RHOT1_expression",
    unit = "RPKM",
    data_tsvs = [summarized_RHOT1_geuvadis_european_expression.out, summarized_RHOT1_geuvadis_african_expression.out],
    group_names = ["European", "African"],
    dosage_threshold = 5
  }

  Int red_blood_cell_distribution_width_idx = phenotype_names.idxs["red_blood_cell_distribution_width"]
  call gwas_tasks.manhattan as supp_fig_25d { input:
    script_dir = script_dir,
    phenotype_name = "red_blood_cell_distribution_width",
    unit = phenotype_names.unit[red_blood_cell_distribution_width_idx],
    chr_lens = files.chr_lens,
    str_gwas_results = files.str_gwas_results[red_blood_cell_distribution_width_idx],
    snp_gwas_results = files.snp_gwas_results[red_blood_cell_distribution_width_idx],
    bounds = {
      "chrom": 17,
      "start": 30287357,
      "end": 30595028
    },
    conditioned_STRs = [30469467],
    conditional_str_results = files.RHOT1_conditioned_str_results,
    conditional_snp_results = files.RHOT1_conditioned_snp_results,
    ext = "png"
  }

  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    call gwas_tasks.reformat_my_str_gwas_table_for_publication { input :
      script_dir = script_dir,
      phenotype = phenotype_names.n[phenotype_idx],
      my_str_gwas = files.str_gwas_results[phenotype_idx],
      flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos,
      str_hg19_pos_bed = files.str_hg19_pos_bed,
      str_hg38_pos_bed = files.str_hg38_pos_bed,
      repeat_units_table = files.repeat_units_table
    }
  }

#  scatter (phenotype_idx in range(length(phenotype_names.n))) {
#    scatter (ethnicity_idx in range(5)) {
#      call gwas_tasks.reformat_my_str_gwas_table_for_publication as reformat_my_ethnic_str_gwas_table_for_publication { input :
#        script_dir = script_dir,
#        phenotype = phenotype_names.n[phenotype_idx],
#        my_str_gwas = files.pheno_to_ethnic_to_str_gwas_results[phenotype_idx][ethnicity_idx],
#        flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos,
#        str_hg19_pos_bed = files.str_hg19_pos_bed,
#        str_hg38_pos_bed = files.str_hg38_pos_bed,
#        repeat_units_table = files.repeat_units_table
#      }
#    }
#  }

  output {
    Float cbl_imperfection_ld_r2 = cbl_imperfection_ld.r2

#    File fig_1b_svg_out = fig_1b.svg
#    File fig_1b_png_out = fig_1b.png
    File fig_1c = bilirubin_overview_manhattan.plot
    File fig_1d = platelet_count_overview_manhattan.plot
    File fig_1e_svg = fig_1ef.barplot_svg
    File fig_1e_png = fig_1ef.barplot_png
    File fig_1f_svg = fig_1ef.heatmap_svg
    File fig_1f_png = fig_1ef.heatmap_png

    File fig_2_png = graph_main_hits.png
    File fig_2_svg = graph_main_hits.svg

    File fig_3_png = concordance_in_other_ethnicities.nonwhite_replication_png
    File fig_3_svg = concordance_in_other_ethnicities.nonwhite_replication_svg

#    File fig_4a_svg = fig_4a.svg
#    File fig_4a_png = fig_4a.png
#    File fig_4b_svg = fig_4b.svg
#    File fig_4b_png = fig_4b.png
    File fig_4ce_out = fig_4ce.plot
    File fig_4d_out = fig_4d.plot
#    File fig_4f_svg = fig_4f.svg
#    File fig_4f_png = fig_4f.png
    File fig_4g_svg = fig_4g.svg
    File fig_4g_png = fig_4g.png
    File fig_4h_svg = fig_4h.svg
    File fig_4h_png = fig_4h.png

    File supp_fig_2_png = supp_fig_2.png
    File supp_fig_2_svg = supp_fig_2.svg

    File stat_statements = post_finemapping.stat_statements
    File supp_fig_3_png = post_finemapping.cs_min_abs_corrs_png
    File supp_fig_3_svg = post_finemapping.cs_min_abs_corrs_svg
    File supp_fig_4_png = post_finemapping.susie_alpha_v_pip_png
    File supp_fig_4_svg = post_finemapping.susie_alpha_v_pip_svg
    File supp_fig_5a_png = post_finemapping.susie_alpha_histogram_png
    File supp_fig_5a_svg = post_finemapping.susie_alpha_histogram_svg
    File supp_fig_5b_png = post_finemapping.finemap_pip_histogram_png
    File supp_fig_5b_svg = post_finemapping.finemap_pip_histogram_svg

    # Supp fig 6 TODO

    File supp_fig_7_png = post_finemapping.susie_cs_finemap_total_pips_png
    File supp_fig_7_svg = post_finemapping.susie_cs_finemap_total_pips_svg
    File supp_fig_8a_png = post_finemapping.finemap_v_susie_consistency_STR_png
    File supp_fig_8a_svg = post_finemapping.finemap_v_susie_consistency_STR_svg
    File supp_fig_8b_png = post_finemapping.finemap_v_susie_consistency_SNP_png
    File supp_fig_8b_svg = post_finemapping.finemap_v_susie_consistency_SNP_svg

    File supp_fig_9a_svg = susie_finemap_venn_diagram.str_svg
    File supp_fig_9a_png = susie_finemap_venn_diagram.str_png
    File supp_fig_9b_svg = susie_finemap_venn_diagram.snp_svg
    File supp_fig_9b_png = susie_finemap_venn_diagram.snp_png

    # used for deciding confidently finemap
    File supp_fig_11ab_png = post_finemapping.susie_best_guess_png
    File supp_fig_11ab_svg = post_finemapping.susie_best_guess_svg
    File supp_fig_12a_png = post_finemapping.finemap_p_thresh_png
    File supp_fig_12a_svg = post_finemapping.finemap_p_thresh_svg
    File supp_fig_12b_png = post_finemapping.finemap_mac_png
    File supp_fig_12b_svg = post_finemapping.finemap_mac_svg
    File supp_fig_12c_png = post_finemapping.finemap_prior_std_derived_png
    File supp_fig_12c_svg = post_finemapping.finemap_prior_std_derived_svg
    File supp_fig_12d_png = post_finemapping.finemap_total_prob_png
    File supp_fig_12d_svg = post_finemapping.finemap_total_prob_svg
    File supp_fig_12e_png = post_finemapping.finemap_conv_tol_png
    File supp_fig_12e_svg = post_finemapping.finemap_conv_tol_svg

    # too conservative
    File supp_fig_12f_png = post_finemapping.finemap_prior_std_low_png
    File supp_fig_12f_svg = post_finemapping.finemap_prior_std_low_svg
    File supp_fig_13a_png = post_finemapping.susie_ratio_png
    File supp_fig_13a_svg = post_finemapping.susie_ratio_svg
    File supp_fig_13b_png = post_finemapping.finemap_ratio_png
    File supp_fig_13b_svg = post_finemapping.finemap_ratio_svg

    File supp_fig_14_png = concordance_in_other_ethnicities.white_replication_png
    File supp_fig_14_svg = concordance_in_other_ethnicities.white_replication_svg

    File supp_fig_16ab_svg = graph_enrichments.regions_svg
    File supp_fig_16ab_png = graph_enrichments.regions_png
    File supp_fig_16cd_svg = graph_enrichments.repeats_svg
    File supp_fig_16cd_png = graph_enrichments.repeats_png

    File supp_fig_18a_svg = supp_fig_18a.svg
    File supp_fig_18a_png = supp_fig_18a.png
    File supp_fig_18b_svg = supp_fig_18b.svg
    File supp_fig_18b_png = supp_fig_18b.png
    File supp_fig_18c_svg = supp_fig_18c.svg
    File supp_fig_18c_png = supp_fig_18c.png
    File supp_fig_18d_out = supp_fig_18d.plot

    File supp_fig_19c_svg = supp_fig_19c.svg
    File supp_fig_19c_png = supp_fig_19c.png
    File supp_fig_19d_svg = supp_fig_19d.svg
    File supp_fig_19d_png = supp_fig_19d.png

    File supp_fig_23b_out = supp_fig_23b.plot
    File supp_fig_25b_svg = supp_fig_25b.svg
    File supp_fig_25b_png = supp_fig_25b.png
    File supp_fig_25d_out = supp_fig_25d.plot

    File supp_table_4 = str_tables_for_paper.singly_finemapped_strs_for_paper
    File singly_finemapped_strs_sorted = str_tables_for_paper.singly_finemapped_strs_sorted
    File supp_table_5 = str_tables_for_paper.confidently_finemapped_strs_for_paper
    File confidently_finemapped_strs_sorted = str_tables_for_paper.confidently_finemapped_strs_sorted
    File supp_table_6 = concordance_in_other_ethnicities.stats

    Array[File] str_gwas_results_for_publication = reformat_my_str_gwas_table_for_publication.unfiltered
    Array[File] filtered_str_gwas_results_for_publication = reformat_my_str_gwas_table_for_publication.filtered
    #Array[Array[File]] ethnic_str_gwas_results_for_publication = reformat_my_ethnic_str_gwas_table_for_publication.unfiltered 
    #Array[Array[File]] filtered_ethnic_str_gwas_results_for_publication = reformat_my_ethnic_str_gwas_table_for_publication.filtered 
  }
}
