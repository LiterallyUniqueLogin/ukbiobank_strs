version 1.0

import "gwas_tasks.wdl"
import "expanse_tasks.wdl"

# figures that use preexisting files
# and so are not part of the main GWAS -> finemapping pipeline
workflow expanse_figures {

  input {
    String script_dir = "."
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample"
    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam"

    File specific_alleles = "association/specific_alleles.tab"

    File CBL_gtex_expression = "misc_data/gtex_yang/CBL_chr11_119206290_GTEX_TPM.tsv"
    File Liver_SLC2A2_exon4_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon4_psi.tsv"
    File Liver_SLC2A2_exon6_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon6_psi.tsv"

    # files that could be computed via WDL but were pregenerated beforehand and include randomness, so need to be reused for consistency
    File unrelated_white_brits_sample_list = "sample_qc/runs/white_brits/no_phenotype/combined_unrelated.sample"
    File unrelated_black_sample_list = "sample_qc/runs/black/no_phenotype/combined_unrelated.sample"
    File unrelated_south_asian_sample_list = "sample_qc/runs/south_asian/no_phenotype/combined_unrelated.sample"
    File unrelated_chinese_sample_list = "sample_qc/runs/chinese/no_phenotype/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_not_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File platelet_count_sample_list = "sample_qc/runs/white_brits/platelet_count/combined_unrelated.sample"
  }

  Array[String] phenotype_names = [
		"alanine_aminotransferase",
		"albumin",
		"alkaline_phosphatase",
		"apolipoprotein_a",
		"apolipoprotein_b",
		"aspartate_aminotransferase",
		"c_reactive_protein",
		"calcium",
		"cholesterol",
		"creatinine",
		"cystatin_c",
		"eosinophil_count",
		"eosinophil_percent",
		"gamma_glutamyltransferase",
		"glucose",
		"glycated_haemoglobin",
		"haematocrit",
		"haemoglobin_concentration",
		"hdl_cholesterol",
		"igf_1",
		"ldl_cholesterol_direct",
		"lymphocyte_count",
		"lymphocyte_percent",
		"mean_corpuscular_haemoglobin",
		"mean_corpuscular_haemoglobin_concentration",
		"mean_corpuscular_volume",
		"mean_platelet_volume",
		"mean_sphered_cell_volume",
		"neutrophil_count",
		"neutrophil_percent",
		"phosphate",
		"platelet_count",
		"platelet_crit",
		"platelet_distribution_width",
		"red_blood_cell_count",
		"red_blood_cell_distribution_width",
		"shbg",
		"total_bilirubin",
		"total_protein",
		"triglycerides",
		"urate",
		"urea",
		"vitamin_d",
		"white_blood_cell_count",
  ]

  # these files were pregenerated but could be generated from the WDL
	scatter (phenotype_name in phenotype_names) {
    String peak_files = "signals/peaks/~{phenotype_name}_250000_5e-8.tab"
  }

  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
    PFiles imputed_snp_p_files = {
      "pgen": "array_imputed/pfile_converted/chr{chrom+1}.pgen",
      "pvar": "array_imputed/pfile_converted/chr{chrom+1}.pvar",
      "psam": "array_imputed/pfile_converted/chr{chrom+1}.psam",
    }
  }

  call expanse_tasks.extract_field as sc_white_brits { input :
    script_dir = script_dir,
    id = 22006
  }

  call gwas_tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = script_dir,
    sc = sc_white_brits.data
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
    fam_file = fam_file,
    sc_pcs = pcs.data,
    sc_assessment_ages = assessment_ages.data
  }

  call gwas_tasks.fig_4a { input :
    script_dir = script_dir,
    all_samples_list = all_samples_list,
    white_brits_sample_list = white_brits_sample_list.data,
    black_sample_list = unrelated_black_sample_list,
    south_asian_sample_list = unrelated_south_asian_sample_list,
    chinese_sample_list = unrelated_chinese_sample_list,
    str_vcf_chr_11 = str_vcfs[10],
    specific_alleles = specific_alleles,
  }

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
    samples_for_phenotype = platelet_count_sample_list,
    is_binary = false
  }

  call gwas_tasks.str_spot_test as CBL_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data,
    transformed_phenotype = transformed_platelet_count.data,
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.locus_plot as fig_4b { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [CBL_assoc.data],
    dosage_fraction_threshold = 0.001
  }

  # use precomputed sample lists
  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_not_SNP_samples { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = unrelated_samples_CBL_hom_not_begin_C_T_snp,
    is_binary = false
  }

  call gwas_tasks.str_spot_test as CBL_hom_not_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data,
    transformed_phenotype = platelet_count_transformed_CBL_hom_not_SNP_samples.data, 
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_SNP_samples { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = unrelated_samples_CBL_hom_begin_C_T_snp,
    is_binary = false
  }
   
  call gwas_tasks.str_spot_test as CBL_hom_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data, 
    transformed_phenotype = platelet_count_transformed_CBL_hom_SNP_samples.data,
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.locus_plot as fig_4f { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [CBL_hom_not_SNP_assoc.data, CBL_hom_SNP_assoc.data],
    group_names = ["homozygous (CGG)n", "homozygous CGGTGG(CGG)m"],
    dosage_threshold = 200,
  }

  call gwas_tasks.summarize_peaks as fig_1ef { input :
    script_dir = script_dir,
    phenotype_names = phenotype_names,
    peak_files = peak_files
  }

  ### plot with data from GTEx Yang

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_CBL_gtex_expression { input :
    script_dir = script_dir,
    individual_tsv = CBL_gtex_expression,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "TPM(expression)"
  }

  call gwas_tasks.locus_plot as fig_4g { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "CBL expression",
    unit = "TPM",
    data_tsvs = [summarized_CBL_gtex_expression.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon4 { input :
    script_dir = script_dir,
    individual_tsv = Liver_SLC2A2_exon4_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_15a { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2 exon 4",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon4.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon6 { input :
    script_dir = script_dir,
    individual_tsv = Liver_SLC2A2_exon6_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_15b { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2 exon 6",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon6.out],
    dosage_threshold = 5
  }

  output {
    File fig_1e_svg_out = fig_1ef.barplot_svg
    File fig_1e_png_out = fig_1ef.barplot_png
    File fig_1f_svg_out = fig_1ef.heatmap_svg
    File fig_1f_png_out = fig_1ef.heatmap_png

    File fig_4a_svg_out = fig_4a.svg
    File fig_4a_png_out = fig_4a.png
    File fig_4b_svg_out = fig_4b.svg
    File fig_4b_png_out = fig_4b.png
    File fig_4f_svg_out = fig_4f.svg
    File fig_4f_png_out = fig_4f.png
    File fig_4g_svg_out = fig_4g.svg
    File fig_4g_png_out = fig_4g.png

    File supp_fig_15a_svg_out = supp_fig_15a.svg
    File supp_fig_15a_png_out = supp_fig_15a.png
    File supp_fig_15b_svg_out = supp_fig_15b.svg
    File supp_fig_15b_png_out = supp_fig_15b.png
  }  
}
