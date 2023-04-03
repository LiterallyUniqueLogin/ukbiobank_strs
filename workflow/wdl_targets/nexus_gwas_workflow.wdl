# reserach_access_platform_workflow.wdl

version 1.0

import "platform_agnostic_workflow.wdl"

workflow rap_main {

  input {
    String script_dir  = "dx://UKB_Test:/imputed_strs_paper"

    String phenotype_name = "platelet_count"
    Int phenotype_id = 30080
    Array[String] categorical_covariate_names = ["platelet_count_device_id"]
    Array[Int] categorical_covariate_ids = [30083]
    #String phenotype_name
    #Int phenotype_id
    #Array[String] categorical_covariate_names = []
    #Array[Int] categorical_covariate_ids = []
    Boolean is_binary = false
    Boolean is_zero_one_neg_nan = false # different binary encoding
    String date_of_most_recent_first_occurrence_update = "2021-04-01"

    File fam_file = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb22418_c1_b0_v2.fam"
    File all_samples_list = "dx://UKB_Test:/Bulk/Imputation/UKB%20imputation%20from%20genotype/ukb22828_c1_b0_v3.sample"
    File withdrawn_sample_list = "dx://UKB_Test:/imputed_strs_paper/sample_qc/withdrawn.sample"
    File kinship = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb_rel.dat"
  }
  
  String sc_data_dir = "dx://UKB_Test:/imputed_strs_paper/main_dataset/extracted_data"

	scatter (categorical_covariate_id in categorical_covariate_ids) {
		String categorical_covariates = "~{sc_data_dir}/~{categorical_covariate_id}.tsv"
	}

  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "dx://UKB_Test:/imputed_strs_paper/imputed_strs/chr~{chrom+1}.vcf.gz",
      "index": "dx://UKB_Test:/imputed_strs_paper/imputed_strs/chr~{chrom+1}.vcf.gz.tbi",
    }
    #PFiles imputed_snp_p_files = {
    #  "pgen": "array_imputed/pfile_converted/chr{chrom+1}.pgen",
    #  "pvar": "array_imputed/pfile_converted/chr{chrom+1}.pvar",
    #  "psam": "array_imputed/pfile_converted/chr{chrom+1}.psam",
    #}
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",

    chr_lens = "dx://UKB_Test:/imputed_strs_paper/misc_data/chr_lens.txt",

    str_vcfs = str_vcfs,
    #imputed_snp_p_files = imputed_snp_p_files,
    specific_alleles = "dx://UKB_Test:/imputed_strs_paper/association/specific_alleles.tab",

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariates,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,

    fam_file = fam_file, # Could instead create a task for downloading this with ukbgene
    all_samples_list = all_samples_list,
    withdrawn_sample_list = withdrawn_sample_list,
    kinship = kinship, # could create a task for downloading this with ukbgene

    sc_white_brits = "~{sc_data_dir}/22006.tsv",
    sc_ethnicity_self_report = "~{sc_data_dir}/21000.tsv",
    sc_sex_aneuploidy = "~{sc_data_dir}/22019.tsv",
    sc_genetic_sex = "~{sc_data_dir}/22001.tsv",
    sc_reported_sex = "~{sc_data_dir}/31.tsv",
    sc_kinship_count = "~{sc_data_dir}/22021.tsv",
    sc_assessment_ages = "~{sc_data_dir}/21003.tsv",
    sc_pcs ="~{sc_data_dir}/22009.tsv",
    sc_year_of_birth = "~{sc_data_dir}/34.tsv",
    sc_month_of_birth = "~{sc_data_dir}/52.tsv",
    sc_date_of_death = "~{sc_data_dir}/40000.tsv",
    sc_phenotype = "~{sc_data_dir}/~{phenotype_id}.tsv",
  }

  output {
    Array[File] out_sample_lists = main.out_sample_lists
  }
}
