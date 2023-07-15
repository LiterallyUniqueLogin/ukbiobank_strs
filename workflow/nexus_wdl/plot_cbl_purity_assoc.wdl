version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "../gwas_wdl/prep_samples_and_phenotype_workflow.wdl"

# same as wgs_tasks.regional_my_str_gwas except calls associaTR instead of the local code
# the output p-values and betas aren't relevant for this - this task doesn't bother to include
# transforming the phenotype or adding covariates
# the reliable bit should be the plotting association statistics
task associaTR_my_str_gwas {
  input {
    VCF str_vcf
    File untransformed_phenotype
    String chrom
    Int pos
    String phenotype_name
  }

  output {
    File data = "out.tab"
  }

  command <<<
    envsetup associaTR \
      out.tab \
      ~{str_vcf.vcf} \
      ~{phenotype_name} \
      ~{untransformed_phenotype} \
      --region ~{chrom}:~{pos}-~{pos} \
      --stats-for-plotting ~{untransformed_phenotype} \
      --plotting-ci-alphas 0.05 \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.5"
    dx_timeout: "2h"
    memory: "8GB"
  }
}

workflow plot_cbl_purity_assoc {

  input {
    String phenotype_name = "platelet_count"
    Int chrom = 11
    Int pos = 119206290

    String script_dir  = "dx://UKB_Test:/imputed_strs_paper"

    VCF wgs_vcf = {
      "vcf": "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz",
      "index": "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz.tbi"
    }

    File fam_file = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb22418_c1_b0_v2.fam"
    File withdrawn_sample_list = "dx://UKB_Test:/imputed_strs_paper/sample_qc/withdrawn.sample"
    File kinship = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb_rel.dat"

    String sc_data_dir = "dx://UKB_Test:/imputed_strs_paper/main_dataset/extracted_data"
    File sc_white_brits = "~{sc_data_dir}/22006.tsv"
    File sc_ethnicity_self_report = "~{sc_data_dir}/21000.tsv"
    File sc_sex_aneuploidy = "~{sc_data_dir}/22019.tsv"
    File sc_genetic_sex = "~{sc_data_dir}/22001.tsv"
    File sc_reported_sex = "~{sc_data_dir}/31.tsv"
    File sc_kinship_count = "~{sc_data_dir}/22021.tsv"
    File sc_assessment_ages = "~{sc_data_dir}/21003.tsv"
    File sc_pcs ="~{sc_data_dir}/22009.tsv"
    File sc_year_of_birth = "~{sc_data_dir}/34.tsv"
    File sc_month_of_birth = "~{sc_data_dir}/52.tsv"
    File sc_date_of_death = "~{sc_data_dir}/40000.tsv"

    Array[File] cbl_purity_sample_lists = ["dx://UKB_Test:/imputed_strs_paper/sample_qc/cbl_pure_white_brits.samples", "dx://UKB_Test:/imputed_strs_paper/sample_qc/cbl_impure_white_brits.samples"]
  }

  call gwas_tasks.phenotype_names

  Int phenotype_idx = phenotype_names.idxs[phenotype_name]
  String unit = phenotype_names.unit[phenotype_idx]

  scatter (sample_list in cbl_purity_sample_lists) {
    call prep_samples_and_phenotype_workflow.prep_samples_and_phenotype { input :
      script_dir = script_dir,
      PRIMUS_command = "run_PRIMUS.pl",

      phenotype_name = phenotype_name,
      categorical_covariate_names = [], # don't actually need these because we're ignoring the regression, just looking at the raw association
      categorical_covariate_scs = [],
      is_binary = false,
      is_zero_one_neg_nan = false,
      date_of_most_recent_first_occurrence_update = "",

      fam_file = fam_file,
      withdrawn_sample_list = withdrawn_sample_list,
      kinship = kinship,

      sc_white_brits = sc_white_brits,
      sc_ethnicity_self_report = sc_ethnicity_self_report,
      sc_sex_aneuploidy = sc_sex_aneuploidy,
      sc_genetic_sex = sc_genetic_sex,
      sc_reported_sex = sc_reported_sex,
      sc_kinship_count = sc_kinship_count,
      sc_assessment_ages = sc_assessment_ages,
      sc_pcs = sc_pcs,
      sc_year_of_birth = sc_year_of_birth,
      sc_month_of_birth = sc_month_of_birth,
      sc_date_of_death = sc_date_of_death,
      sc_phenotype = "~{sc_data_dir}/~{phenotype_names.ID[phenotype_idx]}.tsv",

      subpop_sample_list = sample_list,
    }

    call associaTR_my_str_gwas { input :
      str_vcf = wgs_vcf,
      untransformed_phenotype = prep_samples_and_phenotype.pheno_data[0],
      chrom = "chr~{chrom}",
      pos = pos,
      phenotype_name = phenotype_name,
    }
  }

  call gwas_tasks.locus_plot { input :
    script_dir = script_dir,
    chrom = "chr~{chrom}",
    pos = pos,
    phenotype_name = phenotype_name
    assoc_results = associaTR_my_str_gwas.data,
    dosage_fraction_threshold = 0.001,
    unit = unit,
    total_column_name = "sample_count_per_summed_length",
  }

  output {
    File png = locus_plot.png
    File svg = locus_plot.svg
  }
}
