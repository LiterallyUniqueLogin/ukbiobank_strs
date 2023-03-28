version 1.0

import "platform_agnostic_workflow.wdl"
import "expanse_tasks.wdl"

workflow main {

  input {
    String script_dir  = "."

    File chr_lens = "misc_data/genome/chr_lens.txt"
    File str_loci = "snpstr/str_loci.txt"

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
    String date_of_most_recent_first_occurrence_update = "2021-04-01" # only needed for binary phenotypes

    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam" # Could instead create a task for downloading this with ukbgene
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample" # could instead create a task for downloading this with ukbgene
    File withdrawn_sample_list = "sample_qc/common_filters/remove/withdrawn.sample"
    File kinship = "misc_data/ukbgene/ukb46122_rel_s488282.dat" # could create a task for downloading this with ukbgene
  }

  call expanse_tasks.extract_field as white_brits { input:
    script_dir = script_dir,
    id = 22006
  }

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

  call expanse_tasks.extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }
  
  call expanse_tasks.extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
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

  call expanse_tasks.extract_field as phenotype { input :
    script_dir = script_dir,
    id = phenotype_id
  }

  scatter (categorical_covariate_id in categorical_covariate_ids) {
    call expanse_tasks.extract_field as categorical_covariates { input :
      script_dir = script_dir,
      id = categorical_covariate_id
    }
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
    bgen imputed_snp_bgens = {
      "bgen": "array_imputed/ukb_imp_chr{chrom+1}_v3.bgen",
      "index": "array_imputed/ukb_imp_chr{chrom+1}_v3.bgen.bgi"
    }
    File snp_vars_to_filter_from_finemapping = "finemapping/str_imp_snp_overlaps/chr~{chrom+1}_to_filter.tab'"
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",
    finemap_command = "utilities/finemap/finemap_v1.4_x86_64",

    chr_lens = chr_lens,
    str_loci = str_loci,

    str_vcfs = str_vcfs,
    #imputed_snp_p_files = imputed_snp_p_files,
    imputed_snp_bgens = imputed_snp_bgens,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping, # could write a task for generating this file

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariates.data,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,

    fam_file = fam_file,
    all_samples_list = all_samples_list,
    withdrawn_sample_list = withdrawn_sample_list,
    kinship = kinship, 

    sc_white_brits = white_brits.data,
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
    sc_phenotype = phenotype.data
  }

  output {
    Array[File] out_sample_lists = main.out_sample_lists
  }
}
