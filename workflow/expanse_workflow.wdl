version 1.0

import "platform_agnostic_workflow.wdl"

task extract_field {
  input {
    String script_dir
    File script = "~{script_dir}/main_dataset/decompress_trait.py"
    File ukbconv = "~{script_dir}/ukb_utilities/ukbconv"
    File encoding = "~{script_dir}/ukb_utilities/encoding.ukb"
    Array[File]+ fields_files = ["main_dataset/raw_data/fields46781.ukb", "main_dataset/raw_data/fields46782.ukb"]
    Array[File]+ enc_files = ["main_dataset/raw_data/ukb46781.enc_ukb", "main_dataset/raw_data/ukb46782.enc_ukb"]

    Int id # data field id
  }

  output {
    File data = "~{id}.txt"
  }

  command <<<
    ~{script} \
      ~{id} \
      ~{id} \
      ~{ukbconv} \
      ~{encoding} \
      --fields-files ~{sep=" " fields_files} \
      --enc-files ~{sep=" " enc_files}
  >>>

  runtime {
    dx_timeout: "5h"
  }
}

workflow main {

  input {
    String script_dir  = "."

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

    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam"
    File withdrawn_sample_list = "sample_qc/common_filters/remove/withdrawn.sample"
    File kinship = "misc_data/ukbgene/ukb46122_rel_s488282.dat"
  }

  call extract_field as white_brits { input:
    script_dir = script_dir,
    id = 22006
  }

  call extract_field as ethnicity_self_report { input :
    script_dir = script_dir,
    id = 21000
  }

  call extract_field as sex_aneuploidy { input:
    script_dir = script_dir,
    id = 22019
  }

  call extract_field as genetic_sex { input:
    script_dir = script_dir,
    id = 22001
  }

  call extract_field as reported_sex { input:
    script_dir = script_dir,
    id = 31
  }

  call extract_field as kinship_count { input:
    script_dir = script_dir,
    id = 22021
  }

  call extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }
  
  call extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
  }

  call extract_field as year_of_birth { input :
    script_dir = script_dir,
    id = 34
  }

  call extract_field as month_of_birth { input :
    script_dir = script_dir,
    id = 52
  }

  call extract_field as date_of_death { input :
    script_dir = script_dir,
    id = 40000
  }

  call extract_field as phenotype { input :
    script_dir = script_dir,
    id = phenotype_id
  }

  scatter (categorical_covariate_id in categorical_covariate_ids) {
    call extract_field as categorical_covariates { input :
      script_dir = script_dir,
      id = categorical_covariate_id
    }
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl"

    # TODO chr_lens

    # TODO str_vcfs
    # TODO specific_alleles

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariates.data,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,

    fam_file = fam_file, # Could instead create a task for downloading this with ukbgene
    withdrawn_sample_list = withdrawn_sample_list,
    kinship = kinship, # could create a task for downloading this with ukbgene

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
