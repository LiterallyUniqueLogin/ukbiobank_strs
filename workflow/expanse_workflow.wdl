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
    String script_dir 

    Int project_number
    Int fam_file_sample_num
  }

  call extract_field as ethnicity_self_report { input :
    script_dir = script_dir,
    id = 21000
  }

  call extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,

    # TODO where did this come from?
    white_brits_sample_list = "sample_qc/common_filters/ethnicity/white_brits.sample",
    # Could instead create a task for downloading this with ukbgene
    fam_file = "microarray/ukb~{project_number}_cal_chr1_v2_s~{fam_file_sample_num}.fam",

    sc_ethnicity_self_report = ethnicity_self_report.data,
    sc_assessment_ages = assessment_ages.data
  }
}
