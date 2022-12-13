version 1.0

import "platform_agnostic_workflow.wdl"

#task extract_field {
#  input {
#    String script_dir
#    File script = "~{script_dir}/main_dataset/decompress_trait.py"
#    Array[File]+ enc_files
#    # TODO
#  }
#
#  output {
#    File data
#  }
#
#  command <<<
#    ~{script}  #TODO
#  >>>
#
#  runtime {
#    dx_timeout: "5h"
#  }
#}

workflow main {

  input {
    String script_dir 

    Int project_number
    Int fam_file_sample_num
  }

  call extract_field as assessment_ages { input :
    script_dir = script_dir,
    # TODO
    id = 21003
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,
    # TODO data_dir

    # Could instead create a task for downloading this with ukbgene
    fam_file = "microarray/ukb~{project_number}_cal_chr1_v2_s~{fam_file_sample_num}.fam",
    assessment_ages_file = assessment_ages.data
  }
}
