# reserach_access_platform_workflow.wdl

version 1.0

import "platform_agnostic_workflow.wdl"

task extract_field {

  input {
    String script_dir
    File script = #TODO
    # TODO
  }

  output {
    File data
  }

  # TODO
  command <<<

  >>>
}

workflow main {
  
  input {
    Int project_number
  }

  call extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }

  call platform_agnostic_workflow.main { input:
    script_dir = script_dir,
    #TODO data_dir

    # identical for any chromosome number, so c1 works
    fam_file = "Bulk/Genotype\ Results/Genotype\ calls/ukb~{project_number}_c1_b0_v2.fam"
    assessment_ages_file = assessment_ages.data
  }
}
