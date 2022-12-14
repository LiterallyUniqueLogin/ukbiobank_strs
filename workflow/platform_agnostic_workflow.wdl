# platform agnostic workflow

version 1.0

import "tasks.wdl"

workflow main {

  input {
    String script_dir

    # task for generating these?
    File white_brits_sample_list
    File fam_file 

    # data showcase files
    File sc_assessment_ages
    File sc_ethnicity_self_report
  }

  call tasks.ethnic_sample_lists { input: 
    script_dir = script_dir,
    white_brits_sample_list = white_brits_sample_list,
    sc_ethnicity_self_report = sc_ethnicity_self_report
  }
  
#  scatter (ethnicity in ethnic_sample_lists.ethnicities) {
#    call tasks.qced_sample_list as qced_ethnic_sample_lists { input:
#      script_dir = script_dir,
#      unqced_sample_list = ethnic_sample_lists.sample_lists[ethnicity]
#      # TODO ssample_lists_to_filter =
#    }
#  }
#
#  call tasks.load_shared_covars { input:
#    script_dir = script_dir,
#    fam_file = fam_file,
#    # TODO replace sqc file with PCs file data field 22009
#    sc_assessment_ages = sc_assessment_ages
#  }
#
  output {
    # TODO
    String foo = "bar"
  }
}
