version 1.0

import "tasks.wdl"

workflow main {

  input {
    String script_dir
    String data_dir # locations of pieces of data provided before the run

    File fam_file 

    File assessment_ages_file
  }

  call tasks.ethnic_sample_lists { input: 
    script_dir = script_dir,
    white_brits_sample_list = "~{data_dir}/sample_qc/common_filters/ethnicity/white_brits.sample"
    # TODO data_showcase_ethnicity_self_report = 
  }
  
  scatter (ethnic_sample_list in ethnic_sample_lists.sample_lists) {
    call tasks.qced_sample_list as qced_ethnic_sample_lists { input:
      script_dir = script_dir,
      unqced_sample_list = ethnic_sample_list.right,
      # TODO ssample_lists_to_filter =
    }
  }

  call load_shared_covars { input:
    script_dir = script_dir,
    fam_file = fam_file,
    # TODO replace sqc file with PCs file data field 22009
    assessent_ages_file = assessment_ages_file
  }
}
