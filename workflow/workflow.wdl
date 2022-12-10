version 1.1

import tasks.wdl

workflow main {

  call tasks.ethnic_sample_lists # TODO inputs
  
  scatter (ethnic_sample_list in ethnic_sample_lists.sample_lists) {
    call tasks.qced_sample_list as qced_ethnic_sample_lists { input:
      unqced_sample_list = ethnic_sample_list.right,
      # TODO ssample_lists_to_filter
    }
  }
}
