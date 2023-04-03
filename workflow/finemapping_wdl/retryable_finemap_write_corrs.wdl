version 1.0

import "finemapping_tasks.wdl"

workflow retryable_finemap_write_corrs {

  input {
    String script_dir
    File lds_h5
  }

  call finemapping_tasks.finemap_write_corrs as try_zero { input :
    script_dir = script_dir,
    lds_h5 = lds_h5,
    time = "1h"
  }

  if (!defined(try_zero.all_variants_ld)) {
    call finemapping_tasks.finemap_write_corrs as try_one { input :
      script_dir = script_dir,
      lds_h5 = lds_h5,
      time = "47h30m"
    }
  }

  output {
    File all_variants_ld = select_first([try_zero.all_variants_ld, try_one.all_variants_ld])
  }
}
