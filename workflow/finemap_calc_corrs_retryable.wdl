version 1.0

import "finemapping_tasks.wdl"

workflow finemap_calc_corrs_retryable {

  input {
    String script_dir
    File gts_h5
  }

  call finemapping_tasks.finemap_calc_corrs as try_zero { input :
    script_dir = script_dir,
    gts_h5 = gts_h5,
    time = "1h"
  }

  if (!defined(finemap_calc_corrs.lds_h5)) {
    call finemapping_tasks.finemap_calc_corrs as try_one { input :
      script_dir = script_dir,
      gts_h5 = gts_h5,
      time = "47h30m"
    }
  }

  output {
    File lds_h5 = select_first([try_zero.lds_h5, try_one.lds_h5])
  }
}
