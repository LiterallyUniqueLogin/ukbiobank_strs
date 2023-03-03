version 1.0

import "finemapping_tasks.wdl"

workflow retryable_susie_run {

  input {
    String script_dir
    
    File alpha
    File pheno_residuals_h5
    Int L
    Int max_iter
  }

  call finemapping_tasks.susie_run as try_zero { input :
    script_dir = script_dir,
    time = "1h"
    mem = "16g"
    gts_h5 = gts_h5,
    pheno_residual_h5 = pheno_residual_h5,
    L = L,
    max_iter = max_iter
  }

  if (!defined(try_zero.gts_h5)) {
    call finemapping_tasks.susie_run as try_one { input :
      script_dir = script_dir,
      time = "47h30m"
      mem = "32g"
      gts_h5 = gts_h5,
      pheno_residual_h5 = pheno_residual_h5,
      L = L,
      max_iter = max_iter
    }

    if (!defined(try_one.gts_h5)) {
      call finemapping_tasks.susie_run as try_two { input :
        script_dir = script_dir,
        time = "47h30m"
        mem = "64g"
        gts_h5 = gts_h5,
        pheno_residual_h5 = pheno_residual_h5,
        L = L,
        max_iter = max_iter
      }
    
      if (!defined(try_two.gts_h5)) {
        call finemapping_tasks.susie_run as try_three { input :
          script_dir = script_dir,
          time = "47h30m"
          mem = "128g"
          gts_h5 = gts_h5,
          pheno_residual_h5 = pheno_residual_h5,
          L = L,
          max_iter = max_iter
        }
      }
    }
  }

  output {
    File lbf = select_first([try_three.lbf, try_two.lbf, try_one.lbf, try_zero.lbf])
    File lbf_variable = select_first([try_three.lbf_variable, try_two.lbf_variable, try_one.lbf_variable, try_zero.lbf_variable])
    File sigma2 = select_first([try_three.sigma2, try_two.sigma2, try_one.sigma2, try_zero.sigma2])
    File V = select_first([try_three.V, try_two.V, try_one.V, try_zero.V])
    File converged = select_first([try_three.converged, try_two.converged, try_one.converged, try_zero.converged])
    File lfsr = select_first([try_three.lsfr, try_two.lfsr, try_one.lsfr, try_zero.lsfr])
    File requested_coverage = select_first([try_three.requested_coverage, try_two.requested_coverage, try_one.requested_coverage, try_zero.requested_coverage])
    File alpha = select_first([try_three.alpha, try_two.alpha, try_one.alpha, try_zero.alpha])
    Array[File] CSs = select_first([try_three.CSs, try_two.CSs, try_one.CSs, try_zero.CSs])
  }
}
