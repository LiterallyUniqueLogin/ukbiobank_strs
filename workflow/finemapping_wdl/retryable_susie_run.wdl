version 1.0

import "finemapping_tasks.wdl"

workflow retryable_susie_run {

  input {
    String script_dir
   
    File gts_h5 
    File pheno_residuals_h5
    Int L
    Int max_iter

    Float? tol
    Float? snp_p_over_str_p
    File? varnames_file
    Float? res_var
    Float? prior_var
  }

  call finemapping_tasks.susie_run as try_zero { input :
    script_dir = script_dir,
    time = "1h",
    mem = "16g",
    gts_h5 = gts_h5,
    pheno_residuals_h5 = pheno_residuals_h5,
    L = L,
    max_iter = max_iter,
    tol = tol,
    snp_p_over_str_p = snp_p_over_str_p,
    res_var = res_var,
    prior_var = prior_var,
    varnames_file = varnames_file,
  }

  #if (!defined(try_zero.alpha)) {
  if (!defined(try_zero.susie_output)) {
    call finemapping_tasks.susie_run as try_one { input :
      script_dir = script_dir,
      time = "47h30m",
      mem = "32g",
      gts_h5 = gts_h5,
      pheno_residuals_h5 = pheno_residuals_h5,
      L = L,
      max_iter = max_iter,
      tol = tol,
      snp_p_over_str_p = snp_p_over_str_p,
      res_var = res_var,
      prior_var = prior_var,
      varnames_file = varnames_file,
    }

    #if (!defined(try_one.alpha)) {
    if (!defined(try_one.susie_output)) {
      call finemapping_tasks.susie_run as try_two { input :
        script_dir = script_dir,
        time = "47h30m",
        mem = "64g",
        gts_h5 = gts_h5,
        pheno_residuals_h5 = pheno_residuals_h5,
        L = L,
        max_iter = max_iter,
        tol = tol,
        snp_p_over_str_p = snp_p_over_str_p,
        res_var = res_var,
        prior_var = prior_var,
        varnames_file = varnames_file,
      }
    
      #if (!defined(try_two.alpha)) {
      if (!defined(try_two.susie_output)) {
        call finemapping_tasks.susie_run as try_three { input :
          script_dir = script_dir,
          time = "47h30m",
          mem = "128g",
          gts_h5 = gts_h5,
          pheno_residuals_h5 = pheno_residuals_h5,
          L = L,
          max_iter = max_iter,
          tol = tol,
          snp_p_over_str_p = snp_p_over_str_p,
          res_var = res_var,
          prior_var = prior_var,
          varnames_file = varnames_file,
        }
      }
    }
  }

  output {
    SuSiE_output susie_output = select_first([try_three.susie_output, try_two.susie_output, try_one.susie_output, try_zero.susie_output])
  }
}
