version 1.0

import "finemapping_tasks.wdl"

workflow retryable_susie_run {

  input {
    String script_dir
   
    File gts_h5 
    File pheno_residuals_h5
    Int L
    Int max_iter

    File colnames

    Float? tol
    Float? snp_p_over_str_p
    File? varnames_file
    Float? res_var
    Float? prior_var
    String prefix
  }

  call finemapping_tasks.susie_run as try_zero { input :
    script_dir = script_dir,
    time = "1h",
    mem = "16GB",
    gts_h5 = gts_h5,
    pheno_residuals_h5 = pheno_residuals_h5,
    L = L,
    max_iter = max_iter,
    tol = tol,
    snp_p_over_str_p = snp_p_over_str_p,
    res_var = res_var,
    prior_var = prior_var,
    varnames_file = varnames_file,
    colnames = colnames,
    prefix=prefix,
  }

  if (defined(try_zero.alpha)) {
    SuSiE_output output_zero = object {
      subset: object {
        lbf: select_first([try_zero.lbf]),
        lbf_variable: select_first([try_zero.lbf_variable]),
        sigma2: select_first([try_zero.sigma2]),
        V: select_first([try_zero.V]),
        converged: select_first([try_zero.converged]),
        lfsr: select_first([try_zero.lfsr]),
        requested_coverage: select_first([try_zero.requested_coverage]),
        alpha: select_first([try_zero.alpha]),
        colnames: select_first([try_zero.colnames]),
      },
      CSs: try_zero.CSs
    }
  }

  if (!defined(try_zero.alpha)) {
    call finemapping_tasks.susie_run as try_one { input :
      script_dir = script_dir,
      time = "47h30m",
      mem = "32GB",
      gts_h5 = gts_h5,
      pheno_residuals_h5 = pheno_residuals_h5,
      L = L,
      max_iter = max_iter,
      tol = tol,
      snp_p_over_str_p = snp_p_over_str_p,
      res_var = res_var,
      prior_var = prior_var,
      varnames_file = varnames_file,
      colnames = colnames,
      prefix=prefix,
    }

    if (defined(try_one.alpha)) {
      SuSiE_output output_one = object {
        subset: object {
          lbf: select_first([try_one.lbf]),
          lbf_variable: select_first([try_one.lbf_variable]),
          sigma2: select_first([try_one.sigma2]),
          V: select_first([try_one.V]),
          converged: select_first([try_one.converged]),
          lfsr: select_first([try_one.lfsr]),
          requested_coverage: select_first([try_one.requested_coverage]),
          alpha: select_first([try_one.alpha]),
          colnames: select_first([try_one.colnames]),
        },
        CSs: try_one.CSs
      }
    }

    if (!defined(try_one.alpha)) {
      call finemapping_tasks.susie_run as try_two { input :
        script_dir = script_dir,
        time = "47h30m",
        mem = "64GB",
        gts_h5 = gts_h5,
        pheno_residuals_h5 = pheno_residuals_h5,
        L = L,
        max_iter = max_iter,
        tol = tol,
        snp_p_over_str_p = snp_p_over_str_p,
        res_var = res_var,
        prior_var = prior_var,
        varnames_file = varnames_file,
        colnames = colnames,
        prefix=prefix,
      }

      if (defined(try_two.alpha)) {
        SuSiE_output output_two = object {
          subset: object {
            lbf: select_first([try_two.lbf]),
            lbf_variable: select_first([try_two.lbf_variable]),
            sigma2: select_first([try_two.sigma2]),
            V: select_first([try_two.V]),
            converged: select_first([try_two.converged]),
            lfsr: select_first([try_two.lfsr]),
            requested_coverage: select_first([try_two.requested_coverage]),
            alpha: select_first([try_two.alpha]),
            colnames: select_first([try_two.colnames]),
          },
          CSs: try_two.CSs
        }
      }
   
      if (!defined(try_two.alpha)) {
        call finemapping_tasks.susie_run as try_three { input :
          script_dir = script_dir,
          time = "47h30m",
          mem = "128GB",
          gts_h5 = gts_h5,
          pheno_residuals_h5 = pheno_residuals_h5,
          L = L,
          max_iter = max_iter,
          tol = tol,
          snp_p_over_str_p = snp_p_over_str_p,
          res_var = res_var,
          prior_var = prior_var,
          varnames_file = varnames_file,
          colnames = colnames,
          prefix=prefix,
        }

        SuSiE_output output_three = object {
          subset: object {
            lbf: select_first([try_three.lbf]),
            lbf_variable: select_first([try_three.lbf_variable]),
            sigma2: select_first([try_three.sigma2]),
            V: select_first([try_three.V]),
            converged: select_first([try_three.converged]),
            lfsr: select_first([try_three.lfsr]),
            requested_coverage: select_first([try_three.requested_coverage]),
            alpha: select_first([try_three.alpha]),
            colnames: select_first([try_three.colnames]),
          },
          CSs: try_three.CSs
        }
      }
    }
  }

  output {
    SuSiE_output susie_output = select_first([output_zero, output_one, output_two, output_three])
  }
}
