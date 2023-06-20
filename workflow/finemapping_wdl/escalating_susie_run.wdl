version 1.0

import "retryable_susie_run.wdl"

# increase L and max_retries if they are not sufficient

task should_escalate_susie_run {
  input {
    File converged
    Array[File] CSs
    Int max_CSs
    String assert = "False"
  }

  output {
    Boolean escalate = read_boolean(stdout())
  }

  command <<<
    envsetup python -c "
    with open('~{converged}') as converged:
      did_converge = next(converged).strip()
      if did_converge == 'FALSE':
        if ~{assert}:
          assert False, 'Did not converge'
        print('true')
        exit()
      else:
        assert did_converge == 'TRUE'
    if ~{length(CSs)} < ~{max_CSs}:
        print('false')
        exit()
    assert ~{length(CSs)} == ~{max_CSs}, 'Not the maximum number of CSs'
    CSs = '~{sep=' ' CSs}'.split()
    found_count = 0
    done = False
    redo = False
    for idx, num in enumerate(sorted(str(num) for num in range(1, 1+~{max_CSs}))):
      assert CSs[idx].split('_')[-1] == f'cs{num}.txt', ('CSs in the wrong order', CS)
      with open(CSs[idx]) as cs:
        n_vars = len(next(cs).split())
        next(cs)
        min_ld = float(next(cs).split()[0])
      if min_ld < .2 and n_vars > 50:
        print('false')
        exit()
    if ~{assert}:
      assert False, 'Reached the end and failed'
    print('true')
    exit()
    "
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    memory: "2GB"
  }
}

workflow escalating_susie_run {

  input {
    String script_dir
    
    File gts_h5
    File pheno_residuals_h5

    File colnames

    Float? tol
    Float? snp_p_over_str_p
    File? varnames_file
    Float? res_var
    Float? prior_var
    String prefix
  }

  call retryable_susie_run.retryable_susie_run as try_zero { input :
    script_dir = script_dir,
    gts_h5 = gts_h5,
    pheno_residuals_h5 = pheno_residuals_h5,
    L = 10,
    max_iter = 100,
    tol = tol,
    snp_p_over_str_p = snp_p_over_str_p,
    res_var = res_var,
    prior_var = prior_var,
    varnames_file = varnames_file,
    colnames = colnames,
    prefix=prefix
  }

  call should_escalate_susie_run as escalate_zero { input :
    converged = try_zero.susie_output.subset.converged,
    CSs = try_zero.susie_output.CSs,
    max_CSs = 10
  }

  if (escalate_zero.escalate) {
    call retryable_susie_run.retryable_susie_run as try_one { input :
      script_dir = script_dir,
      gts_h5 = gts_h5,
      pheno_residuals_h5 = pheno_residuals_h5,
      L = 30,
      max_iter = 500,
      tol = tol,
      snp_p_over_str_p = snp_p_over_str_p,
      res_var = res_var,
      prior_var = prior_var,
      varnames_file = varnames_file,
      colnames = colnames,
      prefix=prefix
    }
  
    call should_escalate_susie_run as escalate_one { input :
      converged = try_one.susie_output.subset.converged,
      CSs = try_one.susie_output.CSs,
      max_CSs = 30
    }

    if (escalate_one.escalate) {
      call retryable_susie_run.retryable_susie_run as try_two { input :
        script_dir = script_dir,
        gts_h5 = gts_h5,
        pheno_residuals_h5 = pheno_residuals_h5,
        L = 50,
        max_iter = 500,
        tol = tol,
        snp_p_over_str_p = snp_p_over_str_p,
        res_var = res_var,
        prior_var = prior_var,
        varnames_file = varnames_file,
        colnames = colnames,
        prefix=prefix
     }
  
      call should_escalate_susie_run { input :
        converged = try_two.susie_output.subset.converged,
        CSs = try_two.susie_output.CSs,
        max_CSs = 50,
        assert = "True"
      }
    }
  }

  output {
    SuSiE_output susie_output = select_first([try_two.susie_output, try_one.susie_output, try_zero.susie_output])
  }
}
