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
    escalate = read_boolean(stdout())
  }

  command <<<
    envsetup python -c "
      with open(~{converged}) as converged:
        did_converge = next(converged).strip()
        if did_converge == 'FALSE':
          if ~{assert}:
            assert False
          print("true")
          exit()
        else:
          assert did_converge == 'TRUE'
      if ~{length(CSs)} < ~{max_CSs}:
          print("false")
          exit()
      assert ~{length(CSs)} == ~{max_CSs}:
      CSs = '~{sep=" " CSs}'.split()
      found_count = 0
      done = False
      redo = False
      for num in range(~{maax_CSs}):
        assert CSs[num].split('/')[-1] == f'cs{num}.txt'
        with open(CSs[num]) as cs:
          n_vars = len(next(cs).split())
          next(cs)
          min_ld = float(next(cs).split()[0])
        if min_ld < .2 and n_vars > 50:
          print("false")
          exit()
      if ~{assert}:
        assert False
      print("true")
      exit()
    "
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
  }
}

workflow escalating_susie_run {

  input {
    String script_dir
    
    File alpha
    File pheno_residuals_h5
  }

  call retryable_susie_run.retryable_susie_run as try_zero { input :
    script_dir = script_dir,
    alpha = alpha,
    pheno_residual_h5 = pheno_residual_h5,
    L = 10,
    max_iter = 100
  }

  call should_escalate_susie_run as escalate_zero { input :
    converged = try_zero.converged,
    CSs = try_zero.CSs,
    max_CSs = 10
  }

  if (escalate_zero.escalate)) {
    call retryable_susie_run.retryable_susie_run as try_one { input :
      script_dir = script_dir,
      alpha = alpha,
      pheno_residual_h5 = pheno_residual_h5,
      L = 30,
      max_iter = 500
    }
  
    call should_escalate_susie_run as escalate_one { input :
      converged = try_one.converged,
      CSs = try_one.CSs,
      max_CSs = 30
    }

    if (escalate_one.escalate) {
      call retryable_susie_run.retryable_susie_run as try_two { input :
        script_dir = script_dir,
        alpha = alpha,
        pheno_residual_h5 = pheno_residual_h5,
        L = 50,
        max_iter = 500
      }
  
      call should_escalate_susie_run { input :
        converged = try_two.converged,
        CSs = try_two.CSs,
        max_CSs = 50,
        assert = "True"
      }
    }
  }

  output {
    File lbf = select_first([try_two.lbf, try_one.lbf, try_zero.lbf])
    File lbf_variable = select_first([try_two.lbf_variable, try_one.lbf_variable, try_zero.lbf_variable])
    File sigma2 = select_first([try_two.sigma2, try_one.sigma2, try_zero.sigma2])
    File V = select_first([try_two.V, try_one.V, try_zero.V])
    File converged = select_first([try_two.converged, try_one.converged, try_zero.converged])
    File lfsr = select_first([try_two.lfsr, try_one.lsfr, try_zero.lsfr])
    File requested_coverage = select_first([try_two.requested_coverage, try_one.requested_coverage, try_zero.requested_coverage])
    File alpha = select_first([try_two.alpha, try_one.alpha, try_zero.alpha])
    Array[File] CSs = select_first([try_two.CSs, try_one.CSs, try_zero.CSs])
  }
}
