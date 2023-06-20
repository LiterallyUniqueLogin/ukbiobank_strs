version 1.0

import "finemapping_tasks.wdl"

task should_escalate_finemap_run {
  input {
    Array[File] creds
    Int max_creds
    File log
    String assert = "False"
  }

  output {
    Boolean escalate = read_boolean(stdout())
  }

  command <<<
    envsetup python -c "
    creds = '~{sep=' ' creds}'.split()
    for cred in creds:
      if cred[-2:] == '~{max_creds}':
        if ~{assert}:
          assert False, 'Found a region with maximum allowed number of credible sets'
        else:
          print('true')
          exit()
    print('false')

    import re
    assert re.search('converged after [0-9]* iterations', open('~{log}').read()), 'Did not converge'
    "
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    memory: "2GB"
  }
}
 
workflow escalating_finemap {
  input {
    String script_dir
    String finemap_command
    File master
    File zfile
    File all_variants_ld
    Boolean prior_snps
    Float? prior_std
    Float? prob_conv_sss_tol
    String prefix
    Int cache_breaker = 0
  }

  call finemapping_tasks.finemap_run as try_zero { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    master = master,
    zfile = zfile,
    all_variants_ld = all_variants_ld,
    causal_snps = 20,
    prior_snps = prior_snps,
    prior_std = prior_std,
    prob_conv_sss_tol = prob_conv_sss_tol,
    prefix=prefix,
    cache_breaker = cache_breaker
  }

  call should_escalate_finemap_run as escalate { input :
    creds = try_zero.finemap_output.creds,
    log = try_zero.finemap_output.subset.log_sss,
    max_creds = 20
  }

  if (escalate.escalate) {
    call finemapping_tasks.finemap_run as try_one { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      master = master,
      zfile = zfile,
      all_variants_ld = all_variants_ld,
      causal_snps = 40,
      prior_snps = prior_snps,
      prior_std = prior_std,
      prob_conv_sss_tol = prob_conv_sss_tol,
      prefix=prefix,
      cache_breaker = cache_breaker
    }

    call should_escalate_finemap_run { input :
      creds = try_one.finemap_output.creds,
      log = try_one.finemap_output.subset.log_sss,
      max_creds = 40,
      assert = "True"
    }
  }

  output {
    FINEMAP_output finemap_output = select_first([try_one.finemap_output, try_zero.finemap_output])
    File finemap_input_z = select_first([try_one.finemap_input_z, try_zero.finemap_input_z])
  }
}  
