version 1.0

import "finemapping_tasks.wdl"

workflow retryable_susie_load_gts {

  input {
    String script_dir
    VCF strs
    bgen snps
    File all_samples
    File phenotype_samples
    File pheno_data
    File shared_covars

    File colnames
    String phenotype_name
    region bounds
    Boolean hardcalls = false
  }

  call finemapping_tasks.susie_load_gts as try_zero { input :
    script_dir = script_dir,
    time = "1h",
    strs = strs,
    snps = snps,
    all_samples = all_samples,
    phenotype_samples = phenotype_samples,
    pheno_data = pheno_data,
    shared_covars = shared_covars,

    colnames = colnames,
    phenotype_name = phenotype_name,
    bounds = bounds,
    hardcalls = hardcalls
  }

  if (!defined(try_zero.gts_h5)) {
    call finemapping_tasks.susie_load_gts as try_one { input :
      script_dir = script_dir,
      time = "47h30m",
      strs = strs,
      snps = snps,
      all_samples = all_samples,
      phenotype_samples = phenotype_samples,
      pheno_data = pheno_data,
      shared_covars = shared_covars,
  
      colnames = colnames,
      phenotype_name = phenotype_name,
      bounds = bounds,
      hardcalls = hardcalls
    }
  }

  output {
    File gts_h5 = select_first([try_zero.gts_h5, try_one.gts_h5])
    File pheno_residuals_h5 = select_first([try_zero.pheno_residuals_h5, try_one.pheno_residuals_h5])
    File readme = select_first([try_zero.readme, try_one.readme])
  }
}
