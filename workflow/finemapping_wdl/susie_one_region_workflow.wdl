version 1.0

import "finemapping_tasks.wdl"
import "retryable_susie_load_gts.wdl"
import "escalating_susie_run.wdl"

workflow susie_one_region {

  input {
    String script_dir

    # one per chrom
    VCF str_vcf
    bgen imputed_snp_bgen
    File snp_vars_to_filter_from_finemapping

    File? shared_covars
    File phenotype_samples
    File transformed_phenotype_data

    File? my_str_gwas
    File plink_snp_gwas

    String phenotype_name
    region bounds

    File all_samples_list

    Float? tol
    Float? snp_p_over_str_p
    Float? res_var
    Float? prior_var
    Boolean best_guess = false
    String prefix = ""
  }

  call finemapping_tasks.susie_choose_vars { input :
    script_dir = script_dir,
    str_assoc_results = my_str_gwas,
    snp_assoc_results = plink_snp_gwas,
    variants_to_filter = snp_vars_to_filter_from_finemapping,
    phenotype_name = phenotype_name,
    bounds = bounds,
  }

  call finemapping_tasks.any_vars_for_susie { input :
    colnames = susie_choose_vars.colnames
  }

  if (any_vars_for_susie.b) {
    call retryable_susie_load_gts.retryable_susie_load_gts { input :
      script_dir = script_dir,
      strs = str_vcf,
      snps = imputed_snp_bgen,
      all_samples = all_samples_list,
      phenotype_samples = phenotype_samples,
      pheno_data = transformed_phenotype_data,
      shared_covars = shared_covars,
      colnames = susie_choose_vars.colnames,
      phenotype_name = phenotype_name,
      bounds = bounds,
      best_guess = best_guess
    }

    call escalating_susie_run.escalating_susie_run { input :
      script_dir = script_dir,
      gts_h5 = retryable_susie_load_gts.gts_h5,
      pheno_residuals_h5 = retryable_susie_load_gts.pheno_residuals_h5,
      varnames_file = susie_choose_vars.colnames,
      tol = tol,
      snp_p_over_str_p = snp_p_over_str_p,
      res_var = res_var,
      prior_var = prior_var,
      colnames = susie_choose_vars.colnames,
      prefix=prefix
    }
  }
  if (!any_vars_for_susie.b) {
    SuSiE_output null_output = object {
      subset: object {
        lbf: "/dev/null",
        lbf_variable: "/dev/null",
        sigma2: "/dev/null",
        V: "/dev/null",
        converged: "/dev/null",
        lfsr: "/dev/null",
        requested_coverage: "/dev/null",
        alpha: "/dev/null",
        colnames: "/dev/null"
      },
      CSs: ["/dev/null"]
    }
  }

  output {
    SuSiE_output susie_output = select_first([escalating_susie_run.susie_output, null_output])
  }
}
