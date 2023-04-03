# platform agnostic workflow

version 1.0

import "gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "retryable_finemap_calc_corrs.wdl"
import "retryable_finemap_write_corrs.wdl"
import "retryable_susie_load_gts.wdl"
import "escalating_susie_run.wdl"

# TODO fix chr_lens 21 and 20 the same

workflow finemapping_one_region {

  input {
    String script_dir
    String finemap_command

    File chr_lens
    File str_loci

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    File shared_covars
    File phenotype_samples
    File transformed_phenotype_data

    File my_str_gwas
    File plink_snp_gwas

    String phenotype_name
    Boolean is_binary
    region bounds

    File all_samples_list
  }

  # call FINEMAP
  Int chrom_minus_one = bounds.chrom - 1

  call finemapping_tasks.finemap_write_input_variants { input :
    script_dir = script_dir,
    str_assoc_results = my_str_gwas,
    snp_assoc_results = plink_snp_gwas,
    variants_to_filter = snp_vars_to_filter_from_finemapping[chrom_minus_one],
    phenotype_samples_list = phenotype_samples,
    phenotype = phenotype_name,
    bounds = bounds
  }

  call finemapping_tasks.finemap_load_gts { input :
    script_dir = script_dir,
    strs = str_vcfs[chrom_minus_one],
    snps = imputed_snp_bgens[chrom_minus_one],
    all_samples = all_samples_list,
    phenotype_samples = phenotype_samples, 
    zfile = finemap_write_input_variants.zfile,
    phenotype_name = phenotype_name,
    bounds = bounds
  }

  call retryable_finemap_calc_corrs.retryable_finemap_calc_corrs { input :
    script_dir = script_dir,
    gts_h5 = finemap_load_gts.gts_h5
  }

  call retryable_finemap_write_corrs.retryable_finemap_write_corrs { input :
    script_dir = script_dir,
    lds_h5 = retryable_finemap_calc_corrs.lds_h5
  }

  call finemapping_tasks.finemap_run { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    master = finemap_write_input_variants.master,
    zfile = finemap_write_input_variants.zfile,
    all_variants_ld = retryable_finemap_write_corrs.all_variants_ld
  }

  # TODO need to rerun places where finemapping suggests 20+ causal snps

  # call SuSiE 
 
  call finemapping_tasks.susie_choose_vars { input :
    script_dir = script_dir,
    str_assoc_results = my_str_gwas,
    snp_assoc_results = plink_snp_gwas,
    variants_to_filter = snp_vars_to_filter_from_finemapping[chrom_minus_one],
    phenotype_name = phenotype_name,
    bounds = bounds
  }

  call retryable_susie_load_gts.retryable_susie_load_gts { input :
    script_dir = script_dir,
    strs = str_vcfs[chrom_minus_one],
    snps = imputed_snp_bgens[chrom_minus_one],
    all_samples = all_samples_list,
    phenotype_samples = phenotype_samples,
    pheno_data = transformed_phenotype_data,
    shared_covars = shared_covars,
    colnames = susie_choose_vars.colnames,
    phenotype_name = phenotype_name,
    bounds = bounds,
  }

  call escalating_susie_run.escalating_susie_run { input :
    script_dir = script_dir,
    gts_h5 = retryable_susie_load_gts.gts_h5,
    pheno_residuals_h5 = retryable_susie_load_gts.pheno_residuals_h5,
  }

  output {
    File finemap_snp_file = finemap_run.snp_file
    File finemap_log_sss = finemap_run.log_sss
    File finemap_config = finemap_run.config
    File finemap_creds = finemap_run.creds

    File susie_lbf = escalating_susie_run.lbf
    File susie_lbf_variable = escalating_susie_run.lbf_variable
    File susie_sigma2 = escalating_susie_run.sigma2
    File susie_V = escalating_susie_run.V
    File susie_converged = escalating_susie_run.converged
    File susie_lfsr = escalating_susie_run.lfsr
    File susie_requested_coverage = escalating_susie_run.requested_coverage
    File susie_alpha = escalating_susie_run.alpha
    Array[File] susie_CSs = escalating_susie_run.CSs
  }
}
