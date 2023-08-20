version 1.0

import "finemapping_tasks.wdl"
import "escalating_finemap_run.wdl"
import "retryable_finemap_calc_corrs.wdl"
import "retryable_finemap_write_corrs.wdl"

workflow finemap_one_region {

  input {
    String script_dir
    String finemap_command

    # one per chrom
    VCF str_vcf
    bgen imputed_snp_bgen
    File snp_vars_to_filter_from_finemapping

    File phenotype_samples

    File? my_str_gwas
    File plink_snp_gwas

    String phenotype_name
    region bounds

    File all_samples_list

    # alternative FINEMAP conditions 
    Float? snp_str_ratio
    Float? total_prob
    Int? mac
    File? snp_macs
    Float? inclusion_threshold
    Float? prior_std
    Float? prob_conv_sss_tol
    String prefix = ""
    Int cache_breaker = 0
  }

  call finemapping_tasks.finemap_write_input_variants { input :
    script_dir = script_dir,
    str_assoc_results = my_str_gwas,
    snp_assoc_results = plink_snp_gwas,
    variants_to_filter = snp_vars_to_filter_from_finemapping,
    phenotype_samples_list = phenotype_samples,
    phenotype = phenotype_name,
    bounds = bounds,
    snp_str_ratio = snp_str_ratio,
    total_prob = total_prob,
    mac = mac,
    snp_macs = snp_macs,
    inclusion_threshold = inclusion_threshold,
  }

  call finemapping_tasks.finemap_load_gts { input :
    script_dir = script_dir,
    strs = str_vcf,
    snps = imputed_snp_bgen,
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

  call escalating_finemap_run.escalating_finemap { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    master = finemap_write_input_variants.master,
    zfile = finemap_write_input_variants.zfile,
    all_variants_ld = retryable_finemap_write_corrs.all_variants_ld,
    prior_snps = defined(snp_str_ratio),
    prior_std = prior_std,
    prob_conv_sss_tol = prob_conv_sss_tol,
    prefix=prefix,
    cache_breaker = cache_breaker
  }

  output {
    File finemap_input_z = escalating_finemap.finemap_input_z
    FINEMAP_output finemap_output = escalating_finemap.finemap_output
  }
}
