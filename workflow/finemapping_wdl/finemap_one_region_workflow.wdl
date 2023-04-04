# platform agnostic workflow

version 1.0

import "finemapping_tasks.wdl"
import "retryable_finemap_calc_corrs.wdl"
import "retryable_finemap_write_corrs.wdl"

# TODO fix chr_lens 21 and 20 the same

workflow finemap_one_region {

  input {
    String script_dir
    String finemap_command

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    File shared_covars
    File phenotype_samples

    File my_str_gwas
    File plink_snp_gwas

    String phenotype_name
    region bounds

    File all_samples_list

    # alternative FINEMAP conditions 
    Float? snp_str_ratio
    Float? total_prob
    Int? mac
    Float? inclusion_threshold
    Float? prior_std
    Float? prob_conv_sss_tol
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
    bounds = bounds,
    snp_str_ratio = snp_str_ratio,
    total_prob = total_prob,
    mac = mac,
    inclusion_threshold = inclusion_threshold
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
    all_variants_ld = retryable_finemap_write_corrs.all_variants_ld,
    prior_std = prior_std,
    prob_conv_sss_tol = prob_conv_sss_tol
  }

  # TODO need to rerun places where finemapping suggests 20+ causal snps

  output {
    FINEMAP_output finemap_output = finemap_run.finemap_output
  }
}
