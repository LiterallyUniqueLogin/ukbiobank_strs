# platform agnostic workflow

version 1.0

import "gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "retryable_finemap_calc_corrs.wdl"
import "retryable_finemap_write_corrs.wdl"
import "retryable_susie_load_gts.wdl"
import "escalating_susie_run.wdl"

# TODO fix chr_lens 21 and 20 the same

workflow finemapping {

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

    File all_samples_list
  }

  call finemapping_tasks.generate_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_regions.data)

  # finemap each region
  scatter (region_idx in range(length(finemapping_regions_tsv) - 1)) {
    region bounds = {
      "chrom": read_int(region[0]),
      "start": read_int(region[1]),
      "end": read_int(region[2]),
    }

    call finemapping_one_region.finemapping_one_region { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      chr_lens = chr_lens,
      str_loci = str_loci,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      transformed_phenotype_data = transformed_phenotype_data,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      is_binary = is_binary,
      bounds = bounds,
      all_samples_list = all_samples_list
    }
  }

#  output {
#  }
}
