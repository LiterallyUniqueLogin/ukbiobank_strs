version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"
import "susie_one_region_workflow.wdl"

workflow first_pass_finemapping {

  input {
    String script_dir
    String finemap_command

    File chr_lens

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    String phenotype_name

    File all_samples_list
    File phenotype_samples

    File shared_covars
    File transformed_phenotype_data

    File my_str_gwas
    File plink_snp_gwas
    Array[File] ethnic_my_str_gwass
  }

  call gwas_tasks.generate_finemapping_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_finemapping_regions.data)

  # finemap each region
  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
    region first_pass_bounds = {
      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
    }
    String first_pass_regions = "~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
    Int first_pass_chroms = first_pass_bounds.chrom

    call finemap_one_region_workflow.finemap_one_region as original_finemap_ { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = first_pass_bounds,
      all_samples_list = all_samples_list
    }
    serializable_FINEMAP_output original_finemap_output_ = original_finemap_.finemap_output.subset
    Array[File] original_finemap_creds_ = original_finemap_.finemap_output.creds

		call susie_one_region_workflow.susie_one_region as original_susie_ { input :
			script_dir = script_dir,
			str_vcfs = str_vcfs,
			imputed_snp_bgens = imputed_snp_bgens,
			snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
			shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
			transformed_phenotype_data = transformed_phenotype_data,
			my_str_gwas = my_str_gwas,
			plink_snp_gwas = plink_snp_gwas,
			phenotype_name = phenotype_name,
			bounds = first_pass_bounds,
			all_samples_list = all_samples_list
		}
		serializable_SuSiE_output original_susie_output_ = original_susie_.susie_output.subset
		Array[File] original_susie_CSs_ = original_susie_.susie_output.CSs
	}

	call finemapping_tasks.first_pass_finemapping_df { input :
		script_dir = script_dir,
		phenotype_name = phenotype_name,
		snp_assoc_results = plink_snp_gwas,
		str_assoc_results = my_str_gwas,
		ethnic_str_assoc_results = ethnic_my_str_gwass,
		original_finemap_outputs = original_finemap_output_,
		original_finemap_creds = original_finemap_creds_,
		original_susie_outputs = original_susie_output_,
		original_susie_CSs = original_susie_CSs_,
		regions = first_pass_regions,
		chroms = first_pass_chroms
	}

  output {
    File regions_tsv = generate_finemapping_regions.data
    File regions_readme = generate_finemapping_regions.readme

    Array[serializable_FINEMAP_output] finemap = original_finemap_output_
    Array[Array[File]] finemap_creds = original_finemap_creds_
		Array[serializable_SuSiE_output] susie = original_susie_output_
		Array[Array[File]] susie_CSs = original_susie_CSs_

		File first_pass_df = first_pass_finemapping_df.all_regions_concordance
  }
}
