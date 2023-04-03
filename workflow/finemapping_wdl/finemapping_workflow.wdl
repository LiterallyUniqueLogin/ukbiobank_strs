# platform agnostic workflow

version 1.0

import "../gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"

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
    Array[File] ethnic_my_str_gwass,
    File plink_snp_gwas

    String phenotype_name
    Boolean is_binary

    File all_samples_list
  }

  call gwas_tasks.generate_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_regions.data)

  # finemap each region
  scatter (region_idx in range(length(finemapping_regions_tsv) - 1)) {
    region first_pass_bounds = {
      "chrom": read_int(finemapping_regions_tsv[region_idx+1][0]),
      "start": read_int(finemapping_regions_tsv[region_idx+1][1]),
      "end": read_int(finemapping_regions_tsv[region_idx+1][2]),
    }

    call finemap_one_region_workflow.finemap_one_region as original_finemap { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = first_pass_bounds,
      all_samples_list = all_samples_list
    }

    call susie_one_region_workflow.susie_one_region as original_susie { input :
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
  }

  call first_pass_finemapping_df { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    snp_assoc_results = plink_snp_gwas,
    str_assoc_results = my_str_gwas,
    ethnic_str_assoc_results = ethnic_my_str_gwass,
    original_finemap_outputs = original_finemap.finemap_output,
    original_susie_outputs = original_susie.susie_output
  }

  call generate_followup_regions_tsv { input :
    script_dir = script_dir,
    first_pass_df = first_pass_finemapping_df.all_regions_concordance
  }

  Array[Array[String]] followup_finemapping_regions_tsv = read_tsv(generate_followup_regions_tsv.tsv)

  # finemap each region
  scatter (region_idx in range(length(followup_finemapping_regions_tsv) - 1)) {
    region followup_bounds = {
      "chrom": read_int(finemapping_regions_tsv[region_idx+1][0]),
      "start": read_int(finemapping_regions_tsv[region_idx+1][1]),
      "end": read_int(finemapping_regions_tsv[region_idx+1][2]),
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_ratio { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      snp_str_ratio = 4
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_total_prob { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      total_prob = 4,
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_prior_std_derived { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prior_std = 0.0224
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_prior_std_low { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prior_std = 0.005
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_prob_conv_sss_tol { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      prob_conv_sss_tol = 0.0001
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_mac { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      mac = 100
    }

    call finemap_one_region_workflow.finemap_one_region as finemap_threshold { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      shared_covars = shared_covars,
			phenotype_samples = phenotype_samples,
      my_str_gwas = my_str_gwas,
      plink_snp_gwas = plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      threshold = 0.0005
    }

    call susie_one_region_workflow.susie_one_region as { input :
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
      bounds = followup_bounds,
      all_samples_list = all_samples_list,
      sadf
    }
  }

  output {
    FINEMAP_output finemap_output = finemapping_one_region.finemap_output
    SuSiE_output susie_output = finemapping_one_region.susie_output
#    File finemap_snp_file = finemapping_one_region.finemap_snp_file
#    File finemap_log_sss = finemapping_one_region.finemap_log_sss
#    File finemap_config = finemapping_one_region.finemap_config
#    File finemap_creds = finemapping_one_region.finemap_creds
#
#    File susie_lbf = finemapping_one_region.susie_lbf
#    File susie_lbf_variable = finemapping_one_region.susie_lbf_variable
#    File susie_sigma2 = finemapping_one_region.susie_sigma2
#    File susie_V = finemapping_one_region.susie_V
#    File susie_converged = finemapping_one_region.susie_converged
#    File susie_lfsr = finemapping_one_region.susie_lfsr
#    File susie_requested_coverage = finemapping_one_region.susie_requested_coverage
#    File susie_alpha = finemapping_one_region.susie_alpha
#    Array[File] susie_CSs = finemapping_one_region.susie_CSs
  }
}
