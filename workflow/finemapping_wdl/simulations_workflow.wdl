version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "run_simulation_in_region.wdl"
import "finemap_one_region_workflow.wdl"
import "susie_one_region_workflow.wdl"

workflow simulations {
  input {
    String script_dir
    String finemap_command

		File finemapping_regions
    File all_regions_finemapping_df

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[bgen]+ imputed_snp_bgens
    Array[PFiles] imputed_snp_pfiles
    Array[File] snp_vars_to_filter_from_finemapping

    String phenotype_name
    File phenotype_samples
    File transformed_phenotype_data
    File shared_covars

    File true_my_str_gwas
    File true_plink_snp_gwas

    File all_samples_list

    Array[File] snp_macs

    File flank_start_to_start_and_end_pos
  }

	Array[Array[String]] finemapping_regions_tsv = read_tsv(finemapping_regions)

#  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
#    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
#    region bounds = {
#      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
#      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
#      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
#    }
#
#		Int chrom_minus_one = bounds.chrom - 1
#		call gwas_tasks.imputed_snp_frequencies { input :
#			script_dir = script_dir,
#			imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
#			samples = phenotype_samples,
#			all_samples = all_samples_list,
#			bounds = bounds,
#		}
#	}

  call finemapping_tasks.bin_causal_variants_by_frequency { input :
    script_dir = script_dir,
    all_regions_finemapping_df = all_regions_finemapping_df,
    samples_file = phenotype_samples,
    #mac_files = imputed_snp_frequencies.counts
    mac_files = snp_macs
  }

#  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
  scatter (first_pass_region_idx in range(50)) {
    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
    region bounds = {
      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
    }

#    Int first_pass_region_idx = 215
#    region bounds = {
#      "chrom": finemapping_regions_tsv[216][0],
#      "start": finemapping_regions_tsv[216][1],
#      "end": finemapping_regions_tsv[216][2],
#    }

    Int chrom_minus_one = bounds.chrom - 1

    String region_str = "~{bounds.chrom}_~{bounds.start}_~{bounds.end}"

#    ############# simulation with finemap causal variants only SNPs
#
#    call finemap_one_region_workflow.finemap_one_region as finemap_real_no_strs { input :
#      script_dir = script_dir,
#      finemap_command = finemap_command,
#      str_vcf = str_vcfs[chrom_minus_one],
#      imputed_snp_bgen= imputed_snp_bgens[chrom_minus_one],
#      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
#      phenotype_samples = phenotype_samples,
#      plink_snp_gwas = true_plink_snp_gwas,
#      phenotype_name = phenotype_name,
#      bounds = bounds,
#      all_samples_list = all_samples_list
#    }
#
#    call finemapping_tasks.select_causal_variants_from_finemap_output as causal_from_finemap_no_strs { input :
#      script_dir = script_dir,
#      finemap_output_log = finemap_real_no_strs.finemap_output.subset.log_sss,
#      finemap_input_z = finemap_real_no_strs.finemap_input_z,
#      finemap_output_creds = finemap_real_no_strs.finemap_output.creds,
#      prefix = "from_finemap_no_strs_~{region_str}_"
#    }
#
##    scatter (causal_from_finemap_no_strs_seed in range(5)) {
#      Int causal_from_finemap_no_strs_count = 0
#      Int causal_from_finemap_no_strs_seed = causal_from_finemap_no_strs_count
#      call run_simulation_in_region.run_simulation_in_region as finemap_no_strs_simulation { input :
#        causal_vars = causal_from_finemap_no_strs.vars_and_betas[0],
#        causal_betas = causal_from_finemap_no_strs.vars_and_betas[1],
#        seed = causal_from_finemap_no_strs_seed,
#        script_dir = script_dir,
#        finemap_command = finemap_command,
#        str_vcf = str_vcfs[chrom_minus_one],
#        imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
#        imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
#        snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
#        bounds = bounds,
#        phenotype_samples = phenotype_samples,
#        all_samples_list = all_samples_list,
#        prefix = "from_finemap_no_strs_~{region_str}_rep_~{causal_from_finemap_no_strs_count}_"
#      }
##    }
#
#    ############## simulation with finemap causal variants SNPs and STRs
#
#    call finemap_one_region_workflow.finemap_one_region as finemap_real_snps_strs { input :
#      script_dir = script_dir,
#      finemap_command = finemap_command,
#      str_vcf = str_vcfs[chrom_minus_one],
#      imputed_snp_bgen= imputed_snp_bgens[chrom_minus_one],
#      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
#      phenotype_samples = phenotype_samples,
#      my_str_gwas = true_my_str_gwas,
#      plink_snp_gwas = true_plink_snp_gwas,
#      phenotype_name = phenotype_name,
#      bounds = bounds,
#      all_samples_list = all_samples_list
#    }
#
#    call finemapping_tasks.select_causal_variants_from_finemap_output as causal_from_finemap_snps_strs { input :
#      script_dir = script_dir,
#      finemap_output_log = finemap_real_snps_strs.finemap_output.subset.log_sss,
#      finemap_input_z = finemap_real_snps_strs.finemap_input_z,
#      finemap_output_creds = finemap_real_snps_strs.finemap_output.creds,
#      prefix = "from_finemap_snps_strs_~{region_str}_"
#    }
#
#    call finemapping_tasks.any_causal_strs as any_causal_strs_from_finemap_snps_strs { input :
#      vars_and_betas = causal_from_finemap_snps_strs.vars_and_betas
#    }
#
#    if (any_causal_strs_from_finemap_snps_strs.b) {
##      scatter (causal_from_finemap_snps_strs_count in range(5)) {
#        Int causal_from_finemap_snps_strs_count = 0
#        Int causal_from_finemap_snps_strs_seed = causal_from_finemap_snps_strs_count + 5
#        call run_simulation_in_region.run_simulation_in_region as finemap_snps_strs_simulation { input :
#          causal_vars = causal_from_finemap_snps_strs.vars_and_betas[0],
#          causal_betas = causal_from_finemap_snps_strs.vars_and_betas[1],
#          seed = causal_from_finemap_snps_strs_seed,
#          script_dir = script_dir,
#          finemap_command = finemap_command,
#          str_vcf = str_vcfs[chrom_minus_one],
#          imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
#          imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
#          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
#          bounds = bounds,
#          phenotype_samples = phenotype_samples,
#          all_samples_list = all_samples_list,
#          prefix = "from_finemap_snps_strs_~{region_str}_rep_~{causal_from_finemap_snps_strs_count}_"
#        }
##      }
#    }

    #################### simulation with susie causal variants only SNPs
    call susie_one_region_workflow.susie_one_region as susie_real_no_strs { input :
      script_dir = script_dir,
      str_vcf = str_vcfs[chrom_minus_one],
      imputed_snp_bgen= imputed_snp_bgens[chrom_minus_one],
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
      phenotype_samples = phenotype_samples,
      shared_covars = shared_covars,
      transformed_phenotype_data = transformed_phenotype_data,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

		# just to put betas in a convenient place for select_causal_variants_from_susie_output to access
		# could have been done in a cleaner way
		call finemapping_tasks.finemap_write_input_variants as finemap_input_vars_no_strs { input :
			script_dir = script_dir,
			snp_assoc_results = true_plink_snp_gwas,
			variants_to_filter = snp_vars_to_filter_from_finemapping[chrom_minus_one],
			phenotype_samples_list = phenotype_samples,
			phenotype = phenotype_name,
			bounds = bounds,
		}

    call finemapping_tasks.select_causal_variants_from_susie_output as causal_from_susie_no_strs { input :
      script_dir = script_dir,
      colnames = select_first([susie_real_no_strs.susie_output]).subset.colnames,
      alpha = select_first([susie_real_no_strs.susie_output]).subset.alpha,
      finemap_input_z = finemap_input_vars_no_strs.zfile,
      CSes = select_first([susie_real_no_strs.susie_output]).CSs,
      prefix = "from_susie_no_strs_~{region_str}_"
    }

    # need to exit out if no pure CSs
    if (causal_from_susie_no_strs.vars_and_betas[0][0] != "") { # it should contain nothing, instead contains a string of length zero, this is a bug but oh well.
      scatter (causal_from_susie_no_strs_count in range(3)) {
        String susie_no_strs_method = "susie_no_strs"
        Int susie_no_strs_chrom = bounds.chrom
        String susie_no_strs_region = "~{bounds.start}_~{bounds.end}"
        Int susie_no_strs_replicate = causal_from_susie_no_strs_count
        File susie_no_strs_causal_array = causal_from_susie_no_strs.vars_and_betas_f
#        Int causal_from_susie_no_strs_count = 0
        Int causal_from_susie_no_strs_seed = causal_from_susie_no_strs_count + 10
        call run_simulation_in_region.run_simulation_in_region as susie_no_strs_simulation { input :
          causal_vars = causal_from_susie_no_strs.vars_and_betas[0],
          causal_betas = causal_from_susie_no_strs.vars_and_betas[1],
          seed = causal_from_susie_no_strs_seed,
          script_dir = script_dir,
          finemap_command = finemap_command,
          str_vcf = str_vcfs[chrom_minus_one],
          imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
          imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
          bounds = bounds,
          phenotype_samples = phenotype_samples,
          all_samples_list = all_samples_list,
          prefix = "from_susie_no_strs_~{region_str}_rep_~{causal_from_susie_no_strs_count}_"
        }
        serializable_FINEMAP_output susie_no_strs_finemap_output = susie_no_strs_simulation.finemap_results.subset
        Array[File] susie_no_strs_finemap_creds = susie_no_strs_simulation.finemap_results.creds
        serializable_SuSiE_output susie_no_strs_susie_output = susie_no_strs_simulation.susie_results.subset
        Array[File] susie_no_strs_susie_CSs = susie_no_strs_simulation.susie_results.CSs
      }
    }

    ###################### simulation with susie causal variants SNPs and STRs
    call susie_one_region_workflow.susie_one_region as susie_real_snps_strs { input :
      script_dir = script_dir,
      str_vcf = str_vcfs[chrom_minus_one],
      imputed_snp_bgen= imputed_snp_bgens[chrom_minus_one],
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
      phenotype_samples = phenotype_samples,
      shared_covars = shared_covars,
      transformed_phenotype_data = transformed_phenotype_data,
      my_str_gwas = true_my_str_gwas,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

		# just to put betas in a convenient place for select_causal_variants_from_susie_output to access
		# could have been done in a cleaner way
		call finemapping_tasks.finemap_write_input_variants as finemap_input_vars_snps_strs { input :
			script_dir = script_dir,
			str_assoc_results = true_my_str_gwas,
			snp_assoc_results = true_plink_snp_gwas,
			variants_to_filter = snp_vars_to_filter_from_finemapping[chrom_minus_one],
			phenotype_samples_list = phenotype_samples,
			phenotype = phenotype_name,
			bounds = bounds,
		}

    call finemapping_tasks.select_causal_variants_from_susie_output as causal_from_susie_snps_strs { input :
      script_dir = script_dir,
      colnames = select_first([susie_real_snps_strs.susie_output]).subset.colnames,
      alpha = select_first([susie_real_snps_strs.susie_output]).subset.alpha,
      finemap_input_z = finemap_input_vars_snps_strs.zfile,
      CSes = select_first([susie_real_snps_strs.susie_output]).CSs,
      prefix = "from_susie_snps_strs_~{region_str}_"
    }

    call finemapping_tasks.fix_causal_var_str_coords as fixed_causal_from_susie_snps_strs { input :
      vars_and_betas_f = causal_from_susise_snps_strs.vars_and_betas_f,
      flank_start_to_start_and_end_pos = flank_start_to_start_and_end_pos
    }
    
    call finemapping_tasks.any_causal_strs as any_causal_strs_from_susie_snps_strs { input :
      vars_and_betas = causal_from_susie_snps_strs.vars_and_betas
    }

    if (any_causal_strs_from_susie_snps_strs.b) {
      scatter (causal_from_susie_snps_strs_count in range(3)) {
        String susie_snps_strs_method = "susie_snps_strs"
        Int susie_snps_strs_chrom = bounds.chrom
        String susie_snps_strs_region = "~{bounds.start}_~{bounds.end}"
        Int susie_snps_strs_replicate = causal_from_susie_snps_strs_count
        File susie_snps_strs_causal_array = fixed_causal_from_susie_snps_strs.vars_and_betas_f
#        Int causal_from_susie_snps_strs_count = 0
        Int causal_from_susie_snps_strs_seed = causal_from_susie_snps_strs_count + 15
        call run_simulation_in_region.run_simulation_in_region as susie_snps_strs_simulation { input :
          causal_vars = causal_from_susie_snps_strs.vars_and_betas[0],
          causal_betas = causal_from_susie_snps_strs.vars_and_betas[1],
          seed = causal_from_susie_snps_strs_seed,
          script_dir = script_dir,
          finemap_command = finemap_command,
          str_vcf = str_vcfs[chrom_minus_one],
          imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
          imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
          bounds = bounds,
          phenotype_samples = phenotype_samples,
          all_samples_list = all_samples_list,
          prefix = "from_susie_snps_strs_~{region_str}_rep_~{causal_from_susie_snps_strs_count}_"
        }

        serializable_FINEMAP_output susie_snps_strs_finemap_output = susie_snps_strs_simulation.finemap_results.subset
        Array[File] susie_snps_strs_finemap_creds = susie_snps_strs_simulation.finemap_results.creds
        serializable_SuSiE_output susie_snps_strs_susie_output = susie_snps_strs_simulation.susie_results.subset
        Array[File] susie_snps_strs_susie_CSs = susie_snps_strs_simulation.susie_results.CSs
      }
    }

    ################### simulation with random causal variants

    scatter (causal_from_random_one_var_count in range(3)) {
      String random_one_var_method = "random_one_var"
      Int random_one_var_chrom = bounds.chrom
      String random_one_var_region = "~{bounds.start}_~{bounds.end}"
      Int random_one_var_replicate = causal_from_random_one_var_count
#      Int causal_from_random_one_var_count = 0
      Int causal_from_random_one_var_seed = causal_from_random_one_var_count + 20
      call finemapping_tasks.select_random_causal_variants as one_random_causal { input :
        script_dir = script_dir,
        n_vars_to_choose = 1,
        seed = causal_from_random_one_var_seed,
        samples_file = phenotype_samples,
        snp_macs = snp_macs[first_pass_region_idx],
        bin_weights = [
          bin_causal_variants_by_frequency.data[1][5],
          bin_causal_variants_by_frequency.data[2][5], 
          bin_causal_variants_by_frequency.data[3][5], 
          bin_causal_variants_by_frequency.data[4][5], 
        ],
        bin_effects = bin_causal_variants_by_frequency.bin_effects,
        prefix = "one_random_causal_~{region_str}_rep_~{causal_from_random_one_var_count}_"
      }

      call run_simulation_in_region.run_simulation_in_region as one_random_causal_simulation { input :
        causal_vars = one_random_causal.vars_and_betas[0],
        causal_betas = one_random_causal.vars_and_betas[1],
        seed = causal_from_random_one_var_seed,
        script_dir = script_dir,
        finemap_command = finemap_command,
        str_vcf = str_vcfs[chrom_minus_one],
        imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
        imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
        snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
        bounds = bounds,
        phenotype_samples = phenotype_samples,
        all_samples_list = all_samples_list,
        prefix = "one_random_causal_~{region_str}_rep_~{causal_from_random_one_var_count}_"
      }

      serializable_FINEMAP_output one_random_causal_finemap_output = one_random_causal_simulation.finemap_results.subset
      Array[File] one_random_causal_finemap_creds = one_random_causal_simulation.finemap_results.creds
      serializable_SuSiE_output one_random_causal_susie_output = one_random_causal_simulation.susie_results.subset
      Array[File] one_random_causal_susie_CSs = one_random_causal_simulation.susie_results.CSs
    }

    scatter (causal_from_random_two_var_count in range(3)) {
      String random_two_var_method = "random_two_var"
      Int random_two_var_chrom = bounds.chrom
      String random_two_var_region = "~{bounds.start}_~{bounds.end}"
      Int random_two_var_replicate = causal_from_random_two_var_count
#      Int causal_from_random_two_var_count = 0
      Int causal_from_random_two_var_seed = causal_from_random_two_var_count + 25
      call finemapping_tasks.select_random_causal_variants as two_random_causal { input :
        script_dir = script_dir,
        n_vars_to_choose = 2,
        seed = causal_from_random_two_var_seed,
        samples_file = phenotype_samples,
        snp_macs = snp_macs[first_pass_region_idx],
        bin_weights = [
          bin_causal_variants_by_frequency.data[1][5],
          bin_causal_variants_by_frequency.data[2][5], 
          bin_causal_variants_by_frequency.data[3][5], 
          bin_causal_variants_by_frequency.data[4][5], 
        ],
        bin_effects = bin_causal_variants_by_frequency.bin_effects,
        prefix = "two_random_causal_~{region_str}_rep_~{causal_from_random_two_var_count}_"
      }

      call run_simulation_in_region.run_simulation_in_region as two_random_causal_simulation { input :
        causal_vars = two_random_causal.vars_and_betas[0],
        causal_betas = two_random_causal.vars_and_betas[1],
        seed = causal_from_random_two_var_seed,
        script_dir = script_dir,
        finemap_command = finemap_command,
        str_vcf = str_vcfs[chrom_minus_one],
        imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
        imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
        snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
        bounds = bounds,
        phenotype_samples = phenotype_samples,
        all_samples_list = all_samples_list,
        prefix = "two_random_causal_~{region_str}_rep_~{causal_from_random_two_var_count}_"
      }
      
      serializable_FINEMAP_output two_random_causal_finemap_output = two_random_causal_simulation.finemap_results.subset
      Array[File] two_random_causal_finemap_creds = two_random_causal_simulation.finemap_results.creds
      serializable_SuSiE_output two_random_causal_susie_output = two_random_causal_simulation.susie_results.subset
      Array[File] two_random_causal_susie_CSs = two_random_causal_simulation.susie_results.CSs
    }
    
    scatter (causal_from_random_three_var_count in range(3)) {
      String random_three_var_method = "random_three_var"
      Int random_three_var_chrom = bounds.chrom
      String random_three_var_region = "~{bounds.start}_~{bounds.end}"
      Int random_three_var_replicate = causal_from_random_three_var_count
#      Int causal_from_random_three_var_count = 0
      Int causal_from_random_three_var_seed = causal_from_random_three_var_count + 35
      call finemapping_tasks.select_random_causal_variants as three_random_causal { input :
        script_dir = script_dir,
        n_vars_to_choose = 3,
        seed = causal_from_random_three_var_seed,
        samples_file = phenotype_samples,
        snp_macs = snp_macs[first_pass_region_idx],
        bin_weights = [
          bin_causal_variants_by_frequency.data[1][5],
          bin_causal_variants_by_frequency.data[2][5], 
          bin_causal_variants_by_frequency.data[3][5], 
          bin_causal_variants_by_frequency.data[4][5], 
        ],
        bin_effects = bin_causal_variants_by_frequency.bin_effects,
        prefix = "three_random_causal_~{region_str}_rep_~{causal_from_random_three_var_count}_"
      }

      call run_simulation_in_region.run_simulation_in_region as three_random_causal_simulation { input :
        causal_vars = three_random_causal.vars_and_betas[0],
        causal_betas = three_random_causal.vars_and_betas[1],
        seed = causal_from_random_three_var_seed,
        script_dir = script_dir,
        finemap_command = finemap_command,
        str_vcf = str_vcfs[chrom_minus_one],
        imputed_snp_bgen = imputed_snp_bgens[chrom_minus_one],
        imputed_snp_pfile = imputed_snp_pfiles[chrom_minus_one],
        snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping[chrom_minus_one],
        bounds = bounds,
        phenotype_samples = phenotype_samples,
        all_samples_list = all_samples_list,
        prefix = "three_random_causal_~{region_str}_rep_~{causal_from_random_three_var_count}_"
      }

      serializable_FINEMAP_output three_random_causal_finemap_output = three_random_causal_simulation.finemap_results.subset
      Array[File] three_random_causal_finemap_creds = three_random_causal_simulation.finemap_results.creds
      serializable_SuSiE_output three_random_causal_susie_output = three_random_causal_simulation.susie_results.subset
      Array[File] three_random_causal_susie_CSs = three_random_causal_simulation.susie_results.CSs
    }
  }

  call finemapping_tasks.compile_simulations_df { input :
    script_dir = script_dir,
    methods = flatten(flatten([select_all(susie_no_strs_method), select_all(susie_snps_strs_method), random_one_var_method, random_two_var_method, random_three_var_method])),
    chroms = flatten(flatten([select_all(susie_no_strs_chrom), select_all(susie_snps_strs_chrom), random_one_var_chrom, random_two_var_chrom, random_three_var_chrom])),
    regions = flatten(flatten([select_all(susie_no_strs_region), select_all(susie_snps_strs_region), random_one_var_region, random_two_var_region, random_three_var_region])),
    replicates = flatten(flatten([select_all(susie_no_strs_replicate), select_all(susie_snps_strs_replicate), random_one_var_replicate, random_two_var_replicate, random_three_var_replicate])),
    causal_vars_and_betas = flatten(flatten([select_all(susie_no_strs_causal_array), select_all(susie_snps_strs_causal_array), one_random_causal.vars_and_betas_f, two_random_causal.vars_and_betas_f, three_random_causal.vars_and_betas_f])),
    snp_assocs = flatten(flatten([select_all(susie_no_strs_simulation.snp_assocs), select_all(susie_snps_strs_simulation.snp_assocs), one_random_causal_simulation.snp_assocs, two_random_causal_simulation.snp_assocs, three_random_causal_simulation.snp_assocs])),
    str_assocs = flatten(flatten([select_all(susie_no_strs_simulation.str_assocs), select_all(susie_snps_strs_simulation.str_assocs), one_random_causal_simulation.str_assocs, two_random_causal_simulation.str_assocs, three_random_causal_simulation.str_assocs])),
    finemap_outputs = flatten(flatten([select_all(susie_no_strs_finemap_output), select_all(susie_snps_strs_finemap_output), one_random_causal_finemap_output, two_random_causal_finemap_output, three_random_causal_finemap_output])),
    finemap_creds = flatten(flatten([select_all(susie_no_strs_finemap_creds), select_all(susie_snps_strs_finemap_creds), one_random_causal_finemap_creds, two_random_causal_finemap_creds, three_random_causal_finemap_creds])),
    susie_outputs = flatten(flatten([select_all(susie_no_strs_susie_output), select_all(susie_snps_strs_susie_output), one_random_causal_susie_output, two_random_causal_susie_output, three_random_causal_susie_output])),
    susie_CSs = flatten(flatten([select_all(susie_no_strs_susie_CSs), select_all(susie_snps_strs_susie_CSs), one_random_causal_susie_CSs, two_random_causal_susie_CSs, three_random_causal_susie_CSs])),
  }

  output {
    File bins = bin_causal_variants_by_frequency.data_file
    Array[File] bin_effects = bin_causal_variants_by_frequency.bin_effects

    # Outer array for each is per region, second is per seed
    # for causal vars and betas first row in third array is causal var and second is betas
    # select_alls check for regions that might be filtered (no STRs when STRs were asked for, no pure CSs for SuSiE)

#    File causal_vars_and_betas_from_finemap_no_strs = causal_from_finemap_no_strs.vars_and_betas_f
#    FINEMAP_output finemap_no_strs_simulation_finemap_results = finemap_no_strs_simulation.finemap_results
#    SuSiE_output finemap_no_strs_simulation_susie_results = finemap_no_strs_simulation.susie_results

#    File causal_vars_and_betas_from_finemap_snps_strs = causal_from_finemap_snps_strs.vars_and_betas_f
#    FINEMAP_output finemap_snps_strs_simulation_finemap_results = select_first([finemap_snps_strs_simulation.finemap_results])
#    SuSiE_output finemap_snps_strs_simulation_susie_results = select_first([finemap_snps_strs_simulation.susie_results])

#    File causal_vars_and_betas_from_susie_no_strs = causal_from_susie_no_strs.vars_and_betas_f
#    FINEMAP_output susie_no_strs_simulation_finemap_results = select_first([susie_no_strs_simulation.finemap_results])
#    SuSiE_output susie_no_strs_simulation_susie_results = select_first([susie_no_strs_simulation.susie_results])
#
#    File causal_vars_and_betas_from_susie_snps_strs = causal_from_susie_snps_strs.vars_and_betas_f
#    FINEMAP_output susie_snps_strs_simulation_finemap_results = select_first([susie_snps_strs_simulation.finemap_results])
#    SuSiE_output susie_snps_strs_simulation_susie_results = select_first([susie_snps_strs_simulation.susie_results])
#
#    File causal_vars_and_betas_from_random_one_var = one_random_causal.vars_and_betas_f
#    FINEMAP_output one_random_causal_simulation_finemap_results = one_random_causal_simulation.finemap_results
#    SuSiE_output one_random_causal_simulation_susie_results = one_random_causal_simulation.susie_results
#
#    File causal_vars_and_betas_from_random_two_var = two_random_causal.vars_and_betas_f
#    FINEMAP_output two_random_causal_simulation_finemap_results = two_random_causal_simulation.finemap_results
#    SuSiE_output two_random_causal_simulation_susie_results = two_random_causal_simulation.susie_results
#
#    File causal_vars_and_betas_from_random_three_var = three_random_causal.vars_and_betas_f
#    FINEMAP_output three_random_causal_simulation_finemap_results = three_random_causal_simulation.finemap_results
#    SuSiE_output three_random_causal_simulation_susie_results = three_random_causal_simulation.susie_results

    Array[File] causal_vars_and_betas_from_susie_no_strs = causal_from_susie_no_strs.vars_and_betas_f
    Array[Array[FINEMAP_output]] susie_no_strs_simulation_finemap_results = select_all(susie_no_strs_simulation.finemap_results)
    Array[Array[SuSiE_output]] susie_no_strs_simulation_susie_results = select_all(susie_no_strs_simulation.susie_results)

    Array[File] causal_vars_and_betas_from_susie_snps_strs = causal_from_susie_snps_strs.vars_and_betas_f
    Array[Array[FINEMAP_output]] susie_snps_strs_simulation_finemap_results = select_all(susie_snps_strs_simulation.finemap_results)
    Array[Array[SuSiE_output]] susie_snps_strs_simulation_susie_results = select_all(susie_snps_strs_simulation.susie_results)

    Array[Array[File]] causal_vars_and_betas_from_random_one_var = one_random_causal.vars_and_betas_f
    Array[Array[FINEMAP_output]] one_random_causal_simulation_finemap_results = one_random_causal_simulation.finemap_results
    Array[Array[SuSiE_output]] one_random_causal_simulation_susie_results = one_random_causal_simulation.susie_results

    Array[Array[File]] causal_vars_and_betas_from_random_two_var = two_random_causal.vars_and_betas_f
    Array[Array[FINEMAP_output]] two_random_causal_simulation_finemap_results = two_random_causal_simulation.finemap_results
    Array[Array[SuSiE_output]] two_random_causal_simulation_susie_results = two_random_causal_simulation.susie_results

    Array[Array[File]] causal_vars_and_betas_from_random_three_var = three_random_causal.vars_and_betas_f
    Array[Array[FINEMAP_output]] three_random_causal_simulation_finemap_results = three_random_causal_simulation.finemap_results
    Array[Array[SuSiE_output]] three_random_causal_simulation_susie_results = three_random_causal_simulation.susie_results

    File simulations_df = compile_simulations_df.df
  }
}
