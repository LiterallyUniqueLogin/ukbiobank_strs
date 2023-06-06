version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"
import "susie_one_region_workflow.wdl"

workflow simulations {
  input {
    String script_dir
    String finemap_command

    File chr_lens

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
  }

  Int max_seed = 5

  call gwas_tasks.generate_finemapping_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = true_plink_snp_gwas,
    str_assoc_results = true_my_str_gwas
  }

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_finemapping_regions.data)

  scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
    Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
    region bounds = {
      "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
      "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
      "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
    }

#    region bounds = {
#      "chrom": finemapping_regions_tsv[216][0],
#      "start": finemapping_regions_tsv[216][1],
#      "end": finemapping_regions_tsv[216][2],
#    }
    Int chrom_minus_one = bounds.chrom - 1

    ############# simulation with finemap causal variants only SNPs

    call finemap_one_region_workflow.finemap_one_region as finemap_real_no_strs { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      phenotype_samples = phenotype_samples,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

    call finemapping_tasks.select_causal_variants_from_finemap_output as causal_from_finemap_no_strs { input :
      script_dir = script_dir,
      finemap_output_log = finemap_real_no_strs.finemap_output.subset.log_sss,
      finemap_input_z = finemap_real_no_strs.finemap_input_z,
      finemap_output_creds = finemap_real_no_strs.finemap_output.creds
    }

    scatter (seed in range(max_seed)) {
      call finemapping_tasks.simulate_phenotype as simulate_from_finemap_no_strs { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[chrom_minus_one],
        snp_bgen = imputed_snp_bgens[chrom_minus_one],
        all_samples_list = all_samples_list,
        samples_list = phenotype_samples,
        chrom = bounds.chrom,
        seed = seed,
        causal_vars = causal_from_finemap_no_strs.vars_and_betas[0],
        causal_betas = causal_from_finemap_no_strs.vars_and_betas[1],
      }

      call gwas_tasks.regional_my_str_gwas as str_gwas_from_finemap_no_strs { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[chrom_minus_one],
        transformed_phenotype = simulate_from_finemap_no_strs.phenotype,
        all_samples_list = all_samples_list,
        is_binary = false,
        binary_type = "blah",
        bounds = bounds,
        phenotype_name = "simulation",
        no_details = true,
      }

      call gwas_tasks.prep_plink_input as prep_plink_from_finemap_no_strs { input :
        script_dir = script_dir,
        transformed_phenotype = simulate_from_finemap_no_strs.phenotype,
        pheno_covar_names = "/dev/null",
        is_binary = false,
        binary_type = "blah",
        phenotype_name = "simulation",
      }

      call gwas_tasks.plink_snp_association as plink_snp_from_finemap_no_strs { input :
        script_dir = script_dir,
        plink_command = "plink2",
        imputed_snp_p_file = imputed_snp_pfiles[chrom_minus_one],
        pheno_data = prep_plink_from_finemap_no_strs.data,
        chrom = bounds.chrom,
        phenotype_name = "simulation",
        binary_type = "linear",
        start = bounds.start,
        end = bounds.end
      }

      call finemap_one_region_workflow.finemap_one_region as finemap_from_finemap_no_strs { input :
        script_dir = script_dir,
        finemap_command = finemap_command,
        str_vcfs = str_vcfs,
        imputed_snp_bgens = imputed_snp_bgens,
        snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
        phenotype_samples = phenotype_samples,
        my_str_gwas = str_gwas_from_finemap_no_strs.data,
        plink_snp_gwas = plink_snp_from_finemap_no_strs.data,
        phenotype_name = "simulation",
        bounds = bounds,
        all_samples_list = all_samples_list
      }
    }

    ############## simulation with finemap causal variants SNPs and STRs
    # TODO need to skip this if no STRs

    call finemap_one_region_workflow.finemap_one_region as finemap_real_snps_strs { input :
      script_dir = script_dir,
      finemap_command = finemap_command,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      phenotype_samples = phenotype_samples,
      my_str_gwas = true_my_str_gwas,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

    call finemapping_tasks.select_causal_variants_from_finemap_output as causal_from_finemap_snps_strs { input :
      script_dir = script_dir,
      finemap_output_log = finemap_real_snps_strs.finemap_output.subset.log_sss,
      finemap_input_z = finemap_real_snps_strs.finemap_input_z,
      finemap_output_creds = finemap_real_snps_strs.finemap_output.creds
    }

    call finemapping_tasks.any_causal_strs as any_causal_strs_from_finemap_snps_strs { input :
      vars_and_betas = causal_from_finemap_snps_strs.vars_and_betas
    }

    if (any_causal_strs_from_finemap_snps_strs.b) {
      scatter (seed in range(max_seed)) {
        call finemapping_tasks.simulate_phenotype as simulate_from_finemap_snps_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          snp_bgen = imputed_snp_bgens[chrom_minus_one],
          all_samples_list = all_samples_list,
          samples_list = phenotype_samples,
          chrom = bounds.chrom,
          seed = seed,
          causal_vars = causal_from_finemap_snps_strs.vars_and_betas[0],
          causal_betas = causal_from_finemap_snps_strs.vars_and_betas[1],
        }

        call gwas_tasks.regional_my_str_gwas as str_gwas_from_finemap_snps_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          transformed_phenotype = simulate_from_finemap_snps_strs.phenotype,
          all_samples_list = all_samples_list,
          is_binary = false,
          binary_type = "blah",
          bounds = bounds,
          phenotype_name = "simulation",
          no_details = true,
        }

        call gwas_tasks.prep_plink_input as prep_plink_from_finemap_snps_strs { input :
          script_dir = script_dir,
          transformed_phenotype = simulate_from_finemap_snps_strs.phenotype,
          pheno_covar_names = "/dev/null",
          is_binary = false,
          binary_type = "blah",
          phenotype_name = "simulation",
        }

        call gwas_tasks.plink_snp_association as plink_snp_from_finemap_snps_strs { input :
          script_dir = script_dir,
          plink_command = "plink2",
          imputed_snp_p_file = imputed_snp_pfiles[chrom_minus_one],
          pheno_data = prep_plink_from_finemap_snps_strs.data,
          chrom = bounds.chrom,
          phenotype_name = "simulation",
          binary_type = "linear",
          start = bounds.start,
          end = bounds.end
        }

        call finemap_one_region_workflow.finemap_one_region as finemap_from_finemap_snps_strs { input :
          script_dir = script_dir,
          finemap_command = finemap_command,
          str_vcfs = str_vcfs,
          imputed_snp_bgens = imputed_snp_bgens,
          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
          phenotype_samples = phenotype_samples,
          my_str_gwas = str_gwas_from_finemap_snps_strs.data,
          plink_snp_gwas = plink_snp_from_finemap_snps_strs.data,
          phenotype_name = "simulation",
          bounds = bounds,
          all_samples_list = all_samples_list
        }
      }
    }

    ################### simulation with random finemap causal variants

    #################### simulation with susie causal variants only SNPs
    # TODO need to exit out if no pure CSs
    call susie_one_region_workflow.susie_one_region as susie_real_no_strs { input :
      script_dir = script_dir,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      phenotype_samples = phenotype_samples,
      shared_covars = shared_covars,
      transformed_phenotype_data = transformed_phenotype_data,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

    call finemapping_tasks.select_causal_variants_from_susie_output as causal_from_susie_no_strs { input :
      script_dir = script_dir,
      colnames = susie_real_no_strs.susie_output.subset.colnames,
      finemap_input_z = finemap_real_no_strs.finemap_input_z,
      CSes = susie_real_no_strs.susie_output.CSs
    }

    if (causal_from_susie_no_strs.vars_and_betas[0][0] != "") { # it should contain nothing, instead contains a string of length zero, this is a bug but oh well.
      scatter (seed in range(max_seed)) {
        call finemapping_tasks.simulate_phenotype as simulate_from_susie_no_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          snp_bgen = imputed_snp_bgens[chrom_minus_one],
          all_samples_list = all_samples_list,
          samples_list = phenotype_samples,
          chrom = bounds.chrom,
          seed = seed,
          causal_vars = causal_from_susie_no_strs.vars_and_betas[0],
          causal_betas = causal_from_susie_no_strs.vars_and_betas[1],
        }

        call gwas_tasks.regional_my_str_gwas as str_gwas_from_susie_no_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          transformed_phenotype = simulate_from_susie_no_strs.phenotype,
          all_samples_list = all_samples_list,
          is_binary = false,
          binary_type = "blah",
          bounds = bounds,
          phenotype_name = "simulation",
          no_details = true,
        }

        call gwas_tasks.prep_plink_input as prep_plink_from_susie_no_strs { input :
          script_dir = script_dir,
          transformed_phenotype = simulate_from_susie_no_strs.phenotype,
          pheno_covar_names = "/dev/null",
          is_binary = false,
          binary_type = "blah",
          phenotype_name = "simulation",
        }

        call gwas_tasks.plink_snp_association as plink_snp_from_susie_no_strs { input :
          script_dir = script_dir,
          plink_command = "plink2",
          imputed_snp_p_file = imputed_snp_pfiles[chrom_minus_one],
          pheno_data = prep_plink_from_susie_no_strs.data,
          chrom = bounds.chrom,
          phenotype_name = "simulation",
          binary_type = "linear",
          start = bounds.start,
          end = bounds.end
        }

        call susie_one_region_workflow.susie_one_region as susie_from_susie_no_strs { input :
          script_dir = script_dir,
          str_vcfs = str_vcfs,
          imputed_snp_bgens = imputed_snp_bgens,
          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
          phenotype_samples = phenotype_samples,
          transformed_phenotype_data = simulate_from_susie_no_strs.phenotype,
          my_str_gwas = str_gwas_from_susie_no_strs.data,
          plink_snp_gwas = plink_snp_from_susie_no_strs.data,
          phenotype_name = "simulation",
          bounds = bounds,
          all_samples_list = all_samples_list
        }
      }
    }

    ###################### simulation with susie causal variants SNPs and STRs
    call susie_one_region_workflow.susie_one_region as susie_real_snps_strs { input :
      script_dir = script_dir,
      str_vcfs = str_vcfs,
      imputed_snp_bgens = imputed_snp_bgens,
      snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
      phenotype_samples = phenotype_samples,
      shared_covars = shared_covars,
      transformed_phenotype_data = transformed_phenotype_data,
      my_str_gwas = true_my_str_gwas,
      plink_snp_gwas = true_plink_snp_gwas,
      phenotype_name = phenotype_name,
      bounds = bounds,
      all_samples_list = all_samples_list
    }

    call finemapping_tasks.select_causal_variants_from_susie_output as causal_from_susie_snps_strs { input :
      script_dir = script_dir,
      colnames = susie_real_snps_strs.susie_output.subset.colnames,
      finemap_input_z = finemap_real_snps_strs.finemap_input_z,
      CSes = susie_real_snps_strs.susie_output.CSs
    }
    
    call finemapping_tasks.any_causal_strs as any_causal_strs_from_susie_snps_strs { input :
      vars_and_betas = causal_from_susie_snps_strs.vars_and_betas
    }

    if (any_causal_strs_from_susie_snps_strs.b) {
      scatter (seed in range(max_seed)) {
        call finemapping_tasks.simulate_phenotype as simulate_from_susie_snps_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          snp_bgen = imputed_snp_bgens[chrom_minus_one],
          all_samples_list = all_samples_list,
          samples_list = phenotype_samples,
          chrom = bounds.chrom,
          seed = seed,
          causal_vars = causal_from_susie_snps_strs.vars_and_betas[0],
          causal_betas = causal_from_susie_snps_strs.vars_and_betas[1],
        }

        call gwas_tasks.regional_my_str_gwas as str_gwas_from_susie_snps_strs { input :
          script_dir = script_dir,
          str_vcf = str_vcfs[chrom_minus_one],
          transformed_phenotype = simulate_from_susie_snps_strs.phenotype,
          all_samples_list = all_samples_list,
          is_binary = false,
          binary_type = "blah",
          bounds = bounds,
          phenotype_name = "simulation",
          no_details = true,
        }

        call gwas_tasks.prep_plink_input as prep_plink_from_susie_snps_strs { input :
          script_dir = script_dir,
          transformed_phenotype = simulate_from_susie_snps_strs.phenotype,
          pheno_covar_names = "/dev/null",
          is_binary = false,
          binary_type = "blah",
          phenotype_name = "simulation",
        }

        call gwas_tasks.plink_snp_association as plink_snp_from_susie_snps_strs { input :
          script_dir = script_dir,
          plink_command = "plink2",
          imputed_snp_p_file = imputed_snp_pfiles[chrom_minus_one],
          pheno_data = prep_plink_from_susie_snps_strs.data,
          chrom = bounds.chrom,
          phenotype_name = "simulation",
          binary_type = "linear",
          start = bounds.start,
          end = bounds.end
        }

        call susie_one_region_workflow.susie_one_region as susie_from_susie_snps_strs { input :
          script_dir = script_dir,
          str_vcfs = str_vcfs,
          imputed_snp_bgens = imputed_snp_bgens,
          snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
          phenotype_samples = phenotype_samples,
          transformed_phenotype_data = simulate_from_susie_snps_strs.phenotype,
          my_str_gwas = str_gwas_from_susie_snps_strs.data,
          plink_snp_gwas = plink_snp_from_susie_snps_strs.data,
          phenotype_name = "simulation",
          bounds = bounds,
          all_samples_list = all_samples_list
        }
      }
    }

    ################## simulation with random susie causal variants
  }

  output {
    File regions_tsv = generate_finemapping_regions.data
    File regions_readme = generate_finemapping_regions.readme

    # Outer array for each is per region
    # select_alls check for regions that might be filtered (no STRs when STRs were asked for, no pure CSs for SuSiE)

    # first row: array of causal vars second row: array of betas
    Array[Array[Array[String]]] causal_vars_and_betas_from_finemap_no_strs = causal_from_finemap_no_strs.vars_and_betas
    # array is per seed
    Array[Array[FINEMAP_output]] finemap_results_from_simulation_on_finemap_no_strs = finemap_from_finemap_no_strs.finemap_output
    # first row: array of causal vars second row: array of betas
    Array[Array[Array[String]]] causal_vars_and_betas_from_finemap_snps_strs = causal_from_finemap_snps_strs.vars_and_betas
    # array is per seed
    Array[Array[FINEMAP_output]] finemap_results_from_simulation_on_finemap_snps_strs = select_all(finemap_from_finemap_snps_strs.finemap_output)

    # first row: array of causal vars second row: array of betas
    Array[Array[Array[String]]] causal_vars_and_betas_from_susie_no_strs = causal_from_susie_no_strs.vars_and_betas
    # array is per seed
    Array[Array[SuSiE_output]] susie_results_from_simulation_on_susie_no_strs = select_all(susie_from_susie_no_strs.susie_output)
    # first row: array of causal vars second row: array of betas
    Array[Array[Array[String]]] causal_vars_and_betas_from_susie_snps_strs = causal_from_susie_snps_strs.vars_and_betas
    # array is per seed
    Array[Array[SuSiE_output]] susie_results_from_simulation_on_susie_snps_strs = select_all(susie_from_susie_snps_strs.susie_output)
  }
}
