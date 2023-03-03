# platform agnostic workflow

version 1.0

import "gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_calc_corrs_retryable.wdl"
import "finemap_write_corrs_retryable.wdl"

# TODO fix chr_lens 21 and 20 the same

#struct Phenotype_Input {
#  String name
#  File sc
#  Array[String] categorical_covariate_names
#  Array[File] categorical_covariate_scs
#  Boolean is_binary
#  Boolean is_zero_one_neg_nan
#}

workflow main {

  input {
    String script_dir
    String PRIMUS_command
    String plink_command = "plink2"
    String finemap_command

    File chr_lens
    File str_loci

    Array[VCF]+ str_vcfs
    #Array[PFiles]+ imputed_snp_p_files
    Array[bgen]+ imputed_snp_bgens
    Array[File] snp_vars_to_filter_from_finemapping

    String phenotype_name
    Array[String] categorical_covariate_names
    Array[File] categorical_covariate_scs
    Boolean is_binary
    Boolean is_zero_one_neg_nan
    String date_of_most_recent_first_occurrence_update

    File fam_file # task for generating this?
    File all_samples_list
    File withdrawn_sample_list
    File kinship

    # data showcase files
    File sc_white_brits
    File sc_ethnicity_self_report
    File sc_sex_aneuploidy
    File sc_genetic_sex
    File sc_reported_sex
    File sc_kinship_count
    File sc_assessment_ages
    File sc_pcs
    File sc_year_of_birth
    File sc_month_of_birth
    File sc_date_of_death
    File sc_phenotype
  }

  Array[Int] chroms = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

  call gwas_tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = script_dir,
    sc = sc_white_brits
  }

  call gwas_tasks.ethnic_sample_lists as ethnic_sample_lists_task { input: 
    script_dir = script_dir,
    white_brits_sample_list = white_brits_sample_list.data,
    sc_ethnicity_self_report = sc_ethnicity_self_report
  }

  Array[String] ethnicities = ethnic_sample_lists_task.ethnicities
  Array[String] all_ethnicities = flatten([['white_british'], ethnicities])
  Array[File] ethnic_sample_lists = ethnic_sample_lists_task.sample_lists
  Array[File] all_sample_lists = flatten([
    [white_brits_sample_list.data], ethnic_sample_lists
  ])

  call gwas_tasks.write_sample_list as sex_aneuploidy_sample_list { input:
    script_dir = script_dir,
    sc = sc_sex_aneuploidy
  }

  call gwas_tasks.sex_mismatch_sample_list { input:
    script_dir = script_dir,
    sc_genetic_sex = sc_genetic_sex,
    sc_reported_sex = sc_reported_sex
  }
  
  call gwas_tasks.write_sample_list as low_genotyping_quality_sample_list { input:
    script_dir = script_dir,
    sc = sc_kinship_count,
    value = -1
  }

  scatter (sample_list in all_sample_lists) {
    call gwas_tasks.qced_sample_list as all_qced_sample_lists { input:
      script_dir = script_dir,
      unqced_sample_list = sample_list,
      withdrawn_sample_list = withdrawn_sample_list,
      sex_aneuploidy_sample_list = sex_aneuploidy_sample_list.data,
      sex_mismatch_sample_list = sex_mismatch_sample_list.data,
      low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data
    }
  }

  scatter (sample_list in all_qced_sample_lists.data) {
    call gwas_tasks.unrelated_samples as ethnicity_unrelated_samples { input:
      script_dir = script_dir,
      PRIMUS_command = PRIMUS_command,
      kinship = kinship,
      sample_list = sample_list
    }
  }

  call gwas_tasks.load_shared_covars { input:
    script_dir = script_dir,
    fam_file = fam_file,
    sc_pcs = sc_pcs,
    sc_assessment_ages = sc_assessment_ages
  }

  # get qced unrelated phenotype sample list and transformed phenotype data for each ethnicity
  scatter (sample_list in all_qced_sample_lists.data) {
    if (!is_binary) {
      call gwas_tasks.load_continuous_phenotype { input :
        script_dir = script_dir,
        sc = sc_phenotype,
        qced_sample_list = all_qced_sample_lists.data,
        assessment_ages_npy = load_shared_covars.assessment_ages,
        categorical_covariate_names = categorical_covariate_names,
        categorical_covariate_scs = categorical_covariate_scs
      }
    }
    if (is_binary) {
      call gwas_tasks.load_binary_phenotype { input:
        script_dir = script_dir,
        sc = sc_phenotype,
        qced_sample_list = all_qced_sample_lists.data,
        sc_year_of_birth = sc_year_of_birth,
        sc_month_of_birth = sc_month_of_birth,
        sc_date_of_death = sc_date_of_death,
        date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,
        is_zero_one_neg_nan = is_zero_one_neg_nan
      }
    }
    # regardless of continuous or binary, get the outputs and move on
    File pheno_data = select_first([load_continuous_phenotype.data, load_binary_phenotype.data])
    File covar_names = select_first([load_continuous_phenotype.covar_names, load_binary_phenotype.covar_names])
    File pheno_readme = select_first([load_continuous_phenotype.README, load_binary_phenotype.covar_names])

    call gwas_tasks.write_sample_list_for_phenotype { input:
      script_dir = script_dir,
      pheno_data = pheno_data
    }
  
    if (!is_binary) {
      call gwas_tasks.unrelated_samples as not_binary_phenotype_unrelated_samples { input:
        script_dir = script_dir,
        PRIMUS_command = PRIMUS_command,
        kinship = kinship,
        sample_list = write_sample_list_for_phenotype.data,
      }
    }
    if (is_binary) {
      call gwas_tasks.unrelated_samples as binary_phenotype_unrelated_samples { input:
        script_dir = script_dir,
        PRIMUS_command = PRIMUS_command,
        kinship = kinship,
        sample_list = write_sample_list_for_phenotype.data,
        binary_pheno_data = pheno_data
      }
    }
    # This is the final sample list file for a given phenotype
    File phenotype_unrelated_samples = select_first([not_binary_phenotype_unrelated_samples.data, binary_phenotype_unrelated_samples.data])
    # TODO double check that if we use the unrelated file already generated from this point that we get the same results
  
    call gwas_tasks.transform_trait_values { input :
      script_dir = script_dir,
      pheno_data = pheno_data,
      unrelated_samples_for_phenotype = phenotype_unrelated_samples,
      is_binary = is_binary
    }
 }

   call gwas_tasks.association_regions as str_association_regions { input :
    chr_lens = chr_lens,
    region_len = 10000000
  } 

#  scatter (region in str_association_regions.out_tsv) {
#    call gwas_tasks.regional_my_str_gwas { input :
#      script_dir = script_dir,
#      str_vcf = str_vcfs[read_int(region[0]) - 1],
#      shared_covars = load_shared_covars.shared_covars,
#      untransformed_phenotype = pheno_data[0],
#      transformed_phenotype = transform_trait_values.data[0],
#      all_samples_list = all_samples_list,
#      is_binary = is_binary,
#      binary_type = "linear", # won't be used if not binary
#      bounds = {
#        "chrom": read_int(region[0]),
#        "start": read_int(region[1]),
#        "end": read_int(region[2]),
#      },
#      phenotype_name = phenotype_name,
#    }
#  }
#
#  call gwas_tasks.concatenate_tsvs as my_str_gwas { input :
#    tsvs = regional_my_str_gwas.data
#  }

#  call gwas_tasks.str_spot_test as spot_test_1 { input:
#    script_dir = script_dir,
#    str_vcf = str_vcfs[0],
#    shared_covars = load_shared_covars.shared_covars,
#    untransformed_phenotype = pheno_data,
#    transformed_phenotype = transform_trait_values.data,
#    all_samples_list = all_samples_list,
#    is_binary = is_binary,
#    chrom = 1,
#    pos = 204527033,
#    phenotype_name = "platelet_count"
#  }
#
#  call gwas_tasks.str_spot_test as spot_test_2 { input:
#    script_dir = script_dir,
#    str_vcf = str_vcfs[0],
#    shared_covars = load_shared_covars.shared_covars,
#    untransformed_phenotype = pheno_data,
#    transformed_phenotype = transform_trait_values.data,
#    all_samples_list = all_samples_list,
#    is_binary = is_binary,
#    chrom = 1,
#    pos = 205255038,
#    phenotype_name = "platelet_count"
#  }

  call gwas_tasks.prep_plink_input { input :
    script_dir = script_dir,
    shared_covars = load_shared_covars.shared_covars,
    shared_covar_names = load_shared_covars.covar_names,
    transformed_phenotype = transform_trait_values.data[0],
    pheno_covar_names = covar_names,
    is_binary = is_binary,
    binary_type = "linear", # only used if is_binary
    phenotype_name = phenotype_name
  }

#  scatter (chrom in chroms) {
#    call gwas_tasks.chromosomal_plink_snp_association { input :
#      script_dir = script_dir,
#      plink_command = plink_command,
#      imputed_snp_p_file = imputed_snp_p_files[chrom-1],
#      pheno_data = prep_plink_input.data,
#      chrom = chrom,
#      phenotype_name = phenotype_name,
#      binary_type = if !is_binary then "linear" else "linear_binary",
#    }
#  }
#
#  call gwas_tasks.concatenate_tsvs as plink_snp_association { input :
#    tsvs = chromosomal_plink_snp_association.data
#  }

  # TODO second round for binary
  # TODO do ethnic STR subset differently for binary?

  call gwas_tasks.generate_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas.tsv,
    phenotype = phenotype_name,
    spacing = "250000",
    thresh = "5e-8"
  }

  call gwas_tasks.generate_peaks as summary_manhattan_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas.tsv,
    phenotype = phenotype_name,
    spacing = "20000000",
    thresh = "5e-8"
  }

  # TODO summary manhattan

  call finemapping_tasks.generate_regions { input : 
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas.tsv
  }

  # for future analyses

  call finemapping_tasks.get_strs_in_finemapping_regions { input :
    script_dir = script_dir,
    str_loci = str_loci,
    finemapping_regions_for_pheno = generate_regions.data
  }

  scatter (ethnicity_enumeration in range(length(ethnicities))) {
    scatter (pair in zip(chroms, get_strs_in_finemapping_regions.str_loci)) {
      call gwas_tasks.regional_my_str_gwas as ethnic_regional_my_str_gwas { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[pair.left-1],
        vars_file = pair.right,
        shared_covars = load_shared_covars.shared_covars,
        untransformed_phenotype = pheno_data[ethnicity_enumeration+1],
        transformed_phenotype = transform_trait_values[ethnicity_enumeration+1].data,
        all_samples_list = all_samples_list,
        is_binary = is_binary,
        binary_type = "linear", # won't be used if not binary
        phenotype_name = phenotype_name,
      }
    }
    call gwas_tasks.concatenate_tsvs as ethnic_my_str_gwas { input :
      tsvs = ethnic_regional_my_str_gwas.data
    }
  }

  # continue fine-mapping

  Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_regions.data)

  # finemap each region
  scatter (region_idx in range(length(finemapping_regions_tsv) - 1)) {
    region bounds = {
      "chrom": read_int(region[0]),
      "start": read_int(region[1]),
      "end": read_int(region[2]),
    }

    # first call FINEMAP

    call finemapping_tasks.finemap_write_input_variants { innput :
      str_assoc_results = my_str_gwas.tsv,
      snp_assoc_results = plink_snp_association.tsv,
      variants_to_filter = snp_vars_to_filter_from_finemapping[bounds[0] - 1],
      phenotype_samples_list = phenotype_unrelated_samples[0],
      phenotype = phenotype_name,
      bounds = bounds
    }

    call finemapping_tasks.finemap_load_gts { input :
      script_dir = script_dir,
      strs = str_vcfs[bounds[0] - 1],
      bgen = imputed_snp_bgens[bounds[0] - 1],
      all_samples = all_samples_list,
      phenotype_samples = phenotype_unrelated_samples[0],
      zfile = finemap_write_input_variants.zfile,
      phenotype_name = phenotype_name,
      bounds = bounds
    }

    call finemap_calc_corrs_retryable.finemap_calc_corrs_retryable { input :
      script_dir = script_dir,
      gts_h5 = finemap_load_gts.gts_h5
    }

    call finemap_write_corrs_retryable.finemap_write_corrs_retryable { input :
      script_dir = script_dir,
      lds_h5 = finemap_calc_corrs_retryable.lds_h5
    }

    call finemapping_tasks.finemap_run { input :
      script_dir = script_dir,
      master = finemap_write_input_variants.master,
      zfile = finemap_write_input_variants.zfile,
      all_variants_ld = finemap_write_corrs_retryable.all_variants_ld
    }
  }

#  output {
#    Array[File] out_sample_lists = all_qced_sample_lists.data
#    File assessment_ages = load_shared_covars.assessment_ages
#    File shared_covars = load_shared_covars.shared_covars
#    File shared_covar_names = load_shared_covars.covar_names
#    File pheno_data_out = pheno_data
#    File covar_names_out = covar_names
#    File pheno_readme_out = pheno_readme
#    File unrelated_samples_for_phenotype_out = phenotype_unrelated_samples
#    File transform_trait_values_out = transform_trait_values.data
#    File sp1 = spot_test_1.data
#    File sp2 = spot_test_2.data
#    #File my_str_gwas_out = my_str_gwas.tsv
#    #File plink_snp_gwas_out = plink_snp_association.tsv
#  }
}