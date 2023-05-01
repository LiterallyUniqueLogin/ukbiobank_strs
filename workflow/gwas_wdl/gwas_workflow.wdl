version 1.0

import "gwas_tasks.wdl"

workflow gwas {

  input {
    String script_dir
    String PRIMUS_command
    String plink_command = "plink2"

    File chr_lens
  
    # one per chrom
    Array[VCF]+ str_vcfs
    Array[PFiles]+ imputed_snp_p_files

    File str_loci
    File flank_start_to_start_and_end_pos
    File str_hg19_pos_bed
    File str_hg38_pos_bed
    File repeat_units_table

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

    # If specified, must contain all samples of all ethnicities that you want included
    # (so any samples not included will be omitted)
    # samples that fail QC will still be removed
    # analyses will still be split per ethnicity
    # each ethnicity's sample list will still be shrunk to remove related participants
    File? subpop_sample_list

    # Shortcuts for rerunning the STR paper analyses without redoing the randomness of subsetting
    Array[File]? cached_unrelated_samples_for_phenotype
    File? cached_shared_covars # not sure why this cached version has different samples
  }

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
  Array[String] all_ethnicities_ = flatten([['white_british'], ethnicities])
  Array[File] ethnic_sample_lists = ethnic_sample_lists_task.sample_lists
  Array[File] all_sample_lists_ = flatten([
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

  scatter (sample_list in all_sample_lists_) {
    call gwas_tasks.qced_sample_list as all_qced_sample_lists { input:
      script_dir = script_dir,
      unqced_sample_list = sample_list,
      withdrawn_sample_list = withdrawn_sample_list,
      sex_aneuploidy_sample_list = sex_aneuploidy_sample_list.data,
      sex_mismatch_sample_list = sex_mismatch_sample_list.data,
      low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data,
      subpop_sample_list = subpop_sample_list
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

  File shared_covars_ = select_first([cached_shared_covars, load_shared_covars.shared_covars])

  # get qced unrelated phenotype sample list and transformed phenotype data for each ethnicity
  scatter (sample_list_idx in range(length(all_qced_sample_lists.data))) {
    if (!is_binary) {
      call gwas_tasks.load_continuous_phenotype { input :
        script_dir = script_dir,
        sc = sc_phenotype,
        qced_sample_list = all_qced_sample_lists.data[sample_list_idx],
        assessment_ages_npy = load_shared_covars.assessment_ages,
        categorical_covariate_names = categorical_covariate_names,
        categorical_covariate_scs = categorical_covariate_scs
      }
    }
    if (is_binary) {
      call gwas_tasks.load_binary_phenotype { input:
        script_dir = script_dir,
        sc = sc_phenotype,
        qced_sample_list = all_qced_sample_lists.data[sample_list_idx],
        sc_year_of_birth = sc_year_of_birth,
        sc_month_of_birth = sc_month_of_birth,
        sc_date_of_death = sc_date_of_death,
        date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,
        is_zero_one_neg_nan = is_zero_one_neg_nan
      }
    }
    # regardless of continuous or binary, get the outputs and move on
    File pheno_data_ = select_first([load_continuous_phenotype.data, load_binary_phenotype.data])
    File pheno_covar_names_ = select_first([load_continuous_phenotype.covar_names, load_binary_phenotype.covar_names])
    File pheno_readme_ = select_first([load_continuous_phenotype.README, load_binary_phenotype.covar_names])

    call gwas_tasks.write_sample_list_for_phenotype as write_all_samples_for_phenotype { input:
      script_dir = script_dir,
      pheno_data = pheno_data_
    }
  
    if (!is_binary) {
      call gwas_tasks.unrelated_samples as not_binary_phenotype_unrelated_samples { input:
        script_dir = script_dir,
        PRIMUS_command = PRIMUS_command,
        kinship = kinship,
        sample_list = write_all_samples_for_phenotype.data,
      }
    }
    if (is_binary) {
      call gwas_tasks.unrelated_samples as binary_phenotype_unrelated_samples { input:
        script_dir = script_dir,
        PRIMUS_command = PRIMUS_command,
        kinship = kinship,
        sample_list = write_all_samples_for_phenotype.data,
        binary_pheno_data = pheno_data_
      }
    }
    # This is the final sample list file for a given phenotype

    if (defined(cached_unrelated_samples_for_phenotype)) {
      File samples_for_phenotype_cached = select_first([cached_unrelated_samples_for_phenotype])[sample_list_idx]
    }
    if (!defined(cached_unrelated_samples_for_phenotype)) {
      File samples_for_phenotype_uncached = select_first([
        not_binary_phenotype_unrelated_samples.data,
        binary_phenotype_unrelated_samples.data
      ])
    }
    File samples_for_phenotype_ = select_first([
      samples_for_phenotype_cached,
      samples_for_phenotype_uncached
    ])
  }

  scatter (ethnicity_idx in range(length(pheno_data_))) {
    call gwas_tasks.transform_trait_values { input :
      script_dir = script_dir,
      pheno_data = pheno_data_[ethnicity_idx],
      samples_for_phenotype = samples_for_phenotype_[ethnicity_idx],
      is_binary = is_binary
    }
  }

   call gwas_tasks.association_regions as str_association_regions { input :
    chr_lens = chr_lens,
    region_len = 10000000
  } 

  scatter (association_region in str_association_regions.out_tsv) {

    Int chrom = association_region[0]
    Int start = association_region[1]
    Int end = association_region[2]

    Int chrom_minus_one = chrom - 1

    region bounds = {
        "chrom": chrom,
        "start": start,
        "end": end
    }

    call gwas_tasks.regional_my_str_gwas { input :
      script_dir = script_dir,
      str_vcf = str_vcfs[chrom_minus_one],
      shared_covars = shared_covars_,
      untransformed_phenotype = pheno_data_[0],
      transformed_phenotype = transform_trait_values.data[0],
      all_samples_list = all_samples_list,
      is_binary = is_binary,
      binary_type = "linear", # won't be used if not binary
      bounds = bounds,
      phenotype_name = phenotype_name,
    }

    call gwas_tasks.reformat_my_str_gwas_table_for_publication { input :
      script_dir = script_dir,
      phenotype = phenotype_name,
      my_str_gwas = regional_my_str_gwas.data,
      flank_start_to_start_and_end_pos = flank_start_to_start_and_end_pos,
      str_hg19_pos_bed = str_hg19_pos_bed,
      str_hg38_pos_bed = str_hg38_pos_bed,
      repeat_units_table = repeat_units_table,
    }
  }

  call gwas_tasks.concatenate_tsvs as my_str_gwas_ { input :
    tsvs = regional_my_str_gwas.data
  }

#  call gwas_tasks.concatenate_tsvs as publishable_my_str_gwas_ { input :
#    tsvs = reformat_my_str_gwas_table_for_publication.out
#  }

  call gwas_tasks.prep_plink_input { input :
    script_dir = script_dir,
    shared_covars = shared_covars_,
    shared_covar_names = load_shared_covars.covar_names,
    transformed_phenotype = transform_trait_values.data[0],
    pheno_covar_names = pheno_covar_names_[0],
    is_binary = is_binary,
    binary_type = "linear", # only used if is_binary
    phenotype_name = phenotype_name
  }

  scatter (chrom in range(22)) {
    call gwas_tasks.chromosomal_plink_snp_association { input :
      script_dir = script_dir,
      plink_command = plink_command,
      imputed_snp_p_file = imputed_snp_p_files[chrom],
      pheno_data = prep_plink_input.data,
      chrom = chrom+1,
      phenotype_name = phenotype_name,
      binary_type = if !is_binary then "linear" else "linear_binary",
    }
  }

  call gwas_tasks.concatenate_tsvs as plink_snp_association { input :
    tsvs = chromosomal_plink_snp_association.data
  }

  # TODO second round for binary
  # TODO do ethnic STR subset differently for binary?

  # TODO interactive manhattan

  call gwas_tasks.generate_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas_.tsv,
    phenotype = phenotype_name,
    spacing = "250000",
    thresh = "5e-8"
  }

  call gwas_tasks.generate_peaks as overview_manhattan_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas_.tsv,
    phenotype = phenotype_name,
    spacing = "20000000",
    thresh = "5e-8"
  }

  # TODO overview manhattan

  call gwas_tasks.generate_finemapping_regions { input :
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas_.tsv
  }
  
  call gwas_tasks.get_strs_in_finemapping_regions { input :
    script_dir = script_dir,
    str_loci = str_loci,
    finemapping_regions_for_pheno = generate_finemapping_regions.data
  }
  
  scatter (ethnicity_enumeration in range(length(pheno_data_) -1)) {
    scatter (pair in zip(range(22), get_strs_in_finemapping_regions.str_loci)) {
      Int ethnicity_idx = ethnicity_enumeration + 1
      call gwas_tasks.regional_my_str_gwas as ethnic_regional_my_str_gwas { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[pair.left],
        vars_file = pair.right,
        shared_covars = shared_covars_,
        untransformed_phenotype = pheno_data_[ethnicity_idx],
        transformed_phenotype = transform_trait_values.data[ethnicity_idx],
        all_samples_list = all_samples_list,
        is_binary = is_binary,
        binary_type = "linear", # won't be used if not binary
        phenotype_name = phenotype_name,
      }
    }
    call gwas_tasks.concatenate_tsvs as ethnic_my_str_gwas_ { input :
      tsvs = ethnic_regional_my_str_gwas.data
    }

    #TODO logistic if binary?
  }

  output {
    # arrays are one per the six ethnicities, starting with white brits
    Array[String] all_ethnicities = all_ethnicities_
    Array[File] all_sample_lists = all_sample_lists_ # unqced, and this doesn't take into account the subpop
    Array[File] sample_lists = ethnicity_unrelated_samples.data # unrelated, qced and takes into account the subpop if specified

    File? subpop_sample_list_input = subpop_sample_list

    File shared_covars = shared_covars_
    File shared_covar_names = load_shared_covars.covar_names

    Array[File] all_samples_for_phenotype = write_all_samples_for_phenotype.data # unqced
    Array[File] samples_for_phenotype = samples_for_phenotype_ # unrelated and qced

    Array[File] pheno_data = pheno_data_ # raw
    Array[File] transformed_trait_values = transform_trait_values.data
    Array[File] pheno_covar_names = pheno_covar_names_
    Array[File] pheno_readme = pheno_readme_

    File my_str_gwas = my_str_gwas_.tsv
    #File publishable_my_str_gwas = publishable_my_str_gwas_.tsv
    File plink_snp_gwas = plink_snp_association.tsv
    File peaks = generate_peaks.peaks
    File peaks_readme = generate_peaks.readme
    File finemapping_regions = generate_finemapping_regions.data
    File finemapping_regions_readme = generate_finemapping_regions.readme
	  Array[File] ethnic_my_str_gwas = ethnic_my_str_gwas_.tsv
  }
}
