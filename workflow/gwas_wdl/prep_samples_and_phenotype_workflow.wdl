version 1.0

import "gwas_tasks.wdl"

workflow prep_samples_and_phenotype {

  input {
    String script_dir
    String PRIMUS_command

    String phenotype_name
    Array[String] categorical_covariate_names
    Array[File] categorical_covariate_scs
    Boolean is_binary
    Boolean is_zero_one_neg_nan
    String date_of_most_recent_first_occurrence_update

    File fam_file # task for generating this?
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
  Array[String] all_ethnicities_ = flatten([["white_brits"], ethnicities])
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
      low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data,
      subpop_sample_list = subpop_sample_list
    }
  }

  scatter (pair in zip(all_qced_sample_lists.data, all_ethnicities_)) {
    call gwas_tasks.unrelated_samples as ethnicity_unrelated_samples { input:
      script_dir = script_dir,
      PRIMUS_command = PRIMUS_command,
      kinship = kinship,
      sample_list = pair.left,
      prefix = pair.right
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
        categorical_covariate_scs = categorical_covariate_scs,
        prefix = "~{all_ethnicities_[sample_list_idx]}_original_"
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
        is_zero_one_neg_nan = is_zero_one_neg_nan,
        prefix = "~{all_ethnicities_[sample_list_idx]}_original_"
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
        prefix = all_ethnicities_[sample_list_idx]
      }
    }
    if (is_binary) {
      call gwas_tasks.unrelated_samples as binary_phenotype_unrelated_samples { input:
        script_dir = script_dir,
        PRIMUS_command = PRIMUS_command,
        kinship = kinship,
        sample_list = write_all_samples_for_phenotype.data,
        binary_pheno_data = pheno_data_,
        prefix = all_ethnicities_[sample_list_idx]
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
      is_binary = is_binary,
      prefix = "~{all_ethnicities_[ethnicity_idx]}_pheno"
    }
  }

  output {
    # we can return intermediate sample lists if desired

    # arrays are one per the six ethnicities, starting with white brits
    Array[String] all_ethnicities = all_ethnicities_

    # sample lists 
    # unrelated, qced and takes into account the subpop if specified
    # but not subset to those with the specified phenotype
    # also not guaranteed to be a superset of those with the specified phenotype (due to separate PRIMUS runs)
    # can be rewritten to do that if necessary
    Array[File] sample_lists = ethnicity_unrelated_samples.data 

    File? subpop_sample_list_input = subpop_sample_list

    File shared_covars = shared_covars_
    File shared_covar_names = load_shared_covars.covar_names

    # sample lists subset to those with the specified phenotype
    Array[File] samples_for_phenotype = samples_for_phenotype_ # unrelated, qced and takes into account subpop if specified

    Array[File] pheno_data = pheno_data_ # raw
    Array[File] transformed_trait_values = transform_trait_values.data
    Array[File] pheno_covar_names = pheno_covar_names_
    Array[File] pheno_readme = pheno_readme_
  }
}
