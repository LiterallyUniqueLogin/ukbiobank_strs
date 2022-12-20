# platform agnostic workflow

version 1.0

import "tasks.wdl"

workflow main {

  input {
    String script_dir

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
  }

  call tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = script_dir,
    sc = sc_white_brits
  }

  call tasks.ethnic_sample_lists as ethnic_sample_lists_task { input: 
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

  call tasks.write_sample_list as sex_aneuploidy_sample_list { input:
    script_dir = script_dir,
    sc = sc_sex_aneuploidy
  }

  call tasks.sex_mismatch_sample_list { input:
    script_dir = script_dir,
    sc_genetic_sex = sc_genetic_sex,
    sc_reported_sex = sc_reported_sex
  }
  
  call tasks.write_sample_list as low_genotyping_quality_sample_list { input:
    script_dir = script_dir,
    sc = sc_kinship_count,
    value = -1
  }

  scatter (sample_list in all_sample_lists) {
    call tasks.qced_sample_list as all_qced_sample_lists { input:
      script_dir = script_dir,
      unqced_sample_list = sample_list,
      withdrawn_sample_list = withdrawn_sample_list,
      sex_aneuploidy_sample_list = sex_aneuploidy_sample_list.data,
      sex_mismatch_sample_list = sex_mismatch_sample_list.data,
      low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data
    }
  }

  call tasks.load_shared_covars { input:
    script_dir = script_dir,
    fam_file = fam_file,
    sc_pcs = sc_pcs,
    sc_assessment_ages = sc_assessment_ages
  }

  if (!is_binary) {
    call tasks.load_continuous_phenotype { input :
      script_dir = script_dir,
      sc = sc_phenotype,
      qced_sample_list = all_qced_sample_lists.data[0],
      assessment_ages_npy = load_shared_covars.assessment_ages,
      categorical_covariate_names = categorical_covariate_names,
      categorical_covariate_scs = categorical_covariate_scs
    }
  }
  if (is_binary) {
    call tasks.load_binary_phenotype { input:
      script_dir = script_dir,
      sc = sc_phenotype,
      qced_sample_list = all_qced_sample_lists.data[0],
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

  call tasks.unrelated_samples_for_phenotype { input:
    script_dir = script_dir,
    pheno_data = pheno_data,
    kinship = kinship,
    is_binary = is_binary
  }

  output {
    Array[File] out_sample_lists = all_qced_sample_lists.data
    File assessment_ages = load_shared_covars.assessment_ages
    File shared_covars = load_shared_covars.shared_covars
    File shared_covar_names = load_shared_covars.covar_names
    File pheno_data_out = pheno_data
    File covar_names_out = covar_names
    File pheno_readme_out = pheno_readme
    File unrelated_samples_for_phenotype_out = unrelated_samples_for_phenotype.data
  }
}
