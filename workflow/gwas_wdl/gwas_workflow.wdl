# platform agnostic workflow

version 1.0

import "gwas_tasks.wdl"
import "gwas_given_pheno_data_workflow.wdl"

# TODO fix chr_lens 21 and 20 the same

workflow gwas {

  input {
    String script_dir
    String PRIMUS_command
    String plink_command = "plink2"

    File chr_lens
  
    # one per chrom
    Array[VCF]+ str_vcfs
    Array[PFiles]+ imputed_snp_p_files

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
    File pheno_data_ = select_first([load_continuous_phenotype.data, load_binary_phenotype.data])
    File pheno_covar_names_ = select_first([load_continuous_phenotype.covar_names, load_binary_phenotype.covar_names])
    File pheno_readme_ = select_first([load_continuous_phenotype.README, load_binary_phenotype.covar_names])

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
    File samples_for_phenotype_= select_first([not_binary_phenotype_unrelated_samples.data, binary_phenotype_unrelated_samples.data])
  }

  call gwas_given_pheno_data_workflow.gwas_given_pheno_data { input :
    script_dir = script_dir,
    plink_command = plink_command,

    chr_lens = chr_lens,

    str_vcfs = str_vcfs,
    imputed_snp_p_files = imputed_snp_p_files,

    phenotype_name = phenotype_name,
    is_binary = is_binary,

    all_samples_list = all_samples_list,

    shared_covars = load_shared_covars.shared_covars,
    shared_covar_names = load_shared_covars.covar_names,
    pheno_data = pheno_data_,
    pheno_covar_names = pheno_covar_names_,
    pheno_readme = pheno_readme_,
    samples_for_phenotype = samples_for_phenotype_,
  }

  output {
    # arrays are one per the six ethnicities, starting with white brits
    File shared_covars = load_shared_covars.shared_covars
    File shared_covar_names = load_shared_covars.covar_names
    Array[File] pheno_data = pheno_data_
    Array[File] pheno_covar_names = pheno_covar_names_
    Array[File] samples_for_phenotype = samples_for_phenotype_ # unrelated and qced

    Array[File] transformed_trait_values = gwas_given_pheno_data.transformed_trait_values
    Array[File] transformed_trait_values_readme = gwas_given_pheno_data.transformed_trait_values_readme
    #File my_str_gwas_out = gwas_given_pheno_data.my_str_gwas.tsv
    #File plink_snp_gwas_out = gwas_given_pheno_data.plink_snp_association.tsv
    File peaks = gwas_given_pheno_data.peaks
    File peaks_readme = gwas_given_pheno_data.readme
  }
}
