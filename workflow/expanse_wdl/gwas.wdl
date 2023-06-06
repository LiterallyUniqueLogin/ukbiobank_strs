version 1.0

import "expanse_tasks.wdl"
import "expanse_files.wdl"
import "../gwas_wdl/gwas_workflow.wdl"

workflow gwas {

  input {
    String script_dir  = "."

    Int phenotype_id
    Array[String] categorical_covariate_names
    Array[Int] categorical_covariate_ids
    String phenotype_name
    Boolean is_binary = false
    Boolean is_zero_one_neg_nan = false # different binary encoding

    # If specified, must contain all samples of all ethnicities that you want included
    # (so any samples not included will be omitted)
    # samples that fail QC will still be removed
    # analyses will still be split per ethnicity
    # each ethnicity's sample list will still be shrunk to remove related participants
    File? subpop_sample_list

    Array[File]? cached_unrelated_samples_for_phenotype
    File? cached_shared_covars

    # example inputs:
#    String phenotype_name = "platelet_count"
#    Int phenotype_id = 30080
#    Array[String] categorical_covariate_names = ["platelet_count_device_id"]
#    Array[Int] categorical_covariate_ids = [30083]
  }

  call expanse_files.files

  call expanse_tasks.extract_field as white_brits { input:
    script_dir = script_dir,
    id = 22006
  }

  call expanse_tasks.extract_field as ethnicity_self_report { input :
    script_dir = script_dir,
    id = 21000
  }

  call expanse_tasks.extract_field as sex_aneuploidy { input:
    script_dir = script_dir,
    id = 22019
  }

  call expanse_tasks.extract_field as genetic_sex { input:
    script_dir = script_dir,
    id = 22001
  }

  call expanse_tasks.extract_field as reported_sex { input:
    script_dir = script_dir,
    id = 31
  }

  call expanse_tasks.extract_field as kinship_count { input:
    script_dir = script_dir,
    id = 22021
  }

  call expanse_tasks.extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }
  
  call expanse_tasks.extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
  }

  call expanse_tasks.extract_field as year_of_birth { input :
    script_dir = script_dir,
    id = 34
  }

  call expanse_tasks.extract_field as month_of_birth { input :
    script_dir = script_dir,
    id = 52
  }

  call expanse_tasks.extract_field as date_of_death { input :
    script_dir = script_dir,
    id = 40000
  }

  call expanse_tasks.extract_field as phenotype { input :
    script_dir = script_dir,
    id = phenotype_id
  }

  scatter (categorical_covariate_id in categorical_covariate_ids) {
    call expanse_tasks.extract_field as categorical_covariates { input :
      script_dir = script_dir,
      id = categorical_covariate_id
    }
  }

  call gwas_workflow.gwas { input:
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",

    chr_lens = files.chr_lens,

    str_vcfs = files.str_vcfs,
    imputed_snp_p_files = files.imputed_snp_pfiles,

    str_loci = files.str_loci,
    flank_start_to_start_and_end_pos = files.flank_start_to_start_and_end_pos,
    str_hg19_pos_bed = files.str_hg19_pos_bed,
    str_hg38_pos_bed = files.str_hg38_pos_bed,
    repeat_units_table = files.repeat_units_table,

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariates.data,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = files.date_of_most_recent_first_occurrence_update,

    fam_file = files.fam_file,
    all_samples_list = files.all_samples_list,
    withdrawn_sample_list = files.withdrawn_sample_list,
    kinship = files.kinship, 

    sc_white_brits = white_brits.data,
    sc_ethnicity_self_report = ethnicity_self_report.data,
    sc_sex_aneuploidy = sex_aneuploidy.data,
    sc_genetic_sex = genetic_sex.data,
    sc_reported_sex = reported_sex.data,
    sc_kinship_count = kinship_count.data,
    sc_assessment_ages = assessment_ages.data,
    sc_pcs = pcs.data,
    sc_year_of_birth = year_of_birth.data,
    sc_month_of_birth = month_of_birth.data,
    sc_date_of_death = date_of_death.data,
    sc_phenotype = phenotype.data,

    subpop_sample_list = subpop_sample_list,

    cached_unrelated_samples_for_phenotype = cached_unrelated_samples_for_phenotype,
    cached_shared_covars = cached_shared_covars,
  }

  output {
    Array[File] sample_lists = gwas.sample_lists

    File shared_covars = gwas.shared_covars
    File shared_covar_names = gwas.shared_covar_names

    Array[File] all_samples_for_phenotype = gwas.all_samples_for_phenotype
    Array[File] samples_for_phenotype = gwas.samples_for_phenotype

    Array[File] pheno_data = gwas.pheno_data
    Array[File] transformed_trait_values = gwas.transformed_trait_values
    Array[File] pheno_covar_names = gwas.pheno_covar_names
    Array[File] pheno_readme = gwas.pheno_readme

    File my_str_gwas = gwas.my_str_gwas
    File plink_snp_gwas = gwas.plink_snp_gwas
    File peaks = gwas.peaks
    File peaks_readme = gwas.peaks_readme
    File finemapping_regions = gwas.finemapping_regions
    File finemapping_regions_readme = gwas.finemapping_regions_readme
    Array[File] ethnic_my_str_gwas = gwas.ethnic_my_str_gwas
  }
}
