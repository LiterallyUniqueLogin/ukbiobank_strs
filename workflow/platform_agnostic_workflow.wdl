# platform agnostic workflow

version 1.0

import "tasks.wdl"

task association_regions {
  input {
    File chr_lens
    Int region_len
  }

  output {
    # TODO set these
    Array[Int] chroms
    Array[Int] start_poses
    Array[Int] end_poses
  } 

  command <<<
    python -c "
      import numpy as np
      chr_lens = np.genfromtxt(
        ~{chr_lens},
        usecols=[1],
        skip_header=1,
        dtype=int
      )

      chroms = []
      starts = []
      ends = []
      for chrom in range(1, 23):
        chr_len = chr_lens[chrom-1]
        for start in range(1, chr_len, ~{region_len}):
          if start + region_len - 1 > chr_len:
            end = chr_len
          else:
            end = start + region_len - 1
          chroms.append(chroms)
          starts.append(start)
          ends.append(end)
       TODO
    "
  >>>

  runtime {
    shortTask: true
    dx_timeout: "5m"
  }
}

workflow main {

  input {
    String script_dir
    String PRIMUS_executable

    File chr_lens

    Array[File]+ str_vcfs
    File specific_alleles # could replace this with a different way of providing this info

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
    PRIMUS_executable = PRIMUS_executable
    pheno_data = pheno_data,
    kinship = kinship,
    is_binary = is_binary,
  }

  call tasks.transform_trait_values { input :
    script_dir = script_dir,
    pheno_data = pheno_data,
    unrelated_samples_for_phenotype = unrelated_samples_for_phenotype.data,
    is_binary = is_binary
  }

  call tasks.fig_4a { input :
    script_dir = script_dir,
    str_vcf_chr_11 = str_vcfs[11],
    specific_alleles = specific_alleles
  }

  call association_regions as str_association_regions { input :
    chr_lens = chr_lens,
    region_len = 10000000
  } 

  scatter (region in zip(str_association_regions.chroms, zip(str_association_regions.starts, str_association_regions.ends))) {
    call tasks.regional_my_str_gwas { input :
      script_dir = script_dir,
      str_vcf = str_vcfs[region.left],
      shared_covars = load_shared_covars.shared_covars,
      untransformed_phenotype = pheno_data,
      transformed_phenotype = transform_trait_values.data
      is_binary = is_binary,
      binary_type = "linear", # won't be used if not binary
      chrom = region.left,
      start_pos = region.right.left,
      end_pos = region.right.right
      phenotype_name = phenotype_name
    }
  }

  call tasks.concatenate_csvs as my_str_gwas { input :
    tsvs = region_my_str_gwas.data
  }

  # TODO second round for binary

  output {
    Array[File] out_sample_lists = all_qced_sample_lists.data
    File assessment_ages = load_shared_covars.assessment_ages
    File shared_covars = load_shared_covars.shared_covars
    File shared_covar_names = load_shared_covars.covar_names
    File pheno_data_out = pheno_data
    File covar_names_out = covar_names
    File pheno_readme_out = pheno_readme
    File unrelated_samples_for_phenotype_out = unrelated_samples_for_phenotype.data
    File transform_trait_values_out = transform_trait_values.data
    File fig_4a_svg_out = fig_4a.svg
    File fig_4a_png_out = fig_4a.png
  }
}
