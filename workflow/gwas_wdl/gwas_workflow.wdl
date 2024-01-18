version 1.0

import "gwas_tasks.wdl"
import "prep_samples_and_phenotype_workflow.wdl"

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
    Array[String]? categorical_covariate_names
    Array[File]? categorical_covariate_scs
    Boolean is_binary
    Boolean is_zero_one_neg_nan
    String date_of_most_recent_first_occurrence_update

    File fam_file # task for generating this?
    File withdrawn_sample_list
    File kinship
    File all_samples_list

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

    # should specify sc or the next three
    File? sc_phenotype
    File? premade_pheno_npy # first col ID, second pheno, remaining are covars
    File? premade_pheno_covar_names
    File? premade_pheno_readme

    Boolean transform = true

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

  call prep_samples_and_phenotype_workflow.prep_samples_and_phenotype { input :
    script_dir = script_dir,
    PRIMUS_command = PRIMUS_command,

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariate_scs,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,

    fam_file = fam_file,
    withdrawn_sample_list = withdrawn_sample_list,
    kinship = kinship,

    sc_white_brits = sc_white_brits,
    sc_ethnicity_self_report = sc_ethnicity_self_report,
    sc_sex_aneuploidy = sc_sex_aneuploidy,
    sc_genetic_sex = sc_genetic_sex,
    sc_reported_sex = sc_reported_sex,
    sc_kinship_count = sc_kinship_count,
    sc_assessment_ages = sc_assessment_ages,
    sc_pcs = sc_pcs,
    sc_year_of_birth = sc_year_of_birth,
    sc_month_of_birth = sc_month_of_birth,
    sc_date_of_death = sc_date_of_death,
    sc_phenotype = sc_phenotype,
    premade_pheno_npy = premade_pheno_npy,
    premade_pheno_covar_names = premade_pheno_covar_names,
    premade_pheno_readme = premade_pheno_readme,
    transform = transform,

    subpop_sample_list = subpop_sample_list,

    cached_unrelated_samples_for_phenotype = cached_unrelated_samples_for_phenotype,
    cached_shared_covars = cached_shared_covars
  }

  Array[File] phenos_to_associate = select_first([
    prep_samples_and_phenotype.transformed_trait_values,
    prep_samples_and_phenotype.pheno_data,
  ])


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
      shared_covars = prep_samples_and_phenotype.shared_covars,
      untransformed_phenotype = prep_samples_and_phenotype.pheno_data[0],
      transformed_phenotype = phenos_to_associate[0],
      all_samples_list = all_samples_list,
      is_binary = is_binary,
      binary_type = "linear", # won't be used if not binary
      bounds = bounds,
      phenotype_name = phenotype_name,
    }

#    call gwas_tasks.reformat_my_str_gwas_table_for_publication { input :
#      script_dir = script_dir,
#      phenotype = phenotype_name,
#      my_str_gwas = regional_my_str_gwas.data,
#      flank_start_to_start_and_end_pos = flank_start_to_start_and_end_pos,
#      str_hg19_pos_bed = str_hg19_pos_bed,
#      str_hg38_pos_bed = str_hg38_pos_bed,
#      repeat_units_table = repeat_units_table,
#    }
  }

  call gwas_tasks.concatenate_tsvs as my_str_gwas_ { input :
    tsvs = regional_my_str_gwas.data,
    out = "white_brits_str_gwas"
  }

#  call gwas_tasks.concatenate_tsvs as publishable_my_str_gwas_ { input :
#    tsvs = reformat_my_str_gwas_table_for_publication.out
#  }

  call gwas_tasks.prep_plink_input { input :
    script_dir = script_dir,
    shared_covars = prep_samples_and_phenotype.shared_covars,
    shared_covar_names = prep_samples_and_phenotype.shared_covar_names,
    transformed_phenotype = phenos_to_associate[0],
    pheno_covar_names = prep_samples_and_phenotype.pheno_covar_names[0],
    is_binary = is_binary,
    binary_type = "linear", # only used if is_binary
    phenotype_name = phenotype_name
  }

  scatter (chrom in range(22)) {
    call gwas_tasks.plink_snp_association as chromosomal_plink_snp_association { input :
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
    tsvs = chromosomal_plink_snp_association.data,
    out = "white_brits_snp_gwas"
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

  call gwas_tasks.overview_manhattan as overview_manhattan_ { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    chr_lens = chr_lens,
    str_gwas_results = my_str_gwas_.tsv,
    snp_gwas_results = plink_snp_association.tsv,
    peaks = overview_manhattan_peaks.peaks,
    ext = "png"
  }

  call gwas_tasks.generate_finemapping_regions { input :
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas_.tsv,
    prefix = "~{phenotype_name}_"
  }
  
  call gwas_tasks.get_strs_in_finemapping_regions { input :
    script_dir = script_dir,
    str_loci = str_loci,
    finemapping_regions_for_pheno = generate_finemapping_regions.data
  }
  
  scatter (ethnicity_enumeration in range(5)) {
    Int ethnicity_idx = ethnicity_enumeration + 1
    scatter (pair in zip(range(22), get_strs_in_finemapping_regions.strs_in_finemapping_regions)) {
      call gwas_tasks.regional_my_str_gwas as ethnic_regional_my_str_gwas { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[pair.left],
        vars_file = pair.right,
        shared_covars = prep_samples_and_phenotype.shared_covars,
        untransformed_phenotype = prep_samples_and_phenotype.pheno_data[ethnicity_idx],
        transformed_phenotype = phenos_to_associate[ethnicity_idx],
        all_samples_list = all_samples_list,
        is_binary = is_binary,
        binary_type = "linear", # won't be used if not binary
        phenotype_name = phenotype_name,
      }
    }
    call gwas_tasks.concatenate_tsvs as ethnic_my_str_gwas_ { input :
      tsvs = ethnic_regional_my_str_gwas.data,
      out = "~{prep_samples_and_phenotype.all_ethnicities[ethnicity_idx]}_str_gwas"
    }

    #TODO logistic if binary?
  }

  output {
    # sample lists not subset to those with the specified phenotype
    Array[File] sample_lists = prep_samples_and_phenotype.sample_lists # unrelated, qced and takes into account the subpop if specified

    File? subpop_sample_list_input = subpop_sample_list

    File shared_covars = prep_samples_and_phenotype.shared_covars
    File shared_covar_names = prep_samples_and_phenotype.shared_covar_names

    # sample lists subset to those with the specified phenotype
    Array[File] samples_for_phenotype = prep_samples_and_phenotype.samples_for_phenotype # unrelated, qced and takes into account subpop if specified

    Array[File] pheno_data = prep_samples_and_phenotype.pheno_data # raw
    Array[File]? transformed_trait_values = prep_samples_and_phenotype.transformed_trait_values
    Array[File] pheno_covar_names = prep_samples_and_phenotype.pheno_covar_names
    Array[File] pheno_readme = prep_samples_and_phenotype.pheno_readme

    File my_str_gwas = my_str_gwas_.tsv
    #File publishable_my_str_gwas = publishable_my_str_gwas_.tsv
    File plink_snp_gwas = plink_snp_association.tsv
    File peaks = generate_peaks.peaks
    File peaks_readme = generate_peaks.readme
    File overview_manhattan = overview_manhattan_.plot
    File finemapping_regions = generate_finemapping_regions.data
    File finemapping_regions_readme = generate_finemapping_regions.readme
	  Array[File] ethnic_my_str_gwas = ethnic_my_str_gwas_.tsv
  }
}
