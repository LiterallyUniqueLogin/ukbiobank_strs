version 1.0

import "gwas_tasks.wdl"
import "prep_samples_and_phenotype_workflow.wdl"

workflow gwas {

  input {
    # TODO timing might not be right if n_pcs != 8
    # or if firth

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
    Boolean firth = false
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
    Boolean other_ethnicities = true

    # If specified, must contain all samples of all ethnicities that you want included
    # (so any samples not included will be omitted)
    # samples that fail QC will still be removed
    # analyses will still be split per ethnicity
    # each ethnicity's sample list will still be shrunk to remove related participants
    File? subpop_sample_list

    # Shortcuts for rerunning the STR paper analyses without redoing the randomness of subsetting
    Array[File]? cached_unrelated_samples_for_phenotype
    File? cached_shared_covars # not sure why this cached version has different samples

    Int n_pcs = 40
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
    other_ethnicities = other_ethnicities,

    subpop_sample_list = subpop_sample_list,

    cached_unrelated_samples_for_phenotype = cached_unrelated_samples_for_phenotype,
    cached_shared_covars = cached_shared_covars,
    n_pcs = n_pcs
  }

  call gwas_tasks.write_sample_list_plink_style { input :
    sample_list = prep_samples_and_phenotype.samples_for_phenotype[0]
  }

  scatter (freq_chrom_minus_one in range(22)) {
    Int freq_chrom = freq_chrom_minus_one + 1
    call gwas_tasks.imputed_snp_frequencies as chrom_imputed_snp_frequencies { input :
      imputed_snp_pfiles = imputed_snp_p_files[freq_chrom_minus_one],
      plink_style_sample_file = write_sample_list_plink_style.out,
      out = '~{phenotype_name}_~{freq_chrom}'
    }
  }

  call gwas_tasks.concatenate_tsvs as imputed_snp_frequencies { input :
    tsvs = chrom_imputed_snp_frequencies.afreq,
    out = '~{phenotype_name}_afreq'
  }

  Array[File] phenos_to_associate = select_first([
    prep_samples_and_phenotype.transformed_trait_values,
    prep_samples_and_phenotype.pheno_data,
  ])

  call gwas_tasks.chr_pos_association_regions as continuous_str_association_regions { input :
    chr_pos_file = str_loci,
    region_len = 5000 # one job per max 5k variants
  } 

  scatter (continuous_str_association_region in continuous_str_association_regions.out_tsv) {
    Int continuous_str_chrom = continuous_str_association_region[0]
    Int continuous_str_start = continuous_str_association_region[1]
    Int continuous_str_end = continuous_str_association_region[2]

    Int continuous_str_chrom_minus_one = continuous_str_chrom - 1

    region continuous_str_bounds = {
        "chrom": continuous_str_chrom,
        "start": continuous_str_start,
        "end": continuous_str_end
    }

    call gwas_tasks.regional_my_str_gwas as continuous_regional_str_gwas { input :
      script_dir = script_dir,
      str_vcf = str_vcfs[continuous_str_chrom_minus_one],
      shared_covars = prep_samples_and_phenotype.shared_covars,
      untransformed_phenotype = prep_samples_and_phenotype.pheno_data[0],
      transformed_phenotype = phenos_to_associate[0],
      all_samples_list = all_samples_list,
      is_binary = false,
      bounds = continuous_str_bounds,
      phenotype_name = phenotype_name,
    }
  }

  call gwas_tasks.concatenate_tsvs as continuous_str_gwas { input :
    tsvs = continuous_regional_str_gwas.data,
    out = "white_brits_str_gwas"
  }

  call gwas_tasks.prep_plink_input as continuous_plink_input { input :
    script_dir = script_dir,
    phenotype_name = if transform && !is_binary then "rin_~{phenotype_name}" else phenotype_name,
    pheno_data = phenos_to_associate[0],
    pheno_covar_names = prep_samples_and_phenotype.pheno_covar_names[0],
    shared_covars = prep_samples_and_phenotype.shared_covars,
    shared_covar_names = prep_samples_and_phenotype.shared_covar_names,
    is_binary = false,
  }

  scatter (continuous_snp_chrom_minus_one in range(22)) {
    Int continuous_snp_chrom = continuous_snp_chrom_minus_one + 1

    call gwas_tasks.plink_snp_association as continuous_plink_snp_association { input :
      script_dir = script_dir,
      plink_command = plink_command,
      imputed_snp_p_file = imputed_snp_p_files[continuous_snp_chrom_minus_one],
      pheno_data = continuous_plink_input.data,
      chrom = continuous_snp_chrom,
      phenotype_name = if !is_binary && transform then "rin_~{phenotype_name}" else phenotype_name,
      is_binary = false,
      firth = false,
      time = "48h"
    }
  }

  call gwas_tasks.concatenate_tsvs as continuous_snp_gwas { input :
    tsvs = continuous_plink_snp_association.data,
    out = "white_brits_snp_gwas"
  }

  if (is_binary) {
    call gwas_tasks.subset_assoc_results_to_chr_pos as subsetted_continuous_str_chr_pos { input :
      assoc_results = continuous_str_gwas.tsv,
      p_val_col = 'p_~{phenotype_name}',
      p_val_thresh = 0.1,
      chrom_col = 'chrom',
      pos_col = 'pos'
    }

    call gwas_tasks.chr_pos_association_regions as binary_str_association_regions { input :
      chr_pos_file = subsetted_continuous_str_chr_pos.tab,
      region_len = 2500 # one job per max 2.5k variants
    }

    scatter (binary_str_association_region in binary_str_association_regions.out_tsv) {
      Int binary_str_chrom = binary_str_association_region[0]
      Int binary_str_start = binary_str_association_region[1]
      Int binary_str_end = binary_str_association_region[2]

      Int binary_str_chrom_minus_one = binary_str_chrom - 1

      region binary_str_bounds = {
          "chrom": binary_str_chrom,
          "start": binary_str_start,
          "end": binary_str_end
      }

      call gwas_tasks.regional_my_str_gwas as binary_regional_str_gwas { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[binary_str_chrom_minus_one],
        shared_covars = prep_samples_and_phenotype.shared_covars,
        untransformed_phenotype = prep_samples_and_phenotype.pheno_data[0],
        transformed_phenotype = phenos_to_associate[0],
        all_samples_list = all_samples_list,
        is_binary = true,
        bounds = binary_str_bounds,
        vars_file = subsetted_continuous_str_chr_pos.tab,
        phenotype_name = phenotype_name,
      }
    }

    call gwas_tasks.concatenate_tsvs as binary_str_gwas { input :
      tsvs = binary_regional_str_gwas.data,
      out = "white_brits_str_gwas"
    }

    call gwas_tasks.prep_plink_input as binary_plink_input { input :
      script_dir = script_dir,
      phenotype_name = phenotype_name,
      pheno_data = phenos_to_associate[0],
      pheno_covar_names = prep_samples_and_phenotype.pheno_covar_names[0],
      shared_covars = prep_samples_and_phenotype.shared_covars,
      shared_covar_names = prep_samples_and_phenotype.shared_covar_names,
      is_binary = true,
    }

    call gwas_tasks.subset_assoc_results_to_chr_pos as subsetted_continuous_snp_chr_pos { input :
      assoc_results = continuous_snp_gwas.tsv,
      p_val_col = 'P',
      p_val_thresh = 0.1,
      chrom_col = '#CHROM',
      pos_col = 'POS'
    }

    scatter (binary_snp_chrom_minus_one in range(22)) {
      Int binary_snp_chrom = binary_snp_chrom_minus_one + 1

      call gwas_tasks.plink_snp_association as binary_plink_snp_association { input :
        script_dir = script_dir,
        plink_command = plink_command,
        imputed_snp_p_file = imputed_snp_p_files[binary_snp_chrom_minus_one],
        pheno_data = binary_plink_input.data,
        chrom = binary_snp_chrom,
        phenotype_name = phenotype_name,
        is_binary = true,
        firth = firth,
        time = if !firth then "10h" else "48h",
        vars = subsetted_continuous_snp_chr_pos.tab
      }
    }

    call gwas_tasks.concatenate_tsvs as binary_snp_gwas { input :
      tsvs = binary_plink_snp_association.data,
      out = "white_brits_snp_gwas"
    }
  }

  File str_gwas_ = select_first([binary_str_gwas.tsv, continuous_str_gwas.tsv])
  File plink_snp_gwas_ = select_first([binary_snp_gwas.tsv, continuous_snp_gwas.tsv])

  call gwas_tasks.qq_plot as snp_qq_plot_ { input :
    script_dir = script_dir,
    results_tab = plink_snp_gwas_,
    p_val_col = 'P',
    phenotype_name = phenotype_name,
    variant_type = 'SNP',
    out_name = 'snp_qq_plot',
    null_values = 'NA',
    max_p_val = if !is_binary then 1.5 else 0.1
  }

  call gwas_tasks.qq_plot as str_qq_plot_ { input :
    script_dir = script_dir,
    results_tab = str_gwas_,
    p_val_col = 'p_~{phenotype_name}',
    phenotype_name = phenotype_name,
    variant_type = 'STR',
    out_name = 'str_qq_plot',
    max_p_val = if !is_binary then 1.5 else 0.1
  }

  # TODO interactive manhattan

  call gwas_tasks.generate_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_gwas_,
    str_assoc_results = str_gwas_,
    phenotype = phenotype_name,
    spacing = "250000",
    thresh = "5e-8"
  }

  call gwas_tasks.generate_peaks as overview_manhattan_peaks { input :
    script_dir = script_dir,
    snp_assoc_results = plink_snp_gwas_,
    str_assoc_results = str_gwas_,
    phenotype = phenotype_name,
    spacing = "20000000",
    thresh = "5e-8"
  }

  call gwas_tasks.overview_manhattan as overview_manhattan_ { input :
    script_dir = script_dir,
    phenotype_name = phenotype_name,
    chr_lens = chr_lens,
    str_gwas_results = str_gwas_,
    snp_gwas_results = plink_snp_gwas_,
    peaks = overview_manhattan_peaks.peaks,
    ext = "png"
  }

  call gwas_tasks.generate_finemapping_regions { input :
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_gwas_,
    str_assoc_results = str_gwas_,
    prefix = "~{phenotype_name}_"
  }

  if (other_ethnicities) {
    # TODO do ethnic STR subset differently for binary?
    # TODO do logistic regression if binary

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
          phenotype_name = phenotype_name,
        }
      }
      call gwas_tasks.concatenate_tsvs as ethnic_my_str_gwas_ { input :
        tsvs = ethnic_regional_my_str_gwas.data,
        out = "~{prep_samples_and_phenotype.all_ethnicities[ethnicity_idx]}_str_gwas"
      }
    }
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

    File imputed_snp_freqs = imputed_snp_frequencies.tsv

    File my_str_gwas = str_gwas_
    File plink_snp_gwas = plink_snp_gwas_

    File snp_qq_plot = snp_qq_plot_.plot
    File str_qq_plot = str_qq_plot_.plot
    File peaks = generate_peaks.peaks
    File peaks_readme = generate_peaks.readme
    File overview_manhattan = overview_manhattan_.plot
    File finemapping_regions = generate_finemapping_regions.data
    File finemapping_regions_readme = generate_finemapping_regions.readme
	  Array[File]? ethnic_my_str_gwas = ethnic_my_str_gwas_.tsv
  }
}
