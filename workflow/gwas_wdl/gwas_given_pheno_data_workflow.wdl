 platform agnostic workflow

version 1.0

import "gwas_tasks.wdl"
import "finemapping_tasks.wdl"

# TODO fix chr_lens 21 and 20 the same

workflow gwas_given_pheno_data {

  input {
    String script_dir
    String plink_command = "plink2"

    File chr_lens

    # one per chrom
    Array[VCF]+ str_vcfs
    Array[PFiles]+ imputed_snp_p_files

    String phenotype_name
    Boolean is_binary

    File all_samples_list

    # arrays are one per the six ethnicities, starting with white brits
    File shared_covars
    File shared_covar_names
    Array[File] pheno_data
    Array[File] pheno_covar_names
    Array[File] samples_for_phenotype
  }

  # TODO double check that if we use the unrelated file already generated from this point that we get the same results
  scatter (ethnicity_idx in range(length(pheno_data))) {
    call gwas_tasks.transform_trait_values { input :
      script_dir = script_dir,
      pheno_data = pheno_data[ethnicity_idx],
      samples_for_phenotype = samples_for_phenotype[ethnicity_idx],
      is_binary = is_binary
    }
  }

#   call gwas_tasks.association_regions as str_association_regions { input :
#    chr_lens = chr_lens,
#    region_len = 10000000
#  } 
#
#  scatter (region in str_association_regions.out_tsv) {
#    call gwas_tasks.regional_my_str_gwas { input :
#      script_dir = script_dir,
#      str_vcf = str_vcfs[read_int(region[0]) - 1],
#      shared_covars = shared_covars,
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
#  call gwas_tasks.concatenate_tsvs as my_str_gwas_ { input :
#    tsvs = regional_my_str_gwas.data
#  }

  call gwas_tasks.str_spot_test as spot_test_1 { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[0],
    shared_covars = shared_covars,
    untransformed_phenotype = pheno_data[0],
    transformed_phenotype = transform_trait_values.data[0],
    all_samples_list = all_samples_list,
    is_binary = is_binary,
    chrom = 1,
    pos = 204527033,
    phenotype_name = "platelet_count"
  }

  call gwas_tasks.str_spot_test as spot_test_2 { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[0],
    shared_covars = shared_covars,
    untransformed_phenotype = pheno_data[0],
    transformed_phenotype = transform_trait_values.data[0],
    all_samples_list = all_samples_list,
    is_binary = is_binary,
    chrom = 1,
    pos = 205255038,
    phenotype_name = "platelet_count"
  }

#  call gwas_tasks.prep_plink_input { input :
#    script_dir = script_dir,
#    shared_covars = shared_covars,
#    shared_covar_names = covar_names,
#    transformed_phenotype = transform_trait_values.data[0],
#    pheno_covar_names = covar_names,
#    is_binary = is_binary,
#    binary_type = "linear", # only used if is_binary
#    phenotype_name = phenotype_name
#  }
#
#  scatter (chrom in range(22)) {
#    call gwas_tasks.chromosomal_plink_snp_association { input :
#      script_dir = script_dir,
#      plink_command = plink_command,
#      imputed_snp_p_file = imputed_snp_p_files[chrom],
#      pheno_data = prep_plink_input.data,
#      chrom = chrom+1,
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

#  call gwas_tasks.generate_peaks { input :
#    script_dir = script_dir,
#    snp_assoc_results = plink_snp_association.tsv,
#    str_assoc_results = my_str_gwas_.tsv,
#    phenotype = phenotype_name,
#    spacing = "250000",
#    thresh = "5e-8"
#  }
#
#  call gwas_tasks.generate_peaks as overview_manhattan_peaks { input :
#    script_dir = script_dir,
#    snp_assoc_results = plink_snp_association.tsv,
#    str_assoc_results = my_str_gwas_.tsv,
#    phenotype = phenotype_name,
#    spacing = "20000000",
#    thresh = "5e-8"
#  }

  # TODO overview manhattan

  call finemapping_tasks.generate_regions { input :
    script_dir = script_dir,
    chr_lens = chr_lens,
    phenotype = phenotype_name,
    snp_assoc_results = plink_snp_association.tsv,
    str_assoc_results = my_str_gwas_.tsv
  }
  
  call finemapping_tasks.get_strs_in_finemapping_regions { input :
    script_dir = script_dir,
    str_loci = str_loci,
    finemapping_regions_for_pheno = generate_regions.data
  }
  
  scatter (ethnicity_enumeration in range(length(pheno_data) -1)) {
    scatter (pair in zip(chroms, get_strs_in_finemapping_regions.str_loci)) {
      call gwas_tasks.regional_my_str_gwas as ethnic_regional_my_str_gwas { input :
        script_dir = script_dir,
        str_vcf = str_vcfs[pair.left-1],
        vars_file = pair.right,
        shared_covars = shared_covars,
        untransformed_phenotype = pheno_data[ethnicity_enumeration+1],
        transformed_phenotype = transform_trait_values[ethnicity_enumeration+1].data,
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
    Array[File] transformed_trait_values = transform_trait_values.data
    Array[File] transformed_trait_values_readme = transform_trait_values.readme
    File sp1 = spot_test_1.data
    File sp2 = spot_test_2.data
    #File my_str_gwas = my_str_gwas_.tsv
    #File plink_snp_gwas = plink_snp_association.tsv
    #File peaks = generate_peaks.peaks
    #File peaks_readme = generate_peaks.readme
	  #File ethnic_my_str_gwas = ethnic_my_str_gwas_.tsv
  }
}
