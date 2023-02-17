version 1.0

import "tasks.wdl"
import "expanse_workflow.wdl"

# figures that use preexisting files
# and so are not part of the main GWAS -> finemapping pipeline
workflow expanse_figures {

  input {
    String script_dir = "."
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample"
    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam"

    # files that could be computed via WDL but were pregenerated beforehand
    File white_brits_sample_list = "sample_qc/runs/white_brits/no_phenotype/combined_unrelated.sample", #ethnicity_unrelated_samples.data[0],
    File black_sample_list = "sample_qc/runs/black/no_phenotype/combined_unrelated.sample",
    File south_asian_sample_list = "sample_qc/runs/south_asian/no_phenotype/combined_unrelated.sample",
    File chinese_sample_list = "sample_qc/runs/chinese/no_phenotype/combined_unrelated.sample",
    File unrelated_samples_CBL_hom_not_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/platelet_count/combined.sample"
    File unrelated_samples_CBL_hom_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_begin_C_T_snp/white_brits/platelet_count/combined.sample"
  }

  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
    PFiles imputed_snp_p_files = {
      "pgen": "array_imputed/pfile_converted/chr{chrom+1}.pgen",
      "pvar": "array_imputed/pfile_converted/chr{chrom+1}.pvar",
      "psam": "array_imputed/pfile_converted/chr{chrom+1}.psam",
    }
  }

  call expanse_workflow.extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
  }

  call expanse_workflow.extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }

  call tasks.load_shared_covars { input:
    script_dir = script_dir,
    fam_file = fam_file,
    sc_pcs = pcs.data,
    sc_assessment_ages = assessment_ages.data
  }

  # TODO get location for files
  call tasks.fig_4a { input :
    script_dir = script_dir,
    all_samples_list = all_samples_list,
    white_brits_sample_list = white_brits_sample_list,
    black_sample_list = black_sample_list,
    south_asian_sample_list = south_asian_sample_list,
    chinese_sample_list = chinese_sample_list,
    str_vcf_chr_11 = str_vcfs[10],
    specific_alleles = specific_alleles,
  }

  call expanse_workflow.extract_field as platet_count_sc { input :
    script_dir = script_dir,
    id = 30080,
  }
  
  call expanse_workflow.extract_field as platet_count_covariate_sc { input :
    script_dir = script_dir,
    id = 30083,
  }

  # all, not qced or subset to unrelated, the sample list for this won't be used, only the data
  call tasks.load_continuous_phenotype as platelet_count_all_samples { input:
    script_dir = script_dir,
    sc = platelet_count_sc.data,
    qced_sample_list = all_samples_list, 
    assessment_ages_npy = load_shared_covars.assessment_ages,
    categorical_covariate_names = ["platelet_count_device_id"],
    categorical_covariate_scs = [platelet_count_covariate_sc.data],
  }

  # use precomputed sample lists
  call tasks.transform_trait_values as platelet_count_transformed_CBL_SNP_not_hom_samples { input:
    script_dir = script_dir,
    pheno-data = platelet_count_all_samples.data,
    unrelated_samples_for_phenotype = unrelated_samples_CBL_hom_not_begin_C_T_snp,
    is_binary = False
  }

  call tasks.str_spot_test as CBL_hom_not_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_samples.data 
    transformed_phenotype = platelet_count_transformed_CBL_SNP_not_hom_samples.data, 
    all_samples_list = all_samples_list
    is_binary = false
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call tasks.transform_trait_values as platelet_count_transformed_CBL_SNP_hom_samples { input:
    script_dir = script_dir,
    pheno-data = platelet_count_all_samples.data,
    unrelated_samples_for_phenotype = unrelated_samples_CBL_hom_begin_C_T_snp,
    is_binary = False
  }
   
  call tasks.str_spot_test as CBL_hom_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_samples.data, 
    transformed_phenotype = platelet_count_transformed_CBL_SNP_hom_samples.data,
    all_samples_list = all_samples_list
    is_binary = false
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call tasks.locus_plot as fig_4f { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [CBL_hom_not_SNP_assoc.data, CBL_hom_SNP_assoc.data],
    group_names = ["homozygous (CGG)n", "homozygous CGGTGG(CGG)m"],
    dosage_threshold = 200,
  }

#  call tasks.plot_locus as plot_SLC2A2_locus {
#    phenotype_name = "total_bilirubin"
#    unit = ""
#  }

  output {
    File fig_4a_svg_out = fig_4a.svg
    File fig_4a_png_out = fig_4a.png
  }  
}