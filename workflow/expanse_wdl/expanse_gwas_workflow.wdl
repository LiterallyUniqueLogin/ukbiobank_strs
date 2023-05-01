version 1.0

import "expanse_tasks.wdl"
import "../gwas_wdl/gwas_workflow.wdl"

workflow expanse_gwas {

  input {
    String script_dir  = "."

    File chr_lens = "misc_data/genome/chr_lens.txt"
    File str_loci = "snpstr/str_loci.txt"

    # TODO could generate this set of files
    File flank_start_to_start_and_end_pos = "snpstr/flank_trimmed_vcf/vars.tab"
    File repeat_units_table = "snpstr/repeat_units.tab"
    File str_loci = "snpstr/str_loci.txt"
    File str_hg19_pos_bed = "snpstr/str_loci.bed"
    File str_hg38_pos_bed = "snpstr/str_loci.hg38.bed"

#    String phenotype_name = "platelet_count"
#    Int phenotype_id = 30080
#    Array[String] categorical_covariate_names = ["platelet_count_device_id"]
#    Array[Int] categorical_covariate_ids = [30083]
    String phenotype_name
    Int phenotype_id
    Array[String] categorical_covariate_names
    Array[Int] categorical_covariate_ids
    Boolean is_binary = false
    Boolean is_zero_one_neg_nan = false # different binary encoding
    String date_of_most_recent_first_occurrence_update = "2021-04-01" # only needed for binary phenotypes

    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam" # Could instead create a task for downloading this with ukbgene
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample" # could instead create a task for downloading this with ukbgene
    File withdrawn_sample_list = "sample_qc/common_filters/remove/withdrawn.sample"
    File kinship = "misc_data/ukbgene/ukb46122_rel_s488282.dat" # could create a task for downloading this with ukbgene

    # If specified, must contain all samples of all ethnicities that you want included
    # (so any samples not included will be omitted)
    # samples that fail QC will still be removed
    # analyses will still be split per ethnicity
    # each ethnicity's sample list will still be shrunk to remove related participants
    File? subpop_sample_list

    Array[File]? cached_unrelated_samples_for_phenotype =
    File cached_shared_covars
#    Array[File]? cached_unrelated_samples_for_phenotype = [
#      "sample_qc/runs/white_brits/platelet_count/combined_unrelated.sample",
#      "sample_qc/runs/black/platelet_count/combined_unrelated.sample",
#      "sample_qc/runs/south_asian/platelet_count/combined_unrelated.sample",
#      "sample_qc/runs/chinese/platelet_count/combined_unrelated.sample",
#      "sample_qc/runs/irish/platelet_count/combined_unrelated.sample",
#      "sample_qc/runs/white_other/platelet_count/combined_unrelated.sample",
#    ]
#    File cached_shared_covars = "traits/shared_covars/shared_covars.npy"
  }

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

  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
    PFiles imputed_snp_p_files = {
      "pgen": "array_imputed/pfile_converted/chr~{chrom+1}.pgen",
      "pvar": "array_imputed/pfile_converted/chr~{chrom+1}.pvar",
      "psam": "array_imputed/pfile_converted/chr~{chrom+1}.psam",
    }
    bgen imputed_snp_bgens = {
      "bgen": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen",
      "index": "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.bgi",
      "bgen_reader_metadata" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata",
      "bgen_reader_metadata2" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.metadata2.mmm",
      "bgen_reader_complex_metadata2" : "array_imputed/ukb_imp_chr~{chrom+1}_v3.bgen.complex.metadata2.mmm"
    }
  }

  call gwas_workflow.gwas { input:
    script_dir = script_dir,
    PRIMUS_command = "run_PRIMUS.pl",

    chr_lens = chr_lens,

    str_vcfs = str_vcfs,
    imputed_snp_p_files = imputed_snp_p_files,

    str_loci = str_loci,
    flank_start_to_start_and_end_pos = flank_start_to_start_and_end_pos,
    str_hg19_pos_bed = str_hg19_pos_bed,
    str_hg38_pos_bed = str_hg38_pos_bed,
    repeat_units_table = repeat_units_table,

    phenotype_name = phenotype_name,
    categorical_covariate_names = categorical_covariate_names,
    categorical_covariate_scs = categorical_covariates.data,
    is_binary = is_binary,
    is_zero_one_neg_nan = is_zero_one_neg_nan,
    date_of_most_recent_first_occurrence_update = date_of_most_recent_first_occurrence_update,

    fam_file = fam_file,
    all_samples_list = all_samples_list,
    withdrawn_sample_list = withdrawn_sample_list,
    kinship = kinship, 

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
}
