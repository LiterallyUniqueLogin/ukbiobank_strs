version 1.0

import "../gwas_wdl/gwas_tasks.wdl"

task compare_calls {
  input {
    File script = "dx://UKB_Test:/imputed_strs_paper/wgs/compare_calls.sh"
    File hipstr_strs = "dx://UKB_Test:/imputed_strs_paper/wgs/20230907_old_and_new_WGS_calls.vcf.gz"
    File hipstr_strs_idx = "dx://UKB_Test:/imputed_strs_paper/wgs/20230907_old_and_new_WGS_calls.vcf.gz.tbi"
    File imputed_strs = "dx://UKB_Test:/imputed_strs_paper/wgs/20230907_old_and_new_imputed_singly_finemapped_strs.vcf.gz"
    File imputed_strs_idx = "dx://UKB_Test:/imputed_strs_paper/wgs/20230907_old_and_new_imputed_singly_finemapped_strs.vcf.gz.tbi"
#    File hipstr_strs = "dx://UKB_Test:/TargetedSTR/results/20230831_new_STRs/20230831_new_STRs.merged.sorted.vcf.gz"
#    File hipstr_strs_idx = "dx://UKB_Test:/TargetedSTR/results/20230831_new_STRs/20230831_new_STRs.merged.sorted.vcf.gz.tbi"
#    File imputed_strs = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/20230905_new_singly_finemapped_strs.vcf.gz"
#    File imputed_strs_idx = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/20230905_new_singly_finemapped_strs.vcf.gz.tbi"
#    File hipstr_strs = "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz"
#    File hipstr_strs_idx = "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz.tbi"
#    File imputed_strs = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/finemapped_strs_20230313.vcf.gz"
#    File imputed_strs_idx = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/finemapped_strs_20230313.vcf.gz.tbi"
    File? samples_file
    String prefix
  }

  output {
    File? omitted_samples_hipstr = "~{prefix}-vcf1-omitted-samples.tab"
    File omitted_samples_imputed = "~{prefix}-vcf2-omitted-samples.tab"
    File locuscompare = "~{prefix}-locuscompare.tab"
    File samplecompare = "~{prefix}-samplecompare.tab"
    File overallcompare = "~{prefix}-overall.tab"
  }

  command <<<
    export WGS_STRS=~{hipstr_strs}
    export IMPUTED_STRS=~{imputed_strs}
    export SAMPLES_FILE=~{samples_file}
    export PREFIX=~{prefix}
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "47h30m"
    dx_instance_type: "mem3_ssd1_v2_x16"
  }
}

task double_samples {
  input {
    File in_samples_file
  }

  output {
    File out_samples_file = "doubled_samples.txt"
  }

  command <<<
    awk '{ print $1 "_" $1 }' ~{in_samples_file} > doubled_samples.txt
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
  }
}


workflow compare_calls_w {

  input {
    String script_dir  = "dx://UKB_Test:/imputed_strs_paper"

    File withdrawn_sample_list = "dx://UKB_Test:/imputed_strs_paper/sample_qc/withdrawn.sample"

		String sc_data_dir = "dx://UKB_Test:/imputed_strs_paper/main_dataset/extracted_data"
		File sc_white_brits = "~{sc_data_dir}/22006.tsv"
		File sc_ethnicity_self_report = "~{sc_data_dir}/21000.tsv"
		File sc_sex_aneuploidy = "~{sc_data_dir}/22019.tsv"
		File sc_genetic_sex = "~{sc_data_dir}/22001.tsv"
		File sc_reported_sex = "~{sc_data_dir}/31.tsv"
		File sc_kinship_count = "~{sc_data_dir}/22021.tsv"
  }

	# get QCed but not unrelated samples per ethnicity
  # copie from prep_samples_and_phenotype workflow
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
  Array[String] all_ethnicities_ = flatten([['white_british'], ethnicities])
  Array[File] ethnic_sample_lists = ethnic_sample_lists_task.sample_lists
  Array[File] all_sample_lists_ = flatten([
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

  scatter (sample_list in all_sample_lists_) {
    call gwas_tasks.qced_sample_list as all_qced_sample_lists { input:
      script_dir = script_dir,
      unqced_sample_list = sample_list,
      withdrawn_sample_list = withdrawn_sample_list,
      sex_aneuploidy_sample_list = sex_aneuploidy_sample_list.data,
      sex_mismatch_sample_list = sex_mismatch_sample_list.data,
      low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data,
    }

    call double_samples { input : in_samples_file = all_qced_sample_lists.data } 
  }

	# run the comparison per ethnicity
  scatter (samples_file_idx in range(length(double_samples.out_samples_file))) {
    call compare_calls { input :
      samples_file = double_samples.out_samples_file[samples_file_idx],
      prefix = "results_~{all_ethnicities_[samples_file_idx]}"
    }
  }

  output {
    Array[File?] omitted_samples_hipstr = compare_calls.omitted_samples_hipstr
    Array[File] omitted_samples_imputed = compare_calls.omitted_samples_imputed
    Array[File] locuscompare = compare_calls.locuscompare
    Array[File] samplecompare = compare_calls.samplecompare
    Array[File] overallcompare = compare_calls.overallcompare
  }
}
