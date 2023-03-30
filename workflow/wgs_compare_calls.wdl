version 1.0

task compare_calls {
  input {
    File script = "dx://UKB_Test:/imputed_strs_paper/wgs/compare_calls.sh"
    File hipstr_strs = "dx://UKB_Test:/TargetedSTR/results/100_finemapped_strs/100_finemapped_strs.merged.sorted.vcf.gz"
    File hipstr_strs_idx = "dx://UKB_Test:/TargetedSTR/results/100_finemapped_strs/100_finemapped_strs.merged.sorted.vcf.gz.tbi"
    File imputed_strs = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/finemapped_strs_20230313.vcf.gz"
    File imputed_strs_idx = "dx://UKB_Test:/imputed_strs_paper/imputed_strs/finemapped_strs_20230313.vcf.gz.tbi"
  }

  output {
    File omitted_samples_hipstr = "out-vcf1-omitted-samples.tab"
    File omitted_samples_imputed = "out-vcf2-omitted-samples.tab"
    File locuscompare = "out-locuscompare.tab"
    File samplecompare = "out-samplecompare.tab"
    File overallcompare = "out-overall.tab"
  }

  command <<<
    export MERGED_FILE=~{hipstr_strs}
    export FINEMAPPED_STRS=~{imputed_strs}
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "47h30m"
    dx_instance_type: "mem3_ssd1_v2_x16"
  }
}

workflow compare_calls_w {
  call compare_calls 

  output {
    File omitted_samples_hipstr = compare_calls.omitted_samples_hipstr
    File omitted_samples_imputed = compare_calls.omitted_samples_imputed
    File locuscompare = compare_calls.locuscompare
    File samplecompare = compare_calls.samplecompare
    File overallcompare = compare_calls.overallcompare
  }
}
