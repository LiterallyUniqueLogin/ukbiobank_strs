version 1.0

task check_merged_alleles {
  input {
    File script = "dx://UKB_Test:/imputed_strs_paper/wgs/check_merged_alleles.py"
    File merged_file = "dx://UKB_Test:/TargetedSTR/results/100_finemapped_strs/100_finemapped_strs.merged.sorted.vcf.gz"
    Array[File] batches 
  }

  parameter_meta {
    merged_file : "stream"
    batches : "stream"
  }

  command <<<
    envsetup ~{script} ~{merged_file} ~{sep=" " batches}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "47h30m"
  }

}

workflow check_merged_alleles_w {
  scatter (i in range(40)) {
    String batches = "dx://UKB_Test:/TargetedSTR/results/100_finemapped_strs/100_finemapped_strs_~{i}/targetTR_stage-7/100_finemapped_strs-CHUNK~{i}.filtered.sorted.vcf.gz"
  }

  call check_merged_alleles { input :
    batches = batches
  }
}
