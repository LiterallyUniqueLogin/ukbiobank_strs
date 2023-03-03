version 1.0

# for structs
import "gwas_tasks.wdl"

# any input file with a default relative to the script_dir
# needs to be supplied by the user, it won't be the product of another task
# if input files to tasks can be supplied by another tasks output, 
# there will be a comment specifying
# task input files without comments need to be supplied by the user
# see the expanse workflow for where those are on expanse
# exception: sc (data showcase) tasks are labeled by data field id
# but do need to be supplied by the user

# output files from tasks will be commented with the location
# they reside on expanse
# this isn't necessary for understanding/running the WDL, just useful notes for myself
# for transitioning from snakemake to WDL

# sample_list file format
# first line is 'ID' (case insensitive)
# every successive line is a sample ID

####################### Extracting association signals #########################

task generate_regions {
  input {
    String script_dir
    File script = "~{script_dir}/signals/regions.py"

    File chr_lens # misc_data/genome/chr_lens.txt
    String phenotype
    File str_assoc_results
    File snp_assoc_results
  }

  output {
    File data = "out.tab"
    File readme = "out_README.txt"
  }

  command <<<
    envsetup ~{script} \
      ~{phenotype} \
      ~{chr_lens} \
      ~{str_assoc_results} \
      ~{snp_assoc_results} \
      ~{out.tab}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    mem: "10g"
  }
}

task get_strs_in_finemapping_regions
  input {
    String script_dir
    File script = "~{script_dir}/signals/write_finemapped_strs.py"

    File str_loci
    File finemapping_regions_for_pheno # from generate regions
  }

  output {
    Array[File] str_loci = [
      "out_chr1.tab",
      "out_chr2.tab",
      "out_chr3.tab",
      "out_chr4.tab",
      "out_chr5.tab",
      "out_chr6.tab",
      "out_chr7.tab",
      "out_chr8.tab",
      "out_chr9.tab",
      "out_chr10.tab",
      "out_chr11.tab",
      "out_chr12.tab",
      "out_chr13.tab",
      "out_chr14.tab",
      "out_chr15.tab",
      "out_chr16.tab",
      "out_chr17.tab",
      "out_chr18.tab",
      "out_chr19.tab",
      "out_chr20.tab",
      "out_chr21.tab",
      "out_chr22.tab"
    ] # one per chr
  }

  command <<<
    envsetup ~{script} \
      out \
      ~{finemapping_regions_for_pheno} \
      ~{str_loci}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    mem: "4g"
  }
}

task finemap_write_input_variants {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_write_input_variants.py"
    File sample_utils = "~{script_dir}/finemapping/sample_utils.py"    

    File str_assoc_results
    File snp_assoc_results
    File variants_to_filter #TODO?
    File phenotype_samples_list
    String phenotype
    region bounds
  }

  output {
    File master = "finemap_input.master"
    File zfile = "finemap_input.z"
    File readme = "README.txt"
  }

  command <<<
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      . \
      . \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{variants_to_filter} \
      ~{phenotype_sample_list} \
      ~{phenotype} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
      # TODO snp str ratio
      # TODO total prob
      # TODO inclusion threshold
      # TODO mac
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    mem: "6g"
  }
}

task finemap_load_gts {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemapping_load_gts.py"
    File load_and_filter_genotypes = "~{script_dir}/finemapping/load_and_filter_genotypes.py"
    File python_array_utils = "~{script_dir}/finemapping/python_array_utils.py"
    File sample_utils = "~{script_dir}/finemapping/sample_utils.py"    

    VCF strs
    bgen snps
    File all_samples
    File phenotype_samples
    
    File zfile # from prev finemapping step
    String phenotype_name
    region bounds
  }

  output {
    File gts_h5 = "gts.h5"
    File readme = "readme.txt"
  }

  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{load_and_filter_genotypes} # necessary for dxCompiler to bother to localize this file
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      readme.txt \
      gts.h5 \
      ~{strs.vcf} \
      ~{snps.bgen} \
      ~{zfile} \
      ~{phenotype_name} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
      --varname-header
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "23h30m"
    mem: "6g"
  }
}

task finemap_calc_corrs {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_calc_corrs.py"

    File gts_h5
    String time
  }

  output {
    File? lds_h5 = "lds.h5"
  }

  command <<<
    envsetup ~{script} . ~{gts_h5}
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "~{time}"
    mem: "4g"
  }
}

task finemap_write_corrs {
  input {
    String script_dir
    File script = "~{script_dir}/finemap_write_corrs.py"

    File lds_h5
    String time
  }

  output {
    File? all_variants_ld = "all_variants.ld"
  }

  command <<<
    envsetup ~{script} . ~{lds_h5}
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "~{time}"
    # TODO
    cpus: 4
    mem: "8g"
  }
}

# TODO add FIENMAP to docker
# TODO check that prior snps flag passing makes sense
task finemap_run {
  input {
    String script_dir
    File script = "~{script_dir}/finemap_run.py"
    String finemap_command

    File master
    File zfile
    File all_variants_ld

    Boolean prior_snps = false
    Float? prior_std
    Float? prob_conv_sss_tol
  }

  output {
    #TODO are there other output files
    File snp_file = "finemap_output.snp"
    File config = #TODO
  }

  command <<<
    # need all files in the same directory for FINEMAP to work
    # also necessary for dxCompiler to bother to localize this file
    ln -s ~{master} .
    ln -s ~{zfile} .
    ln -s ~{all_variants_ld} .
    envsetup ~{script} \
      . \
      ~{finemap_command} \
      ~{if prior_snps then "--prior-snps" else if defined(prior_std) then "--prior-std ~{prior_std}" else if defined(prob_conv_sss_tol) then "--prob-conv-sss-tol ~{prob_conv_sss_tol}"
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    mem: "8gb"
  }
}

task susie_choose_vars {
  input {
    String script_dir
    File script = "~{script_dir}/susie_choose_vars.py"
  }

  output {
    File readme = "readme.txt"
    File colnames = "colnames.txt"
  }

  command <<<
    envsetup ~{script} \
      readme.txt \
      colnames.txt \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
  }
}

task susie_load_gts {
  input {
    String script_dir
    File script = "~{script_dir}/"

    String time
  }

  output {

  }

  command <<<
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "~{time}"
    mem: "4g"
  }
}

task susie_run {
  input {
    String script_dir
    File script = "~{script_dir}/"
  }

  output {

  }

  command <<<
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: ""
  }
}

task todo {
  input {
    String script_dir
    File script = "~{script_dir}/"
  }

  output {

  }

  command <<<
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: ""
  }
}
