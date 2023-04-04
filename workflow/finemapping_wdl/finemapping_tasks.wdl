version 1.0

# for structs
import "../gwas_wdl/gwas_tasks.wdl"

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

struct FINEMAP_output {
  File snp_file
  File log_sss
  File config
  Array[File] creds
}

struct SuSiE_output {
  File lbf
  File lbf_variable
  File sigma2
  File V
  File converged
  File lfsr
  File requested_coverage
  File alpha
  Array[File] CSs
}

####################### Extracting association signals #########################

task finemap_write_input_variants {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_write_input_variants.py"
    File python_array_utils = "~{script_dir}/finemapping/python_array_utils.py"
    File sample_utils = "~{script_dir}/finemapping/sample_utils.py"    

    File str_assoc_results
    File snp_assoc_results
    File variants_to_filter
    File phenotype_samples_list
    String phenotype
    region bounds

    Float? snp_str_ratio
    Float? total_prob
    Int? mac
    Float? inclusion_threshold
  }

  output {
    File master = "finemap_input.master"
    File zfile = "finemap_input.z"
    File readme = "README.txt"
  }

  command <<<
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      . \
      . \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{variants_to_filter} \
      ~{phenotype_samples_list} \
      ~{phenotype} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
      ~{if defined(snp_str_ratio) then "--snp-str-ratio ~{snp_str_ratio}" else ""}
      ~{if defined(total_prob) then "--total-prob ~{total_prob}" else ""}
      ~{if defined(mac) then "--mac ~{mac}" else ""}
      ~{if defined(inclusion_threshold) then "--inclusion-threshold ~{inclusion_threshold}" else ""}
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
      ~{all_samples} \
      ~{phenotype_samples} \
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
    dx_timeout: time
    mem: "4g"
  }
}

task finemap_write_corrs {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_write_corrs.py"

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
    dx_timeout: time
    cpu: 4
    mem: "8g"
  }
}

# TODO check that prior snps flag passing makes sense
task finemap_run {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_run.py"
    String finemap_command

    File master
    File zfile
    File all_variants_ld

    Boolean prior_snps = false
    Float? prior_std
    Float? prob_conv_sss_tol
  }

  output {
    FINEMAP_output finemap_output = object {
      snp_file: "finemap_output.snp",
      log_sss: "finemap_output.log_sss",
      config: "finemap_output.config",
      creds: glob("finemap_output.cred*")
    }
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
      ~{if prior_snps then "--prior-snps" else if defined(prior_std) then "--prior-std ~{prior_std}" else if defined(prob_conv_sss_tol) then "--prob-conv-sss-tol ~{prob_conv_sss_tol}" else ""}
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
    File script = "~{script_dir}/finemapping/susie_choose_vars.py"

    File str_assoc_results
    File snp_assoc_results
    File variants_to_filter

    String phenotype_name
    region bounds
    # TODO we could have used a MAF cutoff for susie vars but never did
    # how come? Should this be added as another followup to supp note 3?
  }

  output {
    File readme = "readme.txt"
    File colnames = "colnames.txt"
  }

  command <<<
    envsetup ~{script} \
      readme.txt \
      colnames.txt \
      ~{str_assoc_results} \
      ~{snp_assoc_results} \
      ~{variants_to_filter} \
      ~{phenotype_name} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
  }
}

task susie_load_gts {
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
    File pheno_data
    File shared_covars
    
    File colnames # from prev finemapping step
    String phenotype_name
    region bounds
    Boolean best_guess = false

    String time
  }

  output {
    File? gts_h5 = "gts.h5"
    File? pheno_residuals_h5 = "pheno_residual.h5"
    File? readme = "readme.txt"
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
      ~{colnames} \
      ~{all_samples} \
      ~{phenotype_samples} \
      ~{phenotype_name} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
      --pheno-fname ~{pheno_data} \
      --shared-covars-fname ~{shared_covars} \
      --pheno-out pheno_residual.h5 \
      ~{if best_guess then "--best-guess" else ""}
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: time
    mem: "4g"
  }
}

task susie_run {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/susie_run.r"

    File gts_h5
    File pheno_residuals_h5
    Int L
    Int max_iter

    Float? tol
    Float? snp_p_over_str_p
    File? varnames_file
    Float? res_var
    Float? prior_var

    String mem
    String time
  }

  output {
    SuSiE_output? susie_output = object {
      lbf: "lbf.tab",
      lbf_variable: "lbf_variable.tab",
      sigma2: "sigma2.txt",
      V: "V.tab",
      converged: "converged.txt",
      lfsr: "lfsr.tab",
      requested_coverage: "requested_coverage.txt",
      alpha: "alpha.tab",
      CSs: glob("cs*.txt")
    }
  }

  command <<<
    ACTIVATE=ukb_r envsetup Rscript ~{script} \
      . \
      ~{pheno_residuals_h5} \
      ~{gts_h5} \
      ~{L} \
      ~{max_iter} \
      ~{"--tol " + tol} \
      ~{"--snp-p-over-str-p " + snp_p_over_str_p + " --varnames-file " + varnames_file} \
      ~{"--residual-variance " + res_var} \
      ~{"--scaled-prior-variance " + prior_var}
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: time
    mem: mem
  }
}

task first_pass_finemapping_df {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"

    String phenotype_name
    File snp_assoc_results
    File str_assoc_results
    Array[File] ethnic_str_assoc_results

    Array[FINEMAP_output] original_finemap_outputs
    Array[SuSiE_output] original_susie_outputs
    Array[String] regions
    Array[Int] chroms
  }

  output {
    File all_regions_concordance = "finemapping_all_regions_concordance_~{phenotype_name}.tab"
    File susie_all_regions_min_abs_corrs = "susie_all_regions_min_abs_corrs_~{phenotype_name}.npy"
  }

  command <<<
    REGIONS=~{write_lines(regions)}
    { echo "region" ; cat $REGIONS ; } >> "$REGIONS".headered
    CHROMS=~{write_lines(chroms)}
    { echo "chrom" ; cat $CHROMS ; } >> "$CHROMS".headered
    envsetup ~{script} \
      . \
      True \
      ~{phenotype_name} \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{sep=" " ethnic_str_assoc_results} \
      first_pass
      <(paste ~{write_objects(original_finemap_outputs)} "$CHROMS".headered "$REGIONS".headered)
      <(paste ~{write_objects(original_susie_outputs)} "$CHROMS".headered "$REGIONS".headered)
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    mem: "50g"
  }
}

task generate_followup_regions_tsv {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"

    File first_pass_df
  }

  output {
    File tsv = "followup_regions.tsv"
  }

  command <<<
    envsetup ~{script} . ~{first_pass_df}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
  } 
}

# TODO missing stuff in between

task followup_finemapping_conditions_comparison {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"

    Array[File] all_conditions_tsvs
  }

  output {
    # intermediate tsvs
    File doubly_finemapped_STRs = "doubly_finemapped_STRs.tab"
    File confidently_finemapped_STRs = "confidently_finemapped_STRs.tab"
    File overconfidently_finemapped_STRs = "overconfidently_finemapped_STRs.tab"

    # used for deciding confidently finemapped
    File susie_best_guess_png = "susie_consistency_alpha_best_guess.png"
    File susie_best_guess_svg = "susie_consistency_alpha_best_guess.svg"
    File finemap_conv_tol_png = "finemap_consistency_alpha_conv_tol.png"
    File finemap_conv_tol_svg = "finemap_consistency_alpha_conv_tol.svg"
    File finemap_total_prob_png = "finemap_consistency_alpha_total_prob.png"
    File finemap_total_prob_svg = "finemap_consistency_alpha_total_prob.svg"
    File finemap_prior_std_derived_png = "finemap_consistency_alpha_prior_std_derived.png"
    File finemap_prior_std_derived_svg = "finemap_consistency_alpha_prior_std_derived.svg"
    File finemap_mac_png = "finemap_consistency_alpha_mac.png"
    File finemap_mac_svg = "finemap_consistency_alpha_mac.svg"
    File finemap_p_thresh_png = "finemap_consistency_alpha_p_thresh.png"
    File finemap_p_thresh_svg = "finemap_consistency_alpha_p_thresh.svg"

    # too conservative
    File finemap_ratio_png = "finemap_consistency_alpha_ratio.png"
    File finemap_ratio_svg = "finemap_consistency_alpha_ratio.svg"
    File susie_ratio_png = "susie_consistency_alpha_ratio.png"
    File susie_ratio_svg = "susie_consistency_alpha_ratio.svg"
    File finemap_prior_std_low_png = "finemap_consistency_alpha_prior_std_low.png"
    File finemap_prior_std_low_svg = "finemap_consistency_alpha_prior_std_low.svg"
  }

  command <<<
    envsetup ~{script} . . ~{sep=" " all_conditions_tsvs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    mem: "50g"
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
