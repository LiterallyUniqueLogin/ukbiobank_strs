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
  serializable_FINEMAP_output subset
  Array[File] creds
}

# remove creds
struct serializable_FINEMAP_output {
  File snp_file
  File log_sss
  File config
}

struct SuSiE_output {
  serializable_SuSiE_output subset 
  Array[File] CSs
}

# remove CSs, add colnames
struct serializable_SuSiE_output {
  File lbf
  File lbf_variable
  File sigma2
  File V
  File converged
  File lfsr
  File requested_coverage
  File alpha
  File colnames
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
      serializable_FINEMAP_output subset = object {
        snp_file: "finemap_output.snp",
        log_sss: "finemap_output.log_sss",
        config: "finemap_output.config",
      },
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

    File colnames

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
      serializable_SuSiE_output subset = object {
        lbf: "lbf.tab",
        lbf_variable: "lbf_variable.tab",
        sigma2: "sigma2.txt",
        V: "V.tab",
        converged: "converged.txt",
        lfsr: "lfsr.tab",
        requested_coverage: "requested_coverage.txt",
        alpha: "alpha.tab",
        colnames: colnames
      },
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

####################### Summarize finemapping #########################

task first_pass_finemapping_df {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    String phenotype_name
    File snp_assoc_results
    File str_assoc_results
    Array[File] ethnic_str_assoc_results

    Array[serializable_FINEMAP_output] original_finemap_outputs
    Array[Array[File]] original_finemap_creds
    Array[serializable_SuSiE_output] original_susie_outputs
    Array[Array[File]] original_susie_CSs
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
    original_finemap_creds=~{write_tsv(original_finemap_creds)}
    { echo "creds" ; cat $original_finemap_creds | sed -e 's/\t/,/g' ; } >> "$original_finemap_creds".headered
    original_susie_CSs=~{write_tsv(original_susie_CSs)}
    { echo "CSes" ; cat $original_susie_CSs | sed -e 's/\t/,/g' ; } >> "$original_susie_CSs".headered
    paste ~{write_objects(original_finemap_outputs)} "$CHROMS".headered "$REGIONS".headered "$original_finemap_creds".headered > finemap_files.tsv
    paste ~{write_objects(original_susie_outputs)} "$CHROMS".headered "$REGIONS".headered "$original_susie_CSs".headered > susie_files.tsv
    envsetup ~{script} \
      . \
      df \
      True \
      ~{phenotype_name} \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{sep=" " ethnic_str_assoc_results} \
      first_pass \
      finemap_files.tsv \
      susie_files.tsv
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    mem: "50g"
  }
}

task first_pass_comparison {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    Array[File] first_pass_dfs
    Array[File] susie_all_min_abs_corrs
  }

  output {
    File stat_statements = "stat_statements.txt"
    File cs_min_abs_corrs_png = "cs_min_abs_corrs.png"
    File cs_min_abs_corrs_svg = "cs_min_abs_corrs.svg"
    File susie_alpha_v_pip_png = "susie_alpha_v_pip.png"
    File susie_alpha_v_pip_svg = "susie_alpha_v_pip.svg"
    File susie_alpha_histogram_png = "susie_alpha_histogram.png"
    File susie_alpha_histogram_svg = "susie_alpha_histogram.svg"
    File finemap_pip_histogram_png = "finemap_pip_histogram.png"
    File finemap_pip_histogram_svg = "finemap_pip_histogram.svg"
    File susie_cs_finemap_total_pips_png = "susie_cs_finemap_total_pips.png"
    File susie_cs_finemap_total_pips_svg = "susie_cs_finemap_total_pips.svg"
    File finemap_v_susie_consistency_STR_png = "finemap_v_susie_consistency_STR.png"
    File finemap_v_susie_consistency_STR_svg = "finemap_v_susie_consistency_STR.svg"
    File finemap_v_susie_consistency_SNP_png = "finemap_v_susie_consistency_SNP.png"
    File finemap_v_susie_consistency_SNP_svg = "finemap_v_susie_consistency_SNP.svg"
  }

  command <<<
    envsetup ~{script} \
      . \
      first_pass_comparison \
      --first-pass-dfs ~{sep=" " first_pass_dfs} \
      --susie-all-min-abs-corrs ~{sep=" " susie_all_min_abs_corrs} \
      > stat_statements.txt
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    mem: "50g"
  }
}

task generate_followup_regions_tsv {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    File first_pass_df
  }

  output {
    File tsv = "followup_regions.tsv"
  }

  command <<<
    envsetup ~{script} . followup_regions ~{first_pass_df}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
  } 
}

task followup_finemapping_conditions_df {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    String phenotype_name
    File snp_assoc_results
    File str_assoc_results
    Array[File] ethnic_str_assoc_results

    Array[serializable_FINEMAP_output] original_finemap_outputs
    Array[Array[File]] original_finemap_creds
    Array[serializable_SuSiE_output] original_susie_outputs
    Array[Array[File]] original_susie_CSs

    Array[serializable_FINEMAP_output] total_prob_finemap_outputs
    Array[Array[File]] total_prob_finemap_creds
    Array[serializable_FINEMAP_output] derived_prior_std_finemap_outputs
    Array[Array[File]] derived_prior_std_finemap_creds
    Array[serializable_FINEMAP_output] conv_tol_finemap_outputs
    Array[Array[File]] conv_tol_finemap_creds
    Array[serializable_FINEMAP_output] mac_finemap_outputs
    Array[Array[File]] mac_finemap_creds
    Array[serializable_FINEMAP_output] threshold_finemap_outputs
    Array[Array[File]] threshold_finemap_creds
    Array[serializable_SuSiE_output] best_guess_susie_outputs
    Array[Array[File]] best_guess_susie_CSs

    Array[serializable_FINEMAP_output] low_prior_std_finemap_outputs
    Array[Array[File]] low_prior_std_finemap_creds
    Array[serializable_FINEMAP_output] ratio_finemap_outputs
    Array[Array[File]] ratio_finemap_creds
    Array[serializable_SuSiE_output] ratio_susie_outputs
    Array[Array[File]] ratio_susie_CSs

    Array[String] original_regions
    Array[Int] original_chroms
    Array[String] followup_regions
    Array[Int] followup_chroms
  }

  output {
    File df = "finemapping_followup_concordance_~{phenotype_name}.tab"
  }

  command <<<
    ORIGINAL_REGIONS=~{write_lines(original_regions)}
    { echo "region" ; cat $ORIGINAL_REGIONS ; } >> "$ORIGINAL_REGIONS".headered
    ORIGINAL_CHROMS=~{write_lines(original_chroms)}
    { echo "chrom" ; cat $ORIGINAL_CHROMS ; } >> "$ORIGINAL_CHROMS".headered

    FOLLOWUP_REGIONS=~{write_lines(followup_regions)}
    { echo "region" ; cat $FOLLOWUP_REGIONS ; } >> "$FOLLOWUP_REGIONS".headered
    FOLLOWUP_CHROMS=~{write_lines(followup_chroms)}
    { echo "chrom" ; cat $FOLLOWUP_CHROMS ; } >> "$FOLLOWUP_CHROMS".headered

    original_finemap_creds=~{write_tsv(original_finemap_creds)}
    { echo "creds" ; cat $original_finemap_creds | sed -e 's/\t/,/g' ; } >> "$original_finemap_creds".headered
    original_susie_CSs=~{write_tsv(original_susie_CSs)}
    { echo "CSes" ; cat $original_susie_CSs | sed -e 's/\t/,/g' ; } >> "$original_susie_CSs".headered

    total_prob_finemap_creds=~{write_tsv(total_prob_finemap_creds)}
    { echo "creds" ; cat $total_prob_finemap_creds | sed -e 's/\t/,/g' ; } >> "$total_prob_finemap_creds".headered
    derived_prior_std_finemap_creds=~{write_tsv(derived_prior_std_finemap_creds)}
    { echo "creds" ; cat $derived_prior_std_finemap_creds | sed -e 's/\t/,/g' ; } >> "$derived_prior_std_finemap_creds".headered
    conv_tol_finemap_creds=~{write_tsv(conv_tol_finemap_creds)}
    { echo "creds" ; cat $conv_tol_finemap_creds | sed -e 's/\t/,/g' ; } >> "$conv_tol_finemap_creds".headered
    mac_finemap_creds=~{write_tsv(mac_finemap_creds)}
    { echo "creds" ; cat $mac_finemap_creds | sed -e 's/\t/,/g' ; } >> "$mac_finemap_creds".headered
    threshold_finemap_creds=~{write_tsv(threshold_finemap_creds)}
    { echo "creds" ; cat $threshold_finemap_creds | sed -e 's/\t/,/g' ; } >> "$threshold_finemap_creds".headered
    best_guess_susie_CSs=~{write_tsv(best_guess_susie_CSs)}
    { echo "CSes" ; cat $best_guess_susie_CSs | sed -e 's/\t/,/g' ; } >> "$best_guess_susie_CSs".headered

    low_prior_std_finemap_creds=~{write_tsv(low_prior_std_finemap_creds)}
    { echo "creds" ; cat $low_prior_std_finemap_creds | sed -e 's/\t/,/g' ; } >> "$low_prior_std_finemap_creds".headered
    ratio_finemap_creds=~{write_tsv(ratio_finemap_creds)}
    { echo "creds" ; cat $ratio_finemap_creds | sed -e 's/\t/,/g' ; } >> "$ratio_finemap_creds".headered
    ratio_susie_CSs=~{write_tsv(ratio_susie_CSs)}
    { echo "CSes" ; cat $ratio_susie_CSs | sed -e 's/\t/,/g' ; } >> "$ratio_susie_CSs".headered

    paste ~{write_objects(original_finemap_outputs)} "$ORIGINAL_CHROMS".headered "$ORIGINAL_REGIONS".headered "$original_finemap_creds".headered > original_finemap_files.tsv
    paste ~{write_objects(original_susie_outputs)} "$ORIGINAL_CHROMS".headered "$ORIGINAL_REGIONS".headered "$original_susie_CSs".headered > original_susie_files.tsv
    paste ~{write_objects(total_prob_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$total_prob_finemap_creds".headered > total_prob_finemap_files.tsv
    paste ~{write_objects(derived_prior_std_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$derived_prior_std_finemap_creds".headered > derived_prior_std_finemap_files.tsv
    paste ~{write_objects(conv_tol_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$conv_tol_finemap_creds".headered > conv_tol_finemap_files.tsv
    paste ~{write_objects(mac_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$mac_finemap_creds".headered > mac_finemap_files.tsv
    paste ~{write_objects(threshold_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$threshold_finemap_creds".headered > threshold_finemap_files.tsv
    paste ~{write_objects(best_guess_susie_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$best_guess_susie_CSs".headered > best_guess_susie_files.tsv
    paste ~{write_objects(low_prior_std_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$low_prior_std_finemap_creds".headered > low_prior_std_finemap_files.tsv
    paste ~{write_objects(ratio_finemap_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$ratio_finemap_creds".headered > ratio_finemap_files.tsv
    paste ~{write_objects(ratio_susie_outputs)} "$FOLLOWUP_CHROMS".headered "$FOLLOWUP_REGIONS".headered "$ratio_susie_CSs".headered > ratio_susie_files.tsv
    envsetup ~{script} \
      . \
      df \
      True \
      ~{phenotype_name} \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{sep=" " ethnic_str_assoc_results} \
      followup_conditions \
      original_finemap_files.tsv \
      original_susie_files.tsv \
      total_prob_finemap_files.tsv \
      derived_prior_std_finemap_files.tsv \
      conv_tol_finemap_files.tsv \
      mac_finemap_files.tsv \
      threshold_finemap_files.tsv \
      best_guess_susie_files.tsv \
      low_prior_std_finemap_files.tsv \
      ratio_finemap_files.tsv \
      ratio_susie_files.tsv
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    mem: "50g"
    dx_timeout: "1h"
  }
}

task followup_finemapping_conditions_comparison {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    Array[File] followup_conditions_tsvs
  }

  output {
    # intermediate tsvs
    File doubly_finemapped_STRs = "doubly_finemapped_STRs.tab"
    File confidently_finemapped_STRs = "confidently_finemapped_STRs.tab"
    File overconfidently_finemapped_STRs = "overconfidently_finemapped_STRs.tab"

    # used for deciding confidently finemapped
    File susie_best_guess_png = "susie_consistency_alpha_best_guess.png"
    File susie_best_guess_svg = "susie_consistency_alpha_best_guess.svg"
    File finemap_conv_tol_png = "finemap_consistency_pip_conv_tol.png"
    File finemap_conv_tol_svg = "finemap_consistency_pip_conv_tol.svg"
    File finemap_total_prob_png = "finemap_consistency_pip_total_prob.png"
    File finemap_total_prob_svg = "finemap_consistency_pip_total_prob.svg"
    File finemap_prior_std_derived_png = "finemap_consistency_pip_prior_std_derived.png"
    File finemap_prior_std_derived_svg = "finemap_consistency_pip_prior_std_derived.svg"
    File finemap_mac_png = "finemap_consistency_pip_mac.png"
    File finemap_mac_svg = "finemap_consistency_pip_mac.svg"
    File finemap_p_thresh_png = "finemap_consistency_pip_p_thresh.png"
    File finemap_p_thresh_svg = "finemap_consistency_pip_p_thresh.svg"

    # too conservative
    File finemap_ratio_png = "finemap_consistency_pip_ratio.png"
    File finemap_ratio_svg = "finemap_consistency_pip_ratio.svg"
    File susie_ratio_png = "susie_consistency_alpha_ratio.png"
    File susie_ratio_svg = "susie_consistency_alpha_ratio.svg"
    File finemap_prior_std_low_png = "finemap_consistency_pip_prior_std_low.png"
    File finemap_prior_std_low_svg = "finemap_consistency_pip_prior_std_low.svg"
  }

  command <<<
    envsetup ~{script} . followup_conditions_comparison . ~{sep=" " followup_conditions_tsvs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    mem: "50g"
  }
}

task str_tables_for_paper {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/str_tables_for_paper.py"
    File annotation_utils = "~{script_dir}/post_finemapping/annotation_utils.py"
    File python_array_utils = "~{script_dir}/post_finemapping/python_array_utils.py"

    File str_pos_table
    File repeat_units_table

    #annotations
    File intersects_gene
    File intersects_exon
    File intersects_CDS
    File intersects_five_prime_UTR
    File intersects_three_prime_UTR
    File intersects_UTR

    # only str assocs
    Array[File] phenotype_names # same order as all the assoc files
    Array[File] assocs
    Array[File] black_assocs
    Array[File] south_asian_assocs
    Array[File] chinese_assocs
    Array[File] irish_assocs
    Array[File] white_other_assocs

    Array[File] first_pass_finemapping_dfs
    Array[File] followup_finemapping_dfs
  }

  output {
    File singly_finemapped_strs_for_paper = "singly_finemapped_strs_for_paper.tab"
    File singly_finemapped_strs_sorted = "singly_finemapped_strs_sorted.tab"
    File confidently_finemapped_strs_for_paper = "confidently_finemapped_strs_for_paper.tab"
    File confidently_finemapped_strs_sorted = "confidently_finemapped_strs_sorted.tab"
  }

  command <<<
    ls ~{annotation_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      --outdir . \
      --str-pos-table ~{str_pos_table} \
      --repeat-units-table ~{repeat_units_table} \
      --intersects-gene-annotation ~{intersects_gene} \
      --intersects-exon-annotation ~{intersects_exon} \
      --intersects-CDS-annotation ~{intersects_CDS} \
      --intersects-five-prime-UTR-annotation ~{intersects_five_prime_UTR} \
      --intersects-three-prime-UTR-annotation ~{intersects_three_prime_UTR} \
      --intersects-UTR-annotation ~{intersects_UTR} \
      --assoc-phenotypes ~{sep=" " phenotype_names} \
      --assocs ~{sep=" " assocs} \
      --black-assocs ~{sep=" " black_assocs} \
      --south-asian-assocs ~{sep=" " south_asian_assocs} \
      --chinese-assocs ~{sep=" " chinese_assocs} \
      --irish-assocs ~{sep=" " irish_assocs} \
      --white-other-assocs ~{sep=" " white_other_assocs} \
      --first-pass-finemapping-dfs ~{sep=" " first_pass_finemapping_dfs} \
      --followup-finemapping-dfs ~{sep=" " followup_finemapping_dfs} \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "5h"
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
