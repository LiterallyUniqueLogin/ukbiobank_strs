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

struct serializable_FINEMAP_output {
  File snp_file
  File log_sss
  File config
}

struct FINEMAP_output {
  serializable_FINEMAP_output subset
  Array[File] creds
}

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

struct SuSiE_output {
  serializable_SuSiE_output subset 
  Array[File] CSs
}

####################### Extracting association signals #########################

task finemap_write_input_variants {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_write_input_variants.py"
    File python_array_utils = "~{script_dir}/finemapping/python_array_utils.py"
    File sample_utils = "~{script_dir}/finemapping/sample_utils.py"    

    File? str_assoc_results
    File snp_assoc_results
    File variants_to_filter
    File phenotype_samples_list
    String phenotype
    region bounds

    Float? snp_str_ratio
    Float? total_prob
    Int? mac
    File? snp_macs
    Float? inclusion_threshold
    Boolean is_binary = false
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
      ~{if defined(str_assoc_results) then select_first([str_assoc_results]) else "None"} \
      ~{variants_to_filter} \
      ~{phenotype_samples_list} \
      ~{phenotype} \
      ~{bounds.chrom} \
      ~{bounds.start} \
      ~{bounds.end} \
      ~{if defined(snp_str_ratio) then "--snp-str-ratio ~{snp_str_ratio}" else ""} \
      ~{if defined(total_prob) then "--total-prob ~{total_prob}" else ""} \
      ~{if defined(mac) then "--mac ~{select_first([mac])} ~{select_first([snp_macs])} " else ""} \
      ~{if defined(inclusion_threshold) then "--inclusion-threshold ~{inclusion_threshold}" else ""} \
      ~{if !is_binary then "" else "--is-binary"}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "6GB"
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
    memory: "6GB"
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
    memory: "4GB"
    continueOnReturnCode: [0, 79] # allow for timeouts, this will be handled in the retry workflow
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
    cpus: 4
    memory: "8GB"
    continueOnReturnCode: [0, 79] # allow for timeouts, this will be handled in the retry workflow
  }
}

task finemap_run {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/finemap_run.py"
    String finemap_command

    File master
    File zfile
    File all_variants_ld

    Int causal_snps
    Boolean prior_snps = false
    Float? prior_std
    Float? prob_conv_sss_tol
    String prefix
    Int cache_breaker = 0

    # TODO
    Int another_break = 13
  }

  String snp_file_loc = "~{prefix}finemap_output.snp"
  String log_sss_loc = "~{prefix}finemap_output.log_sss"
  String config_loc = "~{prefix}finemap_output.config"
  String creds_glob = "~{prefix}finemap_output.cred*"
  String finemap_input_z_loc = "~{prefix}finemap_input.z"

  output {
    FINEMAP_output finemap_output = object {
      subset: object {
        snp_file: snp_file_loc,
        log_sss: log_sss_loc,
        config: config_loc,
      },
      creds: glob(creds_glob)
    }
    File finemap_input_z = finemap_input_z_loc
  }

  command <<<
    # need all files in the same directory for FINEMAP to work
    # also necessary for dxCompiler to bother to localize this file
    ln -s ~{master} .
    ln ~{zfile} .
    ln -s ~{all_variants_ld} .
    which -a finemap
    envsetup ~{script} \
      . \
      ~{finemap_command} \
      --n-causal-snps ~{causal_snps} \
      ~{if prior_snps then "--prior-snps" else if defined(prior_std) then "--prior-std ~{prior_std}" else if defined(prob_conv_sss_tol) then "--prob-conv-sss-tol ~{prob_conv_sss_tol}" else ""}
    if [[ "~{prefix}" != "" ]] ; then 
      for file in finemap_input.z finemap_output.cred* finemap_output.config finemap_output.log_sss finemap_output.snp ; do
        ln $file ~{prefix}$file
      done
    fi
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.4"
    dx_timeout: if sub(sub(prefix,
      "gamma_glutamyltransferase_FINEMAP_.*_22_23776282_25846832", "foo"),
      "cystatin_c_FINEMAP_.*_20_19983208_26444535", "foo") == prefix then "4h" else "47h30m"
    memory: "8GB"
  }
}

task susie_choose_vars {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/susie_choose_vars.py"

    File? str_assoc_results
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
      ~{if defined(str_assoc_results) then select_first([str_assoc_results]) else "None"} \
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
    memory: "2GB"
  }
}

task any_vars_for_susie {
  input {
    File colnames
  }

  output {
    Boolean b = read_boolean(stdout())
  }

  command <<<
    if [[ -n "$(cat ~{colnames} | xargs)" ]] ; then
      echo true
    else
      echo false
    fi
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    memory: "2GB"
    shortTask: true
  }
}

task susie_load_gts {
  # TODO won't work properly if there are no shared covars
  # but there are covars in the pheno file
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
    File? shared_covars
    
    File colnames # from prev finemapping step
    String phenotype_name
    region bounds
    Boolean best_guess = false

    String time

    # TODO a random cache breaker
    Int ignored=3
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
      --pheno-out pheno_residual.h5 \
      ~{"--shared-covars-fname " + shared_covars} \
      ~{if best_guess then "--best-guess" else ""}
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: time
    memory: "4GB"
    continueOnReturnCode: [0, 79] # allow for timeouts, this will be handled in the retry workflow
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

    String prefix
  }

  output {
    File? lbf = "~{prefix}lbf.tab"
    File? lbf_variable = "~{prefix}lbf_variable.tab"
    File? sigma2 = "~{prefix}sigma2.txt"
    File? V = "~{prefix}V.tab"
    File? converged = "~{prefix}converged.txt"
    File? lfsr = "~{prefix}lfsr.tab"
    File? requested_coverage = "~{prefix}requested_coverage.txt"
    File? alpha = "~{prefix}alpha.tab"
    File? colnames = "~{prefix}colnames.txt"
    Array[File] CSs = glob("~{prefix}cs*.txt")
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
    if [[ "~{prefix}" != "" ]] ; then 
      for file in lbf.tab lbf_variable.tab sigma2.txt V.tab converged.txt lfsr.tab requested_coverage.txt alpha.tab cs*.txt ; do
        ln $file ~{prefix}$file
      done
    fi
    ln ~{colnames} ~{prefix}colnames.txt
    exit 0 # ignore the previous return code
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: time
    memory: mem
    continueOnReturnCode: [0, 79] # allow for timeouts, this will be handled in the retry workflow
  }
}

####################### Simulations ###################################

task simulate_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/simulate_phenotype.py"
    File load_and_filter_genotypes = "~{script_dir}/finemapping/load_and_filter_genotypes.py"
    File python_array_utils = "~{script_dir}/finemapping/python_array_utils.py"
    File sample_utils = "~{script_dir}/finemapping/sample_utils.py"    

    VCF str_vcf
    bgen snp_bgen
    File all_samples_list
    File samples_list
    Int chrom
    Int seed
    Array[String] causal_vars
    Array[String] causal_betas # should be floats, but cannot cast
  }

  output {
    File phenotype = "out.npy"
  }

  command <<<
    envsetup ~{script} \
      out.npy \
      ~{str_vcf.vcf} \
      ~{snp_bgen.bgen} \
      ~{all_samples_list} \
      ~{samples_list} \
      ~{chrom} \
      ~{seed} \
      --causal-vars ~{sep=" " causal_vars} \
      --causal-betas ~{sep=" " causal_betas} \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "20m"
    memory: "4GB"
    shortTask: true
  }
}

task select_causal_variants_from_finemap_output {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/select_causal_variants_from_finemap_output.py"

    File finemap_output_log
    File finemap_input_z
    Array[File]+ finemap_output_creds

    String prefix
  }

  output {
    Array[Array[String]] vars_and_betas = read_tsv(vars_and_betas_f)
    File vars_and_betas_f = "~{prefix}vars_and_betas.tab"
  }

  command <<<
    envsetup ~{script} \
       out \
       ~{finemap_output_log} \
       ~{finemap_input_z} \
       ~{sep=" " finemap_output_creds} \
      > ~{prefix}vars_and_betas.tab
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "20m"
    memory: "2GB"
    shortTask: true
  }
}

task select_causal_variants_from_susie_output {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/select_causal_variants_from_susie_output.py"

    File colnames
    File alpha
    File finemap_input_z
    Array[File]+ CSes

    String prefix
  }

  output {
    Array[Array[String]] vars_and_betas = read_tsv(vars_and_betas_f)
    File vars_and_betas_f = "~{prefix}vars_and_betas.tab"
  }

  command <<<
    envsetup ~{script} \
       out \
       ~{colnames} \
       ~{alpha} \
       ~{finemap_input_z} \
       ~{sep=" " CSes} \
      > ~{prefix}vars_and_betas.tab
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "20m"
    memory: "2GB"
    shortTask: true
  }
}

task any_causal_strs {
  input {
    Array[Array[String]] vars_and_betas
  }

  output {
    Boolean b = read_boolean(stdout())
  }

  command <<<
    python -c 'print("STR" in "~{sep=" " vars_and_betas[0]}")'
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    memory: "2GB"
    shortTask: true
  }
}

task bin_causal_variants_by_frequency {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/bin_causal_variants_by_frequency.py"
    File all_regions_finemapping_df
    File samples_file
    Array[File] mac_files
  }

  Int n_samples = length(read_lines(samples_file)) - 1

  output {
    Array[Array[String]] data = read_tsv(stdout())
    File data_file = write_tsv(data)
    Array[File] bin_effects = [
      "bin_001gthan0001_effects.txt",
      "bin_01gthan001_effects.txt",
      "bin_1gthan01_effects.txt",
      "bin_5gthan1_effects.txt"
    ]
#    File bin_lthan_0001_effects = "bin_lthan0001.txt"
#    File bin_001_gthan_0001_effects = "bin_001gthan0001_effects.txt"
#    File bin_01_gthan_001_effects = "bin_01gthan001_effects.txt"
#    File bin_1_gthan_01_effects = "bin_1gthan01_effects.txt"
#    File bin_5_gthan_1_effects = "bin_5gthan1_effects.txt"
  }

  command <<<
    envsetup ~{script} \
      ~{all_regions_finemapping_df} \
      ~{n_samples} \
      ~{sep=" " mac_files}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "20GB"
  }
}

task select_random_causal_variants {
  input {
    String script_dir
    File script = "~{script_dir}/finemapping/select_causal_variants_by_frequency.py"

    Int n_vars_to_choose
    Int seed
    File samples_file
    File snp_macs
    Array[Float]+ bin_weights
    Array[File]+ bin_effects

    String prefix
  }
  
  Int n_samples = length(read_lines(samples_file)) - 1

  output {
    Array[Array[String]] vars_and_betas = read_tsv(vars_and_betas_f)
    File vars_and_betas_f = "~{prefix}vars_and_betas.tab"
  }

  command <<<
    envsetup ~{script} \
      ~{n_vars_to_choose} \
      ~{seed} \
      ~{n_samples} \
      ~{snp_macs} \
      ~{sep=" " bin_weights} \
      ~{sep=" " bin_effects} \
      > ~{prefix}vars_and_betas.tab
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "20m"
    memory: "20GB"
    shortTask: true
  }
}

task fix_causal_var_str_coords {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/fix_causal_var_str_coords.py"

    Int chrom
    File in_vars_and_betas_f
    File flank_start_to_start_and_end_pos
  }

  output {
    File vars_and_betas_f = basename(in_vars_and_betas_f)
  }

  command <<<
    envsetup ~{script} \
      ~{chrom} \
      ~{in_vars_and_betas_f} \
      ~{flank_start_to_start_and_end_pos}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "20m"
    memory: "2GB"
    shortTask: true
  }
}

task compile_simulations_df {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/simulations_df.py"
    File script = "~{script_dir}/post_finemapping/simulations_df.py"

    Array[String] methods
    Array[Int] chroms
    Array[String] regions
    Array[Int] replicates
    Array[File] causal_vars_and_betas
    Array[File] snp_assocs
    Array[File] str_assocs
    Array[serializable_FINEMAP_output] finemap_outputs
    Array[Array[File]] finemap_creds
    Array[serializable_SuSiE_output] susie_outputs
    Array[Array[File]] susie_CSs
  }

  output {
    File df = "simulations_df.tab"
  }

  command <<<
    METHODS=~{write_lines(methods)}
    { echo "method" ; cat $METHODS ; } > METHODS.headered
    CHROMS=~{write_lines(chroms)}
    { echo "chrom" ; cat $CHROMS ; } > CHROMS.headered
    REGIONS=~{write_lines(regions)}
    { echo "region" ; cat $REGIONS ; } > REGIONS.headered
    REPLICATES=~{write_lines(replicates)}
    { echo "replicate" ; cat $REPLICATES ; } > REPLICATES.headered
    CAUSAL_VARS_AND_BETAS=~{write_lines(causal_vars_and_betas)}
    { echo "causal_vars_and_betas" ; cat $CAUSAL_VARS_AND_BETAS ; } > CAUSAL_VARS_AND_BETAS.headered
    SNP_ASSOCS=~{write_lines(snp_assocs)}
    { echo "snp_assoc" ; cat $SNP_ASSOCS ; } > SNP_ASSOCS.headered
    STR_ASSOCS=~{write_lines(str_assocs)}
    { echo "str_assoc" ; cat $STR_ASSOCS ; } > STR_ASSOCS.headered
    finemap_creds=~{write_tsv(finemap_creds)}
    { echo "creds" ; cat $finemap_creds | sed -e 's/\t/,/g' ; } > finemap_creds.headered
    susie_CSs=~{write_tsv(susie_CSs)}
    { echo "CSes" ; cat $susie_CSs | sed -e 's/\t/,/g' ; } > susie_CSs.headered
    paste \
      METHODS.headered \
      CHROMS.headered \
      REGIONS.headered \
      REPLICATES.headered \
      CAUSAL_VARS_AND_BETAS.headered \
      SNP_ASSOCS.headered \
      STR_ASSOCS.headered \
      ~{write_objects(finemap_outputs)} \
      finemap_creds.headered \
      ~{write_objects(susie_outputs)} \
      susie_CSs.headered \
      > simulations.tsv
    envsetup ~{script} simulations.tsv
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    memory: "120GB"
  }
}

task analyze_simulations {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/analyze_simulations.py"

    File simulations_df
  }

  output {
    File analysis = "analysis.tab"
  }

  command <<<
    ACTIVATE=new_polars envsetup ~{script} ~{simulations_df}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    memory: "20GB"
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
    Array[File]? ethnic_str_assoc_results

    Array[serializable_FINEMAP_output] original_finemap_outputs
    Array[Array[File]] original_finemap_creds
    Array[serializable_SuSiE_output] original_susie_outputs
    Array[Array[File]] original_susie_CSs
    Array[String] regions
    Array[Int] chroms

    Boolean is_binary = false
  }

  output {
    File all_regions_concordance = "finemapping_all_regions_concordance_~{phenotype_name}.tab"
    File susie_all_regions_min_abs_corrs = "susie_all_regions_min_abs_corrs_~{phenotype_name}.npy"
  }

  command <<<
    REGIONS=~{write_lines(regions)}
    { echo "region" ; cat $REGIONS ; } > REGIONS.headered
    CHROMS=~{write_lines(chroms)}
    { echo "chrom" ; cat $CHROMS ; } > CHROMS.headered
    original_finemap_creds=~{write_tsv(original_finemap_creds)}
    { echo "creds" ; cat $original_finemap_creds | sed -e 's/\t/,/g' ; } > original_finemap_creds.headered
    original_susie_CSs=~{write_tsv(original_susie_CSs)}
    { echo "CSes" ; cat $original_susie_CSs | sed -e 's/\t/,/g' ; } > original_susie_CSs.headered
    paste ~{write_objects(original_finemap_outputs)} CHROMS.headered REGIONS.headered original_finemap_creds.headered > finemap_files.tsv
    paste ~{write_objects(original_susie_outputs)} CHROMS.headered REGIONS.headered original_susie_CSs.headered > susie_files.tsv

    ethnic_str_assoc_results=~{sep=" " ethnic_str_assoc_results} 
    envsetup ~{script} \
      . \
      df \
      True \
      ~{phenotype_name} \
      ~{snp_assoc_results} \
      ~{str_assoc_results} \
      ~{if defined(ethnic_str_assoc_results) then "--ethnic_str_assoc_results $ethnic_str_assoc_results" else ""} \
      first_pass \
      finemap_files.tsv \
      susie_files.tsv \
      ~{if !is_binary then "" else "--is-binary "}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    memory: "50GB"
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
    memory: "50GB"
  }
}

task susie_finemap_venn_diagram {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/venn_diagram.py"

    Array[File] first_pass_dfs
  }

  output {
    File str_png = "all_susie_finemap_venn_STR.png"
    File str_svg = "all_susie_finemap_venn_STR.svg"
    File snp_png = "all_susie_finemap_venn_SNP.png"
    File snp_svg = "all_susie_finemap_venn_SNP.svg"
  }

  command <<<
    envsetup ~{script} \
      all_susie_finemap_venn \
      ~{sep=" " first_pass_dfs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10"
    memory: "8GB"
  }
}

task generate_followup_regions_tsv {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/finemapping_consistency.py"
    File graphing_utils = "~{script_dir}/post_finemapping/graphing_utils.py"

    File first_pass_df
    String phenotype
  }

  output {
    File tsv = "~{phenotype}_followup_finemapping_regions.tsv"
  }

  command <<<
    envsetup ~{script} . followup_regions ~{first_pass_df}
    mv followup_regions.tsv ~{phenotype}_followup_finemapping_regions.tsv
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
    memory: "2GB"
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

    Array[serializable_FINEMAP_output] repeat_finemap_outputs
    Array[Array[File]] repeat_finemap_creds
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
    { echo "region" ; cat $ORIGINAL_REGIONS ; } > ORIGINAL_REGIONS.headered
    ORIGINAL_CHROMS=~{write_lines(original_chroms)}
    { echo "chrom" ; cat $ORIGINAL_CHROMS ; } > ORIGINAL_CHROMS.headered

    FOLLOWUP_REGIONS=~{write_lines(followup_regions)}
    { echo "region" ; cat $FOLLOWUP_REGIONS ; } > FOLLOWUP_REGIONS.headered
    FOLLOWUP_CHROMS=~{write_lines(followup_chroms)}
    { echo "chrom" ; cat $FOLLOWUP_CHROMS ; } > FOLLOWUP_CHROMS.headered

    original_finemap_creds=~{write_tsv(original_finemap_creds)}
    { echo "creds" ; cat $original_finemap_creds | sed -e 's/\t/,/g' ; } > original_finemap_creds.headered
    original_susie_CSs=~{write_tsv(original_susie_CSs)}
    { echo "CSes" ; cat $original_susie_CSs | sed -e 's/\t/,/g' ; } > original_susie_CSs.headered

    repeat_finemap_creds=~{write_tsv(repeat_finemap_creds)}
    { echo "creds" ; cat $repeat_finemap_creds | sed -e 's/\t/,/g' ; } > repeat_finemap_creds.headered
    total_prob_finemap_creds=~{write_tsv(total_prob_finemap_creds)}
    { echo "creds" ; cat $total_prob_finemap_creds | sed -e 's/\t/,/g' ; } > total_prob_finemap_creds.headered
    derived_prior_std_finemap_creds=~{write_tsv(derived_prior_std_finemap_creds)}
    { echo "creds" ; cat $derived_prior_std_finemap_creds | sed -e 's/\t/,/g' ; } > derived_prior_std_finemap_creds.headered
    conv_tol_finemap_creds=~{write_tsv(conv_tol_finemap_creds)}
    { echo "creds" ; cat $conv_tol_finemap_creds | sed -e 's/\t/,/g' ; } > conv_tol_finemap_creds.headered
    mac_finemap_creds=~{write_tsv(mac_finemap_creds)}
    { echo "creds" ; cat $mac_finemap_creds | sed -e 's/\t/,/g' ; } > mac_finemap_creds.headered
    threshold_finemap_creds=~{write_tsv(threshold_finemap_creds)}
    { echo "creds" ; cat $threshold_finemap_creds | sed -e 's/\t/,/g' ; } > threshold_finemap_creds.headered
    best_guess_susie_CSs=~{write_tsv(best_guess_susie_CSs)}
    { echo "CSes" ; cat $best_guess_susie_CSs | sed -e 's/\t/,/g' ; } > best_guess_susie_CSs.headered

    low_prior_std_finemap_creds=~{write_tsv(low_prior_std_finemap_creds)}
    { echo "creds" ; cat $low_prior_std_finemap_creds | sed -e 's/\t/,/g' ; } > low_prior_std_finemap_creds.headered
    ratio_finemap_creds=~{write_tsv(ratio_finemap_creds)}
    { echo "creds" ; cat $ratio_finemap_creds | sed -e 's/\t/,/g' ; } > ratio_finemap_creds.headered
    ratio_susie_CSs=~{write_tsv(ratio_susie_CSs)}
    { echo "CSes" ; cat $ratio_susie_CSs | sed -e 's/\t/,/g' ; } > ratio_susie_CSs.headered

    paste ~{write_objects(original_finemap_outputs)} ORIGINAL_CHROMS.headered ORIGINAL_REGIONS.headered original_finemap_creds.headered > original_finemap_files.tsv
    paste ~{write_objects(original_susie_outputs)} ORIGINAL_CHROMS.headered ORIGINAL_REGIONS.headered original_susie_CSs.headered > original_susie_files.tsv
    paste ~{write_objects(repeat_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered repeat_finemap_creds.headered > repeat_finemap_files.tsv
    paste ~{write_objects(total_prob_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered total_prob_finemap_creds.headered > total_prob_finemap_files.tsv
    paste ~{write_objects(derived_prior_std_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered derived_prior_std_finemap_creds.headered > derived_prior_std_finemap_files.tsv
    paste ~{write_objects(conv_tol_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered conv_tol_finemap_creds.headered > conv_tol_finemap_files.tsv
    paste ~{write_objects(mac_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered mac_finemap_creds.headered > mac_finemap_files.tsv
    paste ~{write_objects(threshold_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered threshold_finemap_creds.headered > threshold_finemap_files.tsv
    paste ~{write_objects(best_guess_susie_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered best_guess_susie_CSs.headered > best_guess_susie_files.tsv
    paste ~{write_objects(low_prior_std_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered low_prior_std_finemap_creds.headered > low_prior_std_finemap_files.tsv
    paste ~{write_objects(ratio_finemap_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered ratio_finemap_creds.headered > ratio_finemap_files.tsv
    paste ~{write_objects(ratio_susie_outputs)} FOLLOWUP_CHROMS.headered FOLLOWUP_REGIONS.headered ratio_susie_CSs.headered > ratio_susie_files.tsv
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
      repeat_finemap_files.tsv \
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
    memory: "50GB"
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
    File finemap_repeat_png = "finemap_consistency_pip_repeat.png"
    File finemap_repeat_svg = "finemap_consistency_pip_repeat.svg"
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
    memory: "50GB"
  }
}

task str_tables_for_paper {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/str_tables_for_paper.py"
    File annotation_utils = "~{script_dir}/post_finemapping/annotation_utils.py"
    File python_array_utils = "~{script_dir}/post_finemapping/python_array_utils.py"

    File flank_start_to_start_and_end_pos
    File str_hg19_pos_bed
    File str_hg38_pos_bed
    File repeat_units_table

    #annotations
    Array[File] intersects_gene
    Array[File] intersects_exon
    Array[File] intersects_CDS
    Array[File] intersects_five_prime_UTR
    Array[File] intersects_three_prime_UTR
    Array[File] intersects_UTR

    # only str assocs
    Array[String] phenotype_names # same order as all the assoc files
    Array[File] assocs
    Array[File] black_assocs
    Array[File] south_asian_assocs
    Array[File] chinese_assocs
    Array[File] irish_assocs
    Array[File] white_other_assocs

    Array[File] first_pass_finemapping_dfs
    Array[File] followup_finemapping_dfs

    Array[File] wgs_comparison_stats
    Array[File] wgs_allele_freqs
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
      --flank-start-to-start-and-end-pos ~{flank_start_to_start_and_end_pos} \
      --str-hg19-pos-bed ~{str_hg19_pos_bed} \
      --str-hg38-pos-bed ~{str_hg38_pos_bed} \
      --repeat-units-table ~{repeat_units_table} \
      --intersects-gene-annotation ~{sep=" " intersects_gene} \
      --intersects-exon-annotation ~{sep=" " intersects_exon} \
      --intersects-CDS-annotation ~{sep=" " intersects_CDS} \
      --intersects-five-prime-UTR-annotation ~{sep=" " intersects_five_prime_UTR} \
      --intersects-three-prime-UTR-annotation ~{sep=" " intersects_three_prime_UTR} \
      --intersects-UTR-annotation ~{sep=" " intersects_UTR} \
      --assoc-phenotypes ~{sep=" " phenotype_names} \
      --assocs ~{sep=" " assocs} \
      --black-assocs ~{sep=" " black_assocs} \
      --south-asian-assocs ~{sep=" " south_asian_assocs} \
      --chinese-assocs ~{sep=" " chinese_assocs} \
      --irish-assocs ~{sep=" " irish_assocs} \
      --white-other-assocs ~{sep=" " white_other_assocs} \
      --first-pass-finemapping-dfs ~{sep=" " first_pass_finemapping_dfs} \
      --followup-finemapping-dfs ~{sep=" " followup_finemapping_dfs} \
      --wgs-comparison-stats ~{sep=" " wgs_comparison_stats} \
      --wgs-allele-freqs ~{sep=" " wgs_allele_freqs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "5h"
    memory: "50GB"
  }
}

task generate_QTL_table {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/gtex_STR2/table.py"

    File confidently_finemapped_table
    File all_QTL_results
    Boolean methylation = false
  }

  output {
    File QTL_table = if (methylation) then "meQTL_STRs.tab" else "qtl_STRs.tab"
  }

  command <<<
    envsetup ~{script} ~{confidently_finemapped_table} ~{all_QTL_results} ~{if(methylation) then "--methylation" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "50GB"
  }
}

task graph_main_hits {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/graph_main_hits_all_assocs.py"
    File annotation_utils = "~{script_dir}/post_finemapping/annotation_utils.py"
    File python_array_utils = "~{script_dir}/post_finemapping/python_array_utils.py"

    File hits_table # singly_finemapped_strs from str_tables_for_paper
    File eQTL_table # post_finemapping/gtex_STR2/blessed_qtl_STRs.tab
    File meQTL_table
    Array[File] closest_gene_annotations
  }

  output {
    File png = "results_graph.png"
    File svg = "results_graph.svg"
  }

  command <<<
    ls ~{annotation_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      results_graph.png \
      ~{hits_table} \
      ~{eQTL_table} \
      ~{meQTL_table} \
      ~{sep=" " closest_gene_annotations}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "2GB"
  }
}

task concordance_in_other_ethnicities {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/replication_comparison.py"

    File confidently_finemapped_STRs_df
    Array[File] first_pass_dfs
  }

  output {
    File nonwhite_replication_png = "replication_by_pval_non_white.png"
    File nonwhite_replication_svg = "replication_by_pval_non_white.svg"
    File white_replication_png = "replication_by_pval_white.png"
    File white_replication_svg = "replication_by_pval_white.svg"
    File stats = stdout()
  }

  command <<<
    envsetup ~{script} \
      . \
      ~{confidently_finemapped_STRs_df} \
      ~{sep=" " first_pass_dfs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "20GB"
  }
}

task generate_enrichments_table {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/enrichments_load.py"
    File annotation_utils = "~{script_dir}/post_finemapping/annotation_utils.py"
    File python_array_utils = "~{script_dir}/post_finemapping/python_array_utils.py"

    File flank_start_to_start_and_end_pos
    File str_loci
    File repeat_units_table
    File eSTR_table
    File gencode
    Array[String] phenotypes
    Array[File] str_assocs
    File confidently_finemapped_STRs_df

    #annotations
    Array[File] intersects_transcript_support_2
    Array[File] intersects_protein_coding_CDS_support_2
    Array[File] intersects_protein_coding_five_prime_UTR_support_2
    Array[File] intersects_protein_coding_three_prime_UTR_support_2
    Array[File] intersects_protein_coding_UTR_support_2
    Array[File] closest_downstream_protein_coding_exon_support_2
    Array[File] closest_upstream_protein_coding_exon_support_2
    Array[File] closest_downstream_protein_coding_gene
    Array[File] closest_upstream_protein_coding_gene
  }

  output {
    File table = "enrichments_df.tab"
  }

  command <<<
    ls ~{annotation_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      --outdir . \
      --flank-start-to-start-and-end-pos ~{flank_start_to_start_and_end_pos} \
      --str-loci ~{str_loci} \
      --repeat-units ~{repeat_units_table} \
      --eSTR-table ~{eSTR_table} \
      --gencode ~{gencode} \
      --phenotypes ~{sep=" " phenotypes} \
      --assocs ~{sep=" " str_assocs} \
      --confidently-finemapped-STRs-df ~{confidently_finemapped_STRs_df} \
      --intersects-protein-coding-CDS-support-2 ~{sep=" " intersects_protein_coding_CDS_support_2} \
      --intersects-protein-coding-UTR-support-2 ~{sep=" " intersects_protein_coding_UTR_support_2} \
      --intersects-protein-coding-five-prime-UTR-support-2 ~{sep=" " intersects_protein_coding_five_prime_UTR_support_2} \
      --intersects-protein-coding-three-prime-UTR-support-2 ~{sep=" " intersects_protein_coding_three_prime_UTR_support_2} \
      --intersects-transcript-support-2 ~{sep=" " intersects_transcript_support_2} \
      --closest-downstream-protein-coding-exon-support-2 ~{sep=" " closest_downstream_protein_coding_exon_support_2} \
      --closest-upstream-protein-coding-exon-support-2 ~{sep=" " closest_upstream_protein_coding_exon_support_2} \
      --closest-downstream-protein-coding-gene ~{sep=" " closest_downstream_protein_coding_gene} \
      --closest-upstream-protein-coding-gene ~{sep=" " closest_upstream_protein_coding_gene}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "3h"
    memory: "50GB"
  }
}

task calc_enrichments {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/enrichments_calc.py"

    File enrichment_df
  }

  output {
    File enrichment_stats = "enrichments.tab"
  }

  command <<<
    envsetup ~{script} \
      . \
      ~{enrichment_df}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "10GB"
  }
}


task graph_enrichments {
  input {
    String script_dir
    File script = "~{script_dir}/post_finemapping/enrichments_plot.py"

    File enrichment_stats
  }

  output {
    File regions_svg = "enrichment_barplots_regions.svg"
    File regions_png = "enrichment_barplots_regions.png"
    File repeats_svg = "enrichment_barplots_repeats.svg"
    File repeats_png = "enrichment_barplots_repeats.png"
  }

  command <<<
    envsetup ~{script} \
      . \
      ~{enrichment_stats}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "10GB"
  }
}

task placeholder {
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
    memory: "2GB"
  }
}
