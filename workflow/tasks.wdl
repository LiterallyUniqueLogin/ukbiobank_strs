version 1.0

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

# TODO: set container for each task

####################### Helper tasks ###########################

task concatenate_tsvs {
  input {
    Array[File]+ tsvs
  }

  output {
    File tsv = "out.tab"
  }

  command <<<
    for file in ~{sep=" " tsvs} ; do 
      head -1 $file
    done | uniq -c | wc -l | # TODO assert 1
    
    head -1 ~{tsvs[0]} > ~{tsv}
    for file in ~{sep=" " tsvs} ; do
      tail -n +2 $file >> ~{tsv}
    done
  >>>

  runtime {
    docker: "ukb_strs"
    dx_timeout: "4h"
  }
}

####################### Loading samples and phenotypes ####################

task write_sample_list {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/write_sample_list.py"

    File sc
    Int? value
  }

  output {
    File data = "data.out"
  }

  command <<<
    ~{script} ~{sc} data.out ~{"--value " +  value}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task ethnic_sample_lists {
  input {
    String script_dir
    File script_ = "~{script_dir}/sample_qc/scripts/ethnicity.py"
    File python_array_utils = "~{script_dir}/sample_qc/scripts/python_array_utils.py"

    File white_brits_sample_list # write_sample_list 22006
    File sc_ethnicity_self_report # 21000
  } 

  output {
    # sample_qc/common_filters/ethnicity/{ethnicity}.sample
    Array[String] ethnicities = [
      "black",
      "south_asian",
      "chinese",
      "irish",
      "white_other",
    ]
    # These can be zipped together to form a map if desired
    Array[File] sample_lists = [
      "black.sample",
      "south_asian.sample",
      "chinese.sample",
      "irish.sample",
      "white_other.sample",
    ]
  }

  command <<<
    ~{script_} . ~{white_brits_sample_list} ~{sc_ethnicity_self_report}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task sex_mismatch_sample_list {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/find_sex_mismatch_list.py"

    File sc_genetic_sex #22001
    File sc_reported_sex #31
  }

  output {
    File data = "out.sample"
  }

  command <<<
    ~{script} ~{sc_genetic_sex} ~{sc_reported_sex} out.sample
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task qced_sample_list {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/combine.py"

    File unqced_sample_list # white brits = write_sample_list 22006 or output from ethnic_sample_lists
    File withdrawn_sample_list 
    File sex_mismatch_sample_list # task above
    File sex_aneuploidy_sample_list # write_sample_list 22019
    File low_genotyping_quality_sample_list # write_sample_list 22021 -1

    File? subpop_sample_list # TODO move this to expanse workflow sample_qc/subpops/{subpop}.txt
  }

  String outfname = "qced.samples"

  output {
     File data = outfname # sample_qc/(subpop_)?runs/({subpop}/)?{ethnicity}/no_phenotype/combined.sample
  }

  command <<<
    ~{script} \
      ~{outfname} \
      discard \
      ~{unqced_sample_list} \
      ~{withdrawn_sample_list} \
      ~{sex_mismatch_sample_list} \
      ~{sex_aneuploidy_sample_list} \
      ~{low_genotyping_quality_sample_list} \
      ~{"--subpop " + subpop_sample_list}
  >>>  

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task load_shared_covars {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_shared_covars.py"

    File fam_file
    File sc_pcs # 22009
    File sc_assessment_ages
  }

  output {
    # all in traits/shared_covars/
    File shared_covars = "shared_covars.npy" 
    File covar_names = "covar_names.txt"
    File assessment_ages = "assessment_ages.npy"
  }

  command <<<
    ~{script} . ~{fam_file} ~{sc_pcs} ~{sc_assessment_ages}
  >>>

  runtime {
    docker: "ukb_strs"
    memory: "10g"

    dx_timeout: "15m"
  }
}

task load_continuous_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_continuous_phenotype_from_main_dataset.py"

    File sc
    File qced_sample_list # from qced_sample_list

    File assessment_ages_npy # from load shared covars
    Array[String] categorical_covariate_names
    Array[File] categorical_covariate_scs
  }

  output {
    File data = "pheno.npy"
    File covar_names = "pheno_covar_names.txt"
    File README = "pheno_README.txt"
  }

  command <<<
    ~{script} \
      ~{sc} \
      '.' \
      ~{qced_sample_list} \
      ~{assessment_ages_npy} \
      --categorical-covariate-names ~{sep=" " categorical_covariate_names} \
      --categorical-covariate-files ~{sep=" " categorical_covariate_scs}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task load_binary_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_binary_phenotype_from_main_dataset.py"

    File sc
    File qced_sample_list # from qced_sample_list

    File sc_year_of_birth # 34
    File sc_month_of_birth # 52
    File sc_date_of_death # 40000
    String date_of_most_recent_first_occurrence_update
    Boolean is_zero_one_neg_nan = false
  }

  output {
    File data = "pheno.npy"
    File covar_names = "pheno_covar_names.txt"
    File README = "pheno_README.txt"
  }

  command <<<
    ~{script} \
      ~{sc} \
      '.' \
      ~{qced_sample_list} \
      ~{sc_year_of_birth} \
      ~{sc_month_of_birth} \
      ~{sc_date_of_death} \
      ~{date_of_most_recent_first_occurrence_update} \
      ~{if is_zero_one_neg_nan then "--zero-one-neg-nan" else ""}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

task unrelated_samples_for_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/unrelated_individuals.py"
    String PRIMUS_command

    File pheno_data # load_xxx_phenotype.data
    File kinship
    Boolean is_binary
  }

  output {
    File data = "out.samples"
  }

  command <<<
    ~{script} out.samples ~{kinship} ~{pheno_data} ~{PRIMUS_command} ~{if is_binary then "--binary-pheno" else ""}
  >>>
  
  runtime {
    docker: "ukb_strs"
    dx_timeout: "24h"
  }
}

task transform_trait_values {
  input {
    String script_dir
    File script = "~{script_dir}/traits/transform_traits.py"

    File pheno_data # from task
    File unrelated_samples_for_phenotype  # from task
    Boolean is_binary
  }

  output {
    File README = "out_README.txt"
    File data = "out.npy"
  }

  command <<<
    ~{script} out ~{pheno_data} ~{unrelated_samples_for_phenotype} ~{if is_binary then "--binary" else ""}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "5m"
  }
}

######## STR calling and QC ###########

# cbl allele dist
task fig_4a {
  input {
    String script_dir = "association/cbl_allele_histo.py'"
    File script = "~{script_dir}/association/cbl_allele_histo.py"
    # TODO sample_utils

    File str_vcf_chr_11
    File specific_alleles
  }

  output {
    File png = "cbl.png"
    File svg = "cbl.svg"
  }

  command <<<
    ~{script} . ~{str_vcf_chr_11} ~{specific_alleles}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "10m"
  }
}

########### Running and plotting associations #############

task str_spot_test {
  input {
    String script_dir
    File script = "~{script_dir}/association/my_regional_gwas.py"
    File weighted_binom_conf = "~{script_dir}/association/weighted_bionm_conf.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"
    File load_and_filter_genotypes = "~{script_dir}/association/load_and_filter_genotypes.py"

    File str_vcf
    File shared_covars # from task
    File untransformed_phenotype # from task
    File transformed_phenotype # from task
    Boolean is_binary
    Int chrom
    Int pos
    String phenotype_name
  }

  output {
    File README = "README.txt"
    File data = "out.tab"
  }

  # TODO remove:
  # project temp
  # sample utils issues
  # trtools
  command <<<
    ~{script} \
      ~{README} \
      strs \
      ~{phenotype_name} \
      --readme \
      --str-vcf ~{str_vcf} \
      ~{if is_binary then "--binary logistic" else ""}
    ~{script} \
      ~{data} \
      strs \
      ~{phenotype_name} \
      --region ~{chrom}:~{pos}-~{pos} \
      --pheno-and-covars ~{transformed_phenotype} \
      --shared-covars ~{shared_covars} \
      --untransformed-phenotypes ~{untransformed_phenotype} \
      --str-vcf ~{str_vcf} \
      ~{if is_binary then "--binary logistic" else ""}
  >>>

  runtime {
    docker: "ukb_strs"
    shortTask: true
    dx_timeout: "20m"
  }
}

# TODO my str gwas README?

task regional_my_str_gwas {
  input {
    String script_dir
    File script = "~{script_dir}/association/my_regional_gwas.py"
    File weighted_binom_conf = "~{script_dir}/association/weighted_bionm_conf.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"
    File load_and_filter_genotypes = "~{script_dir}/association/load_and_filter_genotypes.py"

    File str_vcf
    File shared_covars # from task
    File untransformed_phenotype # from task
    File transformed_phenotype # from task
    Boolean is_binary
    String binary_type # linear or logistic, only needs to be set if is_binary == true
    Int chrom
    Int start_pos
    Int end_pos
    String phenotype_name
  }

  output {
    File data = "out.tab"
  }
  
  # TODO remove:
  # project temp
  # sample utils issues
  # trtools
  command <<<
    ~{script} \
      ~{data} \
      strs \
      ~{phenotype_name} \
      --region ~{chrom}:~{start_pos}-~{end_pos} \
      --pheno-and-covars ~{transformed_phenotype} \
      --shared-covars ~{shared_covars} \
      --untransformed-phenotypes ~{untransformed_phenotype} \
      --str-vcf ~{str_vcf} \
      ~{if is_binary then "--binary " + binary_type else ""}
  >>>

  runtime {
    docker: "ukb_strs"
    dx_timeout: "12h"
    memory: if binary_type == "logistic" then "40g" else "4g"
  }
}

task prep_continuous_plink_input {
  input {
    String script_dir
    File script = "~{script_dir}/prep_plink_input.py "

    File shared_covars # from task
    File shared_covar_names # from task
    File transformed_phenotype # from task
    File pheno_covar_names # from task
    Boolean is_binary
    String binary_type # linear or logistic, only needs to be set if is_binary == true
    String phenotype_name
  }

  output {
    File data = "out.tab"
  }

  command <<<
    ~{script} \
      out.tab \
      ~{phenotype_name} \
      ~{transformed_phenotype} \
      ~{pheno_covar_names} \
      ~{shared_covars} \
      ~{shared_covar_names}
      ~{if is_binary then "--binary " + binary_type else ""}
  >>>

  runtime {
    docker: "ukb_strs"
    dx_timeout: "30m"
  }
}

task chromosomal_plink_snp_association {
  input {
    String script_dir
    File script = "~{script_dir}/plink_association.sh"
    String plink_command

    File imputed_snp_p_file # TODO could generate this with task
    File pheno_data # from task

    Int chrom
    String phenotype_name
    String binary_type # linear or linear_binary or logistic
  }

  output {
    File data = "plink2." + (if binary_type == "linear" then "RIN_" else "") + phenotype_name + ".glm." + (if binary_type == "logistic" then "logistic.hybrid" else "linear") + ".done"
    File log = "plink2.log"
  }

  command <<<
    PHENOTYPE=~{phenotype_name} \
    BINARY_TYPE=~{binary_type} \
    CHROM=~{chrom} \
    OUT_DIR=. \
    PHENO_FILE=~{pheno_data} \
    P_FILE=~{imputed_snp_p_file} \
    PLINK_EXECUTABLE=~{plink_command} \
    ~{script}
  >>>

  runtime {
    docker: "ukb_strs"
    memory: "56g"
    cpu: "28"

    dx_timeout: "24h"
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
    ~{script}
  >>>

  runtime {
    docker: "ukb_strs"
    dx_timeout: ""
  }
}

