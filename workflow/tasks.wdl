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

####################### Common structs #########################

struct VCF {
  File vcf
  File index
}

struct PFiles {
  File pgen
  File psam
  File pvar
}

####################### Helper tasks ###########################

task concatenate_tsvs {
  input {
    Array[File]+ tsvs
  }

  output {
    File tsv = "out.tab"
  }

  command <<<
    envsetup bash -c '
      if (( $(for file in ~{sep=" " tsvs} ; do 
        head -1 $file
      done | uniq -c | wc -l) != 1 )) ; then
        echo "Different headers" 2>&1
        exit 1
      fi
      
      head -1 ~{tsvs[0]} > out.tab
      for file in ~{sep=" " tsvs} ; do
        tail -n +2 $file >> out.tab
      done
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
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
    envsetup ~{script} ~{sc} data.out ~{"--value " +  value}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

task ethnic_sample_lists {
  input {
    String script_dir
    File script_ = "~{script_dir}/sample_qc/scripts/ethnicity.py"

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
    envsetup ~{script_} . ~{white_brits_sample_list} ~{sc_ethnicity_self_report}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
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
    envsetup ~{script} ~{sc_genetic_sex} ~{sc_reported_sex} out.sample
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
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
    envsetup ~{script} \
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
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

task load_shared_covars {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_shared_covars.py"
    File python_array_utils = "~{script_dir}/traits/python_array_utils.py"

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
    envsetup ~{script} . ~{fam_file} ~{sc_pcs} ~{sc_assessment_ages}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    mem: "10g"

    dx_timeout: "15m"
  }
}

task load_continuous_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_continuous_phenotype_from_main_dataset.py"
    File python_array_utils = "~{script_dir}/traits/python_array_utils.py"

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
    envsetup ~{script} \
      ~{sc} \
      'pheno' \
      ~{qced_sample_list} \
      ~{assessment_ages_npy} \
      --categorical-covar-names ~{sep=" " categorical_covariate_names} \
      --categorical-covar-files ~{sep=" " categorical_covariate_scs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
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
    envsetup ~{script} \
      ~{sc} \
      'pheno' \
      ~{qced_sample_list} \
      ~{sc_year_of_birth} \
      ~{sc_month_of_birth} \
      ~{sc_date_of_death} \
      ~{date_of_most_recent_first_occurrence_update} \
      ~{if is_zero_one_neg_nan then "--zero-one-neg-nan" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

task write_sample_list_for_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/traits/write_sample_list_from_pheno_array.py"

    File pheno_data # from load_xxx_phenotype
  }

  output {
    File data = "out.samples"
  }

  command <<<
    envsetup ~{script} ~{pheno_data} out.samples
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

task unrelated_samples {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/unrelated_individuals.py"
    File sample_utils = "~{script_dir}/sample_qc/scripts/sample_utils.py"
    File python_array_utils = "~{script_dir}/sample_qc/scripts/python_array_utils.py"
    String PRIMUS_command

    File kinship
    File sample_list # from task qced_sample_list (for ethnicities) or write_sample_list_for_phenotype (for phenotypes)
    File? binary_pheno_data # task load_binary_phenotype.data
  }

  output {
    File data = "out.samples"
  }

  command <<<
    envsetup ~{script} out.samples ~{kinship} ~{sample_list} ~{PRIMUS_command} ~{"--binary-pheno " + binary_pheno_data}
  >>>
  
  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "24h"
  }
}

task transform_trait_values {
  input {
    String script_dir
    File script = "~{script_dir}/traits/transform_traits.py"
    File python_array_utils = "~{script_dir}/traits/python_array_utils.py"

    File pheno_data # from task
    File unrelated_samples_for_phenotype  # from task
    Boolean is_binary
  }

  output {
    File README = "out_README.txt"
    File data = "out.npy"
  }

  command <<<
    envsetup ~{script} out ~{pheno_data} ~{unrelated_samples_for_phenotype} ~{if is_binary then "--binary" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

######## STR calling and QC ###########

# TODO could do length confusion, pre imputation allele freqs, calc macs

# cbl allele dist
task fig_4a {
  input {
    String script_dir
    File script = "~{script_dir}/association/cbl_allele_histo.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    File all_samples_list
    File white_brits_sample_list # 
    File black_sample_list # task ethnic_samples_lists
    File south_asian_sample_list # task ethnic_samples_lists
    File chinese_sample_list # task ethnic_samples_lists

    VCF str_vcf_chr_11
    File specific_alleles
  }

  output {
    File png = "cbl.png"
    File svg = "cbl.svg"
  }

  command <<<
    envsetup ~{script} . ~{str_vcf_chr_11.vcf} ~{specific_alleles} ~{all_samples_list} \
      ~{white_brits_sample_list} ~{black_sample_list} ~{south_asian_sample_list} ~{chinese_sample_list}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "10m"
  }
}

########### Running and plotting associations #############

task association_regions {
  input {
    File chr_lens
    Int region_len
  }

  output {
    Array[Array[String]] out_tsv = read_tsv(stdout())
  }

  command <<<
    envsetup python -c "
    import numpy as np
    chr_lens = np.genfromtxt(
      '~{chr_lens}',
      usecols=[1],
      skip_header=1,
      dtype=int
    )

    chroms = []
    starts = []
    ends = []
    for chrom in range(1, 23):
      chr_len = chr_lens[chrom-1]
      for start in range(1, chr_len, ~{region_len}):
        if start + ~{region_len} - 1 > chr_len:
          end = chr_len
        else:
          end = start + ~{region_len} - 1
        chroms.append(chroms)
        starts.append(start)
        ends.append(end)
    for chrom, start, end in zip(chroms, starts, ends):
      print(f'{chrom}\t{start}\t{end}')
    "
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
  }
}

task str_spot_test {
  input {
    String script_dir
    File script = "~{script_dir}/association/my_regional_gwas.py"
    File weighted_binom_conf = "~{script_dir}/association/weighted_binom_conf.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"
    File load_and_filter_genotypes = "~{script_dir}/association/load_and_filter_genotypes.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"

    VCF str_vcf
    File shared_covars # from task
    File untransformed_phenotype # from task
    File transformed_phenotype # from task
    File all_samples_list
    Boolean is_binary
    Int chrom
    Int pos
    String phenotype_name
    String temp_dir
  }

  output {
    File README = "README.txt"
    File data = "out.tab"
  }

  # TODO remove:
  # trtools
  command <<<
    envsetup ~{script} \
      README.txt \
      strs \
      ~{phenotype_name} \
      --readme \
      --str-vcf ~{str_vcf.vcf} \
      ~{if is_binary then "--binary logistic" else ""} && \
    envsetup ~{script} \
      out.tab \
      strs \
      ~{phenotype_name} \
      --region ~{chrom}:~{pos}-~{pos} \
      --pheno-and-covars ~{transformed_phenotype} \
      --shared-covars ~{shared_covars} \
      --untransformed-phenotypes ~{untransformed_phenotype} \
      --all-samples-fname ~{all_samples_list} \
      --str-vcf ~{str_vcf.vcf} \
      --temp-dir ~{temp_dir} \
      ~{if is_binary then "--binary logistic" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "20m"
  }
}

# TODO my str gwas README?

task regional_my_str_gwas {
  input {
    String script_dir
    File script = "~{script_dir}/association/my_regional_gwas.py"
    File weighted_binom_conf = "~{script_dir}/association/weighted_binom_conf.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"
    File load_and_filter_genotypes = "~{script_dir}/association/load_and_filter_genotypes.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"

    VCF str_vcf
    File shared_covars # from task
    File untransformed_phenotype # from task
    File transformed_phenotype # from task
    File all_samples_list
    Boolean is_binary
    String binary_type # linear or logistic, only needs to be set if is_binary == true
    Int chrom
    Int start_pos
    Int end_pos
    String phenotype_name
    String temp_dir
  }

  output {
    File data = "out.tab"
  }
  
  # TODO remove:
  # trtools
  command <<<
    envsetup ~{script} \
      out.tab \
      strs \
      ~{phenotype_name} \
      --region ~{chrom}:~{start_pos}-~{end_pos} \
      --pheno-and-covars ~{transformed_phenotype} \
      --shared-covars ~{shared_covars} \
      --untransformed-phenotypes ~{untransformed_phenotype} \
      --all-samples-fname ~{all_samples_list} \
      --str-vcf ~{str_vcf.vcf} \
      --temp-dir ~{temp_dir} \
      ~{if is_binary then "--binary " + binary_type else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "12h"
    mem: if binary_type == "logistic" then "40g" else "4g"
  }
}

# TODO could do my GWAS for SNPs

task prep_plink_input {
  input {
    String script_dir
    File script = "~{script_dir}/association/prep_plink_input.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

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
    envsetup ~{script} \
      out.tab \
      ~{phenotype_name} \
      ~{transformed_phenotype} \
      ~{pheno_covar_names} \
      ~{shared_covars} \
      ~{shared_covar_names}
      ~{if is_binary then "--binary " + binary_type else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
  }
}

task chromosomal_plink_snp_association {
  input {
    String script_dir
    File script = "~{script_dir}/plink_association.sh"
    String plink_command

    PFiles imputed_snp_p_file # TODO could generate this with task
    File pheno_data # from task
    String temp_dir

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
    P_FILE=~{sub(imputed_snp_p_file.pgen, ".pgen$", "")} \
    PLINK_EXECUTABLE=~{plink_command} \
    PROJECT_TEMP=~{temp_dir} \
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    mem: "56g"
    cpus: "28"

    dx_timeout: "24h"
  }
}

## TODO append mfi to plink logistic run
## TODO compare my to plink GWAS

# TODO , also many variants
task manhattan {
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

