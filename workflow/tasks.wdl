version 1.0

# files not generated:
# anything prefixed with data_showcase
# white_brits_sample_list

# sample_list file format:
#  tabular: ID_1 ID_2 missing sex
#  space delim
#  missing is always 0 . sex is 1 (male?) or 2
#  ID1 and 2 and identical
#  ID_1 can be negative for removed samples

# ssample_list file format (= simple sample list):
# first line is 'ID' (case insensitive)
# every successive line is a sample ID

# TODO: set container for each tas

####################### Loading samples and phenotypes ####################

task ethnic_sample_lists {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/ethnicity.py"
    File python_array_utils = "~{script_dir}/sample_qc/scripts/python_array_utils.py"

    File white_brits_sample_list # sample_qc/common_filters/ethnicity_white_brits.sample
    File data_showcase_ethnicity_self_report # 21000
  } 

  output {
    # sample_qc/common_filters/ethnicity/{ethnicity}.sample
    # ethnicity -> File
    Map[String, File] sample_lists = {
      'black': 'black.sample',
      'south_asian': 'south_asian.sample',
      'chinese': 'chinese.sample',
      'irish': 'irish.sample',
      'white_other': 'white_other.sample',
    }

    # I don't think I could scatter over this
    # File black_sample_list = 'black.sample'
    # File south_asian_sample_list = 'south_asian.sample'
    # File chinese_sample_list = 'chinese.sample'
    # File irish_sample_list = 'irish.sample'
    # File white_other_sample_list = 'white_other.sample'
  }

  command <<<
    ~{script} . ~{white_brits_sample_list} ~{data_showcase_ethnicity_self_report}
  >>>

  runtime {
    shortTask: true
    dx_timeout: "5m"
  }
}

task qced_sample_list {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/combine.py"

    File unqced_sample_list # sample_qc/common_filters/ethnicity/{ethnicity}.sample

    # these may also be not-simple sample lists
    # if so, ignore the extra columns
    Array[File]+ ssample_lists_to_filter  # sample_qc/common_filters/remove/*sample

    File? subpop_ssample_list # sample_qc/subpops/{subpop}.txt
  }

  String outfname = 'qced.samples'

  output {
     File qced_sample_list = outfname # sample_qc/(subpop_)?runs/({subpop}/)?{ethnicity}/no_phenotype/combined.sample
  }

  command <<<
    ~{script} 
       ~{outfname}
       discard
       ~{unqced_sample_list}
       ~{sep=' ' ssample_lists_to_filter}
       ~{"--subpop " + subpop_ssample_list}
  >>>  

  runtime {
    shortTask: true
    dx_timeout: "5m"
  }
}

task load_shared_covars {
  input {
    String script_dir
    File script = "~{script_dir}/traits/load_shared_covars.py"

    File fam_file
    File sample_qc_file # TODO replace this with PCs data field 22009
    File assessment_ages_file
  }

  output {
    File shared_covars = "shared_covars.npy"
    File covar_names = "covar_names.txt"
    File assessment_ages = "assessment_ages.npy"
  }

  command <<<
    ~{script} . ~{fam_file} ~{sample_qc_file} ~{assessment_ages_file}
  >>>

  runtime {
    memory: '10g'

    dx_timeout: "15m"
  }
}

#task load_phenotype {
#  input {
#    String script_dir
#    # File script = "~{script_dir}/traits/"
#
#    File qced_sample_list
#    String ethnicity #TODO is this needed if we're passing in files?
#    String phenotype #TODO is this needed if we're passing in files?
#    Int data_field_id
#    # TODO ...
#  }
#
#  output {
#    File data = "phenotype.npy"
#  }
#
#  command <<<
#    ~{script}
#      phenotype
#      ~{qced_sample_list}
#      ~{ethnicity}
#      ~{phenotype}
#  >>>
#
#  runtime {
#
#  }
#}
