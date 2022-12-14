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

# sample_list file format:
#  tabular: ID_1 ID_2 missing sex
#  space delim
#  missing is always 0 . sex is 1 (male?) or 2
#  ID1 and 2 and identical
#  ID_1 can be negative for removed samples

# ssample_list file format (= simple sample list):
# first line is 'ID' (case insensitive)
# every successive line is a sample ID

# TODO: set container for each task

####################### Loading samples and phenotypes ####################

task ethnic_sample_lists {
  input {
    String script_dir
    File script = "~{script_dir}/sample_qc/scripts/ethnicity.py"
    File python_array_utils = "~{script_dir}/sample_qc/scripts/python_array_utils.py"

    File white_brits_sample_list
    File sc_ethnicity_self_report # 21000
  } 

  output {
    # sample_qc/common_filters/ethnicity/{ethnicity}.sample
    # ethnicity -> File
    Array[String] ethnicities = [
      'black',
      'south_asian',
      'chinese',
      'irish',
      'white_other',
    ]
    Map[String, File] sample_lists = {
      'black': 'black.sample',
      'south_asian': 'south_asian.sample',
      'chinese': 'chinese.sample',
      'irish': 'irish.sample',
      'white_other': 'white_other.sample',
    }
  }

  command <<<
    ~{script} . ~{white_brits_sample_list} ~{sc_ethnicity_self_report}
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

    File unqced_sample_list # white_brits user input or output from ethnic_sample_lists

    # these may also be not-simple sample lists
    # if so, ignore the extra columns
    Array[File]+ ssample_lists_to_filter  # TODO move this to expanse workflow sample_qc/common_filters/remove/*sample

    File? subpop_ssample_list # TODO move this to expanse workflow sample_qc/subpops/{subpop}.txt
  }

  String outfname = 'qced.samples'

  output {
     File qced_sample_list = outfname # sample_qc/(subpop_)?runs/({subpop}/)?{ethnicity}/no_phenotype/combined.sample
  }

  command <<<
    ~{script} \
      ~{outfname} \
      discard \
      ~{unqced_sample_list} \
      ~{sep=' ' ssample_lists_to_filter} \
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
    File sc_assessment_ages
  }

  output {
    # all in traits/shared_covars/
    File shared_covars = "shared_covars.npy" 
    File covar_names = "covar_names.txt"
    File assessment_ages = "assessment_ages.npy"
  }

  command <<<
    ~{script} . ~{fam_file} ~{sample_qc_file} ~{sc_assessment_ages}
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
#    ~{script} \
#      phenotype \
#      ~{qced_sample_list} \
#      ~{ethnicity} \
#      ~{phenotype}
#  >>>
#
#  runtime {
#
#  }
#}
