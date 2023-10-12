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

struct bgen {
  File bgen
  File index
  # precomputed indices from the python program bgen_reader
  File? bgen_reader_metadata
  File? bgen_reader_metadata2
  File? bgen_reader_complex_metadata2
}

struct PFiles {
  File pgen
  File psam
  File pvar
}

struct region {
  Int chrom
  Int start
  Int end
}

task phenotype_names {
  output { 
		Array[String] n = [
			"alanine_aminotransferase",
			"albumin",
			"alkaline_phosphatase",
			"apolipoprotein_a",
			"apolipoprotein_b",
			"aspartate_aminotransferase",
			"c_reactive_protein",
			"calcium",
			"cholesterol",
			"creatinine",
			"cystatin_c",
			"eosinophil_count",
			"eosinophil_percent",
			"gamma_glutamyltransferase",
			"glucose",
			"glycated_haemoglobin",
			"haematocrit",
			"haemoglobin_concentration",
			"hdl_cholesterol",
			"igf_1",
			"ldl_cholesterol_direct",
			"lymphocyte_count",
			"lymphocyte_percent",
			"mean_corpuscular_haemoglobin",
			"mean_corpuscular_haemoglobin_concentration",
			"mean_corpuscular_volume",
			"mean_platelet_volume",
			"mean_sphered_cell_volume",
			"neutrophil_count",
			"neutrophil_percent",
			"phosphate",
			"platelet_count",
			"platelet_crit",
			"platelet_distribution_width",
			"red_blood_cell_count",
			"red_blood_cell_distribution_width",
			"shbg",
			"total_bilirubin",
			"total_protein",
			"triglycerides",
			"urate",
			"urea",
			"vitamin_d",
			"white_blood_cell_count",
		]
    Map[String, Int] idxs = {
			"alanine_aminotransferase": 0,
			"albumin": 1,
			"alkaline_phosphatase": 2,
			"apolipoprotein_a": 3,
			"apolipoprotein_b": 4,
			"aspartate_aminotransferase": 5,
			"c_reactive_protein": 6,
			"calcium": 7,
			"cholesterol": 8,
			"creatinine": 9,
			"cystatin_c": 10,
			"eosinophil_count": 11,
			"eosinophil_percent": 12,
			"gamma_glutamyltransferase": 13,
			"glucose": 14,
			"glycated_haemoglobin": 15,
			"haematocrit": 16,
			"haemoglobin_concentration": 17,
			"hdl_cholesterol": 18,
			"igf_1": 19,
			"ldl_cholesterol_direct": 20,
			"lymphocyte_count": 21,
			"lymphocyte_percent": 22,
			"mean_corpuscular_haemoglobin": 23,
			"mean_corpuscular_haemoglobin_concentration": 24,
			"mean_corpuscular_volume": 25,
			"mean_platelet_volume": 26,
			"mean_sphered_cell_volume": 27,
			"neutrophil_count": 28,
			"neutrophil_percent": 29,
			"phosphate": 30,
			"platelet_count": 31,
			"platelet_crit": 32,
			"platelet_distribution_width": 33,
			"red_blood_cell_count": 34,
			"red_blood_cell_distribution_width": 35,
			"shbg": 36,
			"total_bilirubin": 37,
			"total_protein": 38,
			"triglycerides": 39,
			"urate": 40,
			"urea": 41,
			"vitamin_d": 42,
			"white_blood_cell_count": 43,
		}
    Array[String] unit = [
      "U/L",
			"g/L",
			"U/L",
			"g/L",
			"g/L",
			"U/L",
			"mg/L",
			"mmol/L",
			"mmol/L",
			"umol/L",
			"mg/L",
			"10^9 cells/L",
			"%",
			"U/L",
			"mmol/L",
			"mmol/mol",
			"%",
			"g/dL",
			"mmol/L",
			"nmol/L",
			"mmol/L",
			"10^9 cells/L",
			"%",
			"10^-12 g",
			"g/dL",
			"10^-15 L",
			"10^-15 L",
			"10^-15 L",
			"10^9 cells/L",
			"%",
			"mmol/L",
			"10^9 cells/L",
			"%",
			"%",
			"10^12 cells/L",
			"%",
			"nmol/L",
			"umol/L",
			"g/L",
			"mmol/L",
			"umol/L",
			"mmol/L",
			"nmol/L",
			"10^9 cells/L",
    ]
    Array[Int] ID = [
			30620,
			30600,
			30610,
			30630,
			30640,
			30650,
			30710,
			30680,
			30690,
			30700,
			30720,
			30150,
			30210,
			30730,
			30740,
			30750,
			30030,
			30020,
			30760,
			30770,
			30780,
			30120,
			30180,
			30050,
			30060,
			30040,
			30100,
			30270,
			30140,
			30200,
			30810,
			30080,
			30090,
			30110,
			30010,
			30070,
			30830,
			30840,
			30860,
			30870,
			30880,
			30670,
			30890,
			30000,
		]
  }

  command <<< >>>

	runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
		shortTask: true
    memory: "2GB"
  }
}

####################### Helper tasks ###########################

task concatenate_tsvs {
  input {
    # this can be a large enough array that its more than the max number of arguments
    # for a command on the command line. But its not too long for a string variable
    # or an array
    Array[File]+ tsvs
  }

  output {
    File tsv = "out.tab"
  }

  command <<<
    envsetup bash -c '
      tsvs=(~{sep=" " tsvs})
      if (( $(for (( i=0; i<${#tsvs[@]}; i++ )); do
        head -1 "${tsvs[$i]}"
      done | uniq -c | wc -l) != 1 )) ; then
        echo "Different headers" 2>&1
        exit 1
      fi
      
      head -1 ~{tsvs[0]} > out.tab
      for (( i=0; i<${#tsvs[@]}; i++ )); do
        tail -n +2 "${tsvs[$i]}" >> out.tab
      done
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "4h"
    memory: "2GB"
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
    memory: "2GB"
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
    memory: "2GB"
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
    memory: "2GB"
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
    memory: "2GB"
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
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} . ~{fam_file} ~{sc_pcs} ~{sc_assessment_ages}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    memory: "10GB"

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
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
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
    dx_timeout: "1h"
    memory: "2GB"
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
    memory: "2GB"
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
    memory: "2GB"
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
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} out.samples ~{kinship} ~{sample_list} ~{PRIMUS_command} ~{"--binary-pheno " + binary_pheno_data}
  >>>
  
  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "24h"
    memory: "4GB"
  }
}

task transform_trait_values {
  input {
    String script_dir
    File script = "~{script_dir}/traits/transform_traits.py"
    File python_array_utils = "~{script_dir}/traits/python_array_utils.py"

    File pheno_data # from task
    File samples_for_phenotype  # from task
    Boolean is_binary
  }

  output {
    File data = "out.npy"
  }

  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} out ~{pheno_data} ~{samples_for_phenotype} ~{if is_binary then "--binary" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "5m"
    memory: "2GB"
  }
}

######## STR calling and QC, SNP QC ###########

task imputed_snp_frequencies {
  input {
    String script_dir
    File script = "~{script_dir}/association/calc_macs.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    bgen imputed_snp_bgen
    File samples
    File all_samples
    region bounds
  }

  output {
    File counts = "snp_macs.tab"
  }

  command <<<
    envsetup ~{script} \
      snp_macs.tab \
      ~{imputed_snp_bgen.bgen} \
      ~{bounds.start} \
      ~{bounds.end} \
      ~{samples} \
      ~{all_samples}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "12h"
    memory: "8GB" # this is tied to the memory and threads hardcoded in the plink command
  }
}

task STR_frequencies {
  input {
    String script_dir
    File script = "~{script_dir}/side_analyses/freqs/str_allele_freqs.py"
    File sample_utils = "~{script_dir}/side_analyses/freqs/sample_utils.py"
    File python_array_utils = "~{script_dir}/side_analyses/freqs/python_array_utils.py"

    VCF vcf
    File all_samples
    File samples
  }

  output {
    File counts = stdout()
  }

  command <<<
    envsetup ~{script} \
      ~{vcf.vcf} \
      ~{samples} \
      ~{all_samples}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "47h"
    memory: "4GB" 
  }
}

task fraction_of_variation {
  input {
    String script_dir
    File script = "~{script_dir}/side_analyses/freqs/fraction_of_variation.py"

    File snp_counts
    File str_counts
    File snps_to_filter
  }

  output {
    File nums = stdout()
  }

  command <<<
    envsetup ~{script} \
      ~{str_counts} \
      ~{snp_counts} \
      ~{snps_to_filter}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "47h"
    memory: "2GB"
  }
}


# TODO could do length confusion, pre imputation allele freqs, calc macs

task imputed_str_locus_summary {
  input {
    String script_dir
    File script = "~{script_dir}/str_imputed/analyses/loci_stat_inputs.py"

    VCF vcf
    File qced_white_brits
  }

  output {
    File out = "out.tab"
  }

  command <<<
    envsetup ~{script} out.tab ~{vcf.vcf} ~{qced_white_brits}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "24h"
    memory: "5GB"
  }
}



task str_multiallelicness_distro {
  input {
    String script_dir
    File script = "~{script_dir}/str_imputed/analyses/loci_stats.py"

    Float thresh
    Array[File] chrom_locus_summaries # from imputed_str_locus_summary
  }

  output {
    File svg = "out.svg"
    File png = "out.png" 
  }

  command <<<
    envsetup ~{script} out ~{thresh} ~{sep=" " chrom_locus_summaries}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "5m"
    shortTask: true
    memory: "2GB"
  }
}

task cbl_imperfection_ld {
  input {
    String script_dir
    File script = "~{script_dir}/association/CBL_imperfection_LD.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    bgen bgen_chr11 
    File all_samples_file
    File phenotype_samples_file
  }

  output {
    Float r2 = read_float(stdout())
  }

  command <<<
    envsetup ~{script} \
      ~{bgen_chr11.bgen} \
      ~{all_samples_file} \
      ~{phenotype_samples_file}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    memory: "2GB"
  }
}

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
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} . ~{str_vcf_chr_11.vcf} ~{specific_alleles} ~{all_samples_list} \
      ~{white_brits_sample_list} ~{black_sample_list} ~{south_asian_sample_list} ~{chinese_sample_list}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "10m"
    memory: "2GB"
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
        chroms.append(chrom)
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
    memory: "2GB"
  }
}

task prep_conditional_input {
  input {
    String script_dir
    File script = "~{script_dir}/association/prep_conditional_inputs.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"
    File load_and_filter_genotypes = "~{script_dir}/association/load_and_filter_genotypes.py"
    File sample_utils = "~{script_dir}/association/sample_utils.py"

    File all_samples
    Int chrom
    VCF? str_vcf
    Array[Int] strs
    bgen? snp_bgen
    File? snp_mfi
    Array[Int] snps
  }

  output {
    File data = "out.npy"
    File varnames = "out_varnames.txt"
    File readme = "out_README.txt"
  }

  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{load_and_filter_genotypes} # necessary for dxCompiler to bother to localize this file
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      out \
      ~{all_samples} \
      ~{chrom} \
      ~{if defined(str_vcf) then "--str-vcf ~{select_first([str_vcf]).vcf}" else "" } \
      --STRs ~{sep=" " strs} \
      --imputed-SNPs ~{sep=" " snps} \
      ~{"--snp-mfi " + snp_mfi} \
      ~{if defined(snp_bgen) then "--snp-bgen ~{select_first([snp_bgen]).bgen}" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "2GB"
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
  }

  output {
    File README = "README.txt"
    File data = "out.tab"
  }

  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{weighted_binom_conf} # necessary for dxCompiler to bother to localize this file
    ls ~{load_and_filter_genotypes} # necessary for dxCompiler to bother to localize this file
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
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
      --temp-dir . \
      ~{if is_binary then "--binary logistic" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    shortTask: true
    dx_timeout: "20m"
    memory: "10GB"
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
    File? shared_covars # from task
    File? untransformed_phenotype # from task
    File transformed_phenotype # from task
    File? conditional_covars # from task
    File all_samples_list
    Boolean is_binary
    String binary_type # linear or logistic, ignored if not binary
    region? bounds
    File? vars_file
    String phenotype_name
    Boolean no_details = false
  }

  output {
    File data = "out.tab"
  }
  
  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{weighted_binom_conf} # necessary for dxCompiler to bother to localize this file
    ls ~{load_and_filter_genotypes} # necessary for dxCompiler to bother to localize this file
    ls ~{sample_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      out.tab \
      strs \
      ~{phenotype_name} \
      ~{if defined(bounds) then "--region ~{select_first([bounds]).chrom}:~{select_first([bounds]).start}-~{select_first([bounds]).end}" else ""} \
      ~{"--vars-file " + vars_file} \
      --pheno-and-covars ~{transformed_phenotype} \
      ~{"--shared-covars " + shared_covars} \
      ~{if no_details then "" else "--untransformed-phenotypes ~{select_first([untransformed_phenotype])}"} \
      ~{"--conditional-covars " + conditional_covars} \
      --all-samples-fname ~{all_samples_list} \
      --str-vcf ~{str_vcf.vcf} \
      --temp-dir . \
      ~{if is_binary then "--binary " + binary_type else ""} \
      ~{if no_details then "--no-details" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "36h"
    memory: if binary_type == "logistic" then "40GB" else "8GB" # 4g works for all but a few jobs
  }
}

task reformat_my_str_gwas_table_for_publication {
  input {
    String script_dir
    File script = "~{script_dir}/association/reformat_my_str_gwas_for_publication.py"

    String phenotype
    File my_str_gwas
    File flank_start_to_start_and_end_pos
    File str_hg19_pos_bed
    File str_hg38_pos_bed
    File repeat_units_table
  }

  output {
    File unfiltered = "out.tab"
    File filtered = "out.filtered.tab"
  }

  command <<<
    envsetup ~{script} \
      out \
      ~{phenotype} \
      ~{my_str_gwas} \
      ~{flank_start_to_start_and_end_pos} \
      ~{str_hg19_pos_bed} \
      ~{str_hg38_pos_bed} \
      ~{repeat_units_table}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    memory: "20GB"
  }
}

task locus_plot {
  input {
    String script_dir
    File script = "~{script_dir}/association/plot_locus.py"

    String phenotype_name
    # exactly one of the following two must be non empty
    # either results from my STR GWAS, or a tsv specifically formatted to work with this task
    String? chrom # only used if assoc_results
    Int? pos # only used if assoc_results
    Array[File] assoc_results = []
    Array[File] data_tsvs = []
    Array[String] group_names = []

    Int? dosage_threshold
    Float? dosage_fraction_threshold
    String? unit # should not be specified if binary
    Boolean residual = false
    Boolean binary = false
    String total_column_name = "total_subset_dosage_per_summed_gt"
  }

  output {
    File svg = "~{phenotype_name}_~{chrom}_~{pos}.svg"
    File png = "~{phenotype_name}_~{chrom}_~{pos}.png"
  }

  command <<<
    envsetup ~{script} \
      ~{phenotype_name}_~{chrom}_~{pos} \
      '~{phenotype_name}' \
      ~{"--chrom " + chrom} \
      ~{"--pos " + pos} \
      --assoc-results ~{sep=" " assoc_results} \
      --datas-per-length-sum ~{sep=" " data_tsvs} \
      --group-names '~{sep="' '" group_names}' \
      ~{"--dosage-threshold " + dosage_threshold} \
      ~{"--dosage-fraction-threshold " + dosage_fraction_threshold} \
      ~{if defined(unit) then "--unit '~{unit}'" else if binary then "--binary" else ""} \
      ~{if residual then "--residual-phenos" else ""} \
      --total-column-name ~{total_column_name}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30:00"
    memory: "2GB"
  }
}

# for getting Geuvadis and GTEx data into a format that can be passed to locus plot
task summarize_individual_data_for_plotting {
  input {
    String script_dir
    File script = "~{script_dir}/association/summarize_individual_data_for_plotting.py"

    File individual_tsv
    String length_sum_column_name
    String trait_column_name
    String? sep
    String? subset_field
    String? subset_value
    Int? ratio
    Float? offset
  }

  output {
    File out = "out.tsv"
  }

  command <<<
    envsetup ~{script} out.tsv ~{individual_tsv} ~{length_sum_column_name} '~{trait_column_name}' \
      ~{"--sep " + sep} \
      ~{"--subset " + subset_field + " " + subset_value} \
      ~{"--ratio " + ratio} \
      ~{"--offset " + offset} \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
    memory: "2GB"
  }
}

# TODO could do my GWAS for SNPs

task prep_plink_input {
  input {
    String script_dir
    File script = "~{script_dir}/association/prep_plink_input.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    File? shared_covars # from task
    File? shared_covar_names # from task
    File transformed_phenotype # from task
    File pheno_covar_names # from task
    File? conditional_genotypes # npy from task prep conditional input
    File? conditional_covar_names # varnames from task above
    Boolean is_binary
    String binary_type # linear or logistic, only needs to be set if is_binary == true
    String phenotype_name
  }

  output {
    File data = "out.tab"
  }

  command <<<
    ls ~{python_array_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} \
      out.tab \
      ~{phenotype_name} \
      ~{transformed_phenotype} \
      ~{pheno_covar_names} \
      ~{shared_covars} \
      ~{shared_covar_names} \
      ~{"--conditional-genotypes " + conditional_genotypes} \
      ~{"--conditional-covar-names " + conditional_covar_names} \
      ~{if is_binary then "--binary " + binary_type else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "2GB"
  }
}

task plink_snp_association {
  input {
    String script_dir
    File script = "~{script_dir}/association/plink_association.sh"
    String plink_command

    PFiles imputed_snp_p_file # TODO could generate this with task
    File pheno_data # from task

    Int chrom
    Int? start
    Int? end
    String phenotype_name
    String binary_type # linear or linear_binary or logistic
  }

  output {
    File data = "plink2." + (if binary_type == "linear" then "rin_" else "") + phenotype_name + ".glm." + (if binary_type == "logistic" then "logistic.hybrid" else "linear")
    File log = "plink2.log"
  }

  command <<<
    PHENOTYPE=~{phenotype_name} \
    BINARY_TYPE=~{binary_type} \
    CHROM=~{chrom} \
    OUT_DIR=. \
    PHENO_FILE=~{pheno_data} \
    P_FILE=$(echo ~{imputed_snp_p_file.pgen} | sed -e 's/\.pgen$//') \
    PLINK_EXECUTABLE=~{plink_command} \
    PROJECT_TEMP=. \
    ~{"START=" + start} \
    ~{"END=" + end} \
    envsetup ~{script}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    memory: "56GB"
    cpus: 28
    dx_timeout: "24h"
  }
}

## TODO append mfi to plink logistic run
## TODO compare my to plink GWAS

task manhattan {
  input {
    String script_dir
    File script = "~{script_dir}/association/interactive_manhattan_plot.py"
    File region_plots = "~{script_dir}/association/region_plots.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    String phenotype_name
    String unit
    File chr_lens
    File str_gwas_results
    File snp_gwas_results
    region? bounds
    Array[Int] conditioned_STRs = []
    Array[String] conditioned_imputed_SNPs = []
    File? conditional_str_results
    File? conditional_snp_results
    File? peaks
    String ext
    Boolean use_tstat = false
  }

  output {
    File plot = "manhattan.~{ext}"
  }

  command <<<
    envsetup ~{script} \
      "manhattan.~{ext}" \
      ~{phenotype_name} \
      '~{unit}' \
      ~{chr_lens} \
      ~{str_gwas_results} \
      ~{snp_gwas_results} \
      ~{ext} \
      ~{if defined(bounds) then "--chrom ~{select_first([bounds]).chrom}" else ""} \
      ~{if defined(bounds) then "--start ~{select_first([bounds]).start}" else ""} \
      ~{if defined(bounds) then "--end ~{select_first([bounds]).end}" else ""} \
      --conditioned-STRs ~{sep=" " conditioned_STRs} \
      --conditioned-imputed-SNPs ~{sep=" " conditioned_imputed_SNPs} \
      ~{"--conditional-STR-results " + conditional_str_results} \
      ~{"--conditional-imputed-SNP-results " + conditional_snp_results} \
      ~{"--peaks-fname" + peaks} \
      ~{if use_tstat then "--use-tstat " else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "50GB"
  }
}

task overview_manhattan {
  input {
    String script_dir
    File script = "~{script_dir}/association/overview_manhattan_plot.py"
    File region_plots = "~{script_dir}/association/region_plots.py"
    File graphing_utils = "~{script_dir}/association/graphing_utils.py"
    File python_array_utils = "~{script_dir}/association/python_array_utils.py"

    String phenotype_name
    File chr_lens
    File str_gwas_results
    File snp_gwas_results
    File peaks
    String ext
  }

  output {
    File plot = "manhattan.~{ext}"
  }

  command <<<
    envsetup ~{script} \
      ~{phenotype_name} \
      "manhattan.~{ext}" \
      ~{str_gwas_results} \
      ~{snp_gwas_results} \
      ~{peaks} \
      ~{chr_lens} \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "75GB"
  }
}

task compare_bili_to_UKBB {
  input {
    String script_dir
    File script = "~{script_dir}/association/compare_to_UKBB.py"
    
    File pan_ukbb
    File plink_snp
  }

  output {
    File svg = "panukbb_scatter.svg"
    File png = "panukbb_scatter.png"
  }

  command <<<
    envsetup ~{script} \
      . \
      ~{pan_ukbb} \
      ~{plink_snp} 
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "50GB"
    shortTask: true
  }
}

task generate_peaks {
  input {
    String script_dir
    File script = "~{script_dir}/signals/peaks.py"

    File str_assoc_results
    File snp_assoc_results
    String phenotype
    String spacing
    String thresh
  }

  output {
    # originally located in signals/peaks/{phenotype}_{spacing}_{thresh}.tab
    File peaks = "peaks.tab"
    File readme = "readme.txt"
  }

  command <<<
    envsetup ~{script} \
      ~{phenotype} \
      ~{spacing} \
      ~{thresh} \
      ~{str_assoc_results} \
      ~{snp_assoc_results} \
      readme.txt peaks.tab
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "4h"
    memory: "20GB"
  }
}

task summarize_peaks {
  input {
    String script_dir
    File script = "~{script_dir}/signals/peak_summary_images.py"
    File graphing_utils = "~{script_dir}/signals/graphing_utils.py"
    File phenotype_utils = "~{script_dir}/signals/phenotypes.py"

    Array[String] phenotype_names
    Array[File] peak_files
  }

  output {
    # originally located in export_scripts/results/
    File barplot_svg = "peaks_by_pheno.svg"
    File barplot_png = "peaks_by_pheno.png"
    File heatmap_svg = "peak_p_val_heatmap.svg"
    File heatmap_png = "peak_p_val_heatmap.png"
  }

  command <<<
    ls ~{graphing_utils} # necessary for dxCompiler to bother to localize this file
    ls ~{phenotype_utils} # necessary for dxCompiler to bother to localize this file
    envsetup ~{script} . --phenotypes ~{sep=" " phenotype_names} --peak-files ~{sep=" " peak_files}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "1h"
    shortTask: true
    memory: "50GB"
  }
}

####################### Tasks to determine finemapping regions #########################

# these are here because they are used to prioritize STRs for follow up ethnic GWAS

task generate_finemapping_regions {
  input {
    String script_dir
    File script = "~{script_dir}/signals/regions.py"

    File chr_lens # misc_data/genome/chr_lens.txt
    String phenotype
    File str_assoc_results
    File snp_assoc_results
    Boolean remove_skips = false
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
      out.tab \
      ~{if remove_skips then "--remove-skips" else ""}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "30m"
    memory: "10GB"
  }
}

task get_strs_in_finemapping_regions {
  input {
    String script_dir
    File script = "~{script_dir}/signals/write_finemapped_strs.py"

    File str_loci
    File finemapping_regions_for_pheno # from generate regions
  }

  output {
    Array[File] strs_in_finemapping_regions = [
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
    memory: "4GB"
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


