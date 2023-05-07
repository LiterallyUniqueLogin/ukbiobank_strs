version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "../gwas_wdl/prep_samples_and_phenotype_workflow.wdl"

#struct locus {
#  Int chrom
#  Int pos
#}

task any_results_for_phenotype {
  input {
    String phenotype
    File results_tsv
  }

  command <<<
    envsetup python -c '
import polars as pl
df = pl.read_csv("~{results_tsv}", sep="\t")
df = df.filter(
  pl.col("phenotype") == "~{phenotype}"
)
print(df.shape[0] > 0)
'
  >>>

  output {
    Boolean any = read_boolean(stdout())
  }

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
  }
}

#print(
#  df.filter(
#    pl.col("phenotype") == "~{phenotype}"
#  ).select([
#    "chrom",
#    pl.col("start_pos (hg19)").alias("pos")
#  ]).write_csv(sep="\t")
#)

task results_tsv_to_loci {
  input {
    String phenotype
    File results_tsv
  }

  command <<<
    envsetup python -c '
import polars as pl
df = pl.read_csv(
  "~{results_tsv}", sep="\t"
).filter(
  pl.col("phenotype") == "~{phenotype}"
)

print("[", end="")
start = True
for chrom, pos in list(zip(df["chrom"], df["start_pos (hg19)"])):
  if not start:
    print(", ", end="")
  start = False
  print(f"[{chrom}, {pos}]")
print("]", end="")
'
  >>>

  output {
    Array[Array[Int]] loci = read_json(stdout())
    #Array[locus] loci = read_objects(stdout())
  }

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "10m"
    shortTask: true
  }
}

# same as wgs_tasks.regional_my_str_gwas except calls associaTR instead of the local code
# the output p-values and betas aren't relevant for this - this task doesn't bother to include
# transforming the phenotype or adding covariates
# the reliable bit should be the plotting association statistics
task associaTR_my_str_gwas {
  input {
    VCF str_vcf
    File untransformed_phenotype
    Int chrom
    Int pos
    String phenotype_name
  }

  output {
    File data = "out.tab"
  }

  command <<<
    envsetup associaTR \
      out.tab \
      ~{str_vcf.vcf} \
      ~{phenotype_name} \
      ~{untransformed_phenotype} \
      --region ~{chrom}:~{pos}-~{pos} \
      --stats-for-plotting ~{untransformed_phenotype} \
      --plotting-ci-alphas 0.05
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ukb_strs:v1.3"
    dx_timeout: "2h"
    memory: "8GB"
  }
}

workflow plot_wgs_assocs_for_finemapped_loci {

  input {
    String script_dir  = "dx://UKB_Test:/imputed_strs_paper"

    VCF wgs_vcf = {
      "vcf": "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz",
      "index": "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_full_wgs_calls.vcf.gz.tbi"
    }

    File fam_file = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb22418_c1_b0_v2.fam"
    File withdrawn_sample_list = "dx://UKB_Test:/imputed_strs_paper/sample_qc/withdrawn.sample"
    File kinship = "dx://UKB_Test:/Bulk/Genotype%20Results/Genotype%20calls/ukb_rel.dat"

    String sc_data_dir = "dx://UKB_Test:/imputed_strs_paper/main_dataset/extracted_data"
    File sc_white_brits = "~{sc_data_dir}/22006.tsv"
    File sc_ethnicity_self_report = "~{sc_data_dir}/21000.tsv"
    File sc_sex_aneuploidy = "~{sc_data_dir}/22019.tsv"
    File sc_genetic_sex = "~{sc_data_dir}/22001.tsv"
    File sc_reported_sex = "~{sc_data_dir}/31.tsv"
    File sc_kinship_count = "~{sc_data_dir}/22021.tsv"
    File sc_assessment_ages = "~{sc_data_dir}/21003.tsv"
    File sc_pcs ="~{sc_data_dir}/22009.tsv"
    File sc_year_of_birth = "~{sc_data_dir}/34.tsv"
    File sc_month_of_birth = "~{sc_data_dir}/52.tsv"
    File sc_date_of_death = "~{sc_data_dir}/40000.tsv"

    File results_tsv =  "dx://UKB_Test:/imputed_strs_paper/wgs/20230425_singly_finemapped_strs_for_paper.tab"
    File called_wgs_sample_list = "dx://UKB_Test:/imputed_strs_paper/wgs/wgs_samples.txt"
  }

  call gwas_tasks.phenotype_names

  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    String phenotype = phenotype_names.n[phenotype_idx]
    String unit = phenotype_names.unit[phenotype_idx]

    call any_results_for_phenotype { input :
      phenotype = phenotype,
      results_tsv = results_tsv
    }

    if (any_results_for_phenotype.any) {

      call results_tsv_to_loci { input :
        phenotype = phenotype,
        results_tsv = results_tsv
      }

      call prep_samples_and_phenotype_workflow.prep_samples_and_phenotype { input :
        script_dir = script_dir,
        PRIMUS_command = "run_PRIMUS.pl",

        phenotype_name = phenotype,
        categorical_covariate_names = [], # don't actually need these because we're ignoring the regression, just looking at the raw association
        categorical_covariate_scs = [],
        is_binary = false,
        is_zero_one_neg_nan = false,
        date_of_most_recent_first_occurrence_update = "",

        fam_file = fam_file,
        withdrawn_sample_list = withdrawn_sample_list,
        kinship = kinship,

        sc_white_brits = sc_white_brits,
        sc_ethnicity_self_report = sc_ethnicity_self_report,
        sc_sex_aneuploidy = sc_sex_aneuploidy,
        sc_genetic_sex = sc_genetic_sex,
        sc_reported_sex = sc_reported_sex,
        sc_kinship_count = sc_kinship_count,
        sc_assessment_ages = sc_assessment_ages,
        sc_pcs = sc_pcs,
        sc_year_of_birth = sc_year_of_birth,
        sc_month_of_birth = sc_month_of_birth,
        sc_date_of_death = sc_date_of_death,
        sc_phenotype = "~{sc_data_dir}/~{phenotype_names.ID[phenotype_idx]}.tsv",

        subpop_sample_list = called_wgs_sample_list,
      }

      scatter (locus in results_tsv_to_loci.loci) {
        call associaTR_my_str_gwas { input :
          str_vcf = wgs_vcf,
          untransformed_phenotype = prep_samples_and_phenotype.pheno_data[0],
          chrom = locus[0], #locus.chrom,
          pos = locus[1], #locus.pos,
          phenotype_name = phenotype
        }

        call gwas_tasks.locus_plot { input :
          script_dir = script_dir,
          chrom = locus[0], #locus.chrom,
          pos = locus[1], #locus.pos,
          phenotype_name = phenotype,
          assoc_results = [associaTR_my_str_gwas.data],
          dosage_fraction_threshold = 0.001,
          unit = unit,
          total_column_name = "sample_count_per_summed_length"
        }
      } # scatter over loci
    } # any results for pheno
  } # scatter over pheno

  output {
    Array[File] pngs = flatten(select_all(locus_plot.png))
    Array[File] svgs = flatten(select_all(locus_plot.svg))
  }
}
