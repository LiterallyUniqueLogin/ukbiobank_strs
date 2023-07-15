version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "expanse_files.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"

workflow finemapping_files {

  call gwas_tasks.phenotype_names
  call expanse_files.files

  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    File finemapping_region = files.finemapping_regions[phenotype_idx]

    Array[Array[String]] finemapping_regions_tsv = read_tsv(finemapping_region)

    # first pass finemapping results
    if (phenotype_names.n[phenotype_idx] != "ldl_cholesterol_direct" && phenotype_names.n[phenotype_idx] != "total_bilirubin") {
      scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
        Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
        region first_pass_bounds = {
          "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
          "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
          "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
        }

        if (!(
          phenotype_names.n[phenotype_idx] == "urate" &&
          first_pass_bounds.chrom == 4 &&
          first_pass_bounds.start == 8165642 &&
          first_pass_bounds.end == 11717761
        ) && !(
          phenotype_names.n[phenotype_idx] == "total_bilirubin" &&
          first_pass_bounds.chrom == 12 &&
          first_pass_bounds.start == 19976272 &&
          first_pass_bounds.end == 22524428
        ) && !(
          phenotype_names.n[phenotype_idx] == "alkaline_phosphatase" &&
          first_pass_bounds.chrom == 1 &&
          first_pass_bounds.start == 19430673 &&
          first_pass_bounds.end == 24309348
        )) {
          String original_susie_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/susie_results/~{phenotype_names.n[phenotype_idx]}/~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
          serializable_SuSiE_output original_susie_snakemake__ = object {
            lbf: "~{original_susie_dir}/lbf.tab",
            lbf_variable: "~{original_susie_dir}/lbf_variable.tab",
            sigma2: "~{original_susie_dir}/sigma2.txt",
            V: "~{original_susie_dir}/V.tab",
            converged: "~{original_susie_dir}/converged.txt",
            lfsr: "~{original_susie_dir}/lfsr.tab",
            requested_coverage: "~{original_susie_dir}/requested_coverage.txt",
            alpha: "~{original_susie_dir}/alpha.tab",
            colnames: "~{original_susie_dir}/colnames.txt.normal_run",
          }
          Array[String] original_susie_CSs_snakemake__ = read_lines("~{original_susie_dir}/cs_files_list.txt")

#        if (phenotype_names.n[phenotype_idx] == "creatinine") {
#          serializable_FINEMAP_output original_finemap_ = object {
#            # TODO
#          }
#          Array[String] original_finemap_creds_ = read_lines("")
#        }
        }
      }
      Array[serializable_SuSiE_output] original_susie_snakemake_ = select_all(original_susie_snakemake__)
      Array[Array[String]] original_susie_CSs_snakemake_ = select_all(original_susie_CSs_snakemake__)
    }
  }

  output {
    # indexed over snakemake-phenotype, then region
    Array[Array[serializable_SuSiE_output]] original_susie_snakemake = select_all(original_susie_snakemake_)
    # indexed over snakemake-phenotype, then region, then particular CS
    Array[Array[Array[String]]] original_susie_CSs_snakemake = select_all(original_susie_CSs_snakemake_)
  }
}
