version 1.0

import "../platform_wdl/expanse_tasks.wdl"
import "../gwas_wdl/gwas_tasks.wdl"
import "../finemapping_wdl/finemapping_tasks.wdl"
import "../finemapping_wdl/post_finemapping_workflow.wdl"

# also includes tables
workflow expanse_figures {

  input {
    String script_dir = "."
    File all_samples_list = "microarray/ukb46122_hap_chr1_v2_s487314.sample"
    File fam_file = "microarray/ukb46122_cal_chr1_v2_s488176.fam"

    File chr_lens = "misc_data/genome/chr_lens.txt"

    # TODO this group of files could be generated through WDL instead of cached
    File str_pos_table = "snpstr/flank_trimmed_vcf/vars.tab"
    File repeat_units_table = "snpstr/repeat_units.tab"
    File str_pos_table_2 = "snpstr/str_loci.bed"
    File str_hg38_pos_table = "snpstr/str_loci.hg38.bed"
    File str_t2t_pos_table = "snpstr/str_loci.t2tv2.bed"

    File specific_alleles = "association/specific_alleles.tab"

    File CBL_gtex_expression = "misc_data/gtex_yang/CBL_chr11_119206290_GTEX_TPM.tsv"
    File Liver_SLC2A2_exon4_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon4_psi.tsv"
    File Liver_SLC2A2_exon6_psi= "misc_data/gtex_yang/Liver_SLC2A2_exon6_psi.tsv"

    # sample files which include randomness, so need to be reused for consistency
    File unrelated_white_brits_sample_list = "sample_qc/runs/white_brits/no_phenotype/combined_unrelated.sample"
    File unrelated_black_sample_list = "sample_qc/runs/black/no_phenotype/combined_unrelated.sample"
    File unrelated_south_asian_sample_list = "sample_qc/runs/south_asian/no_phenotype/combined_unrelated.sample"
    File unrelated_chinese_sample_list = "sample_qc/runs/chinese/no_phenotype/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_not_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File unrelated_samples_CBL_hom_begin_C_T_snp = "sample_qc/subpop_runs/CBL_hom_begin_C_T_snp/white_brits/platelet_count/combined_unrelated.sample"
    File platelet_count_sample_list = "sample_qc/runs/white_brits/platelet_count/combined_unrelated.sample"

    # other cached inputs are found directly in the workflow because they need scatters or depend on other logic
    # >_<
  }

  # inputs, would like to be in the inputs block but are arrays and so need scatter statements to make
  scatter (chrom in range(22)) {
    VCF str_vcfs = {
      "vcf": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz",
      "index": "str_imputed/runs/first_pass/vcfs/annotated_strs/chr~{chrom+1}.vcf.gz.tbi"
    }
    PFiles imputed_snp_p_files = {
      "pgen": "array_imputed/pfile_converted/chr{chrom+1}.pgen",
      "pvar": "array_imputed/pfile_converted/chr{chrom+1}.pvar",
      "psam": "array_imputed/pfile_converted/chr{chrom+1}.psam",
    }
  }

  # input, but written as a task for reusability
  call gwas_tasks.phenotype_names

  ##### some common task inputs
  call expanse_tasks.extract_field as sc_white_brits { input :
    script_dir = script_dir,
    id = 22006
  }

  call gwas_tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = script_dir,
    sc = sc_white_brits.data
  }

  scatter (phenotype in phenotype_names.n) {
    # GWAS results, cached for efficiency
    File str_gwas_results = "association/results/~{phenotype}/my_str/results.tab"
    File snp_gwas_results = "association/results/~{phenotype}/plink_snp/results.tab"
    scatter (ethnicity in ["black", "south_asian", "chinese", "irish", "white_other"]) {
      File ethnic_str_gwas_resultss = "association/results_finemapped_only/~{ethnicity}//~{phenotype}/my_str/results.tab"
    }
  }

  ##### generate figure 1 (unfinished)
#  # TODO this crashes because of a missing header in the VCF ##command=Hipstr
#  scatter (chrom in range(22)) {
#    call gwas_tasks.imputed_str_locus_summary as imputed_str_locus_summaries { input :
#      script_dir = script_dir,
#      vcf = str_vcfs[chrom],
#      qced_white_brits = unrelated_white_brits_sample_list
#    }
#  }
#
#  call gwas_tasks.str_multiallelicness_distro as fig_1b { input:
#    script_dir = script_dir,
#    thresh = 0.01,
#    chrom_locus_summaries = imputed_str_locus_summaries.out
#  }

#  # cached signal peaks - TODO generate these?
#  scatter (phenotype in phenotype_names.n) {
#    String peak_files = "signals/peaks/~{phenotype}_250000_5e-8.tab"
#  }

  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    call gwas_tasks.generate_peaks { input :
      script_dir = script_dir,
      snp_assoc_results = snp_gwas_results[phenotype_idx],
      str_assoc_results = str_gwas_results[phenotype_idx],
      phenotype = phenotype_names.n[phenotype_idx],
      spacing = "250000",
      thresh = "5e-8"
    }
  }

  call gwas_tasks.summarize_peaks as fig_1ef { input :
    script_dir = script_dir,
    phenotype_names = phenotype_names.n,
    peak_files = generate_peaks.peaks
  }

  ##### generate fine-mapping supplementary figures (3-12, missing 8 and 9)
  # take cached fine-mapping results and summarize them
  scatter (phenotype_idx in range(length(phenotype_names.n))) {
    call gwas_tasks.generate_finemapping_regions { input :
      script_dir = script_dir,
      chr_lens = chr_lens,
      phenotype = phenotype_names.n[phenotype_idx],
      snp_assoc_results = snp_gwas_results[phenotype_idx],
      str_assoc_results = str_gwas_results[phenotype_idx],
      remove_skips = true
    }

    Array[Array[String]] finemapping_regions_tsv = read_tsv(generate_finemapping_regions.data)

    # first pass finemapping results, cached for efficiency
    scatter (first_pass_region_idx in range(length(finemapping_regions_tsv) - 1)) {
      Int first_pass_region_idx_plus_one = first_pass_region_idx + 1
      region first_pass_bounds = {
        "chrom": finemapping_regions_tsv[first_pass_region_idx_plus_one][0],
        "start": finemapping_regions_tsv[first_pass_region_idx_plus_one][1],
        "end": finemapping_regions_tsv[first_pass_region_idx_plus_one][2],
      }
      String first_pass_region_strs = "~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
      Int first_pass_chroms = first_pass_bounds.chrom

      String original_finemap_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/finemap_results/~{phenotype_names.n[phenotype_idx]}/~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
      serializable_FINEMAP_output original_finemap = object {
        snp_file: "~{original_finemap_dir}/finemap_output.snp",
        log_sss: "~{original_finemap_dir}/finemap_output.log_sss",
        config: "~{original_finemap_dir}/finemap_output.config",
      }
      Array[String] original_finemap_creds = read_lines("~{original_finemap_dir}/cred_files_list.txt")

      String original_susie_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/susie_results/~{phenotype_names.n[phenotype_idx]}/~{first_pass_bounds.chrom}_~{first_pass_bounds.start}_~{first_pass_bounds.end}"
      serializable_SuSiE_output original_susie = object {
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
      Array[String] original_susie_CSs = read_lines("~{original_susie_dir}/cs_files_list.txt")
    }

    call finemapping_tasks.first_pass_finemapping_df { input :
      script_dir = script_dir,
      phenotype_name = phenotype_names.n[phenotype_idx],
      snp_assoc_results = snp_gwas_results[phenotype_idx],
      str_assoc_results = str_gwas_results[phenotype_idx],
      ethnic_str_assoc_results = ethnic_str_gwas_resultss[phenotype_idx],
      original_finemap_outputs = original_finemap,
      original_finemap_creds = original_finemap_creds,
      original_susie_outputs = original_susie,
      original_susie_CSs = original_susie_CSs,
      regions = first_pass_region_strs,
      chroms = first_pass_chroms,
    }

    call finemapping_tasks.generate_followup_regions_tsv { input :
      script_dir = script_dir,
      first_pass_df = first_pass_finemapping_df.all_regions_concordance
    }

    Array[Array[String]] followup_finemapping_regions_tsv = read_tsv(generate_followup_regions_tsv.tsv)

    # followup finemapping results, cached for efficiency
    if (length(followup_finemapping_regions_tsv) > 1) {
      scatter (finemap_run in [
        "finemap_results/~{phenotype_names.n[phenotype_idx]}.total_prob_4",
        "finemap_results/~{phenotype_names.n[phenotype_idx]}.prior_std_0.0224",
        "finemap_results/~{phenotype_names.n[phenotype_idx]}.prob_conv_sss_tol_0.0001",
        "finemap_results_mac_100/~{phenotype_names.n[phenotype_idx]}",
        "finemap_results_threshold_0.0005/~{phenotype_names.n[phenotype_idx]}",
        "finemap_results/~{phenotype_names.n[phenotype_idx]}.prior_std_0.005",
        "finemap_results/~{phenotype_names.n[phenotype_idx]}.snp_str_ratio_4",
      ]) {
        scatter (followup_region_idx in range(length(followup_finemapping_regions_tsv) - 1)) {
          Int followup_region_idx_plus_one = followup_region_idx + 1
          region followup_bounds = {
            "chrom": followup_finemapping_regions_tsv[followup_region_idx_plus_one][1],
            "start": sub(sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one][2], "^[^_]*_", ""), "_[^_]*$", ""),
            "end": sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one][2], "^[^_]*_[^_]*_", ""),
          }
          String followup_region_strs = followup_finemapping_regions_tsv[followup_region_idx_plus_one][2]
          Int followup_chroms = followup_bounds.chrom
          
          String finemap_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/~{finemap_run}/~{followup_bounds.chrom}_~{followup_bounds.start}_~{followup_bounds.end}/"
          serializable_FINEMAP_output followup_finemaps = object {
            snp_file: "~{finemap_dir}/finemap_output.snp",
            log_sss: "~{finemap_dir}/finemap_output.log_sss",
            config: "~{finemap_dir}/finemap_output.config",
          }
          Array[String] followup_finemap_creds = read_lines("~{finemap_dir}/cred_files_list.txt")
        }
      } # end scatter over finemapping runs

      scatter (followup_region_idx in range(length(followup_finemapping_regions_tsv) - 1)) {
        Int followup_region_idx_plus_one_ = followup_region_idx + 1
        region followup_bounds_ = {
          "chrom": followup_finemapping_regions_tsv[followup_region_idx_plus_one_][1],
          "start": sub(sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one_][2], "^[^_]*_", ""), "_[^_]*$", ""),
          "end": sub(followup_finemapping_regions_tsv[followup_region_idx_plus_one_][2], "^[^_]*_[^_]*_", ""),
        }

        String susie_best_guess_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/susie_hardcall_results/~{phenotype_names.n[phenotype_idx]}/~{followup_bounds_.chrom}_~{followup_bounds_.start}_~{followup_bounds_.end}/"
        serializable_SuSiE_output best_guess_susies = object {
          lbf: "~{susie_best_guess_dir}/lbf.tab",
          lbf_variable: "~{susie_best_guess_dir}/lbf_variable.tab",
          sigma2: "~{susie_best_guess_dir}/sigma2.txt",
          V: "~{susie_best_guess_dir}/V.tab",
          converged: "~{susie_best_guess_dir}/converged.txt",
          lfsr: "~{susie_best_guess_dir}/lfsr.tab",
          requested_coverage: "~{susie_best_guess_dir}/requested_coverage.txt",
          alpha: "~{susie_best_guess_dir}/alpha.tab",
          colnames: "~{susie_best_guess_dir}/colnames.txt",
        }
        Array[String] best_guess_susie_CSs = read_lines("~{susie_best_guess_dir}/cs_files_list.txt")

        String susie_ratio_dir = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/susie_results/~{phenotype_names.n[phenotype_idx]}_snp_str_ratio_4/~{followup_bounds_.chrom}_~{followup_bounds_.start}_~{followup_bounds_.end}/"
        String original_susie_dir_for_ratio = "/expanse/projects/gymreklab/jmargoli/ukbiobank/finemapping/susie_results/~{phenotype_names.n[phenotype_idx]}/~{followup_bounds_.chrom}_~{followup_bounds_.start}_~{followup_bounds_.end}/"
        serializable_SuSiE_output ratio_susies = object {
          lbf: "~{susie_ratio_dir}/lbf.tab",
          lbf_variable: "~{susie_ratio_dir}/lbf_variable.tab",
          sigma2: "~{susie_ratio_dir}/sigma2.txt",
          V: "~{susie_ratio_dir}/V.tab",
          converged: "~{susie_ratio_dir}/converged.txt",
          lfsr: "~{susie_ratio_dir}/lfsr.tab",
          requested_coverage: "~{susie_ratio_dir}/requested_coverage.txt",
          alpha: "~{susie_ratio_dir}/alpha.tab",
          colnames: if phenotype_names.n[phenotype_idx] != "mean_platelet_volume" then "~{original_susie_dir_for_ratio}/colnames.txt" else "~{susie_ratio_dir}/colnames.txt",
        }
        Array[String] ratio_susie_CSs = read_lines("~{susie_ratio_dir}/cs_files_list.txt")
      }

      call finemapping_tasks.followup_finemapping_conditions_df { input :
        script_dir = script_dir,
        phenotype_name = phenotype_names.n[phenotype_idx],
        snp_assoc_results = snp_gwas_results[phenotype_idx],
        str_assoc_results = str_gwas_results[phenotype_idx],
        ethnic_str_assoc_results = ethnic_str_gwas_resultss[phenotype_idx],
        original_finemap_outputs = original_finemap,
        original_finemap_creds = original_finemap_creds,
        original_susie_outputs = original_susie,
        original_susie_CSs = original_susie_CSs,
        total_prob_finemap_outputs = followup_finemaps[0],
        total_prob_finemap_creds = followup_finemap_creds[0],
        derived_prior_std_finemap_outputs = followup_finemaps[1], 
        derived_prior_std_finemap_creds = followup_finemap_creds[1], 
        conv_tol_finemap_outputs = followup_finemaps[2],
        conv_tol_finemap_creds = followup_finemap_creds[2],
        mac_finemap_outputs = followup_finemaps[3],
        mac_finemap_creds = followup_finemap_creds[3],
        threshold_finemap_outputs = followup_finemaps[4],
        threshold_finemap_creds = followup_finemap_creds[4],
        best_guess_susie_outputs = best_guess_susies,
        best_guess_susie_CSs = best_guess_susie_CSs,
        low_prior_std_finemap_outputs = followup_finemaps[5],
        low_prior_std_finemap_creds = followup_finemap_creds[5],
        ratio_finemap_outputs = followup_finemaps[6],
        ratio_finemap_creds = followup_finemap_creds[6],
        ratio_susie_outputs = ratio_susies,
        ratio_susie_CSs = ratio_susie_CSs,
        original_regions = first_pass_region_strs,
        original_chroms = first_pass_chroms,
        followup_regions = followup_region_strs[0], # any index would do, all identical
        followup_chroms = followup_chroms[0],
      }
    } # end if followup finemapping
  } # end scatter over phenotype

  # cached results of fine-mapping analyses
#  scatter (phenotype in phenotype_names.n) {
#    File first_pass_dfs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/finemapping_all_concordance_~{phenotype}.tab"
#    File susie_all_min_abs_corrs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/susie_all_min_abs_corrs_~{phenotype}.npy"
#    File followup_dfs = "/expanse/projects/gymreklab/jmargoli/ukbiobank/post_finemapping/intermediate_results/finemapping_putatively_causal_concordance_~{phenotype}.tab"
#  }

  call post_finemapping_workflow.post_finemapping { input :
    script_dir = ".",
    first_pass_dfs = first_pass_finemapping_df.all_regions_concordance,
    susie_all_min_abs_corrs = first_pass_finemapping_df.susie_all_regions_min_abs_corrs,
    followup_dfs = select_all(followup_finemapping_conditions_df.df)
  }

  ####### generate figure 4 (missing manhattans)
  call expanse_tasks.extract_field as pcs { input :
    script_dir = script_dir,
    id = 22009
  }

  call expanse_tasks.extract_field as assessment_ages { input :
    script_dir = script_dir,
    id = 21003
  }

  call gwas_tasks.load_shared_covars { input:
    script_dir = script_dir,
    fam_file = fam_file,
    sc_pcs = pcs.data,
    sc_assessment_ages = assessment_ages.data
  }

  call gwas_tasks.fig_4a { input :
    script_dir = script_dir,
    all_samples_list = all_samples_list,
    white_brits_sample_list = white_brits_sample_list.data,
    black_sample_list = unrelated_black_sample_list,
    south_asian_sample_list = unrelated_south_asian_sample_list,
    chinese_sample_list = unrelated_chinese_sample_list,
    str_vcf_chr_11 = str_vcfs[10],
    specific_alleles = specific_alleles,
  }

  call expanse_tasks.extract_field as platelet_count_sc { input :
    script_dir = script_dir,
    id = 30080,
  }
  
  call expanse_tasks.extract_field as platelet_count_covariate_sc { input :
    script_dir = script_dir,
    id = 30083,
  }

  # all, not qced or subset to unrelated, the sample list for this won't be used, only the data
  call gwas_tasks.load_continuous_phenotype as platelet_count_all_white_brits { input:
    script_dir = script_dir,
    sc = platelet_count_sc.data,
    qced_sample_list = white_brits_sample_list.data,
    assessment_ages_npy = load_shared_covars.assessment_ages,
    categorical_covariate_names = ["platelet_count_device_id"],
    categorical_covariate_scs = [platelet_count_covariate_sc.data],
  }

  call gwas_tasks.transform_trait_values as transformed_platelet_count { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = platelet_count_sample_list,
    is_binary = false
  }

  call gwas_tasks.str_spot_test as CBL_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data,
    transformed_phenotype = transformed_platelet_count.data,
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.locus_plot as fig_4b { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [CBL_assoc.data],
    dosage_fraction_threshold = 0.001
  }

  # use precomputed sample lists
  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_not_SNP_samples { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = unrelated_samples_CBL_hom_not_begin_C_T_snp,
    is_binary = false
  }

  call gwas_tasks.str_spot_test as CBL_hom_not_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data,
    transformed_phenotype = platelet_count_transformed_CBL_hom_not_SNP_samples.data, 
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.transform_trait_values as platelet_count_transformed_CBL_hom_SNP_samples { input:
    script_dir = script_dir,
    pheno_data = platelet_count_all_white_brits.data,
    samples_for_phenotype = unrelated_samples_CBL_hom_begin_C_T_snp,
    is_binary = false
  }
   
  call gwas_tasks.str_spot_test as CBL_hom_SNP_assoc { input:
    script_dir = script_dir,
    str_vcf = str_vcfs[10],
    shared_covars = load_shared_covars.shared_covars, 
    untransformed_phenotype = platelet_count_all_white_brits.data, 
    transformed_phenotype = platelet_count_transformed_CBL_hom_SNP_samples.data,
    all_samples_list = all_samples_list,
    is_binary = false,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
  }

  call gwas_tasks.locus_plot as fig_4f { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "platelet_count",
    unit = "10^9 cells/L",
    assoc_results = [CBL_hom_not_SNP_assoc.data, CBL_hom_SNP_assoc.data],
    group_names = ["homozygous (CGG)n", "homozygous CGGTGG(CGG)m"],
    dosage_threshold = 200,
  }

  ### plot with data from GTEx Yang

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_CBL_gtex_expression { input :
    script_dir = script_dir,
    individual_tsv = CBL_gtex_expression,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "TPM(expression)"
  }

  call gwas_tasks.locus_plot as fig_4g { input:
    script_dir = script_dir,
    chrom = 11,
    pos = 119077000,
    phenotype_name = "CBL expression",
    unit = "TPM",
    data_tsvs = [summarized_CBL_gtex_expression.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon4 { input :
    script_dir = script_dir,
    individual_tsv = Liver_SLC2A2_exon4_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_15a { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2 exon 4",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon4.out],
    dosage_threshold = 5
  }

  call gwas_tasks.summarize_individual_data_for_plotting as summarized_SLC2A2_gtex_exon6 { input :
    script_dir = script_dir,
    individual_tsv = Liver_SLC2A2_exon6_psi,
    length_sum_column_name = "Sum_of_allele",
    trait_column_name = "PSI"
  }

  call gwas_tasks.locus_plot as supp_fig_15b { input:
    script_dir = script_dir,
    chrom = 3,
    pos = 17100913,
    phenotype_name = "SCL2A2 exon 6",
    unit = "percent spliced in",
    data_tsvs = [summarized_SLC2A2_gtex_exon6.out],
    dosage_threshold = 5
  }

  call finemapping_tasks.str_tables_for_paper { input :
    script_dir = script_dir,
    str_pos_table = str_pos_table,
    repeat_units_table = repeat_units_table,
    intersects_gene = intersects_gene_annotation,
    intersects_exon = intersects_exon_annotation,
    intersects_CDS = intersects_CDS_annotation,
    intersects_five_prime_UTR = intersects_five_prime_UTR_annotation,
    intersects_three_prime_UTR = intersects_three_prime_UTR_annotation,
    intersects_UTR = intersects_UTR_annotation,
    phenotype_names = phenotype_names.n,
    assocs = str_gwas_results,
    black_assocs = ethnic_str_gwas_resultss[0],
    south_asian_assocs = ethnic_str_gwas_resultss[1],
    chinese_assocs = ethnic_str_gwas_resultss[2],
    irish_assocs = ethnic_str_gwas_resultss[3],
    white_other_assocs = ethnic_str_gwas_resultss[4],
    first_pass_finemapping_dfs = first_pass_finemapping_df.all_regions_concordance,
    followup_finemapping_dfs = select_all(followup_finemapping_conditions_df.df),
  }

  output {
#    File fig_1b_svg_out = fig_1b.svg
#    File fig_1b_png_out = fig_1b.png
    File fig_1e_svg_out = fig_1ef.barplot_svg
    File fig_1e_png_out = fig_1ef.barplot_png
    File fig_1f_svg_out = fig_1ef.heatmap_svg
    File fig_1f_png_out = fig_1ef.heatmap_png

    File fig_4a_svg_out = fig_4a.svg
    File fig_4a_png_out = fig_4a.png
    File fig_4b_svg_out = fig_4b.svg
    File fig_4b_png_out = fig_4b.png
    File fig_4f_svg_out = fig_4f.svg
    File fig_4f_png_out = fig_4f.png
    File fig_4g_svg_out = fig_4g.svg
    File fig_4g_png_out = fig_4g.png

    File stat_statements = post_finemapping.stat_statements
    File supp_fig_3_png = post_finemapping.cs_min_abs_corrs_png
    File supp_fig_3_svg = post_finemapping.cs_min_abs_corrs_svg
    File supp_fig_4_png = post_finemapping.susie_alpha_v_pip_png
    File supp_fig_4_svg = post_finemapping.susie_alpha_v_pip_svg
    File supp_fig_5a_png = post_finemapping.susie_alpha_histogram_png
    File supp_fig_5a_svg = post_finemapping.susie_alpha_histogram_svg
    File supp_fig_5b_png = post_finemapping.finemap_pip_histogram_png
    File supp_fig_5b_svg = post_finemapping.finemap_pip_histogram_svg
    File supp_fig_6_png = post_finemapping.susie_cs_finemap_total_pips_png
    File supp_fig_6_svg = post_finemapping.susie_cs_finemap_total_pips_svg
    File supp_fig_7a_png = post_finemapping.finemap_v_susie_consistency_STR_png
    File supp_fig_7a_svg = post_finemapping.finemap_v_susie_consistency_STR_svg
    File supp_fig_7b_png = post_finemapping.finemap_v_susie_consistency_SNP_png
    File supp_fig_7b_svg = post_finemapping.finemap_v_susie_consistency_SNP_svg

#    File doubly_finemapped_STRs = post_finemapping.doubly_finemapped_STRs
#    File confidently_finemapped_STRs = post_finemapping.confidently_finemapped_STRs
#    File overconfidently_finemapped_STRs = post_finemapping.overconfidently_finemapped_STRs

    # used for deciding confidently finemap
    File supp_fig_10ab_png = post_finemapping.susie_best_guess_png
    File supp_fig_10ab_svg = post_finemapping.susie_best_guess_svg
    File supp_fig_11a_png = post_finemapping.finemap_p_thresh_png
    File supp_fig_11a_svg = post_finemapping.finemap_p_thresh_svg
    File supp_fig_11b_png = post_finemapping.finemap_mac_png
    File supp_fig_11b_svg = post_finemapping.finemap_mac_svg
    File supp_fig_11c_png = post_finemapping.finemap_prior_std_derived_png
    File supp_fig_11c_svg = post_finemapping.finemap_prior_std_derived_svg
    File supp_fig_11d_png = post_finemapping.finemap_total_prob_png
    File supp_fig_11d_svg = post_finemapping.finemap_total_prob_svg
    File supp_fig_11e_png = post_finemapping.finemap_conv_tol_png
    File supp_fig_11e_svg = post_finemapping.finemap_conv_tol_svg

    # too conservative
    File supp_fig_11f_png = post_finemapping.finemap_prior_std_low_png
    File supp_fig_11f_svg = post_finemapping.finemap_prior_std_low_svg
    File supp_fig_12a_png = post_finemapping.susie_ratio_png
    File supp_fig_12a_svg = post_finemapping.susie_ratio_svg
    File supp_fig_12b_png = post_finemapping.finemap_ratio_png
    File supp_fig_12b_svg = post_finemapping.finemap_ratio_svg

    File supp_fig_15a_svg_out = supp_fig_15a.svg
    File supp_fig_15a_png_out = supp_fig_15a.png
    File supp_fig_15b_svg_out = supp_fig_15b.svg
    File supp_fig_15b_png_out = supp_fig_15b.png

   File singly_finemapped_strs_for_paper = str_tables_for_paper.singly_finemapped_strs_for_paper
   File singly_finemapped_strs_sorted = str_tables_for_paper.singly_finemapped_strs_sorted
   File confidently_finemapped_strs_for_paper = str_tables_for_paper.confidently_finemapped_strs_for_paper
   File confidently_finemapped_strs_sorted = str_tables_for_paper.confidently_finemapped_strs_sorted
  }
}
