version 1.0

import "../gwas_wdl/gwas_tasks.wdl"
import "finemapping_tasks.wdl"
import "finemap_one_region_workflow.wdl"
import "susie_one_region_workflow.wdl"

workflow run_simulation_in_region {

  input {
    Array[String] causal_vars
    Array[String] causal_betas
    Int seed

    String finemap_command

    String script_dir

    VCF str_vcf
    bgen imputed_snp_bgen
    PFiles imputed_snp_pfile
    File snp_vars_to_filter_from_finemapping

    region bounds

    File phenotype_samples

    File all_samples_list

    String prefix
  }

  call finemapping_tasks.simulate_phenotype { input :
    script_dir = script_dir,
    str_vcf = str_vcf,
    snp_bgen = imputed_snp_bgen,
    all_samples_list = all_samples_list,
    samples_list = phenotype_samples,
    chrom = bounds.chrom,
    seed = seed,
    causal_vars = causal_vars,
    causal_betas = causal_betas
  }

  call gwas_tasks.regional_my_str_gwas { input :
    script_dir = script_dir,
    str_vcf = str_vcf,
    transformed_phenotype = simulate_phenotype.phenotype,
    all_samples_list = all_samples_list,
    is_binary = false,
    binary_type = "blah",
    bounds = bounds,
    phenotype_name = "simulation",
    no_details = true,
  }

  call gwas_tasks.prep_plink_input { input :
    script_dir = script_dir,
    transformed_phenotype = simulate_phenotype.phenotype,
    pheno_covar_names = "/dev/null",
    is_binary = false,
    binary_type = "blah",
    phenotype_name = "simulation",
  }

  call gwas_tasks.plink_snp_association { input :
    script_dir = script_dir,
    plink_command = "plink2",
    imputed_snp_p_file = imputed_snp_pfile,
    pheno_data = prep_plink_input.data,
    chrom = bounds.chrom,
    phenotype_name = "simulation",
    binary_type = "linear",
    start = bounds.start,
    end = bounds.end
  }

  call finemap_one_region_workflow.finemap_one_region { input :
    script_dir = script_dir,
    finemap_command = finemap_command,
    str_vcf = str_vcf,
    imputed_snp_bgen = imputed_snp_bgen,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    my_str_gwas = regional_my_str_gwas.data,
    plink_snp_gwas = plink_snp_association.data,
    phenotype_name = "simulation",
    bounds = bounds,
    all_samples_list = all_samples_list,
    prefix = "~{prefix}FINEMAP_"
  }

  call susie_one_region_workflow.susie_one_region { input :
    script_dir = script_dir,
    str_vcf = str_vcf,
    imputed_snp_bgen = imputed_snp_bgen,
    snp_vars_to_filter_from_finemapping = snp_vars_to_filter_from_finemapping,
    phenotype_samples = phenotype_samples,
    transformed_phenotype_data = simulate_phenotype.phenotype,
    my_str_gwas = regional_my_str_gwas.data,
    plink_snp_gwas = plink_snp_association.data,
    phenotype_name = "simulation",
    bounds = bounds,
    all_samples_list = all_samples_list,
    prefix = "~{prefix}SuSiE_"
  }

  output {
    File snp_assocs = plink_snp_association.data
    File str_assocs = regional_my_str_gwas.data

    FINEMAP_output finemap_results = finemap_one_region.finemap_output
    SuSiE_output susie_results = susie_one_region.susie_output
  } 
}
