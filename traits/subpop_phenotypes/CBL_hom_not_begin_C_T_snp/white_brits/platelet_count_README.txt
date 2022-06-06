Run date: 2022_05_25
Loading phenotype platelet_count from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/platelet_count_30080.txt 
Subsetting to samples at sample_qc/subpop_runs/CBL_hom_not_begin_C_T_snp/white_brits/no_phenotype/combined.sample
Loading categorical covar with field id 30083 from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/platelet_count_device_id_30083.txt
Choosing phenotype value and age for each participant based on the first visit for which this phenotype had a recorded value. If the participant had this phenotype measured at multiple visits, only the phenotype value and age at the first visit are being used. The age is being loaded from the shared_covars file /expanse/projects/gymreklab/jmargoli/ukbiobank/traits/shared_covars/asssessment_ages.npy . An additional dummy covariate indicating visit number is being added for each visit whose data is used beyond the first.
Dropping 52 samples whose first (non yet dropped) data recording for this phenotype occurred on assessment 2 because that is less than the fraction 0.001 of the total samples for this phenotype which is 213464.
Using categorical covar platelet_count_device_id, default category when all dummy variables are False is AJ38695
All samples have all covariates.
