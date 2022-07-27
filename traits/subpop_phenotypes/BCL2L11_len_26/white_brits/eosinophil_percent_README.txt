Run date: 2022_07_14
Loading phenotype eosinophil_percent from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/eosinophil_percent_30210.txt 
Subsetting to samples at sample_qc/subpop_runs/BCL2L11_len_26/white_brits/no_phenotype/combined.sample
Loading categorical covar with field id 30213 from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/eosinophil_percent_device_id_30213.txt
Choosing phenotype value and age for each participant based on the first visit for which this phenotype had a recorded value. If the participant had this phenotype measured at multiple visits, only the phenotype value and age at the first visit are being used. The age is being loaded from the shared_covars file /expanse/projects/gymreklab/jmargoli/ukbiobank/traits/shared_covars/asssessment_ages.npy . An additional dummy covariate indicating visit number is being added for each visit whose data is used beyond the first.
Dropping 1 samples whose first (non yet dropped) data recording for this phenotype occurred on assessment 1 because that is less than 50 total samples.
Dropping 1 samples whose first (non yet dropped) data recording for this phenotype occurred on assessment 2 because that is less than 50 total samples.
Using categorical covar eosinophil_percent_device_id, default category when all dummy variables are False is AJ38695
All samples have all covariates.
