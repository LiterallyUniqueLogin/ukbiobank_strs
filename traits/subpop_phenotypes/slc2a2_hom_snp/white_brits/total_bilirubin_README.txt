Run date: 2022_05_24
Loading phenotype total_bilirubin from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/total_bilirubin_30840.txt 
Subsetting to samples at sample_qc/subpop_runs/slc2a2_hom_snp/white_brits/no_phenotype/combined.sample
Loading categorical covar with field id 30842 from txt file  /expanse/projects/gymreklab/jmargoli/ukbiobank/main_dataset/extracted_data/total_bilirubin_aliquot_30842.txt
Choosing phenotype value and age for each participant based on the first visit for which this phenotype had a recorded value. If the participant had this phenotype measured at multiple visits, only the phenotype value and age at the first visit are being used. The age is being loaded from the shared_covars file /expanse/projects/gymreklab/jmargoli/ukbiobank/traits/shared_covars/asssessment_ages.npy . An additional dummy covariate indicating visit number is being added for each visit whose data is used beyond the first.
Dropping 15 samples whose first (non yet dropped) data recording for this phenotype occurred on assessment 1 because that is less than 50 total samples.
No qc'ed samples of this ethnicity had value 4 for categorical covariate total_bilirubin_aliquot,30842.
Using categorical covar total_bilirubin_aliquot, default category when all dummy variables are False is 1
All samples have all covariates.
