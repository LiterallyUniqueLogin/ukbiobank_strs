# number of loci removed by plink
cd $UKB/association/results
{ for pheno in $(ls $UKB/finemapping/finemap_results | grep -v height )  ; do grep removed $pheno/plink_snp/chrs/chr*/plink2.log | cut -f2- -d: | cut -f1 -d' ' | datamash sum 1 ; done }

# number of people excluded due to having only measurement in assessment 2
cd $UKB/traits/phenotypes/white_brits
grep 'assessment 2' $({ for pheno in $(ls $UKB/finemapping/finemap_results | grep -v height) ; do echo ${pheno}_README.txt ; done }) | vim -

# number of people excluded due to having aliquot=3
cd $UKB/traits/phenotypes/white_brits
grep -v 'assessment 2' $({ for pheno in $(ls $UKB/finemapping/finemap_results | grep -v height) ; do echo ${pheno}_README.txt ; done }) | grep -v 'Run date' | grep -v Loading | grep -v Subsetting | grep -v 'All samples have' | grep -v 'Dropping 0' | grep -v 'Dropping 26' | grep -vP ':.$' | grep -v 'Using categorical covar' | grep -v 'multiple visits' | cut -f2- -d: | grep 'value 3' | vim -

# number of people excluded due to having aliquot=4
cd $UKB/traits/phenotypes/white_brits
grep -v 'assessment 2' $({ for pheno in $(ls $UKB/finemapping/finemap_results | grep -v height) ; do echo ${pheno}_README.txt ; done }) | grep -v 'Run date' | grep -v Loading | grep -v Subsetting | grep -v 'All samples have' | grep -v 'Dropping 0' | grep -v 'Dropping 26' | grep -vP ':.$' | grep -v 'Using categorical covar' | grep -v 'multiple visits' | cut -f2- -d: | grep 'value 4' | vim -


# check plink error codes
cd $UKB/association/results
{ for phenotype in $(ls $UKB/finemapping/finemap_results | grep -v height) ; do cut -f 15 $phenotype/plink_snp/results.tab ; done } | python -c $'import fileinput\nerrs = set()\nfor line in fileinput.input():\n errs.add(line)\nprint(errs)'
