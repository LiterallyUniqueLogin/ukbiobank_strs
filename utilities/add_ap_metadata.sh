#!/bin/bash
#Beagle currently has a bug where if you ask it for the AP
#Format field but not the GP format field, it will output the correct fields,
#but write out the metadata line for the GP field and not the AP field.
#This utility fixes that problem.
#More specifically: 
#it unzips a vcf.gz file,
#replaces its GP metadata line with the corresponding AP line
#rezips the vcf and reindexes it
#If $1 is foo.vcf.gz, then this will overwrite the files 
#foo.vcf, foo.vcf.gz and foo.vcf.gz.tbi

if [ -z "$1" ] ; then
        echo "Didn't give the first argument - should be the .vcf.gz file without that file extension"
        exit -1
fi

if [ ! -z "$2" ] ; then
        echo "Gave more than one argument, only expecting one argument"
        exit -1
fi

source ~/.bashrc
conda activate bcftools

FILE_NAME=$1

bgzip -df $FILE_NAME.vcf.gz
sed 's/##FORMAT=<ID=GP.*/##FORMAT=<ID=AP1,Number=A,Type=Float,Description="Estimated Allele 1 Probability">\n##FORMAT=<ID=AP2,Number=A,Type=Float,Description="Estimated Allele 2 Probability">/' \
	-i \
	$FILE_NAME.vcf
bgzip -f $FILE_NAME.vcf
tabix -f $FILE_NAME.vcf.gz
conda deactivate
