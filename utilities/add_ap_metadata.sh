#!/bin/bash
#usage:
#add_ap_metadata.sh file.vcf.gz gz
#or
#add_ap_metadata.sh file.vcf no-gz
#Beagle currently has a bug where if you ask it for the AP
#Format field but not the GP format field, it will output the correct fields,
#but write out the metadata line for the GP field and not the AP field.
#This utility fixes that problem.
#More specifically:
#In the gz mode
#it unzips a vcf.gz file,
#replaces its GP metadata line with the corresponding AP line
#rezips the vcf and reindexes it
#If $1 is foo.vcf.gz, then this will overwrite the files 
#foo.vcf, foo.vcf.gz and foo.vcf.gz.tbi
#In the no-gz mode it only does the replacement, not unzipping or rezipping or reindexing

if [ -z "$1" ] ; then
        echo "Error: Didn't give the first argument - should be the .vcf/.vcf.gz file without that file extension"
        exit -1
fi

if [ -z "$2" ] ; then
        echo "Error: Didn't give the second argument - should be either gz or no-gz specifying the file extension"
        exit -1
fi

if [[ "$2" != "gz" ]] && [[ "$2" != "no-gz" ]] ; then
	echo "error: second argument should be gz or no-gz"
	exit -1
fi

if [ ! -z "$3" ] ; then
        echo "Error: Gave more than two arguments, only expecting two arguments"
        exit -1
fi

source ~/.bashrc
conda activate bcftools

FILE_NAME=$1


if [[ "$2" == "gz" ]] ; then
	bgzip -df $file_name.vcf.gz
fi 

sed 's/##FORMAT=<ID=GP.*/##FORMAT=<ID=AP1,Number=A,Type=Float,Description="Estimated Allele 1 Probability">\n##FORMAT=<ID=AP2,Number=A,Type=Float,Description="Estimated Allele 2 Probability">/' \
	-i \
	$FILE_NAME.vcf

if [[ "$2" == "gz" ]] ; then
	bgzip -f $FILE_NAME.vcf
	tabix -f $FILE_NAME.vcf.gz
fi
conda deactivate
