#!/bin/bash

cd $UKB/wdl_cache/associations
for file in *gz ; do 
	population=$(echo "$file" | awk 'BEGIN {FS="_"} {printf $1; if ($1 == "white" || $1 == "south") { print "_" $2 }}')
	phenotype=$(echo "$file" | sed -e 's/_str_gwas_results.*//' | sed -e 's/^\(south_asian\|chinese\|white_british\|white_other\|irish\|black\)_//')
	printf "%s\t" "$(echo "$file" | sed -e 's/\.tab\.gz//')"
	printf "Whole genome sequencing\t"
	printf "NA\t"
	printf "NA\t"
	printf "AssociaTR (TRTools) 5.0.2\t"
	printf "Yes\t"
	printf "Saini et al. SNP + STR reference haplotype panel\t"
	printf "Beagle32 v5.1 (build 25Nov19.28d)\t"
	printf "%s\t" $(( $(zcat $(echo "$file" | sed -e 's/'"$population"'/white_british/') | wc -l) -1))
	printf "Linear regression of trait on repeat length summed over both chromosomes\t"
	printf "Polymorphic short tandem repeats make widespread contributions to blood and serum traits\t"
	printf "TODO adjusted covariates\t"
	printf "$phenotype\t"
	printf "None\t"
	printf "$file\t"
	printf "%s\t" "$(md5sum "$file" | awk '{print $1}')"
	printf "TODO README\t"
	printf "GRCh37\t"
	printf "NA\t"
	printf "UKB $population population\t"
	printf "NA\t"
	printf "Yes\t"
	printf "No\t"
	printf "No\t"
	printf "combined\t"
	printf "1-based\n"
done
