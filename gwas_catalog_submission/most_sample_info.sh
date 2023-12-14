#!/bin/bash

cd $UKB/wdl_cache/associations
for file in *gz ; do 
	population=$(echo "$file" | awk 'BEGIN {FS="_"} {printf $1; if ($1 == "white" || $1 == "south") { print "_" $2 }}')
	phenotype=$(echo "$file" | sed -e 's/_str_gwas_results.*//' | sed -e 's/^\(south_asian\|chinese\|white_british\|white_other\|irish\|black\)_//')
	printf "%s\t" "$(echo "$file" | sed -e 's/\.tab\.gz//')"
	if [[ "$population" == "white_british" ]] ; then
		printf "discovery\t"
	else
		printf "replication\t"
	fi
	printf "TODO number of individuals\t"
	printf "No\t"
	printf "NA\t"
	printf "NA\t"
	printf "NA\t"
	python -c 'print(" ".join(word.capitalize() for word in "'"$population"'".replace("_", " ").split()), end="\t")'
	printf "NA\t"
	printf "NA\t"
	if [[ "$population" == "white_british" ]] ; then
		printf "self reported|genetically determined\t"
	else
		printf "self reported|not in genetically determined White British set\t"
	fi
	printf "United Kingdom\n"
done
