#!/bin/bash

for idx in $(seq 1 22); do 
	correct_checksum=$(head -n $idx checksums.txt | tail -n 1 | cut -f 1 -d' ')
	actual_checksum=$(/projects/ps-gymreklab/jmargoli/ukbiobank/ukb_utilities/ukbmd5 ukb_hap_chr${idx}_v2.bgen | tail -n 1 | cut -d'=' -f3 )
	if [ "$actual_checksum" = "$correct_checksum" ]; then 
		echo "checksum correct for chr $idx" 
	else 
		echo "checksum NOT correct for chr $idx - REDOWNLOAD"
	fi 
done
