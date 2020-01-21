#!/bin/bash

for chr in $(seq 1 22); do
        cut -f3 ../vcf_1_sample/chr${chr}.vcf | sort | uniq -c | sort | \
		awk '{ if ($1 != 1) print $2 }' \
                > chr${chr}.txt
done
