#!/bin/bash

for i in {1..22} ; do
	wget https://s3.amazonaws.com/snp-str-imputation/1000genomes/1kg.snp.str.chr${i}.vcf.gz
	wget https://s3.amazonaws.com/snp-str-imputation/1000genomes/1kg.snp.str.chr${i}.vcf.gz.tbi
done
