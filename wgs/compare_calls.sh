#!/bin/bash
set -euxo pipefail

# requires liftover, compareSTR, bcftools
# uses $MERGED_FILE $FINEMAPPED_STRS

# ---- update merged file ----

# liftover
bcftools query -f'%CHROM\t%POS0\t%POS\t%ID\n' "$MERGED_FILE" -o merged_pos.bed
bcftools query -f'%CHROM\t%POS0\t%INFO/START\t%ID\n' "$MERGED_FILE" -o merged_start.bed
bcftools query -f'%CHROM\t%POS0\t%INFO/END\t%ID\n' "$MERGED_FILE" -o merged_end.bed

chain=hg38ToHg19.over.chain.gz
wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O "$chain"

for file in merged_pos.bed merged_start.bed merged_end.bed ; do 
	liftOver $file "$chain" hg19_"$file" "$file"_unmapped
	bgzip -f hg19_"$file"
	tabix -f hg19_"$file".gz
done

zcat "$MERGED_FILE" | grep -P '^##' > temp.merged.hg19.vcf
printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" >> temp.merged.hg19.vcf
#also make samples match - ID to ID_ID
zcat "$MERGED_FILE" | grep -P '^#CHROM' | cut -f 10- | sed -e 's/\([^\t]\+\)/\1_\1/g' >> temp.merged.hg19.vcf
zcat "$MERGED_FILE" | grep -Pv '^#' | cut -f 1 > merged1.tab
zcat "$MERGED_FILE" | grep -Pv '^#' | cut -f 3- > merged3on.tab
paste merged1.tab <( zcat hg19_merged_pos.bed.gz | cut -f 3 ) merged3on.tab >> temp.merged.hg19.vcf

zcat hg19_merged_start.bed.gz | awk '{print $1 "\t" $2 + 1 "\t" $3}' > hg19_merged_start.annotation
bgzip -f hg19_merged_start.annotation
tabix -f -p bed hg19_merged_start.annotation.gz
bcftools annotate \
	-a hg19_merged_start.annotation.gz -c CHROM,POS,INFO/START \
	temp.merged.hg19.vcf | \
bcftools annotate \
	-a hg19_merged_end.bed.gz -c CHROM,POS,INFO/END,- | \
sed -e 's/ID=chr/ID=/' -e 's/^chr//' > merged.hg19.vcf
# sed matches the different chromosome naming conventions, changing from chr1 to 1

bgzip -f merged.hg19.vcf
tabix -f merged.hg19.vcf.gz

# --- update and subset imputed files
bcftools annotate -h <( echo "##command=HipSTR" ) -o imputed_strs.reheadered.vcf.gz -O z "$FINEMAPPED_STRS"
tabix -f imputed_strs.reheadered.vcf.gz

compareSTR --ignore-phasing --balanced-accuracy --fraction-concordant-len-sum --vcf2-beagle-probabilities \
	--vcf1 merged.hg19.vcf.gz \
	--vcf2 imputed_strs.reheadered.vcf.gz \
	--out out \
	--noplot
