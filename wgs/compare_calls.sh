#!/bin/bash
set -euxo pipefail

# requires liftover, compareSTR, bcftools
# uses $MERGED_FILE $FINEMAPPED_STRS

# ---- update merged file ----

# liftover
bcftools query -f'%CHROM\t%POS0\t%POS0\t%ID\n' "$MERGED_FILE" -o merged_pos.bed
bcftools query -f'%CHROM\t%POS0\t%INFO/START\t%ID\n' "$MERGED_FILE" | \
	awk '{print $1 "\t" $2 "\t" $3 - 1 "\t" $4}' > merged_start.bed
bcftools query -f'%CHROM\t%POS0\t%INFO/END\t%ID\n' "$MERGED_FILE" | \
	awk '{print $1 "\t" $2 "\t" $3 - 1 "\t" $4}' > merged_end.bed

chain=hg38ToHg19.over.chain.gz
wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O "$chain"

for file in merged_pos.bed merged_start.bed merged_end.bed ; do 
	liftOver $file "$chain" "$file"_hg19 "$file"_unmapped
done

zcat "$MERGED_FILE" | cut -f 1 > merged1.tab
zcat "$MERGED_FILE" | cut -f 3- > merged3on.tab
paste merged1.tab <( cut -f 2 merged_pos.bed | awk '{print $1 + 1}' )  merged3on.tab > merged.hg19.vcf.temp

#also make samples match - ID to ID_ID
bcftools \
	-S <(bcftools query -l "$MERGED_FILE" | awk '{print $1 " " $1 "_" $1}' ) \
	-a <(awk '{print $1 "\t" $2 "\t" $3 + 1}' merged_start.bed) -c CHROM,POS,INFO/START \
	merged.hg19.vcf.temp | \
bcftools \
	-a <(awk '{print $1 "\t" $2 "\t" $3 + 1}' merged_end.bed) -c CHROM,POS,INFO/END \
	- | \
sed -e 's/ID=chr/ID=/' -e 's/^chr//' > merged.hg19.vcf
# sed matches the different chromosome naming conventions, changing from chr1 to 1

bgzip merged.hg19.vcf
tabix merged.hg19.vcf.gz

# --- update and subset imputed files
bcftools annotate -h <( echo "##command=HipSTR" ) -o imputed_strs.reheadered.vcf.gz "$FINEMAPPED_STRS"
tabix filtered_imputed_strs.vcf.gz

compareSTR --ignore-phasing --balanced-accuracy --fraction-concordant-len-sum --vcf2-beagle-probabilities \
	--vcf1 merged.hg19.vcf.gz \
	--vcf2 imputed_strs.reheadered.vcf.gz \
	--out out
