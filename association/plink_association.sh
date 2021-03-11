#!/usr/bin/env bash

if [ -z "$PHENOTYPE" ] ; then echo "PHENOTYPE is unset" ; exit 1 ; fi
if [ -z "$CHROM" ] ; then echo "CHROM is unset" ; exit 1 ; fi
echo "PHENOTYPE $PHENOTYPE CHROM $CHROM"
echo "PHENOTYPE $PHENOTYPE CHROM $CHROM" >&2

cd "$UKB"/association/results/"$PHENOTYPE"/plink_snp/chrs/chr"$CHROM" || {
	echo "Failed to move to the output directory" ;
	exit 1 ;
}

mkdir -p ../../logs

"$UKB"/utilities/plink2 \
    --pheno  "$UKB"/association/results/"$PHENOTYPE"/plink_snp/input/rin_phenotype_and_covars.tab \
    --no-psam-pheno \
    --pheno-name rin_"$PHENOTYPE" \
    --covar-name $(cat "$UKB"/traits/shared_covars/covar_names.txt) $(cat "$UKB"/traits/phenotypes/"$PHENOTYPE"_covar_names.txt) \
    --pfile "$UKB"/array_imputed/pfile_converted/chr"$CHROM" \
    --chr "$CHROM" \
    --mac 20 \
    --glm omit-ref pheno-ids hide-covar \
    --ci 0.99999995 \
    --memory 108000 \
    --threads 28 \
    > "$UKB"/association/results/"$PHENOTYPE"/plink_snp/logs/chr"$CHROM".plink.stdout \
    2> "$UKB"/association/results/"$PHENOTYPE"/plink_snp/logs/chr"$CHROM".plink.stderr

mv "$UKB"/association/results/"$PHENOTYPE"/plink_snp/chrs/chr"$CHROM"/plink2.rin_"$PHENOTYPE".glm.linear \
   "$UKB"/association/results/"$PHENOTYPE"/plink_snp/chrs/chr"$CHROM"/plink2.rin_"$PHENOTYPE".glm.linear.done
