#!/usr/bin/env bash

if [ -z "$PHENOTYPE" ] ; then echo "PHENOTYPE is unset" ; exit 1 ; fi
if [ -z "$CHROM" ] ; then echo "CHROM is unset" ; exit 1 ; fi
echo "PHENOTYPE $PHENOTYPE CHROM $CHROM CONDITIONAL $CONDITIONAL START $START END $END"
echo "PHENOTYPE $PHENOTYPE CHROM $CHROM CONDITIONAL $CONDITIONAL START $START END $END" >&2
# conditional should look like 'STR{STRs}__ISNP{ISNPs}__ASNP{ASNPs}'


if [[ -n "$CONDITIONAL" && ( -z "$START" || -z "$END" ) ]] ; then
	echo "Conditional input, but not start or end. Exiting."
	exit 1
fi

if [[ -z "$CONDITIONAL" && ( -n "$START" || -n "$END" ) ]] ; then
	echo "No conditional input, but start or end. Exiting."
	exit 1
fi

if [ -z "$CONDITIONAL" ] ; then
	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp/chrs/chr"$CHROM"
else
	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp_conditional
	OUT_DIR="$OUT_DIR"/chr"$CHROM"_"$START"_"$END"_"$CONDITIONAL"
fi

cd "$OUT_DIR" || {
	echo "Failed to move to the output directory $OUT_DIR" ;
	exit 1 ;
}

if [ -z "$CONDITIONAL" ] ; then
	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/plink_snp/input/rin_phenotype_and_covars.tab 
	mkdir -p ../../logs
	LOG_OUT=../../logs/chr"$CHROM".plink.stdout 
	LOG_ERR=../../logs/chr"$CHROM".plink.stderr
	BED_FILE_COMMAND=" "
else
	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/conditional_inputs/chr"$CHROM"_"$CONDITIONAL"_plink.tab
	LOG_OUT=plink.stdout
	LOG_ERR=plink.stderr
	echo -e chr"$CHROM""\t""$START""\t""$((END + 1))""\n" > region.bed
	BED_FILE_COMMAND=" --extract bed1 region.bed"
fi

"$UKB"/utilities/plink2 \
    --pheno "$PHENO_FILE" \
    --no-psam-pheno \
    --pheno-name rin_"$PHENOTYPE" \
    --covar-name $(head -n 1 "$PHENO_FILE" | cut -f 4-) \
    --pfile "$UKB"/array_imputed/pfile_converted/chr"$CHROM" \
    --chr "$CHROM" \
    $BED_FILE_COMMAND \
    --mac 20 \
    --glm omit-ref pheno-ids hide-covar \
    --ci 0.99999995 \
    --memory 108000 \
    --threads 28 \
    > "$LOG_OUT" \
    2> "$LOG_ERR"

# give filesystem time to register the file plink has written
sleep 60 

mv "$OUT_DIR"/plink2.rin_"$PHENOTYPE".glm.linear \
   "$OUT_DIR"/plink2.rin_"$PHENOTYPE".glm.linear.done
