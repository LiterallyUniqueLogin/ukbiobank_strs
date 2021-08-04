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

if [[ -z "$PROJECT_TEMP" ]] ; then
	echo "PROJECT_TEMP not defined. Exiting."
	exit 1
fi

if [ -z "$CONDITIONAL" ] ; then
	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp/chrs/chr"$CHROM"
else
	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp_conditional
	OUT_DIR="$OUT_DIR"/chr"$CHROM"_"$START"_"$END"_"$CONDITIONAL"
fi

mkdir -p "$OUT_DIR"

cd "$PROJECT_TEMP" || {
	echo "Failed to move to the temp directory $PROJECT_TEMP" ;
	exit 1 ;
}

temp_name=plink_snp/phenotype_"$PHENOTYPE"_chr_"$CHROM"_start_"$START"_end_"$END"_conditional_"$CONDITIONAL"
mkdir -p "$temp_name"
cd "$temp_name" || {
	echo "Failed to move to $temp_name"
	exit 1 ;
}

if [ -z "$CONDITIONAL" ] ; then
	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/plink_snp/input/transformed_phenotype_and_covars.tab 
	mkdir -p "$OUT_DIR"/../../logs
	LOG_OUT="$OUT_DIR"/../../logs/chr"$CHROM".plink.stdout 
	LOG_ERR="$OUT_DIR"/../../logs/chr"$CHROM".plink.stderr
	BED_FILE_COMMAND=" "
else
	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/conditional_inputs/chr"$CHROM"_"$CONDITIONAL"_plink.tab
	LOG_OUT=plink.stdout
	LOG_ERR=plink.stderr
	echo -e chr"$CHROM""\t""$START""\t""$((END + 1))" > region.bed
	BED_FILE_COMMAND=" --extract bed1 region.bed"
fi

if [ -z "$BINARY" ] ;  then
	PHENO_NAME=rin_"$PHENOTYPE"
	COLS=""
else
	PHENO_NAME="$PHENOTYPE"
	COLS="cols=+a1countcc,-tz,-nobs,-test"
fi

"$UKB"/utilities/plink2 \
    --pheno "$PHENO_FILE" \
    --no-psam-pheno \
    --pheno-name "$PHENO_NAME" \
    --covar-name $(head -n 1 "$PHENO_FILE" | cut -f 4-) \
    --pfile "$UKB"/array_imputed/pfile_converted/chr"$CHROM" \
    --chr "$CHROM" \
    $BED_FILE_COMMAND \
    --mac 20 \
    --glm omit-ref pheno-ids hide-covar $COLS \
    --ci 0.99999995 \
    --memory 56000 \
    --threads 28 \
    > "$LOG_OUT" \
    2> "$LOG_ERR"

mv ./* "$OUT_DIR"

mv "$OUT_DIR"/plink2.rin_"$PHENOTYPE".glm.linear \
   "$OUT_DIR"/plink2.rin_"$PHENOTYPE".glm.linear.done
