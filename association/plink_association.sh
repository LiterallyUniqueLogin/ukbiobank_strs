#!/usr/bin/env bash

set -x

if [ -z "$PHENOTYPE" ] ; then echo "PHENOTYPE is unset" ; exit 1 ; fi
if [ -z "$BINARY_TYPE" ] ; then echo "BINARY_TYPE is unset" ; exit 1 ; fi
if [ -z "$CHROM" ] ; then echo "CHROM is unset" ; exit 1 ; fi
if [ -z "$OUT_DIR" ] ; then echo "OUTDIR is unset" ; exit 1 ; fi
if [ -z "$PHENO_FILE" ] ; then echo "PHENO_FILE is unset" ; exit 1 ; fi
if [ -z "$PLINK_EXECUTABLE" ] ; then echo "PLINK_EXECUTABLE is unset" ; exit 1 ; fi
if [ -z "$P_FILE" ] ; then echo "P_FILE is unset" ; exit 1 ; fi
echo "PHENOTYPE $PHENOTYPE BINARY_TYPE $BINARY_TYPE CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END"
echo "PHENOTYPE $PHENOTYPE BINARY_TYPE $BINARY_TYPE CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END" >&2

if [[ ( -n "$START" && -z "$END" ) || ( -z "$START" && -n "$END" ) ]] ; then
	echo "Start but not end or end but not start. Exiting."
	exit 1
fi

if [[ ( $BINARY_TYPE != linear ) && ( $BINARY_TYPE != linear_binary ) && ( $BINARY_TYPE != logistic ) ]] ; then
	echo "BINARY_TYPE not set correctly. Exiting."
	exit 1
fi

if [[ -z "$PROJECT_TEMP" ]] ; then
	echo "PROJECT_TEMP not defined. Exiting."
	exit 1
fi

#if [ -z "$START" ] ; then
#	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp"$SUFFIX"/chrs/chr"$CHROM"
#elif [ -z "$CONDITIONAL" ] ; then
#	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp"$SUFFIX"/batches/chr"$CHROM"_"$START"_"$END"
#else
#	OUT_DIR="$UKB"/association/results/"$PHENOTYPE"/plink_snp"$SUFFIX"_conditional
#	OUT_DIR="$OUT_DIR"/chr"$CHROM"_"$START"_"$END"_"$CONDITIONAL"
#fi

mkdir -p "$OUT_DIR"

temp_name="$PROJECT_TEMP"/plink_snp"$SUFFIX"/phenotype_"$PHENOTYPE"_chr_"$CHROM"_start_"$START"_end_"$END"_conditional_"$CONDITIONAL"
mkdir -p "$temp_name"
cd "$temp_name" || {
	echo "Failed to move to $temp_name"
	exit 1 ;
}

#if [ -z "$CONDITIONAL" ] ; then
#	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/plink_snp"$SUFFIX"/input/transformed_phenotype_and_covars.tab 
#else
#	PHENO_FILE="$UKB"/association/results/"$PHENOTYPE"/conditional_inputs/chr"$CHROM"_"$CONDITIONAL"_plink"$SUFFIX".tab
#fi

if [ -z "$START" ] ; then
	BED_FILE_COMMAND=" "
else
	echo -e chr"$CHROM""\t""$START""\t""$((END + 1))" > region.bed
	BED_FILE_COMMAND=" --extract bed1 region.bed"
fi

if [[ "$BINARY_TYPE" == linear ]] ;  then
	PHENO_NAME=rin_"$PHENOTYPE"
	COLS=""
else
	PHENO_NAME="$PHENOTYPE"
	COLS="cols=+a1countcc,-tz,-nobs,-test"
fi

#$UKB"/utilities/plink2 \
#    --pfile "$UKB"/array_imputed/pfile_converted/chr"$CHROM" \
"$PLINK_EXECUTABLE" \
    --pheno "$PHENO_FILE" \
    --no-psam-pheno \
    --pheno-name "$PHENO_NAME" \
    --covar-name $(head -n 1 "$PHENO_FILE" | cut -f 4-) \
    --pfile "$P_FILE" \
    --chr "$CHROM" \
    $BED_FILE_COMMAND \
    --mac 20 \
    --glm omit-ref pheno-ids hide-covar $COLS \
    --ci 0.99999995 \
    --memory 56000 \
    --threads 28 \
    > plink.stdout \
    2> plink.stderr

mv "$temp_name"/* "$OUT_DIR"

if [[ "$BINARY_TYPE" == linear ]] ;  then
	RIN=rin_
else
	RIN=''
fi

if [[ "$BINARY_TYPE" == logistic ]] ; then
	RUNTYPE="logistic.hybrid"
else
	RUNTYPE="linear"
fi

mv "$OUT_DIR"/plink2."$RIN""$PHENOTYPE".glm."$RUNTYPE" \
   "$OUT_DIR"/plink2."$RIN""$PHENOTYPE".glm."$RUNTYPE".done

