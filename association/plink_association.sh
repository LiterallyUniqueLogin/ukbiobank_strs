#!/usr/bin/env bash

set -x

if [ -z "$PHENOTYPE" ] ; then echo "PHENOTYPE is unset" ; exit 1 ; fi
if [ -z "$BINARY_TYPE" ] ; then echo "BINARY_TYPE is unset" ; exit 1 ; fi
if [ -z "$CHROM" ] ; then echo "CHROM is unset" ; exit 1 ; fi
if [ -z "$OUT_DIR" ] ; then echo "OUTDIR is unset" ; exit 1 ; fi
if [ -z "$PHENO_FILE" ] ; then echo "PHENO_FILE is unset" ; exit 1 ; fi
if [ -z "$PLINK_EXECUTABLE" ] ; then echo "PLINK_EXECUTABLE is unset" ; exit 1 ; fi
if [ -z "$P_FILE" ] ; then echo "P_FILE is unset" ; exit 1 ; fi
if [ -z "$PROJECT_TEMP" ] ; then echo "PROJECT_TEMP is unset" ; exit 1 ; fi
echo "PHENOTYPE $PHENOTYPE BINARY_TYPE $BINARY_TYPE CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END PROJECT_TEMP $PROJECT_TEMP"
echo "PHENOTYPE $PHENOTYPE BINARY_TYPE $BINARY_TYPE CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END PROJECT_TEMP $PROJECT_TEMP" >&2

if [[ ( -n "$START" && -z "$END" ) || ( -z "$START" && -n "$END" ) ]] ; then
	echo "Start but not end or end but not start. Exiting."
	exit 1
fi

if [[ ( $BINARY_TYPE != linear ) && ( $BINARY_TYPE != logistic ) && ( $BINARY_TYPE != firth ) ]] ; then
	echo "BINARY_TYPE not set correctly. Exiting."
	exit 1
fi

start_dir=$(pwd)
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

temp_name="$PROJECT_TEMP"/plink_snp"$SUFFIX"/phenotype_"$PHENOTYPE"_chr_"$CHROM"_start_"$START"_end_"$END"
mkdir -p "$temp_name"
cd "$temp_name" || {
	echo "Failed to move to $temp_name"
	exit 1 ;
}

if [ -z "$START" ] ; then
	BED_FILE_COMMAND=" "
else
	echo -e chr"$CHROM""\t""$START""\t""$((END + 1))" > region.bed
	BED_FILE_COMMAND=" --extract bed1 region.bed"
fi

if [[ "$BINARY_TYPE" == linear ]] ;  then
	COLS="cols=-test,-nobs"
	BINARY_GLM_FLAG=""
else
	COLS="cols=+gcountcc,-nobs,-test"
	#BINARY_GLM_FLAG="single-prec-cc cc-residualize"
	BINARY_GLM_FLAG=""
	if [[ "$BINARY_TYPE" == firth ]] ;  then
		BINARY_GLM_FLAG="$BINARY_GLM_FLAG firth"
	fi
fi


covar_names=$(head -n 1 "$PHENO_FILE" | cut -f 4-)
"$PLINK_EXECUTABLE" \
    --pheno "$PHENO_FILE" \
    --no-psam-pheno \
    --pheno-name "$PHENOTYPE" \
    --pfile "$P_FILE" \
    --chr "$CHROM" \
    $BED_FILE_COMMAND \
    --mac 20 \
    --glm omit-ref pheno-ids hide-covar $COLS $BINARY_GLM_FLAG \
    $(if [[ -n "$covar_names" ]] ; then echo "--covar-name $covar_names" ; else echo "allow-no-covars" ; fi) \
    --ci 0.99999995 \
    --memory 56000 \
    --threads 28 \
    > plink.stdout \
    2> plink.stderr

cd -

mv "$temp_name"/* .

cd "$start_dir"

