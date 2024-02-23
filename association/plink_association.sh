#!/usr/bin/env bash

set -x

if [ -z "$PHENOTYPE" ] ; then echo "PHENOTYPE is unset" ; exit 1 ; fi
if [ -z "$IS_BINARY" ] ; then echo "IS_BINARY is unset" ; exit 1 ; fi
if [ -z "$FIRTH" ] ; then echo "FIRTH is unset" ; exit 1 ; fi
if [ -z "$CHROM" ] ; then echo "CHROM is unset" ; exit 1 ; fi
if [ -z "$OUT_DIR" ] ; then echo "OUTDIR is unset" ; exit 1 ; fi
if [ -z "$PHENO_FILE" ] ; then echo "PHENO_FILE is unset" ; exit 1 ; fi
if [ -z "$PLINK_EXECUTABLE" ] ; then echo "PLINK_EXECUTABLE is unset" ; exit 1 ; fi
if [ -z "$P_FILE" ] ; then echo "P_FILE is unset" ; exit 1 ; fi
if [ -z "$PROJECT_TEMP" ] ; then echo "PROJECT_TEMP is unset" ; exit 1 ; fi
echo "PHENOTYPE $PHENOTYPE IS_BINARY $IS_BINARY FIRTH $FIRTH CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END PROJECT_TEMP $PROJECT_TEMP VARS_FILE $VARS_FILE"
echo "PHENOTYPE $PHENOTYPE IS_BINARY $IS_BINARY FIRTH $FIRTH CHROM $CHROM OUT_DIR $OUT_DIR PHENO_FILE $PHENO_FILE PLINK_EXECUTABLE $PLINK_EXECUTABLE P_FILE $P_FILE START $START END $END PROJECT_TEMP $PROJECT_TEMP VARS_FILE $VARS_FILE" >&2

if [[ ( -n "$START" && -z "$END" ) || ( -z "$START" && -n "$END" ) ]] ; then
	echo "Start but not end or end but not start. Exiting."
	exit 1
fi

if [[ ( $IS_BINARY != true ) && ( $IS_BINARY != false ) ]] ; then
	echo "IS_BINARY not set to either 'true' or 'false'. Exiting."
	exit 1
fi

if [[ ( $FIRTH != true ) && ( $FIRTH != false ) ]] ; then
	echo "FIRTH not set to either 'true' or 'false'. Exiting."
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
	BED_FILE_COMMAND=" --extract-intersect bed1 region.bed"
fi

if [ -z "$VARS_FILE" ] ; then 
	VARS_FILE_COMMAND=" "
else
	tail -n +2 "$VARS_FILE" | awk '{print $1 "\t" $2 "\t" ($2 + 1)}' > vars_file.bed
	VARS_FILE_COMMAND=" --extract-intersect bed1 vars_file.bed"
fi

if [[ "$IS_BINARY" != true ]] ;  then
	COLS="cols=-test,-nobs"
	BINARY_GLM_FLAG=""
else
	COLS="cols=+gcountcc,-nobs,-test"
	#BINARY_GLM_FLAG="single-prec-cc cc-residualize"
	BINARY_GLM_FLAG=""
	if [[ "$FIRTH" == true ]] ;  then
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
    $VARS_FILE_COMMAND \
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

