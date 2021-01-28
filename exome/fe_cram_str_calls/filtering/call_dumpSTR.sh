
source ~/.bashrc ;
conda activate ukb ;
cd "$UKB"/../trtools/repo || { echo "Can't cd properly" 1>&2 ; exit 1 ; }
python -m trtools.dumpSTR.dumpSTR \
        --vcf "$UKB"/exome/fe_cram_str_calls/exome_calls.vcf.gz \
        --vcftype hipstr \
        --out "$UKB"/exome/fe_cram_str_calls/filtering/filtered_exome_calls \
        --hipstr-max-call-flank-indel 0.15 \
        --hipstr-max-call-stutter 0.15 \
        --hipstr-min-supp-reads 10 \
	--hipstr-min-call-allele-bias -2 \
	--hipstr-min-call-strand-bias -2 \
	--filter-hrun \
	--min-locus-hwep 0.00008517887 \
	--use-length \
	--drop-filtered \
	--zip

# hwep at 5e-2/n_STRs = 5e-2/587
# original run was with 630 in the denominator which was inaccurate
