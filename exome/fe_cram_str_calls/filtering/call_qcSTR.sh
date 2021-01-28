
if [ -z "$1" ] ; then echo "Need argument 1: absolute input vcf loc" ; exit 1 ; fi
if [ -z "$2" ] ; then echo "Need argument 2: absolute output loc prefix" ; exit 1 ; fi

source ~/.bashrc ;
conda activate ukb ;
cd "$UKB"/../trtools/repo || { echo "Can't cd properly" 1>&2 ; exit 1 ; }
python -m trtools.qcSTR.qcSTR \
        --vcf "$1"  \
        --vcftype hipstr \
        --out "$2"  \
        --quality per-sample \
        --quality per-call \
        --quality per-locus \
	--outtype png

