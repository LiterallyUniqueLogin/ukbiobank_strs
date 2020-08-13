# $1 should be 'exome' or 'full_genome'

source "$HOME"/.bashrc
cd "$UKB/../trtools/repo" || { echo "Directory $UKB/../trtools/repo doesn't exist" ; exit 1 ; }
conda activate ukb
python -m trtools.statSTR.statSTR \
	--tabfiles $(echo "$UKB"/side_analyses/entropy/"$1"/chr*.tab | sed -e 's/ /,/g') \
	--entropy \
	--out ../../ukbiobank/side_analyses/entropy/"$1"/all_chrs \
	--plot-dists smooth
conda deactivate

