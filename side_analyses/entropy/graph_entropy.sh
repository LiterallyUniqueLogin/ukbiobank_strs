# $1 should be 'combined' or 'combined_het' (could also be 'full_genome' or 'exome' in
# place of 'combined')
# $2 should be entropy or het

source "$HOME"/.bashrc
cd "$UKB/../trtools/repo" || { echo "Directory $UKB/../trtools/repo doesn't exist" ; exit 1 ; }
conda activate ukb
python -m trtools.statSTR.statSTR \
	--tabfiles $(echo "$UKB"/side_analyses/entropy/"$1"/chr*.tab | sed -e 's/ /,/g') \
	--$2 \
	--out ../../ukbiobank/side_analyses/entropy/"$1"/all_chrs \
	--plot-dists smooth
conda deactivate

