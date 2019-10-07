source ~/.bashrc
conda activate bcftools
for i in {1..22} ; do
	INPUT1=$i
	#source ./convert.sh
	#bgzip -c < chr$INPUT1.vcf > chr$INPUT1.vcf.gz
	tabix chr$INPUT1.vcf.gz
done
conda deactivate

