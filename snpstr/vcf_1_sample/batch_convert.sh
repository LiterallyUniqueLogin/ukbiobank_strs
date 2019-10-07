source ~/.bashrc
conda activate bcftools
for i in {1..22} ; do
	#qsub -v "INPUT1=$i" convert.sh &
	bgzip -c < chr$i.vcf > chr$i.vcf.gz
	tabix chr$i.vcf.gz
done
conda deactivate
