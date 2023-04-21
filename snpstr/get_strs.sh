# get str_ids
cut -f 2 Saini_etal_SuppTable2.txt | tail -n +2 | sort | uniq > str_ids.txt
# get bed file
for i in $(seq 1 22) ; do bcftools query -i ID=@str_ids.txt -f '%CHROM\t%INFO/START\t%INFO/END\t%ID\n' info_field/chr${i}.vcf.gz ; done | awk '{print "chr" $1 "\t" $2 -1 "\t" $3 "\t" $4}' > str_loci.bed
# this automatically removes duplicate ids for some reason. All duplicate IDs have the same INFO/START and INFO/END though, so that's fine

# then used https://genome.ucsc.edu/cgi-bin/hgLiftOver on 03/09/2023 for lifting over to hg38 and T2T CHM13v2.0/hs1
# liftover to hg38 using INFO/START and INFO/END produces the same exact coordinates as HipSTR's hg38 reference
# https://github.com/HipSTR-Tool/HipSTR-references/ aside from an off by one representation difference
# this seems to preserve locus order
