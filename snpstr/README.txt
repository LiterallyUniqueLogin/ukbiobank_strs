Saini_etal_SuppTable3.txt is supplementary information assocaited
with SNP STR reference panel paper
( A reference haplotype panel for genome-wide imputation of short tandem repeats,
  Saini, et al,
https://www.nature.com/articles/s41467-018-06694-0 )
I downloaded it from here:
http://gymreklab.com/2018/03/05/snpstr_imputation.html
Opened the file in excel and saved the first sheet 
as a txt file
This table is referred to as supplementary dataset 2 (not supplementary table
3) in the paper link above.

igsr_samples.tsv contains information about the samples in this 
haplotype panel. It was downloaded from the 1000 genomes project
data portal main page (https://www.internationalgenome.org/data-portal/sample)
clicking on the "Download the list" button

eur.samples contains the sample ids of the european descent
samples in the 1000 genomes population. Produced via
awk 'BEGIN { FS = "\t" } {if ($6=="EUR") { print($1); } }' igsr_samples.tsv | \
	sort > eur.sample

