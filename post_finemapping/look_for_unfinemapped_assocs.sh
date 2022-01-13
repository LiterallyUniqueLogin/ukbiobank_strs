# replace 19 and the pos with the correct chrom, pos
{ printf phenotype$'\t' ; head -n 1 association/results/white_blood_cell_count/my_str/results.tab ; for phenotype in $phenotypes ; do printf ${phenotype}$'\t' ; grep -P '^19\t17252977' association/results/${phenotype}/my_str/results.tab ; done } | tabcolumn | vim -
