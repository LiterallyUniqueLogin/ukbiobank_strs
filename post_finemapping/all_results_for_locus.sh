phenotypes=$(ls finemapping/finemap_results/ | grep -v height)
{ printf phenotype$'\t' ; head -n 1 association/results/white_blood_cell_count/my_str/results.tab ; for phenotype in $phenotypes ; do printf ${phenotype}$'\t' ; grep -P '^22\t43385872' association/results/${phenotype}/my_str/results.tab ; done } | tabcolumn | vim -
