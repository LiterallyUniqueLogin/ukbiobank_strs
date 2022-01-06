awk -F"\t" '{if ((NR == 1) || (($9 >= .8) && ($11 >= .8) && ($11 != "NA"))) { print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $9 "\t" $11 "\t" $8 "\t" $15} }' post_finemapping/results/validated/putatively_causal_STRs.tab | grep -P 'p_val|[^,]\t(-|\+),' | grep -P '\-, -$|\+, \+$' | grep -vP '\t\+\t[^\t]*-' | grep -vP '\t-\t[^\t]*\+' | tabcolumn

# has 4 or more replicates, SuSiE >= .3
awk -F"\t" '{if ((NR == 1) || (($9 >= .8) && ($11 >= .3) && ($11 != "NA"))) { print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $9 "\t" $11 "\t" $8 "\t" $15} }' post_finemapping/results/validated/putatively_causal_STRs.tab | grep -P 'p_val|[^,]\t(-|\+),' | awk -F"\t" '{ val=$8;  if ((NR == 1) || (gsub(/-/,"",val) >=4) || (gsub(/+/,"",val) >=4)) { print }}' | grep -vP '\t\+\t[^\t]*-' | grep -vP '\t-\t[^\t]*\+'
