str_num=1
while read -r bedline ; do
	echo "$bedline" | sed -e 's/^chr//' > str_beds/"$str_num".bed
	str_num=$((str_num+1))
done < "$UKB"/side_analyses/exome_strs/snpstr_exome_strs_hg38.bed
