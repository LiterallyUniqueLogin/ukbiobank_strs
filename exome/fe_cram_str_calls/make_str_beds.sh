# beds are named by SNPSTR ID
# (Note: this is NOT equivalent to the STR with the same ID in the HipSTR
# reference panel)
while read -r bedline ; do
	name=$(echo "$bedline" | awk '{print $6}' | sed -e 's%/%_%g')
	if (( $(echo "$name" | wc -w) != 1 )) ; then
		echo "$name has multiple words in it!"
		exit 1
	fi
	echo "$bedline" | sed -e 's/^chr//' > str_beds/"$name".bed
done < "$UKB"/side_analyses/exome_strs/snpstr_exome_strs_38.bed
