for str in $(cat "$UKB"/side_analyses/exome_strs/snpstr_exome_str_cleaned_ids.txt) ; do
	for BATCH in $(seq 1 10) ; do
		check=$(BATCH=$BATCH STR=$str ./check_str_calls.sh)
		id="${str}_batch_$BATCH"
		if [[ -z "$check" ]] ; then
			continue
		else
			if [[ "$check" =~ .*exist.* ]] ; then
				missing="$missing $id"
			elif [[ "$check" =~ .*unzip.* ]] ; then
				unzip="$unzip $id"
			elif [[ "$check" =~ .*file[[:space:]]truncated.* ]] ; then
				truncated="$truncated $id"
			elif [[ "$check" =~ .*line[[:space:]]truncated.* ]] ; then
				line_truncated="$line_truncated $id"
			else
				echo "Found unexpected output from check command"
				exit 1
			fi
		fi
	
		hipstr_err_file=output/"$id".hipstr.stderr
		unset identified
		if [[ -n "$(grep 'many reads' "$hipstr_err_file")" ]] ; then
			many="$many $id"
			identified=yes
		fi
		if [[ -n "$(grep 'open files' "$hipstr_err_file")" ]] ; then
			files="$files $id"
			identified=yes
		fi
		if [[ -n "$(grep few "$hipstr_err_file")" ]] ; then
			few="$few $id"
			identified=yes
		fi
		if [[ -n "$(grep -i killed "$hipstr_err_file")" ]] ; then
			killed="$killed $id"
			identified=yes
		fi
		if [[ -n "$(grep decode "$hipstr_err_file")" ]] ; then
			decode="$decode $id"
			identified=yes
		fi
		if [[ -n "$(grep 'No spanning' "$hipstr_err_file")" ]] ; then
			spanning="$spanning $id"
			identified=yes
		fi
		if [[ -n "$(grep 'Stutter model training failed' "$hipstr_err_file")" ]] ; then
			training="$training $id"
			identified=yes
		fi
		if [[ -z "$identified" ]] ; then
			unidentified="$unidentified $id"
		fi
	done
done

function clean {
	echo "$1" | tr ' ' '\n' | sed -e 's/^$//' | sort
}

if [[ -n "$missing" ]] ; then
	echo "STRs with no VCFs:" 
	clean "$missing"
	echo
else
	echo "No STRs with missing VCFs"
	echo
fi

if [[ -n "$unzip" ]] ; then
	echo "STRs with VCFs that can't be unzipped:" 
	clean "$unzip"
	echo
else
	echo "No STRs with VCFs that can't be unzipped"
	echo
fi

if [[ -n "$truncated" ]] ; then
	echo "STRs with truncated VCFs:" 
	clean "$truncated"
	echo
else
	echo "No STRs with trunacted VCFs"
	echo
fi

if [[ -n "$line_truncated" ]] ; then
	echo "STRs with VCFs ending in truncated lines:" 
	clean "$line_truncated"
	echo
else
	echo "No STRs with VCFs ending in truncated lines"
	echo
fi

echo -e "-----------------------------------\n"

if [[ -n "$files" ]] ; then
	files_not_missing=$(comm -23 <(clean "$files") <(clean "$missing"))
	if [[ -z "$files_not_missing" ]] ; then
		echo "All STRs with too many open files have missing VCFs"
		echo "STRs with too many open files"
		clean "$files"
		echo
	else
		echo "STRs with too many open files and missing VCFs" 
		clean "$(comm -12 <(clean "$files") <(clean "$missing"))"
		echo "STRs with too many open files files and extant VCFs"
		clean "$files_not_missing"
		echo
	fi
else
	echo "No STRs with too many open files"
	echo
fi


if [[ -n "$many" ]] ; then
	many_not_trunc=$(comm -23 <(clean "$many") <(clean "$truncated"))
	if [[ -z "$many_not_trunc" ]] ; then
		echo "All STRs with too many reads have truncated VCFs"
		echo "STRs with too many reads"
		clean "$many"
		echo
	else
		echo "STRs with too many reads and truncated VCFs" 
		clean "$(comm -12 <(clean "$many") <(clean "$truncated"))"
		echo "STRs with too many reads but not truncated VCFs"
		clean "$many_not_trunc"
		echo
	fi
else
	echo "No STRs with too many reads"
	echo
fi

if [[ -n "$few" ]] ; then
	few_not_trunc=$(comm -23 <(clean "$many") <(clean "$truncated"))
	if [[ -z "$few_not_trunc" ]] ; then
		echo "All STRs with too few reads have truncated VCFs"
		echo "STRs with too few reads"
		clean "$few"
		echo
	else
		echo "STRs with too few reads and truncated VCFs" 
		clean "$(comm -12 <(clean "$few") <(clean "$truncated"))"
		echo "STRs with too few reads but not truncated VCFs"
		clean "$few_not_trunc"
		echo
	fi
else
	echo "No STRs with too few reads"
	echo
fi

if [[ -n "$spanning" ]] ; then
	spanning_not_trunc=$(comm -23 <(clean "$spanning") <(clean "$truncated"))
	if [[ -z "$spanning_not_trunc" ]] ; then
		echo "all strs with no spanning reads have truncated vcfs"
		echo "strs with no spanning reads"
		clean "$spanning"
		echo
	else
		echo "strs with no spanning reads and truncated vcfs" 
		clean "$(comm -12 <(clean "$spanning") <(clean "$truncated"))"
		echo "strs with no spanning reads but not truncated vcfs"
		clean "$spanning_not_trunc"
		echo
	fi
else
	echo "no strs with no spanning reads"
	echo
fi

if [[ -n "$training" ]] ; then
	training_not_trunc=$(comm -23 <(clean "$training") <(clean "$truncated"))
	if [[ -z "$training_not_trunc" ]] ; then
		echo "all strs where stutter model training failed have truncated vcfs"
		echo "strs where stutter model training failed"
		clean "$training"
		echo
	else
		echo "strs where stutter model training failed and truncated vcfs" 
		clean "$(comm -12 <(clean "$training") <(clean "$truncated"))"
		echo "strs where stutter model training failed but not truncated vcfs"
		clean "$training_not_trunc"
		echo
	fi
else
	echo "no strs where stutter model training failed"
	echo
fi

if [[ -n "$killed" ]] ; then
	killed_not_unzip=$(comm -23 <(clean "$killed") <(clean "$unzip"))
	if [[ -z "$killed_not_unzip" ]] ; then
		echo "All killed STRs not killed have zip-broken VCFs"
		echo "killed STRs"
		clean "$killed"
		echo
	else
		echo "killed STRs with zip-broken VCFs" 
		clean "$(comm -12 <(clean "$killed") <(clean "$unzip"))"
		echo "killed STRs with zip-working or nonexistant VCFs"
		clean "$killed_not_unzip"
		echo
	fi
else
	echo "No killed STRs"
	echo
fi

if [[ -n "$decode" ]] ; then
	decode_not_unzip=$(comm -23 <(clean "$decode") <(clean "$unzip"))
	if [[ -z "$decode_not_unzip" ]] ; then
		echo "All undecodable STRs have zip-broken VCFs"
		echo "Undecodable STRs"
		clean "$decode"
		echo
	else
		echo "Undecodable STRs with zip-broken VCFs" 
		clean "$(comm -12 <(clean "$decode") <(clean "$unzip"))"
		echo "Undecodable STRs with zip-working or nonexistant VCFs"
		clean "$decode_not_unzip"
		echo
	fi
else
	echo "No STRs with decode errors"
	echo
fi

if [[ -n "$unidentified" ]] ; then
	echo "STRs with unidentified errors"
	clean "$unidentified"
	echo
else
	echo "No STRs with unidentified errors"
	echo
fi

