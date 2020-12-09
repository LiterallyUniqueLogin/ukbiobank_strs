for id in $(cat "$UKB"/side_analyses/exome_strs/snpstr_exome_str_cleaned_ids.txt) ; do
	check=$(BATCH=1 STR=$id	./check_str_calls.sh)
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

	unset identified
	if [[ -n "$(grep many output/"$id"_*err)" ]] ; then
		many="$many $id"
		identified=yes
	fi
	if [[ -n "$(grep few output/"$id"_*err)" ]] ; then
		few="$few $id"
		identified=yes
	fi
	if [[ -n "$(grep killed output/"$id"_*err)" ]] ; then
		killed="$killed $id"
		identified=yes
	fi
	if [[ -n "$(grep decode output/"$id"_*err)" ]] ; then
		decode="$decode $id"
		identified=yes
	fi
	if [[ -n "$(grep 'No spanning' output/"$id"*err)" ]] ; then
		spanning="$spanning $id"
		identified=yes
	fi
	if [[ -z "$identified" ]] ; then
		unidentified="$unidentified $id"
	fi
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
		echo "All STRs with no spanning reads have truncated VCFs"
		echo "STRs with no spanning reads"
		clean "$spanning"
		echo
	else
		echo "STRs with no spanning reads and truncated VCFs" 
		clean "$(comm -12 <(clean "$spanning") <(clean "$truncated"))"
		echo "STRs with no spanning reads but not truncated VCFs"
		clean "$spanning_not_trunc"
		echo
	fi
else
	echo "No STRs with no spanning reads"
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

