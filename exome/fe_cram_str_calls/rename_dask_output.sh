cd output || { echo "Couldn't cd output" ; exit 1 ; } 

for file in call_exome_STRs.* ; do
	str=$(head -n 1 "$file" | awk '{print $4}')
	batch=$(head -n 1 "$file" | awk '{print $2}')
	mv "$file" "${str}_batch_${batch}.dask.$(echo "$file" | cut -f2 -d.)"
done
