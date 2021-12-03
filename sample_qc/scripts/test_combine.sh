#!/bin/bash

if [ -z "$PHEN" ] ; then echo "PHEN is unset" ; exit 1 ; fi
if [ -z "$ETHNICITY" ] ; then echo "ETHNICITY is unset" ; exit 1 ; fi
echo "PHEN $PHEN"
echo "PHEN $PHEN" >&2

TMPDIR=$UKB/temp
RUNDIR=$UKB/sample_qc/runs/"$ETHNICITY"/"$PHEN"
COMMONDIR=$UKB/sample_qc/common_filters

sorted_samples="$TMPDIR/combined.sample"
sort "$RUNDIR"/combined.sample > "$sorted_samples"

for fullfile in "$COMMONDIR"/remove/*sample ; do
	file=$(basename "$fullfile")
	awk '{print $1}' "$fullfile" | sort > "$TMPDIR/$file"
	count=$(comm -12 "$sorted_samples" "$TMPDIR/$file" | grep -cv '^ID')
	if (( count != 0 )) ; then
		echo "Not all IDs from file $file were removed" >&2
		echo "comm output: $(comm -12 "$sorted_samples" "$TMPDIR/$file" | grep -v '^ID')"
		exit 1
	fi
done

fullfile="$COMMONDIR"/ethnicity/"$ETHNICITY".sample
file=$(basename "$fullfile")
awk '{print $1}' "$fullfile" | sort > "$TMPDIR/$file"
count=$(comm -23 "$sorted_samples" "$TMPDIR/$file" | grep -cv '^ID')
if (( count != 0 )) ; then
	echo "Not all IDs were contained in file $file" >&2
	echo "comm output: $(comm -23 "$sorted_samples" "$TMPDIR/$file" | grep -v '^ID')"
	exit 1
fi

echo "Combine was successful"
