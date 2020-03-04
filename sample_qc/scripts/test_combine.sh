#!/bin/bash

TMPDIR=$UKB/temp
RUNDIR=$UKB/sample_qc/runs/height
COMMONDIR=$UKB/sample_qc/common_filters

sorted_samples="$TMPDIR/combined.sample"
sort "$RUNDIR"/combined.sample > "$sorted_samples"

for fullfile in "$COMMONDIR"/remove/* "$RUNDIR"/remove/*; do
	file=$(basename "$fullfile")
	sort "$fullfile" > "$TMPDIR/$file"
	count=$(comm -12 "$sorted_samples" "$TMPDIR/$file" | grep -cv '^ID')
	if (( count != 0 )) ; then
		echo "Not all IDs from file $file were removed" >&2
		echo "comm output: $(comm -12 "$sorted_samples" "$TMPDIR/$file" | grep -v '^ID')"
		exit 1
	fi
done

for fullfile in "$COMMONDIR"/keep/* ; do
	file=$(basename "$fullfile")
	sort "$fullfile" > "$TMPDIR/$file"
	count=$(comm -23 "$sorted_samples" "$TMPDIR/$file" | grep -cv '^ID')
	if (( count != 0 )) ; then
		echo "Not all IDs were contained in file $file" >&2
		echo "comm output: $(comm -23 "$sorted_samples" "$TMPDIR/$file" | grep -v '^ID')"
		exit 1
	fi
done

echo "Success"
