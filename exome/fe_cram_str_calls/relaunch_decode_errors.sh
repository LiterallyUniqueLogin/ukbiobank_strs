for file in $(grep -l decode output/*stderr) ; do
	CUT_SAMPLES=$(head -n 1 "$file" | awk '{print $6}')
	if [ ! -z "$CUT_SAMPLES" ] ; then
		CUT_SAMPLES="$CUT_SAMPLES"-
	fi
	STR=$(echo "$file" | cut -f2 -d/ | sed -e 's/_batch.*//')
	BATCH=$(echo "$file" | cut -f2 -d/ | grep -oP '[0-9]+\.' | cut -f1 -d.)
	CUT_SAMPLES="$CUT_SAMPLES"$(grep encountered "$file" | grep -oP '[0-9]+\.cram' | cut -f1 -d.)
	echo ./relaunch_call.sh "$STR" "$BATCH" "$CUT_SAMPLES"
done
