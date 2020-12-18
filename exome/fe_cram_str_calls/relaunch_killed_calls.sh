for file in $(grep -il 'killed' output/*stderr) ; do
	STR=$(echo "$file" | cut -f2 -d/ | sed -e 's/_batch.*//')
	BATCH=$(echo "$file" | cut -f2 -d/ | grep -oP '[0-9]+\.' | cut -f1 -d.)
	./relaunch_call.sh "$STR" "$BATCH"
done
