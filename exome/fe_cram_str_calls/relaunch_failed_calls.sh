#!/usr/bin/env bash

for file in str_beds/* ; do
	STR=$(echo $file | cut -f2- -d/ | cut -f1 -d.) BATCH=1 ./check_str_calls.sh > /dev/null
	if [[ $? != 0 ]] ; then
		./relaunch_call.sh "$STR"
	fi
done
