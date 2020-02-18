#!/bin/sh

if [ -z "$1" ] ; then
        echo "Didn't give first argument - should be run name"
        exit -1
fi

#Look for any line that isn't the ouput of a time command
#Isn't just a log of the inputs to the script
#And isn't a blank line
for file in $UKB/str_imputed/runs/first_pass/batches/output/*.e* ; do
	if [ ! -z "$(grep -l -v -P '(INPUT)|(RUN_NAME)|(real\t[0-9]+m[0-9.]+s)|(user\t[0-9]+m[0-9.]+s)|(sys\t[0-9]+m[0-9.]+s)|(^$)' $file)" ] ; then
		echo
		echo "Error in file $file :"
		grep -n -v -P '(INPUT)|(real\t[0-9]+m[0-9.]+s)|(user\t[0-9]+m[0-9.]+s)|(sys\t[0-9]+m[0-9.]+s)|(^$)' $file
	fi
done 

