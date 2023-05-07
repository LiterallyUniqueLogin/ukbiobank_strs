#!/bin/bash

# call from $UKB
# make sure to call from an interactive node, this spawns a zillion
# subprocesses
for file in \
	association/*{py,sh} \
	finemapping/*{py,sh} \
	post_finemapping/*{py,sh} \
	sample_qc/scripts/*{py,sh} \
	signals/*{py,sh} \
	traits/*{py,sh} \
	wgs/*{py,sh} \
; do
	{ dx rm -a imputed_strs_paper/"$file" ;
	dx upload --path imputed_strs_paper/"$file" "$file" ; } &
done
