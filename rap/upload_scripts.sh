#!/bin/bash

# call from $UKB
# make sure to call from an interactive node, this spawns a zillion
# subprocesses
for file in \
	association/*py \
	export_scripts/*py \
	finemapping/*py \
	post_finemapping/*py \
	sample_qc/scripts/*py \
	signals/*py \
	traits/*py \
; do
	{ dx rm -a imputed_strs_paper/"$file" ;
	dx upload --path imputed_strs_paper/"$file" "$file" ; } &
done
