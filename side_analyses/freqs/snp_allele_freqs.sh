#!/bin/bash

echo PFILE "$PFILE"
echo PFILE "$PFILE" >&2

$PLINK_COMMAND \
	--pfile "$PFILE" \
	--freq counts cols=+pos \
	--keep "$SAMPLE_FILE" \
	--memory 8000 \
	--threads 4

