#!/bin/sh

#List the participants that are already marked as withdrawn in the 
#downloaded sample file (by sampleID <= 0)
head -n 2 "$UKB"/microarray/*.sample > "$UKB/sample_qc/common_filters/remove/withdrawn_already.sample"
grep -e - "$UKB"/microarray/*.sample >> "$UKB/sample_qc/common_filters/remove/withdrawn_already.sample"
