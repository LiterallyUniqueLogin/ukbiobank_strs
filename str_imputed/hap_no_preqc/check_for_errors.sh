#!/bin/sh

grep -h -b100 -v INPUT $UKB/str_imputed/hap_no_preqc/vcf_batches/output/*.e* | vim -

