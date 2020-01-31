#!/bin/sh

grep -h -b100 -v INPUT $UKB/str_imputed/$1/batches/output/*.e* | vim -

