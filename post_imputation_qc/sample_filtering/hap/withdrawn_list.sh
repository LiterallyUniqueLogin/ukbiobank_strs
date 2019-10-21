#!/bin/sh

head -n 2 $UKB/original/bgen_original/hap/*sample > remove/withdrawn.sample
grep '-' $UKB/original/bgen_original/hap/*sample >> remove/withdrawn.sample
