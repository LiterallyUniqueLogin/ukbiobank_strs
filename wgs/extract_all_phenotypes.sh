#!/bin/bash

# 3 instance traits
for ID in 30000 30010 30020 30030 30040 30050 30060 30070 30080 30090 30100 30110 30120 30140 30150 30180 30200 30210 30270 ; do
  dx run table-exporter \
    --folder /imputed_strs_paper/main_dataset/extracted_data \
    -ioutput="${ID}" \
    -idataset_or_cohort_or_dashboard=/app46122_20220823045256.dataset \
    -ioutput_format=TSV \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ientity=participant \
    -ifield_names=eid \
    $(for i in $(seq 0 2); do echo "-ifield_names=p${ID}_i${i}" ; done) &
done

# 2 instance traits
for ID in 30600 30610 30620 30630 30640 30650 30670 30680 30690 30700 30710 30720 30730 30740 30760 30770 30780 30810 30830 30840 30860 30870 30880 30890 30750 ; do
  dx run table-exporter \
    --folder /imputed_strs_paper/main_dataset/extracted_data \
    -ioutput="${ID}" \
    -idataset_or_cohort_or_dashboard=/app46122_20220823045256.dataset \
    -ioutput_format=TSV \
    -icoding_option=RAW \
    -iheader_style=UKB-FORMAT \
    -ientity=participant \
    -ifield_names=eid \
    $(for i in $(seq 0 1); do echo "-ifield_names=p${ID}_i${i}" ; done) &
done
