#!/bin/bash

export GZIP=-9

{ 
printf '%s\n' "$UKB/wdl_cache/finemapping/finemapping_all_regions_concordance_$1.tab.for_archive"
printf '%s\n' "$UKB/wdl_cache/finemapping/finemapping_followup_concordance_$1.tab.for_archive"
tail -n +2 "$UKB/signals/regions/$1.tab" | cut -f 1-3 | sed -e 's/\t/_/g' | while read region ; do 
	if [[ ( "$region" != "4_8165642_11717761" || "$1" != "urate" ) && ( "$region" != "12_19976272_22524428" || "$1" != "total_bilirubin" ) && ( "$region" != "1_19430673_24309348" || "$1" != "alkaline_phosphatase" ) ]] ; then
		printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_first_pass_${region}_finemap_output."{config,cred*,log_sss,snp}
		printf '%s\n' "$UKB/finemapping/susie_results/${1}/${region}/"{V.tab,alpha.tab,colnames.txt$(if [[ "$1" == mean_platelet_volume ]] ; then echo ".normal_run" ; fi),cs*.txt,lbf.tab,lbf_variable.tab,lfsr.tab,sigma2.txt} | grep -v files_list 
	fi
done
tail -n +2 "$UKB/wdl_cache/finemapping/${1}_followup_finemapping_regions.tsv" | cut -f 3 | while read region ; do 
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_derived_effect_size_prior_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_low_effect_size_prior_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_mac_threshold_100_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_prior_4_signals_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_prior_snps_over_strs_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_pval_threshold_1e4_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_FINEMAP_stricter_stopping_threshold_${region}_finemap_output."{config,cred*,log_sss,snp}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_SuSiE_prior_snps_over_strs_${region}_"{V.tab,alpha.tab,colnames.txt,cs*.txt,lbf.tab,lbf_variable.tab,lfsr.tab,sigma2.txt}
	printf '%s\n' "$UKB/wdl_cache/finemapping/${1}_SuSiE_best_guess_genotypes_${region}_"{V.tab,alpha.tab,colnames.txt,cs*.txt,lbf.tab,lbf_variable.tab,lfsr.tab,sigma2.txt}
done
} | \
	tar -czvf finemapping_"$1".tgz -T - --show-transformed-names \
	--transform='s,log_sss$,log,' \
	--transform='s,'"$1"'_,,' \
	--transform='s,_'"$1"',,' \
	--transform='s,all_regions_concordance,first_pass,' \
	--transform='s,_concordance,,' \
	--transform='s,.for_archive,,' \
	--transform='s,.normal_run,,' \
	--transform='s,.*wdl_cache/finemapping/,,' \
	--transform='s,.*/finemapping/susie_results/'"$1"'/\([0-9_]*\)/,SuSiE_first_pass_\1_,'
