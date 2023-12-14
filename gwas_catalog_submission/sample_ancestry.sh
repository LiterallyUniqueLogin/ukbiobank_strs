#!/bin/bash

cd $UKB/wdl_cache/associations
for file in *gz ; do 
	population=$(echo "$file" | awk 'BEGIN {FS="_"} {printf $1; if ($1 == "white" || $1 == "south") { print "_" $2 }}')
	if [[ "$population" == "black" ]] ; then
		printf "African Unspeicifed\tAfrican and Caribbean\n"
	elif [[ "$population" == "chinese" ]] ; then
		printf "East Asian\tChinese\n"
	elif [[ "$population" == "irish" ]] ; then
		printf "European\tIrish\n"
	elif [[ "$population" == "south_asian" ]] ; then
		printf "South Asian\tIndian, Pakistani and Bangladeshi\n"
	elif [[ "$population" == "white_british" ]] ; then
		printf "European\tWhite British\n"
	elif [[ "$population" == "white_other" ]] ; then
		printf "European\tWhite non-Irish non-British\n"
	else
		echo $((1/0))
	fi
done
