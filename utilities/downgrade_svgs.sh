#!/usr/bin/env bash

# convert v2 svgs to v1 svgs (only fixing href tags)

for file in "$@" ; do
	if ! [[ $file =~ ^.*svg$  ]] ; then
		echo "Expected an svg file. Got $file" 1>&2
		exit 1
	fi
	sed -i "$file" -e 's_xmlns="http://www.w3.org/2000/svg"_\0 xmlns:xlink="http://www.w3.org/1999/xlink"_' -e 's/href/xlink:href/g'
done
