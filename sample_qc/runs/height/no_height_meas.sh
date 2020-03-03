#!/bin/bash

echo "ID" > remove/no_height.sample
cut -f2-3 "$UKB/main_dataset/height.txt" | tail -n+2 | awk '{ if (NF == 1) print $0 }' >> remove/no_height.sample

