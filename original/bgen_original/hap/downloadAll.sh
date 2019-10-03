#!/bin/bash
#Requires the ukbgene tool in this directory
#And .ukbkey in this directory - the credential used for downloading

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;	do 
	./ukbgene hap -c$i 
done

./ukbgene hap c1 -m
