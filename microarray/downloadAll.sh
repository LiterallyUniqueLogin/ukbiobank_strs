#!/bin/bash
#Set $UKB_UTILITIES appropriately
#ukbgene is the downloading tool and 
#k29170.key is the credential used 

for i in $(seq 1 22) ; do 
	$UKB_UTILITIES/ukbgene hap -c$i -a$UKB_UTILITIES/k29170.key
done

$UKB_UTILITIES/ukbgene hap c1 -m -a${UKB_UTILITIES}/k29170.key
