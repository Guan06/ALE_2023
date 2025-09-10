#!/bin/bash

cat raw_vcf.lst |while read file;
do
	cd mincov10_C3_F001/$file
	`less snps.vcf | grep -v '##' |awk -F '\t' '{print $2"\t"$6}' >snps.vcf.qual`
	cd ../../
done
