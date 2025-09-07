#!/bin/bash

cat raw_vcf.lst |while read file;
do
    cd mincov10_C3_F001/$file

    bcftools view --include 'FMT/DP>=20 && (FMT/AO)/(FMT/DP)>=0' snps.raw.vcf  | vt normalize -r reference/ref.fa - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > snps.filt.vcf
    
    snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr -c reference/snpeff.config -dataDir . ref snps.filt.vcf > snps.vcf

    /rds-d6/user/rg684/hpc-work/.conda/envs/snippy/bin/snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf snps.vcf > snps.tab

    /rds-d6/user/rg684/hpc-work/.conda/envs/snippy/bin/snippy-vcf_extract_subs snps.filt.vcf > snps.subs.vcf

    bcftools convert -Oz -o snps.vcf.gz snps.vcf
    
    ## Annotate the 'rm' lines to keep all the raw data
    #rm -f snps.subs.vcf.gz snps.subs.vcf.gz.csi snps.subs.vcf.gz.tbi

    samtools fastq -f 12 -v 20 --threads 10 -c 5 -N -s snps.unmapped_SE.fq.gz -0 snps.unmapped_R0.fq.gz -1 snps.unmapped_R1.fq.gz -2 snps.unmapped_R2.fq.gz snps.bam

    samtools view -h -q 60 snps.bam | samtools sort -l 0 -T /tmp --threads 9 -m 800M > /tmp/snippy.117515.Q60.bam

    samtools index /tmp/snippy.117515.Q60.bam

    /rds-d6/user/rg684/hpc-work/.conda/envs/snippy/bin/snippy-vcf_report --cpus 20 --bam /tmp/snippy.117515.Q60.bam --ref reference/ref.fa --vcf snps.vcf > snps.report.txt
    #rm -f /tmp/snippy.117515.Q60.bam /tmp/snippy.117515.Q60.bam.bai
    #rm snps.bed ref.fa ref.fa.fai snps.consensus.subs.fa snps.vcf.gz.csi
    #rm snps.txt snps.html snps.gff 

    cd ../../
done
