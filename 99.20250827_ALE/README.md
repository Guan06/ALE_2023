
This folder contains all scripts required to process the raw sequencing data through to the final results shown in Figure 2 and the associated EV figures.  

All raw data and intermediate results are stored at:
`/rds-d6/user/rg684/hpc-work/20250827_ALE`.  

__00.data_raw__  
    └─ Symbolic links to the raw sequencing files (already deposited in ENA)
    
__01.data_preprocessing__  
    └─ Fastp v0.23.2 was used to filter the reads and generate clean data for mapping (all html files are uploaded here in this repo)  
    
__02.call_variant_snippy__  
    └─ Variant calling with the Snippy pipeline, followed by VCF re‑filtering (script: `re_filter_raw_vcf.sh`)  
    └─ Sub‑directory “mincov10_C3_F001” – each sample has its own folder  
       containing all intermediate files  
       
__03.merge_vcf__  
    └─ Filtered VCFs from all samples were merged into a single file for downstream analyses
