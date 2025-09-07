
This folder contains the scripts used from processing raw sequencing data to final results presented in Figure 2 and corresponding EV figures.

All the raw data and intermediate results are stored in:
/rds-d6/user/rg684/hpc-work/20250827_ALE.
	--__00.data_raw__
	This folder contains soft link to all raw sequencing data; all data have been uploaded to ENA
	--__01.data_preprocessing__
	Filter data with fastp 0.23.2 to get clean data for mapping
	--__02.call_variant_snippy__
	Call variant with snippy pipeline combined with vcf re-filtering (re_filter_raw_vcf.sh)
		--mincov10_C3_F001
        Each sample has an individual folder with all the intermediate files

	--__03.merge_vcf__
	Merge filtered vcf from all samples into one big file for downstream processing
