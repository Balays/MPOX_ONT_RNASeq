#!/bin/bash


python /mnt/d/ubuntu/programs/SQANTI3/sqanti3_qc.py \
	/mnt/i/data/SARS-CoV2/mapped_v6/TALON_mapv6_pych_cDNA.gtf /mnt/i/data/genomes/NC_045512.2.gtf /mnt/i/data/genomes/NC_045512.2.fasta \
	-o lortia_mapv6 \
	-d /mnt/i/data/SARS-CoV2/mapped_v6/SQUANTI3/lortia_mapv6 \
	--cpus 4 \
	--report both \
	--sites GTAG,GCAG,ATAC,ACAC
