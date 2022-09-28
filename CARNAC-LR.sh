#!/bin/bash

#fastqfrombam.sh .bam.lortia.all .bam.lortia.fastq
#cat .bam.lortia.fastq/*.fastq > all.lortia.fastq
#minimap2 all.lortia.fastq all.lortia.fastq -X > minimap_all.lortia.paf
#ulimit -s unlimited
#paf_to_CARNAC.py minimap_all.lortia.paf all.lortia.fastq CARNAC_all.lortia.txt
#CARNAC-LR -f CARNAC_all.lortia.txt -o CARNAC_all.lortia_out -t 4

####
#cat .bam.lortia.fastq/*.fastq > hpi1-hpi16.lortia.fastq

minimap2 hpi1-hpi16.lortia.fastq hpi1-hpi16.lortia.fastq -X > minimap_hpi1-hpi16.lortia.paf
ulimit -s unlimited
paf_to_CARNAC.py minimap_hpi1-hpi16.lortia.paf hpi1-hpi16.lortia.fastq CARNAC_hpi1-hpi16.lortia.txt
CARNAC-LR -f CARNAC_hpi1-hpi16.lortia.txt -o CARNAC_hpi1-hpi16.lortia_out -t 4


