#!/bin/bash

#/mnt/d/ubuntu/programs/stringtie-2.2.1.Linux_x86_64/stringtie -L -G /mnt/i/data/genomes/NC_045512.2.gtf -p 4 -f 0.0 -t -g 20 -o /mnt/i/data/SARS-CoV2/mapped_v6/stringtie.gtf /mnt/i/data/SARS-CoV2/mapped_v6/.bam/*.bam 
#/mnt/d/ubuntu/programs/stringtie-2.2.1.Linux_x86_64/stringtie -L -p 4 -f 0.0 -t -g 20 -o /mnt/i/data/SARS-CoV2/mapped_v6/stringtie.noref.gtf /mnt/i/data/SARS-CoV2/mapped_v6/.bam/*.bam 

/mnt/d/ubuntu/programs/stringtie-2.2.1.Linux_x86_64/stringtie -L -G /mnt/i/data/genomes/NC_045512.2.gtf -p 4 -f 0.0 -t -g 20 -o /mnt/i/data/SARS-CoV2/mapped_v6/stringtie.pych.gtf mapped_v6/mapv5.filt.pychopped/.bam/*.bam
/mnt/d/ubuntu/programs/stringtie-2.2.1.Linux_x86_64/stringtie -L -p 4 -f 0.0 -t -g 20 -o /mnt/i/data/SARS-CoV2/mapped_v6/stringtie.noref.pych.gtf /mnt/i/data/SARS-CoV2/mapped_v6/stringtie.gtf mapped_v6/mapv5.filt.pychopped/.bam/*.bam
