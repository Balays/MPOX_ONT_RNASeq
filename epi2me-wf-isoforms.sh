#!/bin/bash

OUTPUT=/mnt/i/data/SARS-CoV2/epi2me-wf-isoforms
TEMPDIR=/mnt/i/data/tmp
nextflow run wf-isoforms/ -profile conda --fastq /mnt/i/data/SARS-CoV2/SARS_COV2_fastq --ref_genome /mnt/i/data/genomes/NC_045512.2.fasta --ref_annotation /mnt/i/data/genomes/NC_045512.2.gtf --threads 4 \
--minimap2_opts '-ax splice -un -Y -C5 --MD -g 30000 -G 30000 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -t 4' --out_dir $OUTPUT -w $TEMPDIR -resume