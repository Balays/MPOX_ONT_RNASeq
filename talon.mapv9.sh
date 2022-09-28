#!/bin/bash

###### mapv9

## lortia cDNA
#talon_initialize_database --f /mnt/i/data/genomes/NC_045512.2.gtf --g NC_045512.2 --a NC_045512.2 --idprefix mapv9_pych_cDNA --5p 10 --3p 10 --o TALON_mapv9_LoRTIA_cDNA
talon --f talon_mapv9_cDNA_lortia_os_meta.csv --db TALON_mapv9_LoRTIA_cDNA.db --build NC_045512.2 -t 1 --cov 0.9 --identity 0.8 --o TALON_mapv9_LoRTIA_cDNA
talon_create_GTF --db TALON_mapv9_LoRTIA_cDNA.db -a NC_045512.2 -b NC_045512.2 --o mapv9_LoRTIA_cDNA
agat_convert_sp_gxf2gxf.pl -g mapv9_LoRTIA_cDNA_talon.gtf -o TALON_mapv9_LoRTIA_cDNA.gff
talon_abundance --db TALON_mapv9_LoRTIA_cDNA.db -a NC_045512.2 -b NC_045512.2 --o mapv9_LoRTIA_cDNA


