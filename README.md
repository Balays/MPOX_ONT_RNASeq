# MPOX_ONT_RNASeq
The R workflow that was used to create the results and the figures in *Kakuk et al. In-depth Temporal Transcriptome Profiling of Monkeypox and Host Cells using Nanopore Sequencing. 2022*

After cloning the repository and installing required R packages, the *.bam* files should be downloaded from: 
https://www.ebi.ac.uk/ena/browser/view/PRJEB56841
and should be put into a folder '../mapping/bam' relative to this directory. This can be changed in the 'bamdir' parameter in the worklfow script.
Nothing else is required for the WF to be carried out. Although the scripts can be used to import, process, analyse and generate figures from any *.bam* files, after the appropiate parameters have been set. This mainly includes the metadata and the genome and its annotation.



