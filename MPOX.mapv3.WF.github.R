##### Import libraries and functions
library(hrbrthemes)
library(ggsci)
library(seqinr)
library(Rsamtools)
library(tidyverse)
library(stringi)
library(ggpubr)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)

source('functions/cigar.sum.R')
source('functions/misc.functions.R')
source('functions/ov.from.bam2.R')
source('functions/get.best.aln.R')
source('functions/filter.and.import.bams.R')
source('functions/feature.OV.from.polyC.TR.R')
source('functions/idxstats.bams.R')
source('functions/mapping.eval.R')

bam.flags <- read.delim('functions/bam.flags.tsv')


##### Metadata ####
metadata     <- read.delim('MPOX_metadata.txt')
metadata$sample_id <- gsub('.fastq', '', metadata$file)
metafilt     <- metadata[,c("sample_id", "hpi")]
colnames(metafilt)[1] <- 'sample'
metafilt$hpi  <- factor(metafilt$hpi, levels=c('0h', '1h', '2h' ,'4h', '6h', '12h', '24h', 'MIX'))
metafilt$Time <- as.integer(gsub('h', '', metadata$hpi))
metafilt$libtype <- 'cDNA'
metafilt$libtype[metafilt$hpi == 'MIX'] <- 'dRNA'

## which columns of the metadata table to be included in the downstream analysis?
## The first element of this vector should be sample ID, which should be the same as the .bam file's names.
metacols <- c('sample', 'hpi', 'Time', 'libtype')
#### ####
##

#### Settings ####
## Directory for output tables and figures
outdir <- paste0('output')
dir.create(outdir)


## Miscallenaous
by  <- c('seqnames', 'start', 'end')
palette <- pal_npg()(10)

## Reference genome
genome <- 'ON563414.3'
fasta  <- seqinr::read.fasta(paste(genome, '.fasta', sep = ''))
l_genome <- length(fasta[[1]])

## Annotation
gff          <- as.data.frame(rtracklayer::import.gff('ON563414.3.gff3'))
CDS.df       <- gff[gff$type == 'CDS', c("seqnames", "start", "end", "strand", "type", "Name", "Parent", 'ID')]
CDS.df$ID    <- gsub('gene-MPXV-USA_2022_', '', as.character(unlist(CDS.df$Parent)))
CDS.df$gene  <- as.character(unlist(CDS.df$Parent))

## generate feature table
feature.df <- CDS.df 
feature.colname <- 'ID' 
feature.df <- feature.df[, c(feature.colname, "strand", by, 'gene', 'Name')]

stopifnot(nrow(feature.df) == luniq(feature.df[,feature.colname]))
write.table(feature.df, 'feature.df.tsv', sep = '\t', row.names = F, quote = F)
#### ####
##


######  Import .bam files      ####

## load already saved data or carry out analysis?
load <- F
## filtered or unfiltered bam files?
filtered <- F
##
if (filtered) {
  save.data <- 'MPOX.filt.mapv3.RData'
  if (load) {message('loading ', save.data, '...'); load(save.data)} else {
    bamdir          <- '../mapping/bam.filt'
    source('import.filt.bams.R')
  } 
  } else {
    save.data <- 'MPOX.all.mapv3.RData'
    if (load) {message('loading ', save.data, '...'); load(save.data)} else {
      bamdir          <- '../mapping/bam'
      outfilt         <- '../mapping/bam.filt'
      source('import.all.bams.R')
      
  }
}
#length(unique(bam.filt$qname));length(unique(bam.prep.filt$qname))

#########


#### Subsequent process .bam files?
if (!load) {
    ## give alignment IDs _1 and _2 for secondary alignments!!!
    bam.filt  <- get.best.aln(bam.filt,#[bam.prep$sample == 'dRNA',],
                              bam.flags, best.mapq = F, rm.supplementary = F, rm.secondary = F, keep.chim = F, 
                              give.qname.ID.to.secondary = T)
    bam.filt$qname <- bam.filt$aln_ID
    bam.filt <- bam.filt[,colnames(bam.filt)]
} else {
  bam.filt  <- bam.filt
}
## NOTE::
## in bam.filt, the alignments now have unique identifiers. 
## bam.prep contains the original alignments and IDs
  
  

#### Calculate coverage?
calc.cov <- T
if (calc.cov) {
  window_size <- 100
  window_step <- 100
  whole.aln   <- F
  source('coverage.R')
  colnames(cov.sum)[1] <- 'Hours past infection'
  write_tsv(cov.sum, paste0(outdir, '/viral.coverage.stats.tsv'))
  rm(bam.cov)
}



#### Carry out plottings?
make.plots <- T
if (make.plots) {
    
  #### Plot read lengths
  source('read.length.plot.R')
  rm(ggv, ggh, ggd, ggvs, ggvh, plot.data)
  gc()
  
  
  #### Plot read counts
  source('read.count.plot.R')
  rm(ggcomb, ggrc, ggvr)
  gc()
  
  
  ##### Circos plots
  dir.create(paste0(outdir,'/circos.dRNA'))
  source('circos.dRNA.R')
  dir.create(paste0(outdir,'/circos.cDNA'))
  source('circos.cDNA.bothstrands.R')
  

}



##### Save data
if(!is.na(save.data)) {
  save.image(save.data)
}


##### Filter to viral reads 
rm(bam.all, reads.all)
save.data <- 'MPOX.filt.RData'
if(!is.na(save.data)) {
  save.image(save.data)
}


#### ####
##

