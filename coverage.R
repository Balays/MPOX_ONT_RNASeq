
source('functions/cov.from.bam.R')


### Calculate coverage per bp for mapped regions
bam.cov         <- merge(bam.filt, metafilt, by='sample', all=F)
bam.cov$sample  <- bam.cov$hpi
mapped.cov      <- cov.from.bam(bam.cov, samples = unique(bam.cov$sample), fasta.ref = 'ON563414.3.fasta')

### Summarise coverage per window
win.cov.sum <- window.cov(mapped.cov, 'ON563414.3.fasta', 
                          window_size = window_size, 
                          window_step = window_step
                          )

### Coverage summary
cov.sum     <- mapped.cov %>% group_by(sample) %>% summarise(mean=mean(count), median=median(count), sd=sd(count))

gc()

if (whole.aln) {
  ### Import bamfiles (rm.gaps.in.aln = F !!!)  
  bam.reads <- import.bams(bamfiles, bamnames, write.filtered=write.filtered, force.create=force.create, 
                           rm.gaps.in.aln = F, filtering=filtering, mapq.filt = mapq.filt,
                           flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                           is.lortia = is.lortia)
  
  
  ### Calculate coverage per bp for whole alignments
  aln.cov <- cov.from.bam(bam.reads, samples = unique(bam.reads$sample), fasta.ref = 'ON563414.3.fasta')
  
}

