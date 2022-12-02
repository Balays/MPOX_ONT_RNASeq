
## Alignment import settings
seqnames.tofilt <- genome ## NA ##Filter for alignments that mapped to this contig (viral genome)
flag.tokeep     <- NA # This will not drop alignments in .bam files based on bamflag only
#c(0,16) ## this will leave primary alignments only and filter out supplementary and secondary alignments as well
mapq.filt       <- NA ## Mapping quality threshold (in minimap2 60 is the highest)
write.filtered  <- T ## Write out filtered .bam files
rm.gaps.in.aln  <- T
force.create    <- T ## Overwrite?
filtering       <- F ## filter out supplementary alignments from bamfiles
filter.bams     <- F ## subsequent filtering: get the best alignment for each read (to make qnames unique)
flag.tokeep     <- NA
flag.tocrop     <- NA


## 
is.lortia       <- F

##


pattern   <- '_pass.bam'
bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]

## separate dRNA for now
bamfiles <- bamfiles[grep('dRNA', bamfiles, invert = T)]

nbam <- length(bamfiles)
bamnames <- gsub('.*\\/', '', bamfiles)
bamnames <- gsub(pattern, '', bamnames)

#### Correcting sample names to match with metadata
bamnames  <- gsub('MPOX_inf_',  '', bamnames)
bamnames  <- gsub('MPOX_', '', bamnames)
bamnames  <- gsub('_inf', '', bamnames)

#### Import all reads from the .bam files #####
## Infected
bam.inf  <- import.bams(bamfiles[grep('mock', bamfiles, invert=T) ], # [1], # 
                        bamnames[grep('mock', bamnames, invert=T) ], # [1], # 
                        write.filtered = F, force.create=force.create, #add.primes=F,
                        rm.gaps.in.aln = rm.gaps.in.aln, filtering= F, mapq.filt = mapq.filt,
                        flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = NA, 
                        is.lortia = is.lortia)

## Mock
bam.mock <- import.bams(bamfiles[grep('mock', bamfiles, invert=F) ], 
                        bamnames[grep('mock', bamnames, invert=F) ],
                        #  '../mapping/pychopped_pass_v3/bam/MPOX_mock_A_pass_pychopped.bam', 'mock_A', 
                        write.filtered = F, force.create=force.create, 
                        rm.gaps.in.aln = rm.gaps.in.aln, filtering=F, mapq.filt = mapq.filt,
                        flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = NA,
                        is.lortia = is.lortia)


## dRNA
bamfile.drna <- '../mapping/v3/bam/MPOX_dRNA_inf_pass.bam'; bamname.drna <- 'dRNA'
bam.drna <- import.bams(bamfile.drna, 
                        bamname.drna,
                        write.filtered = F, force.create=force.create, 
                        rm.gaps.in.aln = rm.gaps.in.aln, filtering=F, mapq.filt = mapq.filt,
                        flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = NA,
                        is.lortia = is.lortia)

bamfiles <- c(bamfiles, bamfile.drna); bamnames <- c(bamnames, bamname.drna)
gc()
bam.all <- plyr::rbind.fill(bam.inf, bam.mock, bam.drna)
rm(bam.inf, bam.mock, bam.drna)
gc()
#### ####
##

#### Summarise seqnames and samples ####
all.stats <- idxstats.bams(bamfiles, bamnames = bamnames )

for (i in unique(all.stats$sample)) {
  stats <- all.stats[all.stats$sample == i, ]
  try({ all.stats$unmapped[all.stats$sample == i] <- stats$unmapped[stats$seqnames == '*'] })
}

all.stats <- all.stats[all.stats$seqnames != '*', ]
all.stats$ratio <- all.stats$mapped / all.stats$unmapped

bam.sum <- bam.all %>% group_by(sample, seqnames) %>% summarise(sub_alignment_count=n())
#### ####
##


#### Import only viral reads, and filter out reads with supplementary alignments (potential chimaeras) ####
### Write out filtered .bam files
if(write.filtered) {
  bam.filt <- import.bams(bamfiles[], 
                          bamnames[], 
                          ## for import.bams
                          write.filtered = T, force.create=force.create, 
                          rm.gaps.in.aln = rm.gaps.in.aln, filtering='all.reads.w.supp.ali', 
                          ## for ov.from.bam2
                          mapq.filt = mapq.filt,
                          flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                          is.lortia = is.lortia)
  ##
  target <- gsub('\\.bam', '.filtered.bam', bamfiles)
  target <- c(target, paste0(target, '.bai'))
  dir.create(outfilt)
  destination     <- gsub('/bam/', '/bam.filt/', target)
  destination     <- gsub('.filtered', '', destination)
  file.rename(target, destination)
  
  bamfiles <- destination[grep('.bai', destination, invert = T)] #paste0(stringi::stri_replace_last_regex(bamfiles, '.bam', ''), '.filtered.bam')
} else {
  bamfiles
  bamfiles
  
  bam.filt <- import.bams(bamfiles[], 
                          bamnames[], 
                          ## for import.bams
                          write.filtered = T, force.create=force.create, 
                          rm.gaps.in.aln = rm.gaps.in.aln, filtering='all.reads.w.supp.ali', 
                          ## for ov.from.bam2
                          mapq.filt = mapq.filt,
                          flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                          is.lortia = is.lortia)
}

bam.filt <- bam.filt[!is.na(bam.filt$qname), ]
if(!is.na(seqnames.tofilt)) {stopifnot(all(bam.filt$seqnames == seqnames.tofilt)); bam.filt$seqnames <- as.character(bam.filt$seqnames)}


#### Omit mock from filtered
bam.filt <- bam.filt[!grepl('mock', bam.filt$sample), ]
#### ####
##


#### Summarise seqnames and samples from filtered .bams
filt.stats <- idxstats.bams(bamfiles, bamnames = bamnames )

for (i in unique(filt.stats$sample)) {
  stats <- filt.stats[filt.stats$sample == i, ]
  try({ filt.stats$unmapped[filt.stats$sample == i] <- stats$unmapped[stats$seqnames == '*'] })
}

filt.stats <- filt.stats[filt.stats$seqnames != '*', ]
filt.stats$ratio <- filt.stats$mapped / all.stats$unmapped
##
bam.filt.sum <- bam.filt %>% group_by(sample, seqnames) %>% summarise(sub_alignment_count=n())
#### ####
##







#### Read counts (all)
reads.all          <- bam.all[,] %>% group_by(sample, qname, qwidth) %>% summarise()
reads.all          <- merge(reads.all, metafilt[,metacols], by=metacols[1])
#### Read counts (filtered)
reads.filt      <- bam.filt[,] %>% group_by(sample, qname, qwidth) %>% summarise()
reads.filt      <- merge(reads, metafilt[,metacols], by=metacols[1])
reads.filt$org  <- genome
#### because of potential chimeric reads, viral read names are imported from bam.filt
reads.all$org <- 'unknown/chimeric'
reads.all$org[is.element(reads.all$qname, reads.filt$qname)] <- genome
reads.all$org[!is.na(reads.all$qwidth) & reads.all$org != genome] <- 'host'
plyr::count(reads.all$org)

## check
stopifnot(luniq(reads.filt$qname) == nrow(reads.filt))

##
colnames(read.counts.filt)[5] <- 'viral_nochim_read_count'

## summarise
read.counts.all     <- as.data.frame(reads.all %>% 
                                       group_by(across(any_of(c('seqnames', metacols, 'org')))) %>% 
                                       summarise(read_count=n()) %>%
                                       spread(org, read_count, fill=0))

read.counts.filt    <- as.data.frame(reads.filt %>% 
                                       group_by(across(any_of(c('seqnames', metacols, 'org')))) %>% 
                                       summarise(read_count=n()) %>%
                                       spread(org, read_count, fill=0))
## write read counts
write_tsv(read.counts, paste0(outdir, '/read_counts.all.tsv'))
##
write_tsv(read.counts, paste0(outdir, '/read_counts.filt.tsv'))

