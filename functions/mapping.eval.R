### what is the ratio of the mapped part of the whole read?
calc.alignment.ratio <- function(bam) {
  aln.sum <- bam %>% group_by(sample, qname) %>% summarise(qwidth, sum.aln.width=sum(width), width.ratio=sum.aln.width/qwidth)
  aln.sum <- unique.data.frame(aln.sum)
  bam <- merge(aln.sum, bam, by=c('sample', 'qname', 'qwidth'), all=T)
  return(bam)
}

### how many contigs each read has been mapped onto? (detection of chimeras)
calc.n_seq <- function(bam) {
  bam.uni.seq <- unique.data.frame(bam %>% select(sample, qname, seqnames))
  aln.sum <- bam.uni.seq %>% group_by(sample, qname) %>% summarise(seqnames, n_seq=seq_along(seqnames), max_n=max(n_seq))
  bam <- merge(aln.sum, bam, by=c('sample', 'qname', 'seqnames'), all=T)
  return(bam)
}

### how many organims were each read mapped onto?
calc.n_org <- function(bam) {
  bam.uni.seq <- unique.data.frame(bam %>% select(sample, qname, seqnames, org))
  
  aln.sum     <- bam.uni.seq %>% group_by(sample, qname, org) %>% summarise(n_seq_org=n())
  bam         <- merge(aln.sum, bam, by=c('sample', 'qname', 'org'), all=T)
  
  org.sum     <- aln.sum     %>% group_by(sample, qname, org) %>% summarise(n_org=n())
  org.sum     <- org.sum     %>% group_by(sample, qname) %>% summarise(n_org=sum(n_org))
  bam         <- merge(org.sum, bam, by=c('sample', 'qname'), all=T)
  
  return(bam)
}


### does the reads have supplementary or secondary alignments?
map.eval <- function(bam, bam.flags) {
  
  ## 
  sec.flag   <- 256
  sec.flags  <- c(sec.flag, bam.flags$Decimal + sec.flag)
  
  qname.sec <- bam$qname[is.element(bam$flag, c(sec.flags)) ]
  bam$has.secondary <- F
  bam$has.secondary[is.element(bam$qname, qname.sec) ] <- T
  
  
  supp.flag   <- 2048
  supp.flags  <- c(supp.flag, bam.flags$Decimal + supp.flag)
  
  qname.supp  <- bam$qname[is.element(bam$flag, c(supp.flags)) ]
  bam$has.supplementary <- F
  bam$has.supplementary[is.element(bam$qname, qname.supp) ] <- T
  
  return(bam)
}


