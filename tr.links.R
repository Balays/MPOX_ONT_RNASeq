#### TSS---TES links
sample   <- 'dRNA'
bam.samp <- bam.filt[bam.filt$sample == sample, ]
tr.links <- add.primes(bam.samp)
tr.links <- tr.links %>% group_by(seqnames, prime5, prime3, strand) %>% summarise(count=n())


tr.links$strand <- factor(tr.links$strand, levels = c('+', '-', '*'))

if(is.na(mincount)) {
  mincount  <- plyr::count(tr.links$count)$freq[top]
}
tr.links  <- tr.links[tr.links$count >= mincount, ]


minalpha  <- 0.1
ncolors   <- c(max(tr.links$count[tr.links$strand == '+']),
               max(tr.links$count[tr.links$strand == '-']))


coldf     <- plyr::rbind.fill(data.frame(count=seq(1, ncolors[1]), col= alpha(col[1], seq(minalpha, 1, length.out = ncolors[1])), strand='+' ),
                              data.frame(count=seq(1, ncolors[2]), col= alpha(col[2], seq(minalpha, 1, length.out = ncolors[2])), strand='-' ) )
#scales::show_col(coldf$col)

tr.links  <- merge(tr.links, coldf, by=c('strand', 'count'))

tr.links  <- tr.links[order(tr.links$strand, tr.links$prime5, tr.links$prime3), ]




