library(circlize)
library(viridis)

#### dRNA ####

#### Links settings
plot.links <- T
mincount  <- NA
top       <- 20

#### Circos annotation
ann <- feature.df[,c(by, "strand", feature.colname)]
colnames(ann) <- c('chr', 'start', 'end', 'strand', 'gene')

#### Circos settings
mar <- c(0.01, 0.01)
palette <- 'inferno'
alpha   <- 0.8



#### #####

samples <- unique(as.character(win.cov.sum$sample)) #[8] #[c(4,10)] #,10,16
samples <- samples[7]

col     <- viridis(12, option = palette, alpha = 1)[c(4,8)]
col     <- pal_npg()(10)[c(4,8)]  #((length(samples))
#scales::show_col(col)

coldf   <- data.frame(samples= samples, 
                      col    = col ## "#330A5FCC" ## 
                      )

circos.clear()
pdf(file=paste0(outdir, "/circos.dRNA/circos.dRNA.pdf"), width = 15, height = 15)

col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c(genome), xlim=c(0, l_genome))


i <- 1
#for (i in 1:length(samples)) {
sample <- samples[i]
samp.win.cov.sum <- win.cov.sum[is.element(win.cov.sum$sample, sample), ]

cov <- samp.win.cov.sum[,c("seqnames", "strand", "start.window", "end.window", "window_median.cov")]
colnames(cov) <- c('region', 'strand', 'start', 'end', 'value')
cov <- cov[order(cov$start),]
covplus  <- cov[cov$strand == '+', c('region', 'start', 'end', 'value')]
covminus <- cov[cov$strand == '-', c('region', 'start', 'end', 'value')]
covminus$value <- covminus$value * -1

# coverage plus strand
circos.genomicTrack(data=covplus[,], panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", area=T, col=coldf$col[i], lwd=0.7)
}, track.height=0.2, bg.border=F, track.margin=mar)
# coverage axis
circos.yaxis(labels.cex=0.5, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF" #at=c(100, 500),
             )
# 
#dev.off()

# genome x axis
brk <- seq(0, 2, 0.25)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="bottom", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

# gene labels plus strand
circos.genomicLabels(ann[ann$strand=='+',], labels.column=5, cex=0.5, col=col_text, line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.01, labels_height=0.03, track.margin=c(0.04, 0.005))


# genome track
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.05, track.margin=mar)


# gene labels minus strand
circos.genomicLabels(ann[ann$strand=='-',], labels.column=5, cex=0.5, col=col_text, line_lwd=0.5, line_col="grey80", 
                     side="outside", connection_height=0.01, labels_height=0.03, track.margin=mar)

# genome x axis
brk <- seq(0, 2, 0.25)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="bottom", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

# coverage minus strand
circos.genomicTrack(data=covminus[,], panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", area=T, baseline='top', col=coldf$col[i+1], lwd=0.7)
}, track.height=0.2, bg.border=F, track.margin=mar)
# coverage axis
circos.yaxis(at=c(100, 500), labels.cex=0.5, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF")


if(plot.links) {
  source('tr.links.R')
  
  tss.df <- data.frame(tr.links[,c("seqnames", "strand", "count")], start=tr.links$prime5, end=tr.links$prime5, col=tr.links$col)
  tes.df <- data.frame(tr.links[,c("seqnames", "strand", "count")], start=tr.links$prime3, end=tr.links$prime3, col=tr.links$col)
  
  circos.genomicLink(tss.df[tss.df$strand == '+', c("seqnames", "start", "end")], 
                     tes.df[tes.df$strand == '+', c("seqnames", "start", "end")], 
                     col=tr.links$col[tr.links$strand == '+'], 
                     h.ratio = 0.6,
                     border=NA, directional=1)
  circos.genomicLink(tss.df[tss.df$strand == '-', c("seqnames", "start", "end")], 
                     tes.df[tes.df$strand == '-', c("seqnames", "start", "end")], 
                     col=tr.links$col[tr.links$strand == '-'], 
                     h.ratio = 0.4,
                     border=NA, directional=1)
  
}

dev.off()



