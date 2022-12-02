library(circlize)
library(viridis)

#### cDNA ####

#plot.links <- T

##### Circos annotation
ann <- feature.df[,c(by, "strand", feature.colname)]
colnames(ann) <- c('chr', 'start', 'end', 'strand', 'gene')

#### circos settings
mar <- c(0.01, 0.01)
palette <- 'inferno'
alpha   <- 0.8

#### coverage settings
lognorm <- 10


#### #####

samples <- unique(as.character(win.cov.sum$sample)) #[8] #[c(4,10)] #,10,16
samples  <- samples[grep('MIX', samples, invert = T)]
#samples <- samples[7]

col     <- viridis(length(samples), option = palette, alpha = 1)#[c(4,8)]
## OR
col     <- pal_npg()(10)[-c(10,7,6,1)]  #((length(samples))
col     <- col[c(1,2,5,3,4,6)]
#scales::show_col(col)

coldf   <- data.frame(samples= samples, 
                      col    = col ## "#330A5FCC" ## 
                      )



#### Hpi 1 - Hpi 4 ####
isamples <- 1:3

circos.clear()
pdf(file=paste0(outdir, "/circos.cDNA/circos.cDNA_A.pdf"), width = 15, height = 15)

col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c(genome), xlim=c(0, l_genome))

#for (i in 1:1) { #length(samples)
sample <- samples[isamples] # 1:3
cols   <- c(1:length(sample))
samp.win.cov.sum <- win.cov.sum[is.element(win.cov.sum$sample, sample), ]

message('plotting ', sample, '...')

cov <- samp.win.cov.sum[,c("seqnames", "strand", "start.window", "end.window", "window_median.cov", "sample")]
colnames(cov) <- c('region', 'strand', 'start', 'end', 'value', 'sample')
cov <- cov[order(cov$start),]

if (!is.na(lognorm)) { cov$value <- log(cov$value, base = lognorm); cov$value[cov$value == '-Inf'] <- 0}

cov.both <- cov %>% group_by(region, start, end, sample) %>% summarise(value=sum(value))
cov.both <- cov.both %>% spread(sample, value, fill=0)

# coverage both strands
circos.genomicTrack(data=cov.both[,], stack = F, panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", lwd=0.7, col=coldf$col[isamples])
}, track.height=0.4, bg.border=F, track.margin=mar)

circos.yaxis(labels.cex=0.5, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF") #at=c(100, 500), 
#} 

# genome x axis
brk <- seq(0, 2, 0.25)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="bottom", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

# genome track
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.05, track.margin=mar)


# gene labels both strand
circos.genomicLabels(ann, labels.column=5, cex=0.5, col=col_text, line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.01, labels_height=0.03, track.margin=mar)


dev.off()
#### ####


#### Hpi 6 - Hpi 24 ####
isamples <- 4:6

circos.clear()
pdf(file=paste0(outdir, "/circos.cDNA/circos.cDNA_B.pdf"), width = 15, height = 15)

col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors=c(genome), xlim=c(0, l_genome))

#for (i in 1:1) { #length(samples)
sample <- samples[isamples] # 1:3
cols   <- c(1:length(sample))
samp.win.cov.sum <- win.cov.sum[is.element(win.cov.sum$sample, sample), ]

message('plotting ', sample, '...')

cov <- samp.win.cov.sum[,c("seqnames", "strand", "start.window", "end.window", "window_median.cov", "sample")]
colnames(cov) <- c('region', 'strand', 'start', 'end', 'value', 'sample')
cov <- cov[order(cov$start),]

if (!is.na(lognorm)) { cov$value <- log(cov$value, base = lognorm); cov$value[cov$value == '-Inf'] <- 0}

cov.both <- cov %>% group_by(region, start, end, sample) %>% summarise(value=sum(value))
cov.both <- cov.both %>% spread(sample, value, fill=0)

# coverage both strands
circos.genomicTrack(data=cov.both[,], stack = F, panel.fun=function(region, value, ...) {
  circos.genomicLines(region, value, type="l", lwd=0.7, col=coldf$col[isamples])
}, track.height=0.4, bg.border=F, track.margin=mar)

circos.yaxis(labels.cex=0.5, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF") #at=c(100, 500), 
#} 

# genome x axis
brk <- seq(0, 2, 0.25)*10^5
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="bottom", major.at=brk, labels=round(brk/10^5, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

# genome track
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.05, track.margin=mar)


# gene labels both strand
circos.genomicLabels(ann, labels.column=5, cex=0.5, col=col_text, line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.01, labels_height=0.03, track.margin=mar)


dev.off()
#### ####











