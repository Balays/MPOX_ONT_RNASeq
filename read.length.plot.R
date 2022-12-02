#### Read length plots
require(cowplot)
require(viridisLite)
require(ggsci)
require(hrbrthemes)

#### Color settings
ncolor  <- length(unique(metafilt$hpi[!grepl('mock', metafilt$sample)]))
palette <- viridis::viridis(ncolor, alpha=0.70)
scales::show_col(palette)
vline_col <- pal_aaas()(10)[2]

#### Violin per hpi
ggv <- ggplot(reads.filt[!grepl('mock', reads.filt$sample), ]
              , aes(x = hpi, y=qwidth, fill=hpi)
) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(outlier.size = 0.5, width=0.1) + 
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='red') +
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1) +
  scale_fill_manual(values=palette) +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') 
  # + facet_wrap(~hpi, nrow = 1, scales = 'free_x')

ggsave(paste0(outdir, '/read.lengths.virus.violin.hpi.jpg'), ggv, width = 12, height = 6)

#### Violin per sample
ggv <- ggviolin(reads.filt[!grepl('mock', reads.filt$sample), ]
         , x = 'sample', y='qwidth', fill='hpi', 
         add = 'boxplot', add.params=list(size=0.05)) + 
  scale_fill_manual(values=palette) +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
  facet_wrap(~hpi, nrow = 1, scales = 'free_x')

ggsave(paste0(outdir, '/read.lengths.virus.violin.samples.jpg'), ggv, width = 18, height = 6)

#### Violin -->> SuppFig 1C
ggv <- ggplot(reads.filt[!grepl('mock', reads.filt$sample), ]
              , aes(x = sample, y=qwidth, fill=hpi)
) + 
  geom_violin(trim=FALSE) + 
  #geom_boxplot(outlier.size = 0.5, width=0.2) + 
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='red') +
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1) +
  scale_fill_manual(values=palette) +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
  facet_wrap(~hpi, nrow = 1, scales = 'free_x')

ggsave(paste0(outdir, '/SuppFig1C_read.lengths.virus.violin.samples.v2.jpg'), ggv, width = 18, height = 6)
ggsave(paste0(outdir, '/SuppFig1C_read.lengths.virus.violin.samples.v2.tiff'), ggv, width = 18, height = 6)

#### Density -->> SuppFig 1A
reads.stats <- reads.filt[!grepl('mock', reads.filt$sample), ] %>% group_by(hpi) %>% 
  summarise(min=min(qwidth), max=max(qwidth), mean=mean(qwidth), median=median(qwidth))

ggd <- ggplot(reads.filt[!grepl('mock', reads.filt$sample), ]
              ) + 
  #geom_histogram(aes(x = qwidth)) + 
  geom_density(aes(x = qwidth, fill=hpi)) +
  #scale_color_aaas() +
  #scale_color_manual(values=palette) +
  scale_fill_manual(values=palette) +
  geom_vline(data=reads.stats, aes(xintercept = mean), color=vline_col, size=1, linetype="dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_cartesian(xlim = c(0,5000)) +
  guides(fill = guide_legend(title = "Hours past infection", title.position = "left")) +
  theme_ipsum() + 
  theme(legend.position = 'bottom', 
        axis.title.x= element_blank(), axis.title.y=element_blank(),
        panel.spacing= unit(0.5, units = 'cm'),
        strip.text  = element_blank()) +  
  facet_wrap(~hpi, ncol = 1, scales = 'fixed')
#ggd

ggsave(paste0(outdir, '/SuppFig1A_read.lengths.virus.density.hpi.jpg'),  ggd, width = 12, height = 10)
ggsave(paste0(outdir, '/SuppFig1A_read.lengths.virus.density.hpi.tiff'), ggd, width = 12, height = 10)

#### Histogram -->> SuppFig 1B
bins <- 50
reads.stats <- reads.filt[!grepl('mock', reads.filt$sample), ] %>% group_by(libtype) %>% 
  summarise(min=min(qwidth), max=max(qwidth), mean=mean(qwidth), median=median(qwidth))

ggh <- ggplot(reads.filt[!grepl('mock', reads.filt$sample), ],
              aes(x = qwidth, fill=libtype)
) + 
  geom_histogram(bins = bins, color='black') + 
  #geom_density(aes(x = qwidth, fill=hpi)) +
  #scale_color_aaas() +
  #scale_color_manual(values=palette) +
  scale_fill_manual(values=palette[c(1,length(palette))]) +
  geom_vline(data=reads.stats, aes(xintercept = mean), color=palette[5], size=1, linetype="dashed") +
  scale_y_continuous(labels = scales::number_format()) +
  scale_x_continuous(labels = scales::number_format(), breaks = seq(0, 20000, length.out = bins+1)) +
  coord_cartesian(xlim = c(0,5000)) +
  guides(fill = guide_legend(title = "Library type", title.position = "left")) +
  theme_ipsum() + 
  theme(#legend.position = 'right',
        axis.text.x  = element_text(angle = -60),
        axis.title.x = element_blank(), axis.title.y=element_blank(),
        panel.spacing= unit(0.5, units = 'cm'),
        strip.text  = element_blank()
        ) +  
  facet_wrap(~libtype, ncol = 1, scales = 'fixed')
#ggh

ggsave(paste0(outdir, '/SuppFig1B_read.lengths.virus.histogram.libtype.jpg'),  ggh, width = 12, height = 10)
ggsave(paste0(outdir, '/SuppFig1B_read.lengths.virus.histogram.libtype.tiff'), ggh, width = 12, height = 10)

##### combine
SuppFig1 <- plot_grid(plot_grid(ggd + ggtitle('A'), ggh + ggtitle('B')), ggv + ggtitle('C'), nrow=2)
ggsave(paste0(outdir, '/SuppFig1.jpg' ), SuppFig1, width = 16, height = 12)
ggsave(paste0(outdir, '/SuppFig1.tiff'), SuppFig1, width = 16, height = 12)


####### ALL READS #######

##### Violinplot of host and viral read lengths
plot.data <- plyr::rbind.fill(reads.all[reads.all$org == 'host', ], reads.filt)

### per hpi -->> Fig3
ggvc <- ggplot( plot.data
                , aes(x = hpi, y=qwidth, fill=hpi) ) +
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1) + 
  scale_y_continuous(labels = scales::number_format()) +
  coord_cartesian(ylim = c(0,20000)) +
  scale_fill_viridis_d() +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none')
  # + facet_grid(cols=vars(hpi), scales='free')

ggsave(paste0(outdir, '/read.lengths.violin.all.comb.hpi.jpg' ),  ggvc, width = 12, height = 6)
ggsave(paste0(outdir, '/read.lengths.violin.all.comb.hpi.tiff' ), ggvc, width = 12, height = 6)

### per hpi and organism
ggvh <- ggplot( plot.data
                , aes(x = hpi, y=qwidth, fill=hpi) ) +
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1) +
  scale_y_continuous(labels = scales::number_format()) +
  coord_cartesian(ylim = c(0,20000)) +
  scale_fill_viridis_d() +
  theme_ipsum() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
  facet_grid(rows=vars(org), scales='free')

### per sample and organism
ggvs <- ggplot( plot.data
                  , aes(x = sample, y=qwidth, fill=hpi) ) +
  geom_violin(trim=FALSE) + 
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.1) +
  scale_y_continuous(labels = scales::number_format()) +
  coord_cartesian(ylim = c(0,20000)) +
  scale_fill_viridis_d() +
  theme_ipsum() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none') +
  facet_grid(cols=vars(hpi), rows=vars(org), scales='free')


ggsave(paste0(outdir, '/read.lengths.all.hpi.jpg' ),  ggvh, width = 15, height = 10)
ggsave(paste0(outdir, '/read.lengths.all.samples.jpg' ), ggvs, width = 22, height = 12)

ggsave(paste0(outdir, '/read.lengths.all.samples.tiff'), ggvs, width = 22, height = 12)
ggsave(paste0(outdir, '/read.lengths.all.hpi.tiff'),  ggvh, width = 15, height = 10)

#######  #######

