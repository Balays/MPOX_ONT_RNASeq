library(hrbrthemes)
library(tidyverse)
library(readr)
library(formattable)
library(viridisLite)
library(ggpubr)

read_counts <- read.counts.all[,1:6]
colnames(read_counts)[5:6] <- c("Host read count", "Viral read count nochimaera")
read_counts$`Viral read ratio` <- read_counts$`Viral read count nochimaera` / (read_counts$`Host read count` + read_counts$`Viral read count nochimaera`)
read_counts <- read_counts[order(read_counts$Time, read_counts$sample),]
rownames(read_counts) <- NULL
read_counts$`Viral read ratio` <- round(read_counts$`Viral read ratio`, 4)

write_tsv(read_counts, paste0(outdir, '/Table4.tsv'))

cols <- c(viridisLite::viridis(15, 0.65, option = 'D'),
          viridisLite::viridis(15, 0.65, option = 'B') )

scales::show_col(cols, ncol = 3)


formattable(read_counts %>% select(!any_of(c("Mapped read count", "Unmapped read count"))),
            list(
  `Total read count` = color_bar(cols[5]),
  `Host read count` = color_bar(cols[9]),
  `Viral read count nochimaera` = color_bar(cols[22]),
  `Viral read ratio` = color_bar(cols[24])
    ))
## Export Table2.pdf manually


cols <- c(viridisLite::viridis(15, 1, option = 'D'),
          viridisLite::viridis(15, 1, option = 'B') )

read_counts.gt <- read_counts %>% gather(category, value, -c(1:3))
coeff <- 1
ggrc <- ggplot(read_counts[read_counts$sample != 'dRNA',], 
               aes(x =`Time`)) +
  geom_point(aes(y=`Host read count`), colour=cols[9]) +
  geom_smooth(aes(y=`Host read count`), colour=cols[9]) +
  geom_point(aes(y=`Viral read count nochimaera`), colour=cols[22]) +
  geom_smooth(aes(y=`Viral read count nochimaera`), colour=cols[22]) +
  scale_y_continuous(name ='Viral read count', labels = scales::number_format()
                     , sec.axis = sec_axis(~.*coeff, name='Host read count', labels = scales::number_format())
  ) +
  scale_x_continuous(name = 'Time (hours past infection)')+
  theme_ipsum() +
  theme(axis.title.x       = element_text(size=12, color='black'),
        axis.title.y.right = element_text(size=12, color=cols[9],  hjust = 0.9),
        axis.title.y.left  = element_text(size=12, color=cols[22], hjust = 0.1, vjust = 3),
        plot.margin = unit(c(1,0.25,1,1), units = "cm"))

ggvr <- ggplot(read_counts[read_counts$sample != 'dRNA',],
               aes(x =`Time`)) +
  geom_point(aes(y=`Viral read ratio`), colour=cols[24]) +
  geom_smooth(aes(y=`Viral read ratio`), colour=cols[24]) +
  scale_y_continuous(name ='Viral read ratio', labels = scales::percent_format()
                     #,sec.axis = sec_axis(~.*coeff, name='Total read count')
  ) +
  scale_x_continuous(name = 'Time (hours past infection)')+
  theme_ipsum() +
  theme(axis.title.x       = element_text(size=12, color='black'),
        axis.title.y.left  = element_text(size=12, color=cols[24],  hjust = 0.5),
        plot.margin = unit(c(1,1,1,0.25), units = "cm"))

ggcomb <- cowplot::plot_grid(ggrc, ggvr)

ggsave(paste0(outdir,'/Fig2.tiff'), ggcomb, width = 12, height = 6)
ggsave(paste0(outdir,'/Fig2.jpg'),  ggcomb, width = 12, height = 6)




