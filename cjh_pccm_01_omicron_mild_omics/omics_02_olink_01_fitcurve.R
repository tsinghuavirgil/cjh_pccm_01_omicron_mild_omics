rm(list=ls())
setwd('/xxx/')
getwd()
###################
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
#################################

curve_data <- read.csv('fit_curve_data_01.csv', header=T,row.names=1)
head(curve_data)
#############
final_colbrewer <- c('#99CC33',
                     '#FF0033',
                     '#666699',
                     '#996600',
                     '#CCCCFF',
                     '#CC99CC',
                     '#999966',
                     '#220050', '#b30059','#0091a8','#359023', '#ffa500'
)#
plot(1:12,rep(1, 12), col = final_colbrewer, pch = 16, cex = 2)

brewer_pal<-colorRampPalette(final_colbrewer)
curve_match = length(unique(curve_data$name))

##########
#dot_plot#
##########
curve_01_IL8 <- ggplot(curve_data,aes(x=day,y=IL8, color = name)) +
  theme_classic() + geom_point(size=1,shape=16) +
  scale_color_manual(values = brewer_pal(curve_match)) +
  geom_smooth(method = stats::loess, color = 'black',#lm, glm, mgcv::gam,loess
              fill = 'lightgray',size = 1
              #,weight = 2
  ) +
  xlab('day') + ylab('logIL8') +
  scale_x_continuous(limits = c(0,7),breaks = c(0,1,3,5,7)) +
  theme(legend.position = 'none')#
curve_01_IL8

curve_02_VEGFA <- ggplot(curve_data,aes(x=day,y=VEGFA, color = name)) +
  theme_classic() + geom_point(size=1,shape=16) +
  scale_color_manual(values = brewer_pal(curve_match)) +
  geom_smooth(method = stats::loess, color = 'black',#lm, glm, mgcv::gam,loess
              fill = 'lightgray',size = 1
              #,weight = 2
  ) +
  xlab('day') + ylab('logVEGFA') +
  scale_x_continuous(limits = c(0,7),breaks = c(0,1,3,5,7))+
  theme(legend.position = 'none')
curve_02_VEGFA
##############
cowplot::plot_grid(curve_01_IL8, curve_02_VEGFA
                   ncol = 2)
ggsave(merge_fig_01, file='merge_fig_02.pdf', width = 210, height = 900, units = 'mm')#297

