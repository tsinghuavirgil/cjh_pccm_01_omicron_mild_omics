rm(list = ls())
setwd('/xxx/')
getwd()
library(ggplot2)
################
data<-read.csv('sort_bp_up.csv',header=T)#, sep=', '
data
save_fig_01 <- ggplot(data, aes(x = group, y = Term))+ theme_bw() +
  geom_point(aes(size=-log10(pvalue),color=log2(fold_enrichment))) +
  #scale_color_gradient2(low = "green", midpoint = 3, mid = "blue", high = "red", space="Lab") +
  #scale_color_continuous(low='#77AAAD',high='#DE6449') 
  scale_color_gradient(low='#77AAAD', high='#DE6449')+ #, midpoint = 1, mid = 'white'
  theme(axis.text.x = element_text(angle = 90, hjust = .1, vjust = .5))
print(save_fig_01)

ggsave(save_fig_01, file='sort_01_up_01_bp_dotplot.pdf', width = 8, height = 6)#
