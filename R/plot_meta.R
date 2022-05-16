plot_meta <- function(res,size=9){

  ggplot(res,aes(x=log2FC, y=Study, color=Significance)) +
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
    geom_point(shape=18,size=3) +
    theme_bw() +
    theme(panel.border = element_blank(), axis.line = element_line(),
          plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1,size=size),
          axis.text.y=element_text(size=size),
          axis.text = element_text(size = size,colour = 'black'),
          legend.title = element_text(size = size),
          legend.text = element_text(size = size))+
    xlab('log2(fold difference)')+
    labs(y='Study')
}
