plot.metataxon <- function(Res){
  ggplot(subset(Res$Per_otu,Study=='Meta'),aes(x=log2FoldChange, y=diff_otu, color=Significance)) +
    geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
    geom_errorbarh(aes(xmin=CI_low, xmax=CI_hig), height=0 ) +
    geom_point(size=3)+
    theme_bw() +
    #scale_color_manual(values=c('#E69F00','#0072B2','grey'))  +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    theme(axis.text.x=element_text( hjust=1,size=12,colour = 'black'),
          axis.text.y = element_text(size = 12,colour = 'black'),
          axis.title=element_text(colour='black', size=12),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 11),
    )+
    scale_shape_manual(values=c(17,1,3,4,5,7)) + labs(y='Taxon')

}
