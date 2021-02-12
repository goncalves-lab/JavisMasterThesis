#' Plotting the mutational signature present in a context matrix.
#'
#' @param ctxt Output of the function contextmatrix.

SignaturePlot <- function(ctxt){
  require(tidyverse)
  ggplot(ctxt , aes(x = ctxt, y = value*100 , fill = SBS2 ))+
    geom_col(position="dodge")+facet_grid(~SBS2)+my_theme()+guides(fill = FALSE)+
    expand_limits(y = 0)+
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    ylab("Percentage") + xlab("96-trinucleotide context")+scale_fill_brewer(palette = "Set2")
  }