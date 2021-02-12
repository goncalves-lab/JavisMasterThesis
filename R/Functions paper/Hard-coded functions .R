

#Theme used in plotting.
my_theme<- function (base_size = 11, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.ticks = element_line(colour = "black"),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line	= element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text = element_text(size = rel(0.8), colour = "black"), 
          strip.background = element_blank()
    )
}


#Unused function, plots the signatures directly from the AF table path.
AFtoSignature <- function(AF, species, treatment){
  a <- lapply(AF, readAF, species = species, treatment = treatment)
  b <- bind_rows(a)
  c <- AFtoVranges(b)
  d <- contextmatrix(c, species = species, treatment = treatment) %>% SignaturePlot()
}

#Translates the fastq names to their short version.
translate_names <- function(x){
  sample_names <- c("F7", "F8", "B2", "J4", "L11", "B2", "B7", "I8", "E2", "E8", "B8", "E7", "G8", "E9", "E8", "B12", "D2", "H5", "H12", "G1")
  names(sample_names) <- c("AS-422403", "AS-422405", "AS-422407", "AS-422409", "AS-422411", "AS-422415", "AS-475119", "AS-475121", "AS-475123", "AS-475125", "AS-475127", "AS-475129", "AS-475131", "AS-475133", "AS-475135", "AS-452423", "AS-452425", "AS-452427", "AS-452429", "AS-452431")
  a<-x
  for (i in 1:length(sample_names)){a <- gsub(names(sample_names)[i], sample_names[i], a)}
  return(a)
}
