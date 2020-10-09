library(dndscv)
library(tidyverse)
library(SomaticSignatures)

AFtoVranges <- function(AF){
  VRanges(seqnames = AF$chr,
          ranges = as.character(AF$pos),
          ref = AF$ref,
          alt = AF$alt,
          sampleID = AF$sample,
          treatment = AF$treatment,
          species = AF$species,
          af=AF$af,
          SBS=AF$SBS)
}
AFtoGranges <- function(AF){
  GRanges(seqnames = AF$chr,
          ranges = as.character(AF$pos),
          ref = AF$ref,
          alt = AF$alt,
          sampleID = AF$sample,
          treatment = AF$treatment,
          species = AF$species,
          af=AF$af,
          SBS=AF$SBS)
}

contextmatrix <- function(Vr, species, treatment, normalize = T){
  a <- mutationContext(Vr, ref = get(paste0(species, "_dna")) ) %>% motifMatrix(group = "sampleID", normalize = normalize)
  colnames(a) <- Vr$sampleID %>% unique()
  b <- as.data.frame(a) %>% 
  mutate("SBS2" = rownames(.)) %>% separate(SBS2, into = c("SBS2","ctxt"), sep = " ") %>% 
    pivot_longer(cols = colnames(a), names_to = "sampleID") %>% 
    mutate(species = species, treatment = treatment)
  return(b)
}

meanContext <- function(ctxt){
  a <- ctxt %>% group_by(SBS2, ctxt, species, treatment) %>% summarise("value" = mean(value)) %>% ungroup()
  b <- a %>% group_by(species) %>% pivot_wider(names_from ="treatment", values_from = "value") %>% 
    mutate("subtract" = treated - control) %>% pivot_longer(cols = c(unique(a$treatment), "subtract"), names_to = "treatment")
}

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

SignaturePlot <- function(ctxt){ggplot(ctxt , aes(x = ctxt, y = value*100 , fill = SBS2 ))+geom_col(position="dodge")+facet_grid(~SBS2)+my_theme()+guides(fill = FALSE)+expand_limits(y = 0) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    ylab("Percentage") + xlab("96-trinucleotide context")+scale_fill_brewer(palette = "Set2")}

AFtoSignature <- function(AF, species, treatment){
  a <- lapply(AF, readAF, species = species, treatment = treatment)
  b <- bind_rows(a)
  c <- AFtoVranges(b)
  d <- contextmatrix(c, species = species, treatment = treatment) %>% SignaturePlot()
}
translate_names <- function(x){
  sample_names <- c("F7", "F8", "B2", "J4", "L11", "B2", "B7", "I8", "E2", "E8", "B8", "E7", "G8", "E9", "E8", "B12", "D2", "H5", "H12", "G1")
  names(sample_names) <- c("AS-422403", "AS-422405", "AS-422407", "AS-422409", "AS-422411", "AS-422415", "AS-475119", "AS-475121", "AS-475123", "AS-475125", "AS-475127", "AS-475129", "AS-475131", "AS-475133", "AS-475135", "AS-452423", "AS-452425", "AS-452427", "AS-452429", "AS-452431")
  a<-x
  for (i in 1:length(sample_names)){a <- gsub(names(sample_names)[i], sample_names[i], a)}
  return(a)
}
