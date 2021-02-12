#' Obtaining a Vranges object containing the context of the mutation from the output of AFtoVranges
#'
#' @param Vr Output of the function AFtoVranges.
#' @param species A character string with the name of the species of the samples.
#' @param treatment A character string indicating the treatment received by the samples (Control or Treated).
#' @param normalize Logical value indicating whether to normalize context values to the proportion or not.  

contextmatrix <- function(Vr, species, treatment, normalize = T){
  require("SomaticSignatures")
  require("tidyverse")
  
  a <- mutationContext(Vr, ref = get(paste0(species, "_dna")) ) %>%
    motifMatrix(group = "sampleID", normalize = normalize)
  
  colnames(a) <- Vr$sampleID %>% unique()
  
  b <- as.data.frame(a) %>% 
    mutate("SBS2" = rownames(.)) %>% separate(SBS2, into = c("SBS2","ctxt"), sep = " ") %>% 
    pivot_longer(cols = colnames(a), names_to = "sampleID") %>% 
    mutate(species = species, treatment = treatment)
  return(b)
}