#' Reads the AF table resulting from the Somatic mutation calling pipeline into a table.
#' @param vcf The path to the AF table
#' @param species A character string with the name of the species of the samples.
#' @param treatment A character string indicating the treatment received by the samples (Control or Treated).
#' @param input A character string that can be "vcf", when the input is a vcf file from the Somatic mutation calling pipeline or "AF" when the input is an AF table from the same pipeline. 

readAF <- function(vcf, species, treatment, input = "AF"){
  require(tidyverse)
  name <- gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(hgla|mmus)/(Control|Treated|17618|MEF-(all|1|2|3))(/full_depth|)/(Mutation_calling|Control-1_is_the_normal)/(AF/|)", "", vcf) %>% gsub(".tsv", "", .) %>% gsub("/pass(|2).vcf", "", .)
  if(input == "AF"){
    a<- read.table(vcf, header = FALSE, sep = "\t") %>% dplyr::rename(., "chr"=V1, "pos"=V2, "SBS"=V3,"af"=V4) %>% 
    mutate("sampleID" = name, "species" = species, "treatment" = treatment) %>% 
    separate(SBS, into = c("ref", "alt"), sep = ">", remove = F) %>% mutate(SBS = gsub("A>C", "T>G",.$SBS) %>% gsub("A>G", "T>C", .) %>% gsub("A>T", "T>A", .) %>% gsub("G>C", "C>G", .) %>% gsub("G>T", "C>A", .) %>% gsub("G>A", "C>T", .))
    return(a)
  }
  if(input == "vcf"){
  a <- readPassVcf(vcf) %>% dplyr::select(c(chr,pos,ref,alt,tumorAF)) %>% mutate("af" = as.numeric(tumorAF)) %>% 
    mutate(sampleID = name, species = species, treatment = treatment) %>% 
  unite("SBS", c("ref", "alt"), sep = ">", remove = F) %>% mutate(SBS = gsub("A>C", "T>G",.$SBS) %>% gsub("A>G", "T>C", .) %>% gsub("A>T", "T>A", .) %>% gsub("G>C", "C>G", .) %>% gsub("G>T", "C>A", .) %>% gsub("G>A", "C>T", .))
  return(a)
  }
}