library(dndscv)
library(tidyverse)
library(SomaticSignatures)

AFtoDnds <- function(AF, species, treatment){
  if(class(AF) == "data.frame")
    a <- AF
  else
    a <- lapply(AF, readAF, species = species, treatment = treatment)
  a <- bind_rows(a)
  b <- a[,c("sampleID", "chr", "pos", "ref", "alt", "af")]
  c <- dndscv(b, refdb = paste0("/icgc/dkfzlsdf/analysis/B210/Javi/",species, "_refcds.rda"))
  return(c)
}

AFtoGeneNames <- function(AF, dnds = NULL, species, treatment){
  if(is.null(dnds))
    a <- AFtoDnds(AF, species = species, treatment = treatment)
  else
    a <- dnds
  d<- a$annotmuts %>% unite(data = ., col = "chr_pos", c("chr", "pos"), sep = "_")
  e <- readAF(AF, species, treatment) %>% unite(data = . , col = "chr_pos", c("chr", "pos"), sep = "_")
  e <- merge(d, e, by = c("chr_pos", "sampleID"))
  return(e)
}
