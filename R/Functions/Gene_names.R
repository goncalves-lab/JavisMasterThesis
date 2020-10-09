library(dndscv)
library(tidyverse)
library(SomaticSignatures)


SlicedtoGeneNames <- function(Sliced, AF, species, treatment = NA){
  a<- read.table(Sliced,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="") %>%
    mutate("species" = species, "treatment" = treatment)
  a$sampleID <- gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(mmus|hgla)/(Treated|Control|17618|MEF-1|MEF-all|MEF-2|MEF-3)(/full_depth|)/Mutation_calling/genes/", "", Sliced) %>% 
    gsub("_sliced.bed", "", .)
  b <- a %>% separate(V6, into = c("gene_id", "other"), sep = ";")
  b$gene_id <- gsub("gene_id ", "", b$gene_id) %>% gsub("\"", "", .) %>% gsub(";", "", .) %>% gsub(" ", "", .)
  b$gene_name <- b[,(grepl("gene_name", b))] %>% gsub("%gene_name ", "", .) %>% 
    gsub("\"", "", .) %>% gsub(";", "", .) %>% gsub(" ", "", .) %>% gsub("%gene_sourceensembl", "", .) %>% casefold(upper = T)
  c <- b %>% dplyr::rename("chr" = V1, "pos" = V3, "frame" = V5, "Strand" = V4, "region" = V2) %>% mutate("V7" = 0) %>% dplyr::select(-c(V7, other))
  d<- c %>% unite(data = ., col = "chr_pos", c("chr", "pos"), sep = "_")
  e <- readAF(AF, species, treatment) %>% unite(data = . , col = "chr_pos", c("chr", "pos"), sep = "_")
  e <- merge(d, e, by = c("chr_pos", "sampleID", "species", "treatment"))
  return(e)
}

CompareGeneNamesMartincorena <- function(GeneNames, genes = "Martin"){
  name <- colnames(GeneNames)[grep("gene", colnames(GeneNames))] %>% gsub("gene_id", NA, .) %>% unique() %>% na.omit()
  if(genes == "Martin")
    Martincorena_genes <- read.table("/icgc/dkfzlsdf/analysis/B210/Javi/dndscv/Genes - Martincorena et al 2017.tsv", sep = " ", header=F) %>%
    t() %>% as.data.frame() 
  if(genes == "dndscv")
    Martincorena_genes <- read.table("/icgc/dkfzlsdf/analysis/B210/Javi/dndscv/known_cancergenes_dndscv.tsv", sep = " ", header=T)
  colnames(Martincorena_genes) <- name
  a <- inner_join(Martincorena_genes, GeneNames, by = name)
  return(a)
}
plotGeneNamesPoints <- function(GeneNames, labels = F, thr = 0.1){
  name <- colnames(GeneNames)[grep("gene_name", colnames(GeneNames))]
  mut <- GeneNames %>% mutate("random" = runif(nrow(GeneNames), max = 3), "random2" = runif(nrow(GeneNames), max =3))
  if (labels == F)
    a <- ggplot(mut, 
         aes(x=random, y = random2, color = get(name)))+ggforce::geom_circle(aes(x0 = random, y0 = random2, r = af, fill=get(name), alpha = 1, color= get(name)), inherit.aes = F)+facet_wrap(~sampleID)+my_theme()
  if(labels == T)
    a <- ggplot(mut, aes(x=random, y = random2, color = get(name)))+ggforce::geom_circle(aes(x0 = random, y0 = random2, r = af, fill=get(name), alpha = 1, color = get(name)), inherit.aes = F)+
                           ggrepel::geom_label_repel(data = mut %>% subset(af > thr), aes(x=random, y = random2, color = get(name), size = 0.2, label=gene_name))+
                          facet_wrap(~sampleID)+my_theme()+theme(legend.position = 0)
  return(a)
}
plotGeneNamesPoints_grid <- function(GeneNames, labels = F, thr = 0.1){
  name <- colnames(GeneNames)[grep("gene_name", colnames(GeneNames))]
  mut <- GeneNames %>% mutate("random" = runif(nrow(GeneNames), max = 3), "random2" = runif(nrow(GeneNames), max =3)) %>% 
    mutate( "letterz" = sub("^([[:alpha:]]*).*", "\\1", sampleID))
  lt = vector()
  for(i in mut$letterz){lt <- c(lt, grep(i, LETTERS))}
  mut$letters <- lt
  if (labels == F)
    a <- ggplot(mut, 
                aes(x=2*random + 2*parse_number(sampleID), y = 2*random2 + 2*letters, color = get(name)))+
    ggforce::geom_circle(aes(x0 = 2*random + 2*parse_number(sampleID), y0 = 2*random2 + 2*letters, r = af, fill=get(name), alpha = 1, 
                         color= get(name)), inherit.aes = F)+facet_wrap(~sampleID)+my_theme()
  if(labels == T)
    a <- ggplot(mut, aes(x=2*random + 2*parse_number(sampleID), y = 2*random2 + 2*letters, color = get(name)))+ggforce::geom_circle(aes(x0=2*random + 2*parse_number(sampleID), y0 = 2*random2 + 2*letters, r = af, fill=get(name), alpha = 1, color = get(name)), inherit.aes = F)+
    ggrepel::geom_label_repel(data = mut %>% subset(af > thr), aes(x=2*random + 2*parse_number(sampleID), y = 2*random2 + 2*letters, color = get(name), size = 0.2, label=gene_name))+
    facet_wrap(~sampleID)+my_theme()+theme(legend.position = 0)
  return(a)
}



CommonGeneNames <- function(x, y){
  genes <- merge(x, y, by = c("gene_name", "chr_pos", "species", "treatment")) %>% .$chr_pos
  a <- GeneNames_b %>% subset(sampleID == unique(x$sampleID) | sampleID == unique(y$sampleID)) %>% 
    subset(chr_pos %in% genes) %>% subset(region == "gene")
}

cosineDist <- function(x){
  as.matrix(x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}


my_palette <- list(`Martin`  = c("a" = "#978BBC", "b" =  "#F5B577", "c" =  "#6A589E","d" = "#F1963F", "e" = "#59C2C7", "f" = "#F06DB3", "g" = "#F2913C"),
                   `Mannheim` = c("blue1" = "#5CC3FA", "blue2" =  "#7AD5F7", "red" =  "#E47354", "red2" = "#A52D1E","yellow" = "#F2E079", "gold" = "#E5A74E"),
                   `Sunset` = c("#eeaf61" , "#fb9062" , "#ee5d6c" , "#ce4993" ,"#AA4C8A", "#6a0d83"),
                   `Sunset2` = c("#ffa4ab", "#e95a71"),
                   `DMBA` = c("#AA4C8A", "#999999", "#777777", "#555555", "#333333")
                   )

my_pal <- function(palette = "Mannheim") {
  pal <- my_palette[[palette]]
  colorRampPalette(pal)
}
scale_color_my <- function(palette = "Mannheim") {
  pal <- my_pal(palette = palette)
  discrete_scale("colour", "My", palette = pal)
}
scale_fill_my <- function(palette = "Mannheim") {
  pal <- my_pal(palette = palette)
  discrete_scale("fill", "My", palette = pal)
}
