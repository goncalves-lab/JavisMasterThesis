ReadIntervals <- function(x){
  read.delim(x, header = F, sep = "\t")%>% 
    mutate("Mb3" = c(rep(1:(nrow(.)/3000), each = 3000), rep(as.integer(nrow(.)/3000 + 1), (nrow(.) - as.integer(nrow(.)/3000)*3000))))%>% 
    group_by(Mb3) %>% summarise("number" = sum(V7))
}
ReadIntervalsMatrix <- function(x, length = 1000){
  print(paste("Creating a matrix with number of mutations in intervals of length", length, "base pairs"))
  a <- length/1000
  name <- gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(hgla|mmus)/(Treated|Control|17618)/(full_depth/|)Mutation_calling/1kb_intervals/", "", x) %>% gsub(".tsv", "", .)
  if(a == 1)
    b <- read.delim(x, header = F, sep = "\t")%>% mutate("bin" = 1:nrow(.), "number" = V7) %>% dplyr::select(number) %>% as.matrix() %>% t()
  else
    b <- read.delim(x, header = F, sep = "\t")%>%
    mutate("bin" = c(rep(1:(nrow(.)/a), each = a), rep(as.integer(nrow(.)/a + 1), (nrow(.) - as.integer(nrow(.)/a)*a))))%>% 
    group_by(bin) %>% summarise("number" = sum(V7)) %>% dplyr::select(-bin) %>% as.matrix() %>% t()
  rownames(b) <- name
  return(b)
}
