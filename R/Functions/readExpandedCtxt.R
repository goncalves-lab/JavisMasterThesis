#Reads output of siglasso function siglasso:vcf2spec resulting in a matrix where the rows represent different samples and the columns proportion of SNVs that are from a specific substitution and context.
#Context lengths can be different and indicate the number of nucleotides considered in one direction. i.e AA.AA is a length of 2. 
#An example of input for the function is /icgc/dkfzlsdf/analysis/B210/Javi/ctxt_hgla2.tsv


readExpandedCtxt <- function(x, length = 1, df= FALSE){
  a <- read.delim(x, header=FALSE) %>% 
    mutate(V1 = str_split(V1, "")) %>% 
    separate(V1, into = c("kk", paste0("a", 1:length),"ref",paste0("c", 1:length), "k")) %>% dplyr::select(-c(kk, k)) %>% 
    dplyr::rename("alt"= V2, "sampleID" = V3) %>% unite("a", c(paste0("a", 1:length)), sep="")%>% 
    unite("c", c(paste0("c", 1:length)), sep="") %>% unite("ctxt", c(a,c), sep = "") %>% RefAltToSBS6() %>% dplyr::select(-c("ref","alt"))%>% 
    group_by(sampleID, SBS, ctxt) %>% mutate(value = n()) %>% 
    group_by(sampleID) %>% mutate(value = value/n()) %>% unique() %>% .[order(.$SBS, .$ctxt),]
  b <- a %>% pivot_wider(names_from = sampleID, values_from = "value", values_fill = 0)
  n <- 4^(2*length)
  ctxt <- c()
  for(i in 1:(2*length)){ctxt <-paste0(ctxt, rep(c("A", "C","G","T"), each = n/(4^i)))}
  data <- data.frame( "SBS" = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = n), 
                    "ctxt" = ctxt)
  kk <- merge(data, b, by = c("SBS", "ctxt"), all.x = T) 
  kk[is.na(kk)] <- 0
  c <- kk %>% dplyr::select(-c("ctxt", "SBS")) %>% t() %>% as.matrix()
  if(df == TRUE)
    return(a)
  else
    return(c)
}


readExpandedCtxt_old <- function(x, length = 1){
  if(length == 1){
    a <- read.delim(x, header=FALSE) %>% 
      mutate(V1 = str_split(V1, "")) %>% 
      separate(V1, into = c("kk", "a","ref","c", "k")) %>% dplyr::select(-c(kk, k)) %>% 
      dplyr::rename("alt"= V2, "sampleID" = V3) %>% unite("ctxt", c(a,c), sep = ".") %>% RefAltToSBS6() %>% dplyr::select(-c("ref","alt"))%>% 
      group_by(sampleID, SBS, ctxt) %>% mutate(value = n()) %>% 
      group_by(sampleID) %>% mutate(value = value/n()) %>% unique() %>% .[order(.$SBS, .$ctxt),]
    b <- a %>% pivot_wider(names_from = sampleID, values_from = "value", values_fill = 0)
    df <- data.frame( "SBS" = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 16), 
                    "ctxt" = paste0(rep(c("A", "C","G","T"), each = 4),".",
                                    rep(c("A", "C","G","T"), 4)))
    kk <- merge(df, b, by = c("SBS", "ctxt"), all.x = T) 
    kk[is.na(kk)] <- 0
    c <- kk %>% dplyr::select(-c("ctxt", "SBS")) %>% t() %>% as.matrix()
    return(c)}
  if(length == 2){
    a <- read.delim(x, header=FALSE) %>% 
    mutate(V1 = str_split(V1, "")) %>% 
    separate(V1, into = c("kk", "a","b","ref","c","d", "k")) %>% dplyr::select(-c(kk, k)) %>% 
    dplyr::rename("alt"= V2, "sampleID" = V3) %>% unite("a", c(a,b), sep="")%>% 
    unite("c", c(c,d), sep="") %>% unite("ctxt", c(a,c), sep = ".") %>% RefAltToSBS6() %>% dplyr::select(-c("ref","alt"))%>% 
    group_by(sampleID, SBS, ctxt) %>% mutate(value = n()) %>% 
    group_by(sampleID) %>% mutate(value = value/n()) %>% unique() %>% .[order(.$SBS, .$ctxt),]
    b <- a %>% pivot_wider(names_from = sampleID, values_from = "value", values_fill = 0)
    df <- data.frame( "SBS" = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 256), 
                    "ctxt" = paste0(rep(c("A", "C","G","T"), each = 64),rep(c("A", "C","G","T"), each = 16),".",
                                    rep(c("A", "C","G","T"), each = 4), rep(c("A", "C","G","T"), 64)))
    kk <- merge(df, b, by = c("SBS", "ctxt"), all.x = T) 
    kk[is.na(kk)] <- 0
    c <- kk %>% dplyr::select(-c("ctxt", "SBS")) %>% t() %>% as.matrix()
    return(c)}
  if(length == 3){
    a <- read.delim(x, header=FALSE) %>% 
    mutate(V1 = str_split(V1, "")) %>% 
    separate(V1, into = c("kk","w", "a","b","ref","c","d","x", "k")) %>% dplyr::select(-c(kk, k)) %>% 
    dplyr::rename("alt"= V2, "sampleID" = V3) %>% unite("a", c(w,a,b), sep="")%>% 
    unite("c", c(c,d,x), sep="") %>% unite("ctxt", c(a,c), sep = ".") %>% RefAltToSBS6() %>% dplyr::select(-c("ref","alt"))%>% 
    group_by(sampleID, SBS, ctxt) %>% mutate(value = n()) %>% 
    group_by(sampleID) %>% mutate(value = value/n()) %>% unique() %>% .[order(.$SBS, .$ctxt),]
    b <- a %>% pivot_wider(names_from = sampleID, values_from = "value", values_fill = 0)
    df <- data.frame( "SBS" = rep(c("C>A","C>G","C>T","T>A","T>C","T>G"), each = 4096), 
                      "ctxt" = paste0(rep(c("A", "C","G","T"), each = 1024), rep(c("A", "C","G","T"), each = 256), rep(c("A", "C","G","T"), each = 64),".",rep(c("A", "C","G","T"), each = 16),
                                      rep(c("A", "C","G","T"), each = 4), rep(c("A", "C","G","T"), 1024) ))
    kk <- merge(df, b, by = c("SBS", "ctxt"), all.x = T) 
    kk[is.na(kk)] <- 0
    c <- kk %>% dplyr::select(-c("ctxt", "SBS")) %>% t() %>% as.matrix()
    return(c)}
}
