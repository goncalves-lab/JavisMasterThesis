SignatureAFthreshold <- function(AF = AF, animal = animal, int = int, width = width, CosOrTA = "TA", quantiles = F){
  if(quantiles == TRUE){
    int1 <- quantile(AF%>% subset(species==animal) %>% .$af, int)
    int2 <- quantile(AF%>% subset(species==animal) %>% .$af, int+width)
  }
  if(quantiles == FALSE){
    int1 <- int-width
    int2 <- int
  }
  check_t <-  AF %>% subset(species==animal & treatment == "treated"  & af > int1 & af < int2 ) %>% nrow()
  check_c <-  AF %>% subset(species==animal & treatment == "control"  & af > int1 & af < int2 ) %>% nrow()
  if(check_t == 0 | check_c == 0){return(0)}
  a <- bind_rows(AF %>% subset(species==animal & treatment == "treated"  & af > int1 & af < int2 )%>% AFtoVranges()%>%
              contextmatrix(species = animal, treatment = "treated"), AF %>%
              subset(species==animal & treatment == "control"   & af > int1 & af < int2)%>%
              AFtoVranges()%>%
              contextmatrix(species = animal, treatment = "control"))%>% meanContext()
  treatmentEffect <-  a %>% subset(treatment == "subtract") %>% .$value %>% as.numeric()
  DMBA_sig <- read_delim("/icgc/dkfzlsdf/analysis/B210/Javi/DMBA-signature", delim = "\t") %>% .$DMBA %>% as.numeric() %>% pmax(0)
  #DMBA_sig <- c(rep(0,48), DMBA_sig[c(49:64)], rep(0,32))
  DMBA_sig <- c(rep(0,48), 0.01,0,0.04,0,0.02,0.03,0.15,0.04,0.02,0.01,0.08,0.02,0.01,0,0.05,0.01, rep(0,32))
  cosSim <- MutationalPatterns::cos_sim(treatmentEffect, DMBA_sig)
  #dist(rbind(treatmentEffect, DMBA_sig), method = "euclidean")
  TA <- a %>% subset(treatment == "subtract") %>% group_by(SBS2) %>% summarise("SBS"=sum(value)) %>%
    .[(.$SBS2=="TA"),2] %>% as.numeric()
  if(CosOrTA == "cos"){return(cosSim)}
  if(CosOrTA == "TA"){return(TA)}
}

WhereIsTheSignature <- function(AF = AF_b, width = 0.1,  CosOrTA = "TA", quantiles = F){
  length <- 1/width
  rep <- width*1:(length-1)
  v <- c()
  for (i in c(rep)){
    v <- c(v, SignatureAFthreshold( AF, "hgla", i, width,  CosOrTA, quantiles))
  }
  w <- c()
  for (i in c(rep)){
    w <- c(w, SignatureAFthreshold(AF = AF, "mmus", i, width,  CosOrTA, quantiles))
  }
  data.frame("quant" = c(rep), "value" = c(w,v),
             "species" = rep(c("mmus", "hgla"), each = length(rep)))%>%
    mutate("int" = quantile(AF %>%subset(species == animal)%>%.$af, quant))
}
