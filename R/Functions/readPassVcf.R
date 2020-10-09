#' Reads the vcf resulting from the Somatic mutation calling pipeline into a table.
#' @param vcf The path to the vcf file
#' @param flag logical value indicating whether to add the filtering flag of mutect or not
#'
#'
readPassVcf <- function(vcf, flag = F){
  a <- read.table(vcf, sep = "\t", header = F)
  if (flag == F)
    b <- a %>% dplyr::select(-c(V3,V6,V7))
  else
    b <- a %>% dplyr::select(-c(V3,V6))
  c <- b %>% unite("all", paste0("V",10:ncol(a)), sep = ",")
  d <- c
  d$all <- gsub("\\.,", "", c$all) %>% gsub(",\\.", "", .) %>% gsub("0/1", "%0/1", .)
  e <- d %>% separate(all, into = c("normal", "tumor"), sep = "%") %>% 
    separate(V8, into =c("AC", "AN", "DP", "ECNT", "NLOD", "N_ART_LOD", "POP_AF", "P_CONTAM", "P_GERMLINE", "SF", "TLOD") , sep = ";")%>% separate(tumor, into =paste0("tumor", c("GT", "AD", "AF","DP","F1R2","F2R1","MBQ","MFRL","MMQ","MPOS")) , sep = ":") %>% separate(normal, into =paste0("normal", c("GT", "AD", "AF","DP","F1R2","F2R1","MBQ","MFRL","MMQ","MPOS")) , sep = ":")
  f <- e %>% dplyr::select(-V9) %>% 
    mutate(AC = readr::parse_number(AC), AN = readr::parse_number(AN), DP = readr::parse_number(DP),
           ECNT = readr::parse_number(ECNT), NLOD = readr::parse_number(NLOD), N_ART_LOD = readr::parse_number(N_ART_LOD), 
           POP_AF = readr::parse_number(POP_AF), P_CONTAM = readr::parse_number(P_CONTAM), P_GERMLINE = readr::parse_number(P_GERMLINE),
           SF = readr::parse_number(SF), TLOD = readr::parse_number(TLOD))
  f <- f %>% dplyr::rename("chr" = V1, "pos" = V2, "ref" = V4, "alt" = V5)
  return(f)
}
