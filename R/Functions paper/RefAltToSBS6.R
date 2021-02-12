#' Converting two columns of a dataframe that contain the ref and alt bases into a column containing the 6 types of single base substitutions
#'
#' @param dataframe A dataframe containing two columns ref and alt


RefAltToSBS6 <- function(dataframe){
  dataframe %>% unite("SBS", c("ref", "alt"), sep = ">", remove = F) %>% 
    mutate(SBS = gsub("A>C", "T>G",.$SBS) %>% gsub("A>G", "T>C", .) %>% gsub("A>T", "T>A", .) %>% gsub("G>C", "C>G", .) %>% gsub("G>T", "C>A", .) %>% 
             gsub("G>A", "C>T", .))
}