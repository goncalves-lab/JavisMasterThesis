#' Obtaining a Dataframe containing the average context per group (treatment).
#'
#' @param ctxt Output of the function contextmatrix..


meanContext <- function(ctxt){
  require(tidyverse)
  a <- ctxt %>% group_by(SBS2, ctxt, species, treatment) %>%
    summarise("value" = mean(value)) %>% ungroup()
  
  b <- a %>% group_by(species) %>% pivot_wider(names_from ="treatment", values_from = "value") %>% 
    mutate("subtract" = treated - control) %>% 
    pivot_longer(cols = c(unique(a$treatment), "subtract"), names_to = "treatment")
}
