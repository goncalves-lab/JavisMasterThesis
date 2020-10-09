EuclideanPhysicaldistance <- function(ctxt){
  Physical_matrix <- matrix(rep(0,28), ncol = 2)
  for(i in 1:nrow(ctxt)){
  str <- str_split(row.names(ctxt)[i], pattern = "_") %>% .[[1]] %>% .[1]
  x <- str %>% str_split(pattern = "") %>% .[[1]] %>% .[1]
  x<- grep(x, LETTERS)
  y <- parse_number(str)
  Physical_matrix[i,1] = as.integer(x)
  Physical_matrix[i,2] = as.integer(y)
}
rownames(Physical_matrix) <- rownames(ctxt)

eumat_phy <- as.matrix(dist(Physical_matrix, method = "euclidean"))
eumat_immediate <- as.matrix(dist(ctxt, method = "euclidean"))
phy_and_sig <- function(x){
  mat_phy <- as.dist(eumat_phy[grepl(x, rownames(eumat_phy)),grepl(x, rownames(eumat_phy)) ])
  mat_eu <- as.dist(eumat_immediate[grepl(x, rownames(eumat_immediate)),grepl(x, rownames(eumat_immediate)) ])
  k <- vector()
  for(i in 1:length(names(mat_phy))){
    k <- c(k, paste(names(mat_phy)[i], names(mat_phy)[i+1]), paste(names(mat_phy)[i], names(mat_phy)[i+2]), paste(names(mat_phy)[i], names(mat_phy)[i+3])
           , paste(names(mat_phy)[i], names(mat_phy)[i+3]),  paste(names(mat_phy)[i], names(mat_phy)[i+4]),  paste(names(mat_phy)[i], names(mat_phy)[i+5]))
    k <- k[!grepl("NA", k)] %>% unique()
  }
  data.frame("physical" = as.vector(mat_phy),"euclidean" = as.vector(mat_eu), "sample" = k) %>%
    .[(.$physical != 0),]}
hgla_treated_phy <- phy_and_sig("hgla_treated")
hgla_control_phy <- phy_and_sig("hgla_control")
mmus_treated_phy <- phy_and_sig("mmus_treated")
mmus_control_phy <- phy_and_sig("mmus_control")

all_phy <- rbind(hgla_treated_phy, hgla_control_phy, mmus_treated_phy, mmus_control_phy) %>% separate(sample, into = c("sample1", "sample2"), sep = " ") %>% mutate("sample2" = gsub("(hgla|mmus)_(control|treated)", "", sample2) %>%  gsub("_", "", .)) %>%
  separate(sample1, into = c("sample1", "species", "treatment"))
}


lm_eqn <- function(df){
    m <- lm(physical ~ euclidean, df);
    eq <- paste0("y = ", format(unname(coef(m)[1]), digits = 2), "+", format(unname(coef(m)[2]), digits = 2), "x ,  R^2 =", format(summary(m)$r.squared, digits = 3))
    return(eq)
}
