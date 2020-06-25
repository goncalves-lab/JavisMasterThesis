  library(SomaticSignatures)
  library(tidyverse)
  library(patchwork)
  library(ggtext)
  #All samples here need to be first selected for pass mutations with egrep "PASS|##|#CHROM" SAMPLE.vcf > SAMPLE.pass.vcf
  
  path_treated <- "/icgc/dkfzlsdf/analysis/B210/Javi/vcf_analysis/Treated_pass_vcfs/"
  path_untreated <- "/icgc/dkfzlsdf/analysis/B210/Javi/vcf_analysis/Untreated_pass_vcfs/"
  
  #Treated_vcfs <- list.files(path = path_treated , pattern = "Newpass", full.names = FALSE)
  Treated_vcfs <- c("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AS-452423/pass.vcf", "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AS-452425/pass.vcf", "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AS-452427/pass.vcf", "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AS-452429/pass.vcf")
  #Treated_vcfs <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/18431_T/", pattern = "pass.vcf", full.names = TRUE )
  #setwd(path_treated)
  Vranges_T_2 <- list()
  Vranges_T_1 <- list()
  Vranges_T <- list()
  for (i in Treated_vcfs){ 
    Vranges_T_1[[i]] <- readVcfAsVRanges(i)
    Vranges_T_2[[i]] <- Vranges_T_1[[i]][Vranges_T_1[[i]]@sampleNames == "TUMOR"]
    Vranges_T[[i]] <- Vranges_T_2[[i]][nchar(Vranges_T_2[[i]]@ref) == 1 & nchar(Vranges_T_2[[i]]@alt) == 1]
    print(sum(!is.na(Vranges_T[[i]]@ref)))
    Vranges_T[[i]]$sample <- i
    print(Vranges_T[[i]][1,]$sample)
  }
  
  #Untreated_vcfs <- list.files(path = path_untreated, pattern = "pass.vcf", full.names = FALSE)
  Untreated_vcfs <- ("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AS-452431/pass.vcf")
  #Untreated_vcfs <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/18431_C/", pattern = "pass.vcf", full.names = TRUE )
  #setwd(path_untreated)
  Vranges_U_2 <- list()
  Vranges_U_1 <- list()
  Vranges_U <- list()
  for (i in Untreated_vcfs){
    Vranges_U_1[[i]] <- readVcfAsVRanges(i)
    Vranges_U_2[[i]] <- Vranges_U_1[[i]][Vranges_U_1[[i]]@sampleNames == "TUMOR"]
    Vranges_U[[i]] <- Vranges_U_2[[i]][nchar(Vranges_U_2[[i]]@ref) == 1 & nchar(Vranges_U_2[[i]]@alt) == 1 ]
    print(sum(!is.na(Vranges_U[[i]]@ref)))
    Vranges_U[[i]]$sample <- i
    print(Vranges_U[[i]][1,]$sample)
  }
  
  
  #dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Heterocephalus_glaber_female.HetGla_female_1.0.dna_rm.toplevel_short_ids.fa")
  dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Mus_musculus.GRCm38.dna_rm.toplevel_short_IDs.fa")
  ctxt_T <-list()
  ctxt_mat_T <- matrix(0, ncol = length(Treated_vcfs), nrow = 96)
  ctxt_mat_T_mean <- matrix (0, ncol = 1, nrow = 96)
  for (i in 1:length(Vranges_T)) {
      ctxt_T[[i]] <- mutationContext(Vranges_T[[i]], ref = dna)
      ctxt_mat_T[,i] <- motifMatrix(ctxt_T[[i]], group = "sample", normalize = TRUE)
      print(ctxt_mat_T[1,i])
  }
  for (i in 1:nrow(ctxt_mat_T)) {
        ctxt_mat_T_mean[i,] = mean(ctxt_mat_T[i,])
  }  
      
  ctxt_U <-list()
  ctxt_mat_U <- matrix(0, ncol = length(Untreated_vcfs), nrow = 96)
  ctxt_mat_U_mean <- matrix(0, ncol = 1, nrow = 96)
  for (i in 1:length(Vranges_U)) {
    ctxt_U[[i]] <- mutationContext(Vranges_U[[i]], ref = dna)
    ctxt_mat_U[,i] <- motifMatrix(ctxt_U[[i]], group = "sample", normalize = TRUE)
    print(ctxt_mat_T[1,i])
  }
  for (i in 1:nrow(ctxt_mat_U)) {
    ctxt_mat_U_mean[i,] = mean(ctxt_mat_U[i,])
  } 
  
  SBS_96 <- rep(c("A.A", "A.C","A.G","A.T","C.A", "C.C","C.G","C.T","G.A", "G.C","G.G","G.T","T.A", "T.C","T.G","T.T"), 6)
  SBS_6 <- c(rep("C>A",16), rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
  result = vector()
  for (i in 1:length(Treated_vcfs)) {
    result = c(result, sum((!is.na(Vranges_T[[i]]@ref))))
    print(result)
    meanSNV_Treated = sum(result)/length(Treated_vcfs)
    print(meanSNV_Treated)
  }
  
  Treated_matrix_df <- data.frame(SBS = SBS_96, value = ctxt_mat_T_mean[,1], ref = SBS_6, condition = "treated")
  Treated_matrix_df_N <- data.frame(SBS = SBS_96, value = ctxt_mat_T_mean[,1]*meanSNV_Treated, ref = SBS_6, condition = "treated")
  
  theme_BM <- function (base_size = 11, base_family = "") 
  {
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
      theme(axis.ticks = element_line(colour = "black"),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.backgroun = element_blank(), 
            axis.line	= element_line(colour = "black", size = 0.5),
            axis.ticks.x = element_line(colour = "black", size = 0.5),
            axis.text = element_text(size = rel(0.8), colour = "black"), 
            strip.background = element_blank()
      )
  }
  
  Treated <- ggplot(Treated_matrix_df, aes(x= SBS, y= value, fill = ref))+ geom_col()+ 
    facet_grid(.~ ref)+
    guides(fill = FALSE)+ 
    theme_BM()  + 
    expand_limits(y = 0) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    ylab("Relative contribution") + xlab("96-trinucleotide context")+
    scale_fill_brewer(palette="Set2")
  
  result2 = vector()
  for (i in 1:length(Untreated_vcfs)) {
    result2 = c(result2, sum((!is.na(Vranges_U[[i]]@ref))))
    print(result2)
    meanSNV_Untreated = sum(result2)/length(Untreated_vcfs)
    print(meanSNV_Untreated)
  }
  
  Control_matrix_df <- data.frame(SBS = SBS_96, value = ctxt_mat_U_mean[,1], ref = SBS_6, condition = "Untreated")
  
  Control_matrix_df_N <- data.frame(SBS = SBS_96, value = ctxt_mat_U_mean[,1]*meanSNV_Untreated, ref = SBS_6, condition = "Untreated")
  
  Untreated <- ggplot(Control_matrix_df, aes(x= SBS, y= value, fill = ref))+ geom_col()+ 
    facet_grid(.~ ref)+
    guides(fill = FALSE)+ 
    theme_BM()  + 
    expand_limits(y = 0) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    ylab("Relative contribution") + xlab("96-trinucleotide context")+
    scale_fill_brewer(palette="Pastel2")
  
  matrix_substract = ctxt_mat_T_mean- ctxt_mat_U_mean
  matrix_substract_N = ctxt_mat_T_mean*meanSNV_Treated- ctxt_mat_U_mean*meanSNV_Untreated
  substract_matrix_df <- data.frame(SBS = SBS_96, value = matrix_substract[,1], ref = SBS_6, condition = "Treated - Untreated")
  substract_matrix_df_N <- data.frame(SBS = SBS_96, value = matrix_substract_N[,1], ref = SBS_6, condition = "Treated - Untreated")
  
  Substract <- ggplot(substract_matrix_df, aes(x= SBS, y= value, fill = ref))+ geom_col()+ 
    facet_grid(.~ ref)+
    guides(fill = FALSE)+ 
    theme_BM()  + 
    expand_limits(y = 0) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    ylab("Relative contribution") + xlab("96-trinucleotide context")+
    scale_fill_brewer(palette="Set2")
  
  Treated / Untreated / Substract
  
  
  #We already have the substracted mutational profile, now, the barplot with the two main SBS.
  df_list_T <- list()
  sumaT.N = vector()
  sumaC.N = vector()
  for (i in 1:length(Treated_vcfs)) {
    df_list_T[[i]] <- data.frame(sample = "a" , SBS = c("T>N", "C>N"), value = c(sum(Vranges_T[[i]]@ref == "A" | Vranges_T[[i]]@ref == "T"), sum(Vranges_T[[i]]@ref == "C" | Vranges_T[[i]]@ref == "G")))
      df_list_T[[i]]$sample <- Treated_vcfs[[i]]
    print(df_list_T[[i]])
    sumaT.N = c(sumaT.N, df_list_T[[i]][1,]$value )%>% print()
    meanT.N = mean(sumaT.N) %>% print()
    sdT.N = sd(sumaT.N)%>% print()
    sumaC.N = c(sumaC.N, df_list_T[[i]][2,]$value )%>% print()
    meanC.N = mean(sumaC.N)%>% print()
    sdC.N = sd(sumaC.N)%>% print()
  }  
  df_list_U <- list()
  sumaT.N_U = vector()
  sumaC.N_U = vector()
  for (i in 1:length(Untreated_vcfs)) {
    df_list_U[[i]] <- data.frame(sample = "a" , SBS = c("T>N", "C>N"), value = c(sum(Vranges_U[[i]]@ref == "A" | Vranges_U[[i]]@ref == "T"), sum(Vranges_U[[i]]@ref == "C" | Vranges_U[[i]]@ref == "G")))
    df_list_U[[i]]$sample <- Untreated_vcfs[[i]]
    print(df_list_U[[i]])
    sumaT.N_U = c(sumaT.N_U, df_list_U[[i]][1,]$value )%>% print()
    meanT.N_U = mean(sumaT.N_U) %>% print()
    sdT.N_U = sd(sumaT.N_U)%>% print()
    sumaC.N_U = c(sumaC.N_U, df_list_U[[i]][2,]$value )%>% print()
    meanC.N_U = mean(sumaC.N_U)%>% print()
    sdC.N_U = sd(sumaC.N_U)%>% print()
  }  
  Treated_df = data.frame( condition = "Treated", SBS = c("T>N", "C>N"), mean = c(meanT.N, meanC.N), sd = c(sdT.N, sdC.N))
  
  Control_df = data.frame(condition = "Untreated",SBS = c("T>N", "C>N"), mean = c(meanT.N_U, meanC.N_U), sd = c(sdT.N_U, sdC.N_U))
  
  data <- rbind(Control_df, Treated_df)
  ymin = c((Control_df[1,]$mean - Control_df[1,]$sd), (Control_df[1,]$mean + Control_df[2,]$mean - Control_df[2,]$sd),(Treated_df[1,]$mean - Treated_df[1,]$sd), (Treated_df[1,]$mean + Treated_df[2,]$mean - Treated_df[2,]$sd))
  ymax = c((Control_df[1,]$mean + Control_df[1,]$sd), (Control_df[1,]$mean + Control_df[2,]$mean + Control_df[2,]$sd),(Treated_df[1,]$mean + Treated_df[1,]$sd), (Treated_df[1,]$mean + Treated_df[2,]$mean + Treated_df[2,]$sd))
  ggplot(data, aes(fill=SBS, y=mean, x=condition))+
    geom_bar(position = "stack", stat = "identity")+
    geom_errorbar(aes(ymin= ymin , ymax = ymax), width = .1)
  sample_names_True <- c("TRUE", "FALSE")
  names(sample_names_True) <- c("Untreated", "Treated")
  data_per_sample <- rbind(bind_rows(df_list_T),bind_rows(df_list_U))
  data_per_sample[1:12,] %>% ggplot(aes(x=sample, y = value, fill = SBS))+geom_col(position = "stack")+facet_grid(~(sample=="AS-422415.pass.vcf"))+theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text.x = element_blank())+xlab("")+ylab("Number of mutations")

  
 persample_mat <- as.data.frame(ctxt_mat_T) %>% mutate(SBS = SBS_6) %>% gather(sample, value , V1, V2, V3, V4)
  
  
  
  
  
  
  