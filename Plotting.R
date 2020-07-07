library("vcfR")
library("ggplot2")
library("pinfsc50")
library(data.table)
library(tidyverse)
library("ggtext")
samples_list <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Final_analysis/AF", pattern = ".tsv", full.names = T)
#samples_list <- c(list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/Treated", full.names = TRUE), list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/Control", full.names = TRUE))
path <- "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Final_analysis"
intervals <- list()
scaffold <- list()
empty_scaffold <- list()
AF <- list()
empty_AF <- list()
scaffold_AF <- list()
for (i in samples_list){
  AF[[i]] <- read.table(i)
  AF[[i]]$sample <- gsub( "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/(Control|Treated)/", "", i) %>% gsub(".tsv", "", . )
  AF[[i]]$state <- gsub("AS-4224[03,05,07,09,11]", "Treated" , AF[[i]]$sample) %>% gsub("AS-422415|AS-47513[1,3,5]", "Untreated" , .)
  #empty_AF[[i]] <- data.frame(Pos = AF[[i]]$V1, AF = AF[[i]]$V4 )
  #scaffold_AF[[i]] <- rbind(empty_scaffold[[i]], empty_AF[[i]]) %>% .[order(.$Pos),]
  #scaffold_AF[[i]]$sample <- i
}
theme_nothing <- theme(panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sample_new_names <- c("hgla_DT25_skin_T3_F7", "hgla_DT25_skin_T3_F8", "hgla_DT25_skin_T3_B2", "hgla_DT25_skin_T3_J4", "hgla_DT25_skin_T3_M11", "hgla_AC25_skin_C3_B2")
#sample_labels <- paste(gsub(".tsv", "", samples_list),":", lapply(AF, nrow), "mutations") %>% set_names(. , samples_list)
sample_labels <- paste(sample_new_names,":", lapply(AF, nrow), "mutations") %>% set_names(. , gsub( "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/(Control|Treated)/", "", samples_list) %>% gsub(".tsv", "", . ))


#Alternative to all the loops
samples_AF <- c(list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/Treated", full.names = TRUE), list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/all/AF/Control", full.names = TRUE))
AF_list <- lapply(samples_AF, function(x){
  dat <- read.table(x, header = FALSE, sep = "\t") %>% dplyr::rename(., "CHROM"=V1, "POS"=V2, "SBS"=V3,"AF"=V4) %>% mutate(sample = gsub("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF/", "", x)) 
})
AF_df <- bind_rows(AF_list) %>% ggplot()+geom_density(aes(x=AF, fill=sample))+facet_wrap(~sample)+coord_cartesian(xlim=c(0,1))+theme(legend.position = 0)
AF_df
AF_df_experiment <- list()
AF_df_experiment_sigs <- list()
AF_df_experiment_sigs_binded <- list()
AF_df_experiment_sigs_sliced <- list()
for (f in c(0.05,0.1,0.3,0.5,1)){
  AF_df_experiment[[as.character(f)]] <- lapply(AF_list, function(x){dat <- subset(x, AF < f)})
  AF_df_experiment_sigs[[as.character(f)]] <- lapply(AF_df_experiment[[as.character(f)]], function(x){ dat <- data.frame("C>G"=sum(x$SBS == "C>G"| x$SBS == "G>C"), "C>T"=sum(x$SBS == "C>T"| x$SBS == "G>A"), "C>A"=sum(x$SBS == "C>A"| x$SBS == "G>A"), "T>G"=sum(x$SBS == "T>G"| x$SBS == "A>C"), "T>C"=sum(x$SBS == "T>C"| x$SBS == "A>G"), "T>A"=sum(x$SBS == "T>A"| x$SBS == "A>T"), "sample"=x$sample)})
  AF_df_experiment_sigs_sliced[[as.character(f)]] <- lapply(AF_df_experiment_sigs[[as.character(f)]], function(x){slice(x, 1)})
  AF_df_experiment_sigs_binded[[as.character(f)]] <- bind_rows(AF_df_experiment_sigs_sliced[[as.character(f)]])
  AF_df_experiment_sigs_binded[[as.character(f)]]$thr <- as.character(f)
}
AF_df_experiment_sigs_binded_binded <- bind_rows(AF_df_experiment_sigs_binded) %>%  gather(key = "SBS", "value", C.G, C.T, C.A, T.G, T.C, T.A) 
ggplot(AF_df_experiment_sigs_binded_binded, aes(x= thr, y=value, fill = SBS))+geom_col(position = "fill")+facet_wrap(~sample)



plot_intervals <- bind_rows(intervals)
intervals_only_mut <- plot_intervals %>% subset(V1 == "JHQ*")
ggplot(intervals_only_mut, aes(x=V1, y = V7, fill = sample))+geom_point()
ggplot(plot_intervals, aes(x=V1, y=V7, fill = sample))+geom_point()+facet_wrap(~sample, labeller = labeller(sample = sample_labels))+
  theme_nothing+theme(axis.text.x = element_blank(), axis.ticks=element_blank())+coord_cartesian(ylim = c(0,12))+ggtitle("Number of mutations per genomic bin")+xlab("1 kb genomic bin")+ylab("Number of mutations")
ggsave("plot_intervals.png")

plot_scaffold <- bind_rows(scaffold) %>% mutate(V1 = V1 %/% 1000) %>% ggplot(aes(x=V1, y = V4))+geom_col()+facet_wrap(~sample)+theme_Pdro+coord_cartesian(ylim = c(0,10))

plot_AF_scaffols <- bind_rows(scaffold_AF) %>% ggplot(aes(x=Pos, y=AF, fill = (sample == "AS-422415.tsv")))+geom_col()+theme_Pdro+coord_cartesian(ylim = c(0,5))
number_mutations <- data.frame(sample = as.vector(samples_list), n = as.vector(as.numeric(lapply(AF, nrow))) )
ggplot(number_mutations, aes(x=sample, y=n, fill=sample))+geom_col()+theme_Pdro
ggsave("Number_mutations.png")
AF_binded <- bind_rows(AF) %>% dplyr::rename("AF" = V4, "SBS" = V3, "CHROM" = V1, "POS" = V2)
AF_binded %>% subset(grepl("AS-422", sample)) %>% ggplot(aes(x=AF, fill = (grepl("AS-422415", sample))))+geom_density()+coord_cartesian(xlim= c(0,1), ylim = c(0,30))+facet_wrap(~ sample,  labeller = labeller(sample = sample_labels))+theme(legend.position = 0, panel.grid = element_blank(), panel.background = element_blank())+xlab("Variant allelic frequency")

AF_binded %>% subset(sample == "AS-422403.tsv"| sample == "AS-422405.tsv"| sample == "AS-422407.tsv"| sample == "AS-422409.tsv"| sample == "AS-422411.tsv") %>% ggplot(aes(x=V4, fill = "s"))+geom_density()+facet_wrap(~sample,  labeller = labeller(sample = sample_labels))+theme(legend.position = 0)+coord_cartesian(ylim = c(0,9))
AF_binded %>% subset(sample == "AS-422415.tsv"| sample == "AS-475131.tsv"| sample == "AS-475133.tsv"| sample == "AS-475135.tsv") %>% ggplot(aes(x=V4, fill = (sample != "A")))+geom_density()+facet_wrap(~sample,  labeller = labeller(sample = sample_labels))+theme(legend.position = 0)+coord_cartesian(ylim = c(0,9))
ggsave("VAF_density.png")
#TRying to get the log2 FC
bind_cols(scaffold) %>% mutate("AS-422403" = V4, "AS-422405" =V41, "AS-422407" = V42, "AS-422409" = V43, "AS-422411" = V44, "AS-422415" = V45) %>%  select(V1, "AS-422403", "AS-422405", "AS-422407", "AS-422409", "AS-422411", "AS-422415") %>% mutate(log2 = log2(mean(c("AS-422403", "AS-422405", "AS-422407", "AS-422409", "AS-422411"))/"AS-422415"))


#The VAF experiment

AF_0.05 <- list()
df_0.05 <- list()
AF_0.1 <- list()
df_0.1 <- list()
AF_1 <- list()
df_1 <- list()
AF_0.4 <- list()
df_0.4 <- list()
AF_0.5 <- list()
df_0.5 <- list()
for (f in samples_list){
    AF_0.05[[f]] <- AF[[f]] %>% subset(V4<0.05)
    df_0.05[[f]] <- data.frame("C>G"=sum(AF_0.05[[f]]$V3 == "C>G"| AF_0.05[[f]]$V3 == "G>C"), "C>T"=sum(AF_0.05[[f]]$V3 == "C>T"| AF_0.05[[f]]$V3 == "G>A"), "C>A"=sum(AF_0.05[[f]]$V3 == "C>A"| AF_0.05[[f]]$V3 == "G>A"), "T>G"=sum(AF_0.05[[f]]$V3 == "T>G"| AF_0.05[[f]]$V3 == "A>C"), "T>C"=sum(AF_0.05[[f]]$V3 == "T>C"| AF_0.05[[f]]$V3 == "A>G"), "T>A"=sum(AF_0.05[[f]]$V3 == "T>A"| AF_0.05[[f]]$V3 == "A>T"))
    df_0.05[[f]]$thr <- 0.05 
    df_0.05[[f]]$sample <- f
    AF_0.1[[f]] <- AF[[f]] %>% subset(V4<0.1)
    df_0.1[[f]] <- data.frame("C>G"=sum(AF_0.1[[f]]$V3 == "C>G"| AF_0.1[[f]]$V3 == "G>C"), "C>T"=sum(AF_0.1[[f]]$V3 == "C>T"| AF_0.1[[f]]$V3 == "G>A"), "C>A"=sum(AF_0.1[[f]]$V3 == "C>A"| AF_0.1[[f]]$V3 == "G>A"), "T>G"=sum(AF_0.1[[f]]$V3 == "T>G"| AF_0.1[[f]]$V3 == "A>C"), "T>C"=sum(AF_0.1[[f]]$V3 == "T>C"| AF_0.1[[f]]$V3 == "A>G"), "T>A"=sum(AF_0.1[[f]]$V3 == "T>A"| AF_0.1[[f]]$V3 == "A>T"))
    df_0.1[[f]]$thr <- 0.1
    df_0.1[[f]]$sample <- f
    AF_1[[f]] <- AF[[f]] %>% subset(V4<1)
    df_1[[f]] <- data.frame("C>G"=sum(AF_1[[f]]$V3 == "C>G"| AF_1[[f]]$V3 == "G>C"), "C>T"=sum(AF_1[[f]]$V3 == "C>T"| AF_1[[f]]$V3 == "G>A"), "C>A"=sum(AF_1[[f]]$V3 == "C>A"| AF_1[[f]]$V3 == "G>A"), "T>G"=sum(AF_1[[f]]$V3 == "T>G"| AF_1[[f]]$V3 == "A>C"), "T>C"=sum(AF_1[[f]]$V3 == "T>C"| AF_1[[f]]$V3 == "A>G"), "T>A"=sum(AF_1[[f]]$V3 == "T>A"| AF_1[[f]]$V3 == "A>T"))
    df_1[[f]]$thr <- 1
    df_1[[f]]$sample <- f
    AF_0.5[[f]] <- AF[[f]] %>% subset(V4<0.5)
    df_0.5[[f]] <- data.frame("C>G"=sum(AF_0.5[[f]]$V3 == "C>G"| AF_0.5[[f]]$V3 == "G>C"), "C>T"=sum(AF_0.5[[f]]$V3 == "C>T"| AF_0.5[[f]]$V3 == "G>A"), "C>A"=sum(AF_0.5[[f]]$V3 == "C>A"| AF_0.5[[f]]$V3 == "G>A"), "T>G"=sum(AF_0.5[[f]]$V3 == "T>G"| AF_0.5[[f]]$V3 == "A>C"), "T>C"=sum(AF_0.5[[f]]$V3 == "T>C"| AF_0.5[[f]]$V3 == "A>G"), "T>A"=sum(AF_0.5[[f]]$V3 == "T>A"| AF_0.5[[f]]$V3 == "A>T"))
    df_0.5[[f]]$thr <- 0.5
    df_0.5[[f]]$sample <- f
    AF_0.4[[f]] <- AF[[f]] %>% subset(V4<0.4)
    df_0.4[[f]] <- data.frame("C>G"=sum(AF_0.4[[f]]$V3 == "C>G"| AF_0.4[[f]]$V3 == "G>C"), "C>T"=sum(AF_0.4[[f]]$V3 == "C>T"| AF_0.4[[f]]$V3 == "G>A"), "C>A"=sum(AF_0.4[[f]]$V3 == "C>A"| AF_0.4[[f]]$V3 == "G>A"), "T>G"=sum(AF_0.4[[f]]$V3 == "T>G"| AF_0.4[[f]]$V3 == "A>C"), "T>C"=sum(AF_0.4[[f]]$V3 == "T>C"| AF_0.4[[f]]$V3 == "A>G"), "T>A"=sum(AF_0.4[[f]]$V3 == "T>A"| AF_0.4[[f]]$V3 == "A>T"))
    df_0.4[[f]]$thr <- 0.4
    df_0.4[[f]]$sample <- f
}
df_0.05_binded <- bind_rows(df_0.05)
df_0.1_binded <- bind_rows(df_0.1)
df_1_binded <- bind_rows(df_1)
df_0.5_binded <- bind_rows(df_0.5)
df_0.4_binded <- bind_rows(df_0.4)
df_all_AF_binded <- bind_rows(df_0.05_binded , df_0.1_binded, df_1_binded, df_0.5_binded, df_0.4_binded)
df_all_AF_binded$thr <- as.character(df_all_AF_binded$thr)
df_all_AF_binded %>% 
  #subset(sample == "AS-422403.tsv"| sample == "AS-422405.tsv"| sample == "AS-422407.tsv"| sample == "AS-422409.tsv"| sample == "AS-422411.tsv" | sample == "AS-422415.tsv") %>% 
  gather(key = "SBS", "value", C.G, C.T, C.A, T.G, T.C, T.A)  %>% ggplot(aes(x= thr, y=value, fill = SBS))+geom_col(position = "fill")+facet_wrap(~sample)+xlab("VAF threshold")+ylab("")

DP <- list.files("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Final_analysis/DP/", full.names = TRUE)
DP_files <- lapply(DP, read.table)
names(DP_files) <- samples_list[6:9]
for(i in names(DP_files)){DP_files[[i]]$sample <- i}
DP_binded <- bind_rows(DP_files)
ggplot(DP_binded)+geom_density(aes(x=V5, fill = sample))+facet_wrap(~sample)+theme_minimal()
