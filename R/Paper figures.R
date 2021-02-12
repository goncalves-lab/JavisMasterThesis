################################################
#                                              #
# Script with the figures needed               #
# for the paper. Includes the following plot:  #
# VAF distribution,                            #
# signature similarity,                        # 
# signature plot in interval                   #
# signature fitting                            # 
################################################

set.seed(123)

##1. Functions
#All the functions and packages used for the plots.

#Packages required
packages <- c("tidyverse", "SomaticSignatures", "MutationalPatterns", "siglasso")

sapply(packages, require, character.only = TRUE)

#Functions
sapply(list.files("Functions paper/", full.names = T), source)

#Analyses and plotting


### The first analyses are about the numbers. Number of mutations, mutational burden, variant allelic frequency distribution.
### Therefore, I used the subsampled version for those analyses. 
### Using a table that is the final output of the Somatic variant calling pipeline (github.com/JaviBotey/MasterThesis/Somatic Variant Calling Pipeline/Snakefile)

hgla_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Treated/100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) %>% 
  .[-grep("AS-422", .)]
hgla_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Control/100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) %>%
  .[-grep("AS-422", .)]
mmus_treated <-  list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Treated_100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE)
mmus_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Control_100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) 

### The genome of both species in the form of a FaFile is needed for obtaining the mutational context.
hgla_dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Heterocephalus_glaber_female.HetGla_female_1.0.dna_rm.toplevel_short_ids.fa")
mmus_dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Mus_musculus.GRCm38.dna_rm.toplevel_short_IDs.fa")

### The same way of reading is used for all the samples.

AF<- list()
VR <- list()
CTXT <- list()

for (i in c( "hgla_treated", "hgla_control", "mmus_control", "mmus_treated")){
  
  #Get species and treatment
  species <- str_match(i, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", i)
  
  #Transform AF to table
  AF[[i]] <- lapply(get(i), readAF, species = species, treatment = treatment, input = "AF") #The function readAF to read the table
  
  #AF to Vranges for context extraction
  VR[[i]] <- lapply(AF[[i]], AFtoVranges)
  
  #Get the context
  CTXT[[i]] <- lapply(VR[[i]], contextmatrix, species = species, treatment = treatment) # And a context matrix is obtained for each sample
}

AF_b <- bind_rows(AF) # All samples are bound into a data set of AF, contexts
CTXT_b <- bind_rows(CTXT)
CTXT_b_mean <- meanContext(CTXT_b) # The mean context for each species and treatment is obtained with this function and a subtracted average context

#Figure . VAF distribution 
pdf("AF distribution.pdf")

ggplot(AF_b, aes(x=af, fill=sampleID, apha = .1))+
  geom_density()+my_theme()+scale_fill_my("Sunset")+facet_wrap(treatment~species)

dev.off()

#Figure . Signature similarity. The box

where<- WhereIsTheSignature(AF = AF_b,
                            width = 0.01, CosOrTA = "cos", quantiles = T)

where <- where %>% group_by(species) %>% mutate(int0 = lag(int, default = 0))

pdf("Signature similarity box.pdf")

ggplot(where, aes(x=int, fill=pmax(value,0)))+
  geom_rect(aes(xmin = int0, xmax = int, ymin = 0, ymax = 1))+
  facet_wrap(~species, ncol = 1 )+
  scale_fill_gradient(low = my_palette[["Sunset"]][1], high = my_palette[["Sunset"]][6])+
  my_theme()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))

dev.off()

#Figure . Signature spectrum plot for specific AF ranges.

#Setting the lower and upper AF thresholds for hgla
h_thr_l<- 0.2
h_thr_u<- 0.25

pdf("Signature spectrum NMR.pdf")
#Plotting
bind_rows(AF_b %>%
            subset(species=="hgla" & treatment == "treated"  & af > h_thr_l & af < h_thr_u) %>%
            AFtoVranges()%>%
            contextmatrix(species = "hgla", treatment = "treated"), AF_b %>%
            subset(species=="hgla" & treatment == "control"   & af > h_thr_l & af < h_thr_u)%>%
            AFtoVranges()%>%
            contextmatrix(species = "hgla", treatment = "control"))%>% meanContext() %>% 
  mutate("value" = pmax(value, 0)) %>% subset(treatment == "subtract")%>% SignaturePlot()+
  facet_grid(treatment~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim = c(0,1.5))

dev.off()
#Setting the lower and upper AF thresholds for mmus
m_thr_l<- 0.12
m_thr_u<- 0.2

pdf("Signature spectrum mouse.pdf")
bind_rows(AF_b %>% subset(species=="mmus" & treatment == "treated"  & af > m_thr_l & af < m_thr_u)%>% AFtoVranges()%>%
            contextmatrix(species = "mmus", treatment = "treated"), AF_b %>%
            subset(species=="mmus" & treatment == "control"   & af > m_thr_l & af < m_thr_u)%>%
            AFtoVranges()%>%
            contextmatrix(species = "mmus", treatment = "control"))%>% meanContext() %>%
  mutate("value" = pmax(value, 0)) %>% subset(treatment == "subtract")%>% SignaturePlot()+
  facet_grid(treatment~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim = c(0,1.5))
dev.off()

#Full depth for mutational signatures
#Moving on now to the full depth version of the sequencing, more suitable for identifying signatures and mutated genes


hgla_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Treated/full_depth/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE)
hgla_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Control/full_depth/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) 
mmus_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF", pattern = "G1", full.names = TRUE)
mmus_treated <-  list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF", pattern = "B12|D2|H12|H5", full.names = TRUE)

#The same as before is followed
### The same way of reading is used for all the samples.

AF<- list()
VR <- list()
CTXT <- list()

for (i in c( "hgla_treated", "hgla_control", "mmus_control", "mmus_treated")){
  
  #Get species and treatment
  species <- str_match(i, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", i)
  
  #Transform AF to table
  AF[[i]] <- lapply(get(i), readAF, species = species, treatment = treatment, input = "AF") #The function readAF to read the table
  
  #AF to Vranges for context extraction
  VR[[i]] <- lapply(AF[[i]], AFtoVranges)
  
  #Get the context
  CTXT[[i]] <- lapply(VR[[i]], contextmatrix, species = species, treatment = treatment) # And a context matrix is obtained for each sample
}

AF_b <- bind_rows(AF) # All samples are bound into a data set of AF, contexts
CTXT_b <- bind_rows(CTXT)
CTXT_b_mean <- meanContext(CTXT_b) # The mean context for each species and treatment is obtained with this function and a subtracted average context




#Reading now the mutations in MEF, to extract the mutational signature

mmus_MEF_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/Control1/AF", pattern = "a.tsv",  full.names = TRUE, recursive = T)
mmus_MEF_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/Control1/AF", pattern = "o.tsv",  full.names = TRUE, recursive = T)

AF<- list()
for (i in c("mmus_MEF_control", "mmus_MEF_treated")){
  species <- str_match(i, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", i)
  AF[[i]] <- lapply(get(i), readAF, species = species, treatment = treatment, input = "AF")
}
AF_b_MEF <- bind_rows(AF)

#Extracting the mutational signatures found in MEF 
ctxt_to_extract <- AF_b_MEF %>% 
  slice_max(af, prop = 0.01)%>% AFtoVranges() %>% contextmatrix(species = "mmus", treatment = "treated") %>% 
  mutate("sampleID" = translate_names(sampleID))%>%
  mutate("sampleID" = paste(species, treatment, sampleID)) %>% dplyr::select(-c(species, treatment)) %>%  
  pivot_wider(names_from = "sampleID", values_from = "value") %>% dplyr::select(-c(SBS2, ctxt))
sigs <- SomaticSignatures::identifySignatures(as.matrix(ctxt_to_extract), 3, nmfDecomposition, method = "lee") #Selected 3 signatures to extract, could do also 2,4 or 5

#Plotting the mutational signatures (To check if the signature corresponds to the expected)
SBS_96 <- rep(c("A.A", "A.C","A.G","A.T","C.A", "C.C","C.G","C.T","G.A", "G.C","G.G","G.T","T.A", "T.C","T.G","T.T"), 6)
SBS_6 <- c(rep("C>A",16), rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
sigs_df <- sigs@signatures %>% as.data.frame() %>% mutate("SBS" = SBS_6, "ctxt" = SBS_96) %>% pivot_longer(cols = colnames(as.data.frame(sigs@signatures)))
ggplot(sigs_df, aes(x=ctxt, y=value*100, fill=SBS))+geom_col()+facet_grid(name~SBS)+my_theme()+scale_fill_my("Sunset")+
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = 0)

#Fitting the signatures to the skin samples  

ctxt_to_fit <- CTXT_b  %>%  mutate("sampleID" = paste(species, treatment, sampleID)) %>% 
  dplyr::select(-c(species, treatment)) %>%  pivot_wider(names_from = "sampleID", values_from = "value") %>% 
  dplyr::select(-c(SBS2, ctxt))
DMBA_cosmic <- cbind(as.matrix(sigs@signature), siglasso::cosmic_v3_exo[,c(7:10,43)]) %>%
  as.matrix() #Adding the COSMIC signatures to fit them witht the extracted signatures

#Figure . Signature fitting to mouse and NMR skin samples

pdf("Fitted signatures.pdf")
MutationalPatterns::fit_to_signatures(mut_matrix = as.matrix(ctxt_to_fit),
                                      signatures = DMBA_cosmic) %>%
  .$contribution %>% t() %>% as.data.frame() %>% mutate("sampleID" = rownames(.)) %>%  
  pivot_longer(names_to = "sig", cols = colnames(.)[-length(colnames(.))]) %>% 
  mutate("sig" = gsub("a|b|c|d", "", sig))%>%
  ggplot( aes(x=grepl("treated", sampleID), y=100*value, fill=sig))+
  geom_col(width = .8, position = "fill")+my_theme()+
  facet_wrap(~grepl("hgla", sampleID), nrow = 1)+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_my("Sunset")
dev.off()
