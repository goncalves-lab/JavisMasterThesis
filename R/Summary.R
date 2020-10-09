################################################
#                                              #
# Script with all the analysis performed for   #
# my master thesis. Includes:                  #
# Visualizing the number of variants,          #
# VAF distribution,                            #
# mutational spectrum,                         # 
# signature extraction,                        #
# signature fitting,                           # 
# nonsynonymous mutant clones.                 # 
################################################

sapply(list.files('/icgc/dkfzlsdf/analysis/B210/Javi/code/SkinCancerEvolution/wes/R/Functions', full.names = T), source)

hgla_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Treated/full_depth/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE)
hgla_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Control/full_depth/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) 
mmus_treated <-  list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Treated_100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE)
mmus_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Control_100X/Mutation_calling/AF", pattern = ".tsv", full.names = TRUE) 
mmus_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF", pattern = "G1", full.names = TRUE)
mmus_treated <-  list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF", pattern = "B12|D2|H12|H5", full.names = TRUE)



hgla_treated_vcf <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Treated/full_depth/Mutation_calling", pattern = "pass.vcf", full.names = TRUE, recursive = T)
hgla_control_vcf <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Control/full_depth/Mutation_calling", pattern = "pass.vcf", full.names = TRUE, recursive = T) 
mmus_vcf <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling", pattern = "pass.vcf", full.names = TRUE, recursive = T) 
mmus_treated_vcf <- mmus_vcf[1:4]
mmus_control_vcf <- mmus_vcf[5]

mmus_MEF_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/Control-1_is_the_normal/AF", pattern = "a.tsv",  full.names = TRUE, recursive = T)
mmus_MEF_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/Control-1_is_the_normal/AF", pattern = "o.tsv",  full.names = TRUE, recursive = T)
mmus_MEF_vcf <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/PON", pattern = "pass2.vcf",  full.names = TRUE, recursive = T)
mmus_MEF_treated_vcf <- mmus_MEF_vcf[c(1,2,3)]
mmus_MEF_control_vcf <- mmus_MEF_vcf[c(3,5)]

hgla_F_treated <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/3o_is_the_normal/AF", pattern = "a.tsv",  full.names = TRUE, recursive = T)
hgla_F_control <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/3o_is_the_normal/AF", pattern = "o.tsv",  full.names = TRUE, recursive = T)
hgla_F_vcf <- list.files(path = "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/3o_is_the_normal", pattern = "pass.vcf",  full.names = TRUE, recursive = T)
hgla_F_control_vcf <- hgla_F_vcf[grepl("o/pass", hgla_F_vcf)]
hgla_F_treated_vcf <- hgla_F_vcf[grepl("a/pass", hgla_F_vcf)]

hgla_dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Heterocephalus_glaber_female.HetGla_female_1.0.dna_rm.toplevel_short_ids.fa")
mmus_dna <- FaFile("/icgc/dkfzlsdf/analysis/B210/references_genome/Mus_musculus.GRCm38.dna_rm.toplevel_short_IDs.fa")


AF<- list()
VR <- list()
CTXT <- list()
mutations <- list()
for (i in c( "hgla_treated", "hgla_control", "mmus_control", "mmus_treated")){
  species <- str_match(i, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", i)
  AF[[i]] <- lapply(get(i), readAF, species = species, treatment = treatment, input = "AF")
  mutations[[i]] <- lapply(AF[[i]], function(x){ 
      a <- x %>% unite("sampleID_treatment", c("sampleID", "treatment"), sep = "_")
      dplyr::select(a, c("sampleID_treatment", "chr", "pos", "ref", "alt"))}
      )
    VR[[i]] <- lapply(AF[[i]], function(x){
      a <- AFtoVranges(x) 
      b <- mutationDistance(a)
      a@elementMetadata@listData[["distance"]] <- b@elementMetadata@listData[["distance"]]
      return(a)
      })
    CTXT[[i]] <- lapply(VR[[i]], contextmatrix, species = species, treatment = treatment)
}

AF_b <- bind_rows(AF)
CTXT_b <- bind_rows(CTXT)
CTXT_b_mean <- meanContext(CTXT_b)

#Removing control positions from the naked mole rat fibroblasts

control_positions_F <- AF_b %>% unite(c(chr, pos),sep = "_", col =  "chr_pos") %>% 
  subset(treatment == "F_control")
treated_F_no_ctrl <- anti_join(AF_b %>% unite(c(chr, pos),sep = "_", col =  "chr_pos") %>% 
            subset(treatment == "F_treated"), control_positions_F, by = "chr_pos") %>% 
  separate("chr_pos", c("chr", "pos"), sep = "_")

AFtoVranges(treated_F_no_ctrl) %>% contextmatrix(species = "hgla", treatment = "treated") %>% 
  SignaturePlot()+facet_grid(sampleID~SBS2)+scale_fill_my("Sunset")

#Removing control positions from mouse fibroblasts
control_positions_MEF <- AF_b %>% unite(c(chr, pos),sep = "_", col =  "chr_pos") %>% 
  subset(treatment == "MEF_control")
treated_MEF_no_ctrl <- anti_join(AF_b %>% unite(c(chr, pos),sep = "_", col =  "chr_pos") %>% 
                                 subset(treatment == "MEF_treated"), control_positions_MEF, by = "chr_pos") %>% 
  separate("chr_pos", c("chr", "pos"), sep = "_")

AFtoVranges(treated_MEF_no_ctrl) %>% contextmatrix(species = "mmus", treatment = "treated") %>% 
  SignaturePlot()+facet_grid(sampleID~SBS2)+scale_fill_my("Sunset")


#Different SBS proportions in different AF?

ggplot(AF_b %>% subset(species=="hgla"), aes(x=af, fill=SBS))+geom_histogram(binwidth = .2, position = "fill")+
  facet_wrap(~treatment, ncol=1)+my_theme()+scale_fill_my("Sunset")+coord_cartesian(xlim=c(0.2,1))

AF_b %>% subset(species=="hgla" & af > .3 & af < .6 & !grepl("AS-42", sampleID))%>% AFtoVranges()%>%
  contextmatrix(species = "hgla", treatment = "treated")%>%SignaturePlot()+facet_grid(sampleID~ SBS2)
bind_rows(AF_b %>% subset(species=="mmus" & treatment == "treated" & !grepl("AS-42", sampleID) & af > .17 & af < .175 )%>% AFtoVranges()%>%
  contextmatrix(species = "mmus", treatment = "treated"),AF_b %>%
  subset(species=="mmus" & treatment == "control"  & !grepl("AS-42", sampleID) & af > .17 & af < .175)%>%
  AFtoVranges()%>%
  contextmatrix(species = "mmus", treatment = "control"))%>% meanContext()%>% 
  mutate("value" = pmax(value, 0))%>% subset(treatment == "subtract")%>%
  SignaturePlot()+
  facet_grid(treatment~SBS2)+my_theme()+scale_fill_my("Sunset")+
  theme(axis.text.x = element_text(angle=90, size = 5))

Background_ctxt_0.2_hgla <- subset(AF_b, species == "hgla" & treatment == "control"  & !grepl("AS-42", sampleID) & af > .2) %>% AFtoVranges() %>% contextmatrix(species = "hgla", treatment = "control") %>%  subset(treatment == "control") %>% group_by(species, SBS2, ctxt) %>% summarise(value = mean(value), .groups = "keep")
CTXT_b_background_removed_0.2_hgla <- subset(AF_b, species == "hgla" & treatment == "treated" & af > .2)%>% AFtoVranges() %>% contextmatrix(species = "hgla", treatment = "treated")%>% inner_join(subset(Background_ctxt_0.2_hgla), by= c("SBS2", "ctxt", "species")) %>% mutate("value" = value.x - value.y)
SignaturePlot(CTXT_b_background_removed_0.2_hgla %>% subset(species == "hgla" & !grepl("AS-42", sampleID)))+facet_grid(translate_names(gsub("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Treated/100X/Mutation_calling/AF/", "",sampleID))~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim=c(0,6))+theme(axis.text.x = element_text(size=5))

ggplot(AF_b %>% mutate("af_ints" = cut(af, 20)) %>% group_by(af_ints, treatment, species) %>%
         summarise("A_T" = count(SBS=="T>A")/n()),
          aes(x=af_ints, y = A_T, color = paste(treatment, species)))+geom_point()+stat_summary( geom="line")

where<- WhereIsTheSignature(AF = AF_b,
                              width = 0.01, CosOrTA = "cos", quantiles = T)
where <- where %>% group_by(species) %>% mutate(int0 = lag(int, default = 0))
ggplot(where, aes(x=int, fill=pmax(value,0)))+
  geom_rect(aes(xmin = int0, xmax = int, ymin = 0, ymax = 1))+
  facet_wrap(~species, ncol = 1 )+
  scale_fill_gradient(low = my_palette[["Sunset"]][1], high = my_palette[["Sunset"]][6])+
  my_theme()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
ggplot(where, aes(x=quant, y=100*value, color = species))+
  geom_line()

#Figure1 - Number of Somatic Variants

VCF <- list()
hgla_exome <- 62083075
mmus_exome <- 49370117
genome <- 2.7e9
exome <- c(hgla_exome, mmus_exome)
for (i in c("hgla_treated_vcf", "hgla_control_vcf", "mmus_treated_vcf", "mmus_control_vcf")){
  species <- str_match(i, "mmus|hgla")
  treatment <- str_match(i, "treated|control")
  VCF[[i]] <- lapply(get(i), function(x) {readPassVcf(x) %>% mutate("sampleID" = x)})
}
VCF_b <- bind_rows(VCF)
VCF_b_type <- mutate(VCF_b, "Indel" = (nchar(as.character(alt))+nchar(as.character(ref))) == 2)
VCF_b_type %>% 
  group_by(species, treatment, sampleID, Indel) %>% summarise("n" = n()) %>% 
  group_by(species, treatment, sampleID) %>% mutate(total = sum(n)) %>%
  group_by(species, treatment) %>% mutate("sd"=sd(total)) %>% 
  group_by(species, treatment, Indel) %>% 
  summarise("mean"=mean(n), "sd" = sd) %>% unique() %>%group_by(species, treatment) %>%
  mutate("total" = sum(mean)) %>% 
  ggplot()+geom_col(aes(fill = Indel, x=treatment, y=mean))+
  geom_errorbar(aes(x=treatment, ymin=total-sd, ymax= total +sd), 
                width = 0.3, size = .2)+
  facet_wrap(~species)+my_theme()+scale_fill_my("Sunset2")

AF_b %>% group_by(species,treatment, sampleID) %>% summarise(n=n()) %>% 
  group_by(species,treatment) %>% mutate(mean=mean(n), sd=sd(n)) %>% 
  subset( grepl("F_",treatment)==F) %>% subset(species=="hgla" & !grepl("AS-42", sampleID))%>%
  ggplot(aes(x=paste(treatment, species), y=(n/hgla_exome)*1000000))+
  geom_boxplot()+geom_point()+my_theme()+
  ylab("Mutations per Megabase")+xlab("")+coord_cartesian(ylim=c(0,500))
#Figure2 - VAF of Somatic Variants
ggplot(AF_b %>% subset(grepl("F_", treatment) == F)  , aes(x=paste(species, treatment), y=af, fill=treatment))+
  geom_violin()+geom_boxplot(width=.1, outlier.shape = 1, alpha = .1)+my_theme()+scale_fill_my("Sunset2")
ggplot(AF_b %>% subset(chr != "7" & !grepl("19|X", chr)), aes(x=af, fill=sampleID, apha = .1))+
  geom_density()+my_theme()+scale_fill_my("Sunset")+facet_wrap(treatment~species)

#Figure3 - 
 #Intervals
#Figure 4 - T > A numbers

CTXT_b %>% subset(!grepl("MEF", treatment)) %>% subset(SBS2 == "TA") %>% 
  group_by(species, treatment, ctxt) %>% summarise("sd"=sd(value), "mean" = mean(value)) %>% 
  ggplot()+geom_col(aes(fill = grepl("C.G", ctxt), x=treatment, y=mean*100))+
  facet_wrap(~species)+my_theme()+scale_fill_my("Sunset2")
#Figure 5 - Signatures
SignaturePlot(CTXT_b_mean %>% mutate("value" = pmax(value, 0)) %>% subset(species == "hgla" & treatment != "subtract"))+
  facet_grid(treatment~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim = c(0,8))
SignaturePlot(CTXT_b_mean %>% mutate("value" = pmax(value, 0)) %>% subset(species == "hgla" & treatment == "subtract"))+
  facet_grid(treatment~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim = c(0,2))
SignaturePlot(CTXT_b_background_removed %>% subset(species == "hgla"))+facet_grid(sampleID~SBS2)



#Plot per sample
SignaturePlot(CTXT_b)+facet_grid(sampleID~SBS2)+scale_fill_my("Sunset")

#Strand bias
ggplot(AF_b, aes(x=paste0(ref, ">", alt), fill = SBS))+geom_bar()+facet_grid(species+treatment~SBS, scales = "free")+my_theme()+theme(axis.text.x = element_text(angle = 90))

#Background mutational signature
Background_ctxt <- CTXT_b %>%  subset(treatment == "control") %>% group_by(species, SBS2, ctxt) %>% summarise(value = mean(value), .groups = "keep")
SignaturePlot(Background_ctxt)+facet_grid(species~SBS2)+scale_fill_brewer(palette = "Pastel2")


#Removing the background from each sample
CTXT_b_background_removed <- CTXT_b %>% subset(treatment == "treated") %>% inner_join(subset(Background_ctxt), by= c("SBS2", "ctxt", "species")) %>% mutate("value" = value.x - value.y)
SignaturePlot(CTXT_b_background_removed %>% subset(species == "hgla"))+facet_grid(sampleID~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim=c(-3,3))
SignaturePlot(CTXT_b_background_removed %>% subset(species == "mmus"))+facet_grid(sampleID~SBS2)+scale_fill_my("Sunset")+coord_cartesian(ylim=c(0,3))

#Extracting signatures from each sample

ctxt_to_extract <- subset(AF_b, species == "hgla" & grepl("F", treatment)) %>% 
  slice_max(af, prop = 0.01)%>% AFtoVranges() %>% contextmatrix(species = "hgla", treatment = "treated") %>% 
  mutate("sampleID" = translate_names(gsub("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/3o_is_the_normal/AF/", "",sampleID)))%>%
  mutate("sampleID" = paste(species, treatment, sampleID)) %>% dplyr::select(-c(species, treatment)) %>%  
  pivot_wider(names_from = "sampleID", values_from = "value") %>% dplyr::select(-c(SBS2, ctxt))
sigs <- MutationalPatterns::extract_signatures(as.matrix(ctxt_to_extract), rank = 3)
sigs <- SomaticSignatures::identifySignatures(as.matrix(ctxt_to_extract)+0.000001, 3, nmfDecomposition, method = "lee")
nmf <- nmf(as.matrix(ctxt_to_extract), 3, nrun = 10, maxIter = 200, method = "lee")
SBS_96 <- rep(c("A.A", "A.C","A.G","A.T","C.A", "C.C","C.G","C.T","G.A", "G.C","G.G","G.T","T.A", "T.C","T.G","T.T"), 6)
SBS_6 <- c(rep("C>A",16), rep("C>G",16),rep("C>T",16),rep("T>A",16),rep("T>C",16),rep("T>G",16))
sigs_df <- sigs@signatures %>% as.data.frame() %>% mutate("SBS" = SBS_6, "ctxt" = SBS_96) %>% pivot_longer(cols = colnames(as.data.frame(sigs@signatures)))
sigs_contr_df <- sigs@samples %>% as.data.frame() %>% mutate("sigs" = paste0("S",1:3)) %>% 
  pivot_longer(cols=colnames(as.data.frame(sigs@samples)),names_to = "sampleID" )
ggplot(sigs_contr_df, aes(x=sampleID, y= value, fill = sigs))+
  geom_col()+my_theme()+scale_fill_my("Sunset")+
  theme(axis.text.x = element_text(angle=90))

ggplot(sigs_df, aes(x=ctxt, y=value*100, fill=SBS))+geom_col()+facet_grid(name~SBS)+my_theme()+scale_fill_my("Sunset")+
  theme(axis.text.x = element_text(angle = 90, size = 5), legend.position = 0)
plotSampleMap(sigs)
plotSamples(sigs)+my_theme()+scale_fill_my("Sunset")+
  geom_col(aes(x=sample, y=value*100, fill = signature))+theme(axis.text.x = element_text(angle = 90))
colnames(sigs@signatures) <- c("B1","DMBA1", "DMBA2", "B2", "B3")

#Fitting them to the mouse samples
DMBA_sig <- read_delim("/icgc/dkfzlsdf/analysis/B210/Javi/DMBA-signature", delim = "\t")
ctxt_to_fit <- CTXT_b %>% subset(grepl("F", treatment) == F) %>%  mutate("sampleID" = paste(species, treatment, sampleID)) %>% dplyr::select(-c(species, treatment)) %>%  pivot_wider(names_from = "sampleID", values_from = "value") %>% dplyr::select(-c(SBS2, ctxt))
DMBA_cosmic <- cbind(as.matrix(DMBA_sig), siglasso::cosmic_v3_exo[,c(7:10,43)]) %>% as.matrix()
siglasso <- siglasso::siglasso(ctxt_to_fit, signature = as.matrix(DMBA_cosmic), gamma = 100)
siglasso_df <- siglasso %>% t() %>%  as.data.frame() %>%
  mutate("sampleID" = rownames(.)) %>%
  pivot_longer(rownames(siglasso), names_to = "sig")
ggplot(siglasso_df, aes(x=sampleID, y=100*value, fill=sig))+geom_col()+my_theme()+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_my("Sunset")

MutationalPatterns::fit_to_signatures(mut_matrix = as.matrix(ctxt_to_fit),
                                      signatures = DMBA_cosmic) %>%
   .$contribution %>% t() %>% as.data.frame() %>% mutate("sampleID" = rownames(.)) %>%  
  pivot_longer(names_to = "sig", cols = colnames(.)[-length(colnames(.))]) %>% 
  mutate("sig" = gsub("a|b|c|d", "", sig))%>%
  ggplot( aes(x=grepl("treated", sampleID), y=100*value, fill=sig))+
  geom_col(width = .8, position = "fill")+my_theme()+
  facet_wrap(~grepl("hgla", sampleID), nrow = 1)+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_my("Sunset")
  

Fitted_sig <- MutationalPatterns::fit_to_signatures_bootstrapped(as.matrix(ctxt_to_fit),(as.matrix(DMBA_sig))) 
Fitted_sig %>% as.data.frame() %>%
  tibble::rownames_to_column("exp") %>%
  tidyr::gather(key = "sig", value = "contri", -exp) %>%
  dplyr::mutate(
    sampleID = gsub("_[^_]+$", "", exp),
    sampleID = factor(sampleID, levels = unique(sampleID)),
    sig = factor(sig, levels = unique(sig))) %>% group_by(sig, sampleID) %>% 
  summarise("mean" = mean(contri), "sem" = sd(contri)/sqrt(n())) %>% subset(sig == "DMBA") %>%  ggplot(aes(x=sampleID, y=mean, fill=sig))+
  geom_col(position="stack")+my_theme()+scale_fill_my("Sunset")+
  geom_errorbar(aes(x=sampleID, ymin=mean-sem, ymax= (mean+sem)),width = 0.3, size = .2)
MutationalPatterns::plot_correlation_bootstrap(Fitted_sig)
MutationalPatterns::plot_contribution_heatmap(MutationalPatterns::fit_to_signatures(as.matrix(ctxt_to_fit),(as.matrix(DMBA_sig))) %>% .$contribution)
Fitted_sig_strict <- MutationalPatterns::fit_to_signatures_strict(as.matrix(ctxt_to_fit),(as.matrix(DMBA_sig))) 


MutationalPatterns::calculate_lesion_segregation(
  lapply(mmus_treated, function(x){AFtoGranges(readAF(x, species = "hgla", treatment = "treated"))}), 
  sample_names = gsub("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF/", "", mmus_treated)
  ) 

MutationalPatterns::calculate_lesion_segregation(
lapply(mmus_treated, function(x){a <- AFtoGranges(readAF(x, species = "hgla", treatment = "treated"))
  GenomeInfoDb::genome(a) = 'mm10'
  return(a)}),
sample_names = gsub("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/17618/Mutation_calling/AF/", "", mmus_treated), split_by_type = T, ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
)




VR[["mmus_MEF_treated"]][[1]] %>% ggplot(aes(x=distance, fill = SBS))+geom_histogram()+facet_wrap(SBS~sampleID)+coord_cartesian(xlim=c(-1000000,1500000))

plotRainfall(VR[["mmus_MEF_treated"]][[1]])


#########

GeneNames <- list()
for (i in c("hgla_treated", "hgla_control", "mmus_treated", "mmus_control")){
  Sliced <- get(i) %>% gsub("AF", "genes", .) %>% gsub(".tsv", "_sliced.bed", .)
  species <- str_match(i, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", i)
  for (j in 1:length(Sliced)){
  GeneNames[[i]][[j]] <- SlicedtoGeneNames(Sliced[j], AF= get(i)[j] ,species = species, treatment = treatment)
}}
GeneNames_b <- bind_rows(GeneNames)
GeneNames_b_Martin <- CompareGeneNamesMartincorena(GeneNames_b, genes = "dndscv")
GeneNames_b_Martin %>% separate(chr_pos, into = c("chr", "pos"), sep = "_") %>% subset(species == "hgla" & treatment == "treated")%>% AFtoVranges()%>% 
  contextmatrix(species = "hgla", treatment = "treated") %>% SignaturePlot()+facet_grid(sampleID~SBS2)

ggplot(GeneNames_b_Martin, aes(x=sampleID, fill=SBS))+geom_bar(position="stack")+facet_wrap(species~treatment, scales = "free")



plotGeneNamesPoints(GeneNames_b_Martin %>% subset(species == "hgla") %>% subset(region == "gene"), labels = T, thr = 0.15
                    )+theme(axis.text = element_blank(), axis.ticks = element_blank())+
  xlab("")+ylab("")+facet_wrap(treatment~sampleID)+
  theme(legend.position = 0)+scale_color_Martin()+scale_fill_Martin()

GeneNames_b_Martin_n_mut <- GeneNames_b_Martin %>% subset(region == "gene") %>% group_by(gene_id, chr_pos) %>% mutate(n=n()) %>% ungroup()
GeneNames_b_n_mut <- GeneNames_b %>% subset(region == "gene") %>% group_by(gene_id, chr_pos) %>% mutate(n=n()) %>% ungroup()
sampleplot <- function(x){plotGeneNamesPoints(GeneNames_b_n_mut %>%  subset(n > 1 & sampleID == x & species == "hgla" & treatment == "treated"), labels = F)+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_blank(), strip.text = element_blank())+
  xlab("")+ylab("")+
  scale_color_my("Martin")+scale_fill_my("Martin")}
B8 <- sampleplot("B8")
B7 <- sampleplot("B7")
cowplot::plot_grid(B7,B8)

plotGeneNamesPoints_grid(GeneNames_b_n_mut %>% subset(SBS == "T>A" & af > 0.3) %>% subset(species == "mmus" & treatment == "treated"), labels = T, thr = 0.3)+
  theme(axis.ticks = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_blank(), legend.position = 0)+facet_wrap(treatment~species)+
  xlab("")+ylab("")+
  scale_color_my("Martin")+scale_fill_my("Martin")



###Clones figure
GeneNames_b %>% subset(af > .2 & SBS == "T>A" & region == "gene" & gene_name != "" & treatment == "treated" & grepl("B7|B8|E2|E7|D2|B12|H12|H5", sampleID)) %>% group_by(gene_name) %>% mutate("n" = length(unique(sampleID))) %>%
ungroup() %>% slice_max(order_by = n, n = 100 ) %>%
plotGeneNamesPoints()+facet_wrap(species+treatment ~ sampleID, nrow = 2)+
scale_fill_my("Sunset")+scale_color_my("Sunset")+
theme(panel.border = element_rect(fill=NA), axis.ticks = element_blank(),
axis.text = element_blank(), legend.title = element_blank(), strip.text = element_blank(),
axis.ticks.x=element_blank(), axis.line = element_blank())+xlab("")+ylab("")

GeneNames_b_Martin %>% subset(af > .2 & region == "gene" & gene_name != "" & treatment == "treated" & grepl("B7|B8|E2|E7|D2|B12|H12|H5", sampleID)) %>% group_by(gene_name) %>% mutate("n" = length(unique(sampleID))) %>%
  ungroup() %>% slice_max(order_by = n, n = 1000 ) %>%
  plotGeneNamesPoints()+facet_wrap(species+treatment ~ sampleID, nrow = 2)+
  scale_fill_my("Sunset")+scale_color_my("Sunset")+
  theme(panel.border = element_rect(fill=NA), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.title = element_blank(), strip.text = element_blank(),
        axis.ticks.x=element_blank(), axis.line = element_blank())+xlab("")+ylab("")




#####
GeneNames_b %>% subset(region == "gene" & sampleID == "B7" | region == "gene" & sampleID == "B8") %>% group_by(gene_id, chr_pos) %>% mutate(n=n()) %>% ungroup() %>%
  subset(n > 1) %>% plotGeneNamesPoints(., labels = T)+facet_wrap(treatment~sampleID)+scale_color_Martin()


genes_23_25 <- merge(GeneNames[["hgla_treated"]][[4]] %>% subset(region == "gene"), GeneNames[["hgla_treated"]][[6]] %>% subset(region == "gene"), by = c("gene_name", "chr_pos", "species", "treatment")) %>% .$gene_name
plotGeneNamesPoints(GeneNames_b %>% subset(sampleID == "AS-475125" | sampleID == "AS-475129") %>% 
subset(gene_name %in% genes_23_25) %>% subset(region == "gene") %>% subset(SBS == "T>A"), labels = T, thr = 0.4)+
  theme(axis.text = element_blank(), axis.ticks = element_blank())+xlab("")+ylab("")

hgla_T_2 <- list()
for(i in c(1,3)){
  hgla_T_2[[i]] <- CommonGeneNames(GeneNames[["hgla_control"]][[2]], GeneNames[["hgla_control"]][[i]]) %>% subset(sampleID != unique(GeneNames[["hgla_control"]][[2]]$sampleID))
}
hgla_T_2_b <- bind_rows(hgla_T_2)
plotGeneNamesPoints(hgla_T_2_b, labels = T)+facet_wrap(~sampleID, ncol = 2)

hgla_1_2 <- CommonGeneNames(GeneNames[["hgla_treated"]][[1]], GeneNames[["hgla_treated"]][[2]])

plotGeneNamesPoints(jj,  labels = T)

######Intervals

intervals <- lapply(c("hgla_F_treated", "hgla_F_control"), function(x){
  species <- str_match(x, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", x)
  intervals <- list()
  for(i in get(x)){
  file <- gsub("AF", "1kb_intervals", i)
  name <- gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(hgla|mmus)/(Treated|Control|17618)/(full_depth/|)Mutation_calling/1kb_intervals/", "", file) %>% gsub(".tsv", "", .)
  intervals[[name]] <- ReadIntervals(file) %>% mutate("sampleID" = name, "species" = species, "treatment" = treatment)
  }
  k <- bind_rows(intervals)
  return(k)
  })
intervals_b <- bind_rows(intervals)
ggplot(intervals_b , aes(x=Mb3, y=number, color = treatment))+geom_jitter()+geom_col(position="dodge")+facet_wrap(species~treatment)+my_theme()+theme(legend.position = 0)+scale_color_my("Sunset") 

intervals_t <- intervals_b %>% unite("sampleID", c("sampleID", "species", "treatment"), sep = "_", remove = T)%>% pivot_wider(names_from = "sampleID", values_from = "number") %>% dplyr::select(-Mb3) %>% as.matrix() %>% t()
intervals_t_hgla <- intervals_t[grepl("hgla", rownames(intervals_t)),] %>% .[,(1:874)]
intervals_t_mmus <- intervals_t[grepl("mmus", rownames(intervals_t)),]

cosmat_intervals_hgla <- as.matrix(cosineDist(as.matrix(intervals_t_hgla)))
MutationalPatterns::plot_cosine_heatmap(cosmat_intervals_hgla)

cosmat_intervals_mmus <- as.matrix(cosineDist(as.matrix(intervals_t_mmus)))
MutationalPatterns::plot_cosine_heatmap(cosmat_intervals_mmus)

matrix_mut_ctxt_hgla <- cbind(intervals_t_hgla, ctxt_expanded_hgla_5)
cosmat_mut_ctxt_hgla <- as.matrix(cosineDist(as.matrix(matrix_mut_ctxt_hgla)))
MutationalPatterns::plot_cosine_heatmap(cosmat_mut_ctxt_hgla)


bin_matrices <- lapply(c("hgla_treated", "hgla_control", "mmus_treated", "mmus_control"), function(x){
  species <- str_match(x, "mmus|hgla")
  treatment <- gsub("(hgla|mmus)_", "", x)
  bin_matrix <- list()
  for(i in get(x)){
    file <- gsub("AF", "1kb_intervals", i)
    name <- gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(hgla|mmus)/(Treated|Control|17618)/(full_depth/|)Mutation_calling/1kb_intervals/", "", file) %>% gsub(".tsv", "", .)
    bin_matrix[[paste0(name, "_", species, "_", treatment)]] <- ReadIntervalsMatrix(file)
    bin_matrix[[paste0(name, "_", species, "_", treatment)]] <- bin_matrix[[paste0(name, "_", species, "_", treatment)]]/rowSums(bin_matrix[[paste0(name, "_", species, "_", treatment)]])
  }
  return(bin_matrix)
})
bin_matrices_hgla <- sapply(c(bin_matrices[[1]],bin_matrices[[2]] ), rbind) %>% t()
bin_matrices_mmus <- sapply(c(bin_matrices[[3]],bin_matrices[[4]] ), rbind) %>% t()

cosmat_1kb_bin_hgla <- as.matrix(cosineDist(as.matrix(bin_matrices_hgla)))
MutationalPatterns::plot_cosine_heatmap(cosmat_1kb_bin_hgla)

cosmat_1kb_bin_mmus <- as.matrix(cosineDist(as.matrix(bin_matrices_mmus)))
MutationalPatterns::plot_cosine_heatmap(cosmat_1kb_bin_mmus)

ctxt_5_1kb_bin_hgla <- cbind(bin_matrices_hgla, ctxt_expanded_hgla_5)
cosmat_ctxt_5_1kb_bin_hgla <- as.matrix(cosineDist(as.matrix(ctxt_5_1kb_bin_hgla)))
MutationalPatterns::plot_cosine_heatmap(cosmat_ctxt_5_1kb_bin_hgla)


ctxt_5_1kb_bin_mmus <- cbind(bin_matrices_mmus, ctxt_expanded_mmus_5)
cosmat_ctxt_5_1kb_bin_mmus <- as.matrix(cosineDist(as.matrix(ctxt_5_1kb_bin_mmus)))
MutationalPatterns::plot_cosine_heatmap(cosmat_ctxt_5_1kb_bin_mmus)

umap_9samples <- umap.defaults
umap_9samples$n_neighbors <- 9
umap_ctxt_5_1kb_bin_hgla <- umap::umap(config = umap_9samples, cosineDist(as.matrix(ctxt_5_1kb_bin_hgla)))
cbind(as.data.frame(umap_ctxt_5_1kb_bin_hgla$layout), samples = row.names(ctxt_5_1kb_bin_hgla)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())


#########
setwd("/icgc/dkfzlsdf/analysis/B210/Javi")
ctxt_df <- lapply(c("ctxt_hgla5.tsv", "ctxt_mmus5.tsv"),  readExpandedCtxt, length = 5, df = TRUE) %>% bind_rows()%>% 
  mutate("ctxt2" =str_split(ctxt, "")) %>% 
  separate(ctxt2, into = c("kk", paste0("L", 1:5),paste0("R", 1:5), "k")) %>% dplyr::select(-c(kk,k)) %>% 
  separate(sampleID, into = c("sample","species", "treatment"))%>% 
  pivot_longer(cols = c(paste0("L", 1:5),paste0("R", 1:5)), names_to = "pos",values_to = "NT")%>%  
  group_by(species, treatment, SBS) %>% mutate(n_mut = n()/10) %>% 
  group_by(species,treatment,SBS,NT,pos) %>% mutate(prop = n()/n_mut)
ctxt_df_short <- ctxt_df %>% dplyr::select(-c(sample, value, ctxt, n_mut)) %>% unique()
ggplot(ctxt_df_short %>% subset(SBS == "T>A")%>% 
         unique(), aes(x=pos, fill=NT))+
  geom_col(aes(y=prop*100), position=position_dodge())+
  geom_text(aes(size = prop, y=(prop*100+2), color=NT, label = NT),
            position = position_dodge(width=1))+
  facet_wrap(species~treatment)+my_theme()+
  coord_cartesian(ylim = c(0,50))+theme(legend.position = 0, axis.text.x = element_blank(), 
                                           )+xlab("")+ylab("")+
  scale_fill_my("Sunset")+scale_color_my("Sunset")




ctxt_immediate_hgla <- readExpandedCtxt("ctxt_hgla.tsv", length = 1) %>% .[(order(row.names(.))),]
ctxt_immediate_mmus <- readExpandedCtxt("ctxt_mmus.tsv") %>% .[(order(row.names(.))),]
ctxt_immediate_MEF <- readExpandedCtxt("ctxt_MEF1.tsv") %>% .[(order(row.names(.))),]
ctxt_immediate <- rbind(ctxt_immediate_hgla, ctxt_immediate_mmus, ctxt_immediate_MEF)


cosine <- ProjectionBasedClustering::tSNE(as.matrix(ctxt_immediate), OutputDimension=2,Algorithm='tsne_cpp',
                                           method="euclidean",Whitening=FALSE, Iterations=1000,PlotIt=F,k = 4)

cbind(as.data.frame(cosine$ModelObject$Y), samples = row.names(ctxt_immediate)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())



ctxt_expanded_hgla <- readExpandedCtxt("ctxt_hgla2.tsv", 2) %>% .[(order(row.names(.))),]
ctxt_expanded_mmus <- readExpandedCtxt("ctxt_mmus2.tsv", 2) %>% .[(order(row.names(.))),]
ctxt_expanded_MEF <- readExpandedCtxt("ctxt_MEF2.tsv", 2) %>% .[(order(row.names(.))),]
ctxt_expanded <- rbind(ctxt_expanded_hgla, ctxt_expanded_mmus, ctxt_expanded_MEF)


cosine2 <- ProjectionBasedClustering::tSNE(as.matrix(ctxt_expanded), OutputDimension=2,Algorithm='tsne_cpp',
                                          method="euclidean",Whitening=FALSE, Iterations=1000,PlotIt=F,k = 4, check_duplicates = F)

cbind(as.data.frame(cosine2$ModelObject$Y), samples = row.names(ctxt_expanded)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())

ctxt_expanded_hgla_3 <- readExpandedCtxt("ctxt_hgla3.tsv", length = 3) %>% .[(order(row.names(.))),]
ctxt_expanded_mmus_3 <- readExpandedCtxt("ctxt_mmus3.tsv", length = 3) %>% .[(order(row.names(.))),]
ctxt_expanded_MEF_3<- readExpandedCtxt("ctxt_MEF3.tsv", length = 3) %>% .[(order(row.names(.))),]
ctxt_expanded_3 <- rbind(ctxt_expanded_hgla_3, ctxt_expanded_mmus_3, ctxt_expanded_MEF_3)    
cosine3 <- ProjectionBasedClustering::tSNE(as.matrix(ctxt_expanded_3), OutputDimension=2,Algorithm='tsne_cpp',
                                           method="euclidean",Whitening=FALSE, Iterations=1000,PlotIt=F,k = 4, check_duplicates = FALSE)

cbind(as.data.frame(cosine3$ModelObject$Y), samples = row.names(ctxt_expanded_3)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())

ctxt_expanded_hgla_5 <- readExpandedCtxt("ctxt_hgla5.tsv", length = 5) %>% .[(order(row.names(.))),]
ctxt_expanded_mmus_5 <- readExpandedCtxt("ctxt_mmus5.tsv", length = 5) %>% .[(order(row.names(.))),]
ctxt_expanded_MEF_5 <- readExpandedCtxt("ctxt_MEF5.tsv", length = 5) %>% .[(order(row.names(.))),]
ctxt_expanded_5 <- rbind(ctxt_expanded_hgla_5, ctxt_expanded_mmus_5, ctxt_expanded_MEF_5)    
cosine5 <- ProjectionBasedClustering::tSNE(as.matrix(ctxt_expanded_5), OutputDimension=2,Algorithm='tsne_cpp',
                                           method="euclidean",Whitening=FALSE, Iterations=1000,PlotIt=F,k = 4, check_duplicates = FALSE)

cbind(as.data.frame(cosine5$ModelObject$Y), samples = row.names(ctxt_expanded_5)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())




cosmat <- as.matrix(cosineDist(as.matrix(ctxt_immediate)))
MutationalPatterns::plot_cosine_heatmap(cosmat)

cosmat2 <- as.matrix(cosineDist(as.matrix(ctxt_expanded)))
MutationalPatterns::plot_cosine_heatmap(cosmat2)


cosmat3 <- as.matrix(cosineDist(as.matrix(ctxt_expanded_3)))
MutationalPatterns::plot_cosine_heatmap(cosmat3)

cosmat5 <- as.matrix(cosineDist(as.matrix(ctxt_expanded_5)))
MutationalPatterns::plot_cosine_heatmap(cosmat5)
MutationalPatterns::cluster_signatures(cosmat5)


umap <- umap::umap(cosineDist(as.matrix(ctxt_immediate)))
cbind(as.data.frame(umap$layout), samples = row.names(ctxt_immediate)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())

umap2 <- umap::umap(config = umap.defaults_euclidean, cosineDist(as.matrix(ctxt_expanded)))
cbind(as.data.frame(umap2$layout), samples = row.names(ctxt_expanded)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())

umap3 <- umap::umap(config = umap.defaults_euclidean, cosineDist(as.matrix(ctxt_expanded_3)))
cbind(as.data.frame(umap3$layout), samples = row.names(ctxt_expanded_3)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())

umap5 <- umap::umap(cosineDist(as.matrix(ctxt_expanded_5)))
cbind(as.data.frame(umap5$layout), samples = row.names(ctxt_expanded_5)) %>% separate(samples, into = c("sample", "species", "treatment"), sep = "_") %>% ggplot(aes(x=V1, y=V2, color = paste(species, treatment)))+geom_point()+geom_text(aes(label=sample))+my_theme()+theme(legend.title = element_blank())



Eu_phy <- EuclideanPhysicaldistance(ctxt_5_1kb_bin_hgla)
ggplot(Eu_phy, aes(x= physical, y=euclidean, color = treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+facet_wrap(~species)+my_theme()+coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))
summary(lm(data = subset(Eu_phy, species=="mmus"), formula= physical~euclidean))$r.squared
Eu_phy%>% group_by(species,treatment) %>% group_modify(~as.data.frame(lm_eqn(.x))) %>% merge(., Eu_phy, by=c("species","treatment")) %>% 
  ggplot(aes(x=physical, y=euclidean, color=treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+ geom_smooth(method = lm, se = FALSE)+geom_text(aes(y = as.numeric(gsub("treated", 0.08 ,treatment) %>% gsub("control", 0.078, .)), x = 5, label = `lm_eqn(.x)`))+facet_wrap(~species)+my_theme()+
  coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))

Eu_phy_2 <- EuclideanPhysicaldistance(ctxt_expanded)
ggplot(Eu_phy_2, aes(x= physical, y=euclidean, color = treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+facet_wrap(~species)+my_theme()+coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))
summary(lm(data = subset(Eu_phy_2, species=="hgla"), formula= physical~euclidean))$r.squared
Eu_phy_2%>% group_by(species,treatment) %>% group_modify(~as.data.frame(lm_eqn(.x))) %>% merge(., Eu_phy_2, by=c("species","treatment")) %>% 
  ggplot(aes(x=physical, y=euclidean, color=treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+ geom_smooth(method = lm, se = FALSE)+geom_text(aes(y = as.numeric(gsub("treated", 0.08 ,treatment) %>% gsub("control", 0.078, .)), x = 5, label = `lm_eqn(.x)`))+facet_wrap(~species)+my_theme()+
  coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))


Eu_phy_3 <- EuclideanPhysicaldistance(ctxt_expanded_3)
ggplot(Eu_phy_3, aes(x= physical, y=euclidean, color = treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+facet_wrap(~species)+my_theme()+coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))
summary(lm(data = subset(Eu_phy_3, species=="hgla"), formula= physical~euclidean))$r.squared
Eu_phy_3 %>% group_by(species,treatment) %>% group_modify(~as.data.frame(lm_eqn(.x))) %>% merge(., Eu_phy_3, by=c("species","treatment")) %>% 
  ggplot(aes(x=physical, y=euclidean, color=treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+ geom_smooth(method = lm, se = FALSE)+geom_text(aes(y = as.numeric(gsub("treated", 0.08 ,treatment) %>% gsub("control", 0.078, .)), x = 5, label = `lm_eqn(.x)`))+facet_wrap(~species)+my_theme()+
  coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))



Eu_phy_5 <- EuclideanPhysicaldistance(ctxt_expanded_5)
Eu_phy_5 %>% group_by(species,treatment) %>% group_modify(~as.data.frame(lm_eqn(.x))) %>% merge(., Eu_phy_5, by=c("species","treatment")) %>% 
ggplot(aes(x=physical, y=euclidean, color=treatment))+geom_jitter()+geom_text(aes(label=paste(sample1, sample2)))+ geom_smooth(method = lm, se = FALSE)+geom_text(aes(y = as.numeric(gsub("treated", 0.08 ,treatment) %>% gsub("control", 0.078, .)), x = 5, label = `lm_eqn(.x)`))+facet_wrap(~species)+my_theme()+
  coord_cartesian(ylim = c(0.04,0.08), xlim = c(0,12))
summary(lm(data = subset(Eu_phy_5, species=="hgla"), formula= physical~euclidean))$r.squared

#dndsCV

dnds_hgla <- dndscv(mutations = bind_rows(mutations[["hgla_treated"]], mutations[["hgla_control"]]), refdb = "hgla_refcds.rda", cv = NULL)
ggplot(dnds$annotmuts, aes(x=sampleID, fill=impact))+geom_bar()+scale_fill_my("Sunset")+facet_wrap(~grepl("control", sampleID), scales = "free_x")
dnds_mmus <- dndscv(mutations = bind_rows(mutations[["mmus_treated"]], mutations[["mmus_control"]]), refdb = "mmus_refcds.rda", cv = NULL)
ggplot(dnds_mmus$annotmuts %>% subset(impact != "Synonymous"), aes(x=sampleID, fill=impact))+geom_bar()+scale_fill_my("Sunset")+facet_wrap(~grepl("control", sampleID), scales = "free_x")


RefCDS_mmus <- load("mmus_refcds.rda")
RefCDS_mmus <- buildcodon(RefCDS_mmus)
dnds_codon <- codondnds(dnds_mmus, RefCDS_mmus)
RefCDS_hgla <- load("hgla_refcds.rda")
RefCDS_hgla <- buildcodon(RefCDS_hgla)


annot_mut_af_hgla <- inner_join(dnds_hgla$annotmuts %>% dplyr::rename("gene_name" = gene) %>% 
                                  separate("sampleID", into = c("sampleID", "treatment"), sep = "_"), AF_b %>% subset(species == "hgla"), by = c("chr","pos", "sampleID", "treatment")) 
annot_mut_af_mmus <- inner_join(dnds_mmus$annotmuts %>% dplyr::rename("gene_name" = gene) %>% 
                                  separate("sampleID", into = c("sampleID", "treatment"), sep = "_"), AF_b %>% subset(species == "mmus"), by = c("chr","pos", "sampleID", "treatment")) 

plotGeneNamesPoints_grid(subset(annot_mut_af_hgla, impact != "Synonymous"))+theme(legend.position = 0)+facet_wrap(~treatment)
plotGeneNamesPoints_grid(subset(annot_mut_af_mmus, impact != "Synonymous"))+theme(legend.position = 0)+facet_wrap(~treatment)

#vep 

vep <- list()
for (i in c("hgla_control_vcf", "hgla_treated_vcf", "mmus_control_vcf", "mmus_treated_vcf")){
  Sliced <- get(i) %>% gsub("calling/", "calling/vep/", .) %>% gsub("/pass.vcf", ".vcf", .)
  species <- str_match(i, "mmus|hgla")
  treatment <-  str_match(i, "treated|control")
  for (j in 1:length(Sliced)){
    sample <- Sliced[j] %>%  gsub("/icgc/dkfzlsdf/analysis/B210/Javi/(hgla|mmus)/(Treated|Control|17618)/(full_depth/|)Mutation_calling/vep/", "", .) %>% gsub(".vcf", "", .) %>% translate_names()
    vep[[i]][[j]] <- read.table(Sliced[j], sep = "\t") %>% unique() %>% mutate("species" = species, "treatment" = treatment, "sampleID" = sample) %>% 
      separate(V1, c("chr", "pos", "ref"), sep = "_") %>% separate(ref, c("ref", "alt"), sep = "/") %>% dplyr::select(-c(V2,V3,V6,paste0("V",8:13))) %>%
      unite("chr_pos", c(chr,pos), sep = "_") %>% separate(V14, into = c("impact"), sep = ";") %>%group_by(chr_pos, V4, V5) %>% 
      mutate(V7 = paste(V7, collapse = ",") %>% gsub("_variant","",.), "impact"=paste(impact, collapse = ",") ) %>%
      separate(V7, into =c("consequence", NA), sep = ",", remove = F) %>% mutate("consequence" = gsub("(downstream|upstream)_gene|intergenic|non_coding_transcript_exon", "noncoding", consequence) %>% gsub("coding_sequence|incomplete_terminal_codon|stop_lost|protein_altering|start_lost|inframe_insertion|inframe_deletion", NA, .) %>% gsub("splice_(donor|region|acceptor)|intron", "splice_site", .)) %>%
      dplyr::rename("gene_id"=V4, "transcript_id"=V5) %>% dplyr::select(-V7)%>% unique() 
  }}
vep_b <- bind_rows(vep)

vep_GeneNames <- inner_join(GeneNames_b, vep_b, by = c("sampleID", "species", "chr_pos", "gene_id", "ref", "alt", "treatment"))
vep_GeneNames_Martin <- inner_join(GeneNames_b_Martin, vep_b, by = c("sampleID", "species", "chr_pos", "gene_id", "ref", "alt", "treatment"))


vep_GeneNames_Martin %>% subset(grepl("missense|frameshift|splice|stop", consequence)) %>% 
  subset(af > .2 & region == "gene" & gene_name != "" & treatment == "treated" & grepl("B7|B8|E2|E7|D2|B12|H12|H5", sampleID)) %>% group_by(gene_name) %>% mutate("n" = length(unique(sampleID))) %>%
  ungroup() %>% slice_max(order_by = n, n = 50 ) %>%
  plotGeneNamesPoints()+facet_wrap(species+treatment ~ sampleID, nrow = 2)+
  scale_fill_my("Sunset")+scale_color_my("Sunset")+
  theme(panel.border = element_rect(fill=NA), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.title = element_blank(), strip.text = element_blank(),
        axis.ticks.x=element_blank(), axis.line = element_blank())+xlab("")+ylab("")

vep_GeneNames_Martin %>% subset(grepl("missense|frameshift|splice|stop", consequence)) %>%
  subset( region == "gene" & gene_name != "" & treatment == "treated") %>% group_by(gene_name) %>% mutate("n" = length(unique(sampleID))) %>%
  ungroup()%>%
  ggplot(aes(x=gene_name, y=paste(species, sampleID), fill = species))+geom_tile()+my_theme()+
  scale_fill_my("Sunset2")+theme(legend.position = 0, axis.text.x = element_text(angle = 90, size = 8),
                                 panel.grid.minor = element_line(color = "black"))
                                
                                 
#########Spectrums in just reference

mouse_treated <- lapply(
  c("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-2-1a/pass2.vcf",
    "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-1-4a/pass2.vcf",
    "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-3-4a/pass2.vcf"),
  function(x) {readAF(x, species = "mmus", treatment = "treated") %>% subset(nchar(ref)==1 & nchar(alt) ==1) %>%
    AFtoVranges()%>% contextmatrix(species = "mmus", treatment = "treated")})%>% bind_rows()

mouse_control <- lapply(
  c("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-2-2o/pass2.vcf",
    "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-3-5o/pass2.vcf",
    "/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/just_reference/MEF-1-3o/pass2.vcf"),
  function(x){readAF(x, species = "mmus", treatment = "control") %>% subset(nchar(ref)==1 & nchar(alt) ==1) %>%
    AFtoVranges()%>% contextmatrix(species = "mmus", treatment = "control")})%>% bind_rows()


#######Fibroblastsss

ispr07 <- lapply(
  list.files("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/Mutation_calling", full.names = T, recursive = T, pattern = "glaber_scaffolds_aab_Mutect2_output_filtered.vcf$"),
  function(x) readPassVcf(x)%>%mutate("sampleID" = gsub("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/ispr07/Mutation_calling/", "", x)%>%gsub("/Mutect2/glaber_scaffolds_aab_Mutect2_output_filtered.vcf", "", .)
                                        ))%>%bind_rows()


ggplot(ispr07, aes(x=as.numeric(pmin(tumorAF, 1)), fill=sampleID))+geom_density(position="stack")

ggplot(RefAltToSBS6(ispr07%>% subset(nchar(ref)==1&nchar(alt)==1)), aes(x=as.numeric(pmin(tumorAF, 1)), fill=SBS))+
  geom_histogram(binwidth = .1, position = "fill")+
  facet_wrap(~grepl("a", sampleID), ncol=1)+my_theme()+scale_fill_my("Sunset")

MEF <- lapply(
  list.files("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/PON", full.names = T, recursive = T, pattern = "pass2.vcf"),
  function(x) readPassVcf(x)%>%mutate("sampleID" = gsub("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/MEF-all/PON/", "", x)%>%gsub("pass2.vcf", "", .)
  ))%>%bind_rows()

ggplot(RefAltToSBS6(MEF), aes(x=as.numeric(pmin(tumorAF, 1)), fill=SBS))+
  geom_histogram(binwidth = .1, position = "fill")+geom_text()+
  facet_wrap(~grepl("a", sampleID), ncol=1)+my_theme()+scale_fill_my("Sunset")
