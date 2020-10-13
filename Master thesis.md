**UNIVERSITÄTSMEDIZIN MANNHEIM MEDICAL FACULTY MANNHEIM OF THE
UNIVERSITY OF HEIDELBERG**

**Master of Science in Translational Medical Research**

**Mutational Landscape of a chemical skin carcinogenesis model of Mouse
and Naked Mole-rat**

Author: Francisco Javier Botey Bataller

From: Valencia, Spain

Thesis supervisor and tutor: Angela Gonçalves Filimon

Second supervisor: Luca Penso Dolfin

Second tutor: Michaela Frye

Thesis submitted October 2020

**Table of contents**
=====================

[Table of contents 2](#table-of-contents)

[Abstract 3](#_Toc53043204)

[1 Introduction 4](#_Toc49248808)

[1.1 Cancer, a genetic disease 4](#cancer-a-genetic-disease)

[1.2 Analyzing the mutational landscape
5](#analysing-the-mutational-landscape)

[1.2.1 Mutational processes 5](#mutational-processes)

[1.2.2 Somatic mutation calling 5](#somatic-mutation-calling)

[1.2.4 Mutational signatures 6](#mutational-signatures)

[1.3 Two-Stage carcinogenic model of skin cancer
9](#two-stage-carcinogenic-model-of-skin-cancer)

[1.4 The Naked mole-rat: a model for cancer resistance
9](#the-naked-mole-rat-a-model-for-cancer-resistance)

[2 Aims and Objectives 11](#aims-and-objectives)

[3. Methods and Materials 12](#methods-and-materials)

[3.1 DMBA/TPA treatment 12](#dmbatpa-treatment)

[3.1.1 In vitro: on cultured fibroblasts
12](#in-vitro-on-cultured-fibroblasts)

[3.1.2 In vivo: on mouse and Naked mole-rat
12](#in-vivo-on-mouse-and-naked-mole-rat)

[3.2 Whole-exome and Whole-genome sequencing
13](#whole-exome-and-whole-genome-sequencing)

[3.3 Somatic Variant calling 13](#somatic-variant-calling)

[3.4 Somatic variant analysis and signature extraction
13](#somatic-variant-analysis-and-signature-extraction)

[4 Results 14](#results)

[4.1 Mutational burden 14](#mutational-burden)

[4.2 Mutational spectrum 15](#mutational-spectrum)

[4.3 Mutational Signature in mouse embryonic fibroblasts
20](#mutational-signature-in-mouse-embryonic-fibroblasts)

[4.4 Fitting the signature found in MEF to the skin samples
22](#fitting-the-signature-found-in-mef-to-the-skin-samples)

[4.5 Mutational signature in variant allelic frequency clusters
23](#mutational-signature-in-variant-allelic-frequency-clusters)

[4.6 Nonsynonymous mutations in cancer driver genes
24](#nonsynonymous-mutations-in-cancer-driver-genes)

[5. Discussion 26](#discussion)

[Bibliography 29](#_Toc53043229)

<span id="_Toc49248808" class="anchor"></span>**1 Introduction**

**1.1 Cancer, a genetic disease**
---------------------------------

Cancer is one of the leading causes of death worldwide. Data from the
Global Cancer Observatory indicates that, in 2018, 18 million people
were diagnosed with cancer and it caused 9.6 million deaths. Among the
different aetiologies of the disease, lung cancer represents the one
with the highest incidence and mortality and non-melanoma skin cancer
ranks the fifth in incidence, representing 5,8% of the cases diagnosed
and the twenty-third in mortality. Included in the latter group are
squamous cell carcinomas (SCC), and cutaneous squamous cell carcinoma
(cSCC) the central aetiology of this thesis (Bray et al., 2018).

Abiding by the hallmarks of cancer, all cancer types are defined by six
principles: sustained proliferative signalling, cell death resistance,
evasion of growth suppressors, induction of angiogenesis, enabled
replicative immortality and metastasis. This is the case of SCC, in
which processes such as cell death or cell cycle are affected, causing
an aberrant cell proliferation. Triggering these disruptions is genomic
instability. Copy number alterations and mutations located in genes
coding for proteins. Examples of the latter in SCC are the *Notch*
signalling pathway, the *Ras* family of gtpases or the tumor suppressor
gene *TP53* (Gao et al., 2014; Hanahan & Weinberg, 2011)*.*

Cells can acquire one of the hallmarks, cell proliferation, by acquiring
somatic mutations. Cell proliferation is mainly regulated by the binding
of ligands on cell-surface receptors, which activate downstream pathways
and modulate cell metabolism. Cancerous cells can sustain proliferation
and evade growth suppression affecting these mechanisms. They can do so
by increasing the concentration of ligands, synthesizing more receptors,
constitutively activating the downstream pathways or disabling
suppression pathways. Which can be consequences of changes at the DNA
level. Cancerous cells acquire somatic mutations that alter normal
proliferation pathways. This can give them a selective advantage over
non-cancerous cells, allowing them to clonally expand (Hanahan &
Weinberg, 2011).

Somatic mutations are not only present in disease, they are also common
in tissues exposed to potential mutagens, such as skin. UV light and
other mutational processes can alter the DNA sequence of exposed cells,
causing the formation of clones that contain mutations in cancer driver
genes (Martincorena et al., 2015). Examples of genes reported as cancer
drivers are: *NOTCH1*, *TP53, HRAS* or *FAT1,* each one having a role in
both homeostasis and cancer*. NOTCH1* is part of the Notch signaling
pathway. It is regulated by ligand binding and has a role in
development, modulating cell growth and survival. Genes related to this
pathway have been found altered in different cancer aetiologies such as
SCC or leukemia (Wang et al., 2011). TP53 is a tumor suppressor gene,
its mediation in many different cellular mechanisms is usually disrupted
in cancer(Pitolli et al., 2019). Adding to this list is *HRAS*, a gene
coding a protein in the Ras family of small GTPases, mutated in some
cancer etiologies and linked to tumor progression in tumor models of SCC
(Lowry et al., 2016). *FAT1,* the last in the list, is a gene expressed
in mammals epithelial cells, which is involved in the downregulation of
𝛃-Catenin, a protein linked to Wnt signaling , and mutated in SCC
(Morris et al., 2013). Mutations in the coding sequences of these genes
are commonly found in cancers and have been found to alter
proliferation. It is also known that they are also mutated in normal
tissues and therefore thought that it is in the accumulation of several
of these mutations when malignancies arise. Studying, therefore, how
this mutations are produced in the first place and why some tissues
progress into tumors and others do not although having somatic mutations
in cancer driver genes can be determinant for achieving an earlier
detection of cancer and a better treatment.

**1.2 Analysing the mutational landscape**
------------------------------------------

###  **1.2.1 Mutational processes**

The Darwinian model for cancer evolution assumes somatic mutations
happen randomly in the genome. They hit genes and alter cell function.
If this alteration grants the cell a selective advantage, such as a
mutation in a driver gene, it will clonally expand, forming a clone of
cells with aberrant proliferation, a tissue that is not controlled by
tissue homeostasis and can progress into a malignancy.

Somatic mutations do happen spontaneously, as was mentioned before. And
understanding the mutational processes that cause them and how is key to
understanding malign neoplasia (Martincorena et al., 2015). Mutational
processes known are diverse. In the case of SCC, somatic mutations can
be a result of UV light exposure, (Brandt & Moore, 2019) in lung cancer,
tobacco compounds can be causing the alterations, as can DNA repair
defects be the cause in others. In order to identify the causing agents
of the somatic mutations it is vital to understand the process by which
each one is causing the mutations.

UV-radiation produces dimeric pyrimidine lesions that alter DNA
metabolism and can be repaired by the nucleotide excision repair (NER)
system, either coupled to transcription or not. If, instead, lesions are
kept they can lead to the formation of a malignancy. It has been
described that they mainly give rise to Cytosine (C) to Thymine (T) (or
Guanine (G) to Adenine (A)) substitutions and can therefore be
identified in skin (Choi et al., 2006; Jans et al., 2005). The
mutational effect of tobacco is not as specific as the one caused by
UV-radiation. Tobacco contains 60 different carcinogens that act in
different ways, overall causing DNA adducts that, if are not repaired,
produce mainly G-T and G-A substitutions (Hecht, 2003). Another known
mutational process is defective mismatch repair (MMR) mechanism.
Mutations altering the structure of proteins acting in this process are
found in some cancers and the characteristics of these mutations are
well-described in the literature (Drost et al., 2017).

The use of next generation sequencing makes it possible to identify
accurately somatic mutations. But characterizing processes that caused
these mutations is not trivial. Many processes can be involved and
acting at the same time and only few of them are characterized. In order
to characterize the mutational processes, the mechanism by which they
act is summarized in a mutational signature.

###  **1.2.2 Somatic mutation calling**

Next generation sequencing (NGS) techniques have opened new avenues for
the study of germline and somatic mutations, with new techniques and
algorithms improving performance in the detection of low frequency
somatic mutations. Deep sequencing technologies result in a higher
number of reads per base (a higher sequencing coverage) and are
especially relevant for variant calling. Sequencing reads from a deep
sequencing mapped to the reference genome are the input used in this
type of calling. Two samples are needed: the sample in which somatic
mutations will be called , in the literature and from now on in this
thesis the “tumor” sample, and a sample from which germline variants
will be extracted, the “normal” sample. From this point, different
variant calling algorithms emerge to call mutations on the tumor sample
(Xu, 2018).

Relevant for this thesis are *Mutect2* and *Strelka because* they were
used throughout this work. Some algorithms are better suited for an
analysis where different clones may be mixed in the same samples, as is
the case of tumor samples. The two algorithms mentioned are suited as
they do not assume diploidy in the sample, considering there may be
subclones present in the same sample. Other than that, *Strelka*
considers joint allelic frequencies detected in tumor and normal
samples, and computes the probability of the allelic frequencies found
with the number of reads, to then compute the probability of the
position being a somatic variant (Saunders et al., 2012). *Mutect2* uses
a haplotype-based approach, also used in *Haplotypecaller* (Poplin et
al., 2017) for germline variants, in which “active regions” are
identified in the genome, where reads are assembled and then variants
called. This is performed by position. Log odds ratio (LOD) are computed
for the tumor read positions, if the LOD is above a pre-set threshold
these positions are considered potential variants and the LOD for the
normal base in the normal sample is then taken into account to decide if
it is a somatic variant. After this, filtering is performed on the
output (Xu, 2018). Once somatic mutations are called, a lot of
information can be extracted from them about the cells harbouring them:
how many of the sequenced cells have that mutation, when did that
mutation occur, how “mutated” are the sequenced cells or which process
was responsible for these mutations.

###  **1.2.4 Mutational signatures**

Mutational processes do not act evenly across the genome, they tend to
prefer certain substitutions and patterns, leaving a fingerprint, a
mutational signature. If the signature for a mutational process is
characterized, it can then be inferred from the mutational spectrum of a
sample. The most common way of characterizing a mutational signature is
by dissecting it into a vector containing the six different
substitutions possible (grouping them with their complementary
substitution), and the context of these mutations, the combination of
possible bases in the 5’ and 3’ base to the substitution, adding to 96
components. Using this technique, it has been possible to gather a
repertoire of signatures acting in human cancers and the processes
responsible for them (Alexandrov et al., 2013; Alexandrov et al., 2020).
Moreover, some mutational processes also have a bias for the strand in
which they cause the mutations. This, the strand bias, is another
feature that can help identify the causing agent that is a consequence
of transcription and replication-coupled repair mechanisms (Robinson et
al., 2020). The expanded context in which the mutational process acts
can also be specific, meaning, not only the first 5’ and 3’ base to the
mutations, but up to the 10th base. Being this a consequence of the fact
that some mutagens, like colibactin, a toxic metabolite of human
microbiome, damage the DNA by causing adenine adducts that are not
formed in adjacent bases, but rather in bases that are up to 3
nucleotides away (Pleguezuelos-Manzano et al., 2020). Another feature of
a mutational process is the location of the mutation in the genome.
Genomic structures are not even across the genome, chromatin,
heterochromatin and chromosomal structures are not repaired with the
same frequency and machinery, which results in different mutation rates
and signatures along the genome, also specific for the mutational
process (Zhang et al., 2020).

Mutational processes are defined by the mutational signature they print
on the genome. These are numerically analysed as the probability of
leaving each of the single base substitutions (SBS) in a specific
mutational context. Multiplying these values by the exposure of the
process, the number of mutations it causes on the genome, gives us the
resulting mutational spectrum that is found in the genome.
Mathematically, the mutational spectrum found in *g* samples can be
represented as the matrix *M,* a matrix of order 96 x g, where values
represent number of mutations, rows are the 96 SBS and columns each of
the samples. Each column in *M* is a vector *m* with 96 components
(*m*<sub>*g*</sub><sup>1</sup>, *m*<sub>*g*</sub><sup>2</sup>,
…<sub>,</sub> *m*<sub>*g*</sub><sup>96</sup>). The mutational spectrum
is the result of different mutational processes acting with different
intensities (exposures). Therefore, *m* is the result of multiplying the
vectors *s* and *e. s* contains the signatures of a mutational process
*n* , representing the probability of each SBS
(*s*<sub>*g*</sub><sup>1</sup>, *s*<sub>*g*</sub><sup>2</sup>, …,
*s*<sub>*g*</sub><sup>96</sup>). And *e* the exposures of the process
*n*, the number of mutations it will cause on the genome g,
*e*<sub>*g*</sub><sup>*n*</sup> (Alexandrov et al., 2012) . Vectors from
different mutated samples can be grouped in a matrix *M* of order 96 x
g.

$$M = \\begin{pmatrix}
m\_{1}^{1} & m\_{2}^{1} & \\cdots & m\_{g}^{1} \\\\
m\_{1}^{2} & m\_{2}^{2} & \\cdots & m\_{g}^{2} \\\\
 \\vdots & \\vdots & \\  & \\vdots \\\\
m\_{1}^{96} & m\_{2}^{96} & \\cdots & m\_{g}^{96} \\\\
\\end{pmatrix}$$

All acting mutational processes’ signatures in a matrix S of order 96 x
n.

$S = \\begin{pmatrix}
s\_{1}^{1} & s\_{2}^{1} & \\cdots & s\_{n}^{1} \\\\
s\_{1}^{2} & s\_{2}^{2} & \\cdots & s\_{n}^{2} \\\\
 \\vdots & \\vdots & \\  & \\vdots \\\\
s\_{1}^{96} & s\_{2}^{96} & \\cdots & s\_{n}^{96} \\\\
\\end{pmatrix}$

And the exposures of the acting mutational processes in a matrix *E* of
order n x g.

$$E = \\begin{pmatrix}
e\_{1}^{1} & e\_{2}^{1} & \\cdots & e\_{g}^{1} \\\\
e\_{1}^{2} & e\_{2}^{2} & \\cdots & e\_{g}^{2} \\\\
 \\vdots & \\vdots & \\  & \\vdots \\\\
e\_{1}^{n} & e\_{2}^{n} & \\cdots & e\_{g}^{n} \\\\
\\end{pmatrix}$$

Considering that multiplying the signature vector
*p*<sub>*k*</sub><sup>96</sup> and the exposure
*e*<sub>*n*</sub><sup>*k*</sup>, the mutational spectrum
*m*<sub>*n*</sub><sup>96</sup> of a sample is the result, the same can
be applied to the matrices, resulting in Equation 1.

  *M* ≈ *S* × *E*     (*E**q*. 1)

Using this framework, two main analyses can be performed, I will refer
to them as signature extraction and signature fitting. Signature
extraction is estimating which signatures (and therefore which
mutational processes) are present in a set of mutated genomes. It is
thought as a problem where two latent variables (matrices *S* and *E*)
need to be estimated .The most common approach for this is the use of
the Poisson nonnegative matrix factorization (NMF), to extract matrices
*S* and *E* from the mutational spectrum of matrix *M.* NMF minimizes
the distance between *M* and *S x E.* And the probability of a mutation
is considered an independent Bernoulli trial whose mean is equal to the
variance, that is, a Poisson distributed probability. Many are the
variations made to this approach, but here I will center my efforts in
the most widely used and “gold standard” by Alexandrov et., al 2012.
This approach is based on the implementation by Brunet et al., 2004 of
the multiplicative update algorithm (D. D. Lee & Seung, 1999). First,
the number of components in M is reduced by removing mutations with a
proportion below 1%. Then, the reduced sample is resampled with Monte
Carlo bootstrapping and the multiplicative update algorithm is applied
to the resample, aiming to reduce the distance between *M* and *S x E*.
These two steps are iterated, and an average *S* matrix is obtained
(Alexandrov et al., 2013). Other approaches differ with this one in the
NMF step, where different models are applied (Baez-Ortega & Gori, 2019;
Lyu et al., 2020).

With this approach, a repertoire of somatic mutational signatures in
cancer has been generated as part of the ‘Catalogue Of Somatic Mutations
In Cancer’ (COSMIC), with the goal of understanding the causes of the
mutations in the different cancer aetiologies and being able to estimate
the contribution in each of them in new samples and therefore estimate
the causing agents (Alexandrov et al., 2020). Having a repertoire of
signatures extracted from different aetiologies, the causes of the
somatic mutations in a given genome can be inferred by estimating which
signatures are present and to which extent are they acting on the
genome. This is signature fitting. It is mainly performed by estimating
non-negative linear combinations of the known mutational signatures that
explain the observed catalogue. Approaches are based on multiple lineal
regression algorithms minimizing the distance between the linear
combination and the observed catalogue (Omichessan et al., 2019).
Considering the framework in Equation 1, signature fitting is estimating
for a given *M* the parameter *E* knowing the signatures present in *S*.
Being *S* a set o mutational signatures that can be found in the genome,
i.e signatures from the COSMIC catalogue. An algorithm developed for
this purpose is sigLASSO. sigLASSO estimates *E* minimizing the distance
between *M* and *S x E.* It does so assuming that the number of
signatures present in a sample needs to be small, taking into account
the mutational burden of each sample and allowing the use of priors.
Fitting fewer signatures is achieved by applying LASSO L1
regularization, which penalizes the weights of the regression, in this
case *e*, applying a factor λ that leads to some of the weights being 0.
The use of priors also penalizes some of the signatures, those that are
not expected in the sample, such as UV light exposure signature in
tissues that are internal. And it raises others, making it possible, for
instance, to add information about the patient, such as smokers, in
which the signature of smoking is more likely to be found (Li et al.,
2020).

The ability to infer the causing agents from the mutational signatures
found in a sample and the common framework that has been established the
past years in the literature, making it possible to build many different
computational methods (Omichessan et al., 2019) to extract the most
information out of them makes their study a valuable approach to
understand the nature and origin of the somatic mutations found in
healthy tissues and cancers.

**1.3 Two-Stage carcinogenic model of skin cancer**
---------------------------------------------------

In vivo animal models are widely used for their suitability in cancer
research. In SCC, genetic models and chemical models are the more
common. The first ones are based on disrupting genes that will cause a
cancer, and the latter on applying mutagens that will give raise to
lesions that can proliferate into in situ and invasive carcinomas. If
the goal is studying the early stages of SCC, how primary lesions are
produced and how they progress into malignancies, chemical models are
better suited (Amberg et al., 2015). Other mutagens are used to study
other cancer aetiologies, such as DEN in liver cancer, a model that
allows the study of the clonal expansion of mutations caused by the
carcinogen (Aitken et al., 2020). Likewise, a two-stage model is used in
SCC, the DMBA-TPA model for skin carcinogenesis (Nassar et al., 2015).

One of the more widely used models in the study of skin cancer is the
treatment with 9,10-dimethy-1,2-benzanthracene (DMBA) and
tetradecanoyl-phorbol acetate TPA. DMBA, when applied to the skin,
causes primary lesions in the form of a benign papilloma, which, after
TPA treatment can advance into SCC. DMBA is a compound that can be
formed in the non-complete combustion of organic compounds such as coil
(Rajapaksa, 2010). When it is metabolized by the cell, by the enzyme
microsomal epoxide hydrolase (mEH), secondary metabolites are formed,
among them DMBA-3,4-dio which is thought to be the major responsible for
DMBA genotoxicity (Miyata et al., 2001). This results in the formation
of both deoxyadenosine and deoxyguanosine adducts, which mainly end up
causing T to A (T&gt;A) transversions and G to T transversions (Ross &
Nesnow, 1999). The study of the mutational signatures found in mouse
models treated with DMBA/TPA is already described. Its main substitution
and context is the T&gt;A transversion in a CG context. Which is thought
to be a consequence of mutagens acting on this context and base excision
repair (BER) being less efficient in CAG (and CTG) trinucleotides (Cai
et al., 2010; Nassar et al., 2015).

Overall, the use of this model is a valid approach to investigate the
initiation and progression of SCC, being the mutational signature of
DMBA well described and having been found genes mutated in human SCC
(Nassar et al., 2015).

**1.4 The Naked mole-rat: a model for cancer resistance**
---------------------------------------------------------

The naked mole-rat (NMR) is an African subterranean rodent that is
surprisingly long-lived. While other rodents of similar size tend to
live around 5 years, this species lives on average over 30 years and is
cancer resistant. The hallmarks of these unique properties are yet
unknown. Discovering them would be key for the development of early
treatments for cancer and other diseases associated with aging.

The theories behind the cancer resistance are many. An enhanced
mechanism of contact inhibition (Seluanov et al., 2009), a difference in
hyaluronan, an extracellular matrix component (Tian et al., 2013), a
reduction of the Ras signaling pathway (Zhao et al., 2020), or a change
in the immune system (Hilton et al., 2019) among others. Contact
inhibition is the process by which cells stop proliferating when there’s
cell to cell contact. This mechanism is altered in cancer cells, as they
will keep proliferating regardless of the contact to other cells.
Gorbunova et al hypothesized in 2009 that NMR cells had a mechanism of
“Early contact Inhibition” (ECI), by which NMR fibroblasts had an
enhanced contact inhibition in vitro compared to other species, which
would confer them a resistance to abnormal proliferation, accumulation
of malignant clones via clonal expansion and, therefore, cancer
(Seluanov et al., 2009). They also hypothesized in 2013 that ECI could
be mediated by their extracellular matrix (ECM). NMR’s ECM contains a
larger hyaluronan than other species and a lower expression of
hyaluronidases, enzymes that digest hyaluronan, as well as differences
in the hyaluronan synthesis machinery. This causes a more viscous ECM
when culturing NMR fibroblasts and could mediate the contact between
cells (Tian et al., 2013). The same group (Gorvunova et al) reported in
2020 that the basal levels of expression of components in the Ras
signalling pathway is lower in NMR compared to other species. And that
the overexpression of these factors can counteract the cancer resistant
properties of the NMR and transform their cells into cancer cells. This
could, therefore, be a prove of a cancer resistance mechanism of the NMR
(Zhao et al., 2020). Another recent study points at differences in cell
type composition of the immune system, although they do not link it with
cancer (Hilton et al., 2019).

Other reports state that the NMR longevity and cancer resistance is
sustained at a molecular level. Genome stability, DNA repair levels,
which are comparable in terms of transcription with those of human
beings (MacRae et al., 2015), or a well-regulated RNA splicing (B. P.
Lee et al., 2020) are a few of the advantages that are thought help the
NMR live longer (Petruseva et al., 2017). These molecular differences
would result in an animal with a better protected DNA, less prone to
malignant genomic changes and, maybe, shielded against external
mutagens. Models reported in the used of carcinogenic models in the NMR
are scarce, almost non-existent.

**2 Aims and Objectives**
=========================

The aim of this project is examining the genomic differences that arise
from using a chemical carcinogen on the skin of both mouse and naked
mole-rat. With this, we aim to answer the following questions:

-   How many mutations are found in naked mole-rat and mouse treated
    skin?

-   Is the characteristic DMBA mutational signature found in naked
    mole-rat and mouse treated skin?

-   Is DMBA affecting mouse embryonic fibroblasts in a similar manner it
    is affecting naked mole-rat and mouse skin?

-   Are common cancer driver genes affected in naked mole-rat and mouse
    treated skin?

**3. Methods and Materials**
============================

**3.1 DMBA/TPA treatment** 
--------------------------

### **3.1.1 In vitro: on cultured fibroblasts**

For the treatment of fibroblasts, C57BL/6 mouse embryonic fibroblasts
(MEFs) (ATCC, SCRC-1008) were used. MEFs were cultured at 37°C, 5%
CO<sub>2</sub>, in Dulbecco´s Modificed Eagles´s Medium (DMEM) (Gibco,
Life Technologies, 41966-029) supplemented with 15% fetal bovine serum
(Biomol, FBS-01-0500) and 1% PEST (Sigma-Aldrich, P0781).

Sub-confluent (50%) fibroblasts were treated with 50 μM DMBA
(Sigma-Aldrich, D3254) or 0.05% DMSO for 3 hours and single sorted using
Aria I flow sorter, 100 μm nozzle, 20 psi. Single-sorted cells were
manually inspected 3-24 hours after plating to exclude potential well
with more than one cell. After 9 days the cells were all transferred to
a 6-well well/clone) and for a total of 14 days (MEFs).

Genomic DNA was extracted from 200,000 fibroblasts from each expanded
clone using the AllPrep micro DNA/RNA kit (Qiagen, 80284) according to
the manufacturer’s instructions. DNA quantity was determined using the
Qubit dsDNA HS Assay Kit (Life Technologies) and gDNA integrity was
assessed using a Genomic DNA ScreenTape (Agilent Technologies) (DIN
average 7.5 +/- 0.24).

### **3.1.2 In vivo: on mouse and Naked mole-rat** 

For the skin carcinogenic model, 7-week-old female C57Bl6/J Janvier mice
and 6-month-old NMR were used. A single dose of DMBA (400 nmol/100 μl
acetone; Sigma-Aldrich) was applied to the shaved dorsal skin of the
animal, followed by three times a week application of TPA (10 nmol/100
μl acetone). Animals were sacrificed when papilloma reached 1 mm and
skin biopsies were obtained, labelling the biopsies as their position in
a 20x20 mm grid. Liver samples were also obtained from the same
individual and an untreated C57Bl6/J mouse or NMR served as a control.
All samples were kept in deep-well 96-well plates and frozen at -80°C.

In every biopsy, genomic DNA was extracted using the AllPrep DNA/RNA
Micro Kit (Qiagen) and DNeasy Blood and Tissue Kit (Qiagen). AllPrep
DNA/RNA RLT buffer was added before homogenizing the sample, adding
proteinase K, incubation and DNA extraction. DNA quality was assessed
with Agilent Tapestation.

**3.2 Whole-exome and Whole-genome sequencing**
-----------------------------------------------

For whole exome sequencing, libraries were prepared following the
instructions of the SureSelectXT Low Inp ut HS and XT fractionation
mouse all exon target enrichment system (Agilent Technologies) or a
custom-designed NMR all exon capture array. 15-plex pools of exome
libraries were sequenced on two lanes Illumina NovaSeq6000 S2 flow cell
using a 100 bp paired-end read protocol.

For whole genome sequencing, WGS libraries were generated from 500 ng
high molecular weight genomic DNA using the NEBNext Ultra II FS DNA
Library Prep kit (NEB), according to the manufacturer’s instructions.
gDNA was fragmented for 15 min at 37C to target \~ 250 bp.
Quantification of the final libraries was carried out using the Qubit
dsDNA HS Assay Kit (Life Technologies), and the cDNA integrity was
assessed using D1000 ScreenTapes (Agilent Technologies). Libraries were
pooled to 10nM in 6-plex and sequenced on a HiSeq X Ten (Illumina) to
produce paired-end 150-bp reads. Each pool of 6 libraries was sequenced
over three (mouse) or four (NMR) lanes (minimum of 18X).

**3.3 Somatic Variant calling**
-------------------------------

From both whole genome and whole exome sequences, raw reads were stored
in a fastq file. Reads are first trimmed with *Trim Galore* and mapped
to the reference genome using *bwa mem*, different lanes from the same
sample are then merged and duplicate reads removed with gatk
*MarkDuplicates.* Depending on the comparison to be made, mapped reads
are subsampled, removing a percentage of reads to normalize all samples
to the minimal sequencing depth (number of reads per base) or not. After
this, somatic mutations are called using gatk *Mutect2* tumor-normal
calling using the mapped reads from the sample were somatic variants are
called (fibroblasts or skin samples) as “tumor” sample and the sample
containing germline variants (one of the untreated cultured fibroblasts
or the matched liver) as the “normal” sample. Somatic variants are then
filtered with gatk *FilterMutectCalls,* and the result is a variant
calling format (vcf) file.

**3.4 Somatic variant analysis and signature extraction**
---------------------------------------------------------

Somatic single nucleotide variants in a sample are considered based on
their mutational context and type of substitution, their mutational
spectrum. In order to infer which mutational signatures are present in a
sample non-negative matrix factorization (NMF) is used. Signatures are
extracted using the data analysis software and programming language R,
in which is used the implementation of NMF in the package *NMF* using D.
D. Lee & Seung, 1999 method. All other analyses and visualizations of
the variant calling data is also performed using R and different
packages implemented there. All code used for this thesis is available
in <https://github.com/JaviBotey/MasterThesis>.

**4 Results**
=============

The skin of a naked mole-rat and a mouse was treated with the mutagen
DMBA and the proliferation agent TPA. The somatic mutations present in
the skin of the treated animal and a control are relevant for the study
of clonal evolution of potentially cancerous lesions, it is relevant to
study the burden of these mutations, the process that originated them
and which gene functions do they affect. To start, the mutational burden
of treated and control skin is presented:

 **4.1 Mutational burden**
-------------------------

The general numbers of the sequencing can affect the downstream analysis
and are summarized on table 1.

| **Table 1**. Summary of samples sequenced from the two species (Naked Mole-Rat) and two conditions (DMBA/TPA treated and control).                                                                                                                                                       |                      |             |                      |             |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------|-------------|----------------------|-------------|
|                                                                                                                                                                                                                                                                                          | **Naked mole-rat**   | **Mouse**   |                      |             |
|                                                                                                                                                                                                                                                                                          | **DMBA/TPA-Treated** | **Control** | **DMBA/TPA-Treated** | **Control** |
| Number of biopsies                                                                                                                                                                                                                                                                       | 7                    | 4           | 5                    | 1           |
| Mean coverage                                                                                                                                                                                                                                                                            | ≈130,06X             | ≈225,08X    | ≈241,32X             | ≈241,24X    |
| On-target ratio                                                                                                                                                                                                                                                                          | 0,79                 | 0,8         | 0,77                 | 0,8         |
| Duplication rate                                                                                                                                                                                                                                                                         | 0,67                 | 0,6         | 0,51                 | 0,5         |
| Mean coverage represents average number of reads per base. On target ratio is calculated as the number of reads that cover the target sequences, the exome, divided by the total number of reads. Duplication rate represents the proportion of reads that are identified as duplicates. |                      |             |                      |             |

Among the four groups of the study, treated and control, naked mole-rat
(NMR) and mouse, one individual was used for each group, from which a
variable number of biopsies were obtained. The number of biopsies makes
statistical testing possible between treated and control in NMR and not
in mouse. All biopsies were deeply sequenced, as can be seen in the mean
coverage in table 1. Sequencing depth was normalised in cases where it
was necessary and full depth was used otherwise. The on-target ratio
represents the proportion of reads covering the exome and it is very
similar among all groups, where the exome coverage is adequate (table
1). The duplication rate indicates the proportion of reads identified as
duplicates and removed for downstream analysis and rates are variable
among groups, with treated NMR having the highest duplication rate and
lowest mean coverage (table 1).

<img src="media/image3.png" style="width:2.23958in;height:1.9125in" /><img src="media/image4.png" style="width:2.61667in;height:1.91736in" />

**Figure 1.** Mutation rate per group. Number of somatic single base
substitutions (SBS) called per megabase in naked mole-rat and mouse skin
both treated and control. Each dot represents a different biopsy.
Samples were normalized to the same sequencing depth (100X) before
somatic mutation calling was performed.

<img src="media/image5.png" style="width:5.53542in;height:2.93611in" />Mutational
burdens are represented in figure 1, calculated as the number of single
base substitutions (SBS) per megabase. If they are compared for each
group, it is found that the biggest difference due to treatment is found
in mouse, where the number of mutations called is 1.5 times bigger in
the treated in average. In NMR, however, changes are not so striking.
Also relevant to point is the difference in average coverage per group,
a factor that can affect the number of mutations called per sample. To
account for this factor the mutational burden is compared in Figure 1
subsampling all the samples to the same sequencing depth and calling
somatic mutations in those samples.

**Figure 2.** Variant Allelic Frequency (VAF) distribution per group.
Density plot representing how are VAF values distributed. Each color in
the plot represents a different biopsy. VAF are calculated based on the
proportion of reads that contain the mutation. Samples were normalized
to the same sequencing depth (100X) before somatic mutation calling was
performed.

Somatic mutations called are also found in different proportion of
sequenced cells, that is, in different allelic frequencies. The
distribution of allelic frequencies is therefore, also considered, and
showed for each group, subsampled to the same sequencing depth, in
Figure 2. Variant allelic frequencies (VAF) represent the proportion of
reads in which the variant is present. Average VAF are in a similar
range in all groups, around 0.1, representing that, on average, each
somatic variant is found in 1 out of 10 reads. Comparing them
species-wise and condition-wise differences in the distributions are
nonsignificant.

**4.2 Mutational spectrum**
---------------------------

SBS types and context is a relevant feature in the study of the effect
of a DMBA treatment, as its bias for T&gt;A transversions in a C.G, G.G
and T.G context is well characterized (McCreery et al., 2015). Focusing
on the SBS type proportions is, therefore, relevant among groups.

<img src="media/image6.png" style="width:4.38765in;height:2.64782in" />

**Figure 3.** Proportion of Thymine to Adenine (T&gt;A) transversions in
control and treated Naked mole-rat and mouse. Percentage of single base
substitutions (SBS) that are T&gt;A and proportion of those that are in
a cytosine guanine context (C.G). Wilcoxon rank sum test was performed
comparing the proportion of T&gt;A substitutions and T&gt;A
substitutions in a C.G context between the naked mole-rat (NMR) treated
skin samples (n=6) and the control skin samples (n=3). The difference in
T&gt;A proportion was nonsignificant, while the proportion of C.G T&gt;A
was significantly different between the treatment and the control (\*
p-value &lt; 0.05).

In Figure 3, the percentage of T&gt;A is displayed for the somatic
mutations in the two species and conditions. An increase of T&gt;A
proportion is seen in the treated biopsies for both species. But the
number of control samples (three in NMR and one in mouse) only allows
for statistical testing in NMR, in which no significant difference is
found. It is significant, nevertheless, the change of T&gt;A
transversion in a C.G context in both conditions of NMR.

<img src="media/image7.png" style="width:5.67926in;height:1.68344in" />

**Figure 4.** Proportion of each single base substitution (SBS) type in
naked mole-rat (NMR) and mouse treated and control skin samples.
Percentages of all six SBS types are displayed in different colours for
the different conditions: treated and control skin.

In Figure 4, the percentage of all six SBS types is displayed for both
species and conditions. Percentage changes between conditions are found,
but when testing for significance in NMR, no significant changes are
found. C&gt;A and T&gt;A both increase their proportion with the
treatment in both species. And T&gt;C proportion decreases in both
species. The other three SBS types diverge in the sense of their change
in both species.

<img src="media/image8.png" style="width:7.48264in;height:6.85764in" />

**Figure 5.** Somatic mutational spectrum of the Naked Mole-Rat.
Percentage of somatic Single Base Substitutions (SBS) for each
substitution type and context (5’ and 3’ bases). Each row represents the
average proportion of SBS type and context for each condition. DMBA
(n=6), the samples from the individual treated with DMBA/TPA. Control
(n=3), the samples from the individual treated with acetone. And the
subtraction of both conditions indicates the mutation types and contexts
enriched in the DMBA treatment. Colours indicate SBS types, and each SBS
type includes the complementary substitution (i.e C&gt;T includes C&gt;T
and G&gt;A).

<img src="media/image9.png" style="width:7.47917in;height:6.85694in" />

**Figure 6.** Somatic mutational spectrum of the mouse skin. Percentage
of somatic Single Base Substitutions (SBS) for each substitution type
and context (5’ and 3’ bases). Each row represents the average
proportion of SBS type and context for each condition. DMBA (n=4), the
samples from the individual treated with DMBA/TPA. Control (n=1), the
samples from the individual treated with acetone. And the subtraction of
both conditions indicates the mutation types and contexts enriched in
the DMBA treatment. Colours indicate SBS types, and each SBS type
includes the complementary substitution (i.e C&gt;T includes C&gt;T and
G&gt;A).

In Figure 5 and 6 the proportion of mutations of each substitution type
and context is shown for both species. And this is shown for the
mutations in the treated individual, in the control individual and the
difference between the treated and the control, showing which
substitutions and contexts have an increased presence in the treated
skin.

Figure 5 shows that the main substitution types in both treated and
control NMR skin are C to T and T to G SBS. And that the increased
mutation types in the treated skin are found in all substitutions and a
clear pattern is difficult to be extracted. However, according to figure
3, it is known that there is a significant increase in the presence of
the main SBS type and context found in this type of treatment
(CTG&gt;CAG).

Figure 6 shows that the main SBS and context in mouse skin, both control
and treated, are C to T and T&gt;A SBS. And that the increased mutation
types are mainly T&gt;A transversions in a CG context or a GG context,
with a two percent increase.

<img src="media/image10.png" style="width:6.68573in;height:3.79167in" />

**Figure 7.** Expanded context of T&gt;A substitutions. Values indicate
percentage of T&gt;A transversions that contain each nucleotide in each
position. Positions relative to the substitution in indicate the
direction (5’ or 3’) and the distance in number of nucleotides (from 1
until 5 nucleotides away from the substitution). Colors and displayed
labels indicate the nucleotide (C, T, A, G). Bars indicate percentages.

In figure 7, the proportion of bases in each position relative to a
T&gt;A transversion is shown. Considering that four different bases can
occupy a position, each has a 25% chance of being there. And looking at
figure 7, almost all positions and bases are in the 20%-30% range.
Different to this are the positions immediate to the transversion, where
C is predominant in 5’ and G in 3’ in all groups. This prevalent
presence of a C.G context is increased in the treatment condition for
both species, in which almost a half of transversions have a C in the 5’
base and almost a half of transversions have a G in the 3’ base.

**4.3 Mutational Signature in mouse embryonic fibroblasts** 
-----------------------------------------------------------

In order to obtain the mutational context and substitution type left by
the DMBA treatment, mouse Embryonic fibroblasts were used.

| **Table 2**. Summary of fibroblast samples whole genome sequenced from mouse embryonic fibroblasts in two conditions (DMBA treated and control)                                                                                                 |                  |             |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|-------------|
|                                                                                                                                                                                                                                                 | **DMBA-Treated** | **Control** |
| Number of biopsies                                                                                                                                                                                                                              | 3                | 3           |
| Mean sequencing depth                                                                                                                                                                                                                           | ≈10,18X          | ≈10,03X     |
| Number of SNVs per sample                                                                                                                                                                                                                       | 12060            | 21157       |
| Number of Mutations per sample                                                                                                                                                                                                                  | 13046            | 23021       |
| Mean coverage represents average number of reads per base. Number of SNVs per sample and Number of Mutations per sample indicate the number of somatic SNVs called in skin samples and total number of somatic mutations called in each sample. |                  |             |

Overall sequencing depth of all samples only had slight changes. Number
of somatic mutations called did show differences in the conditions and
the same holds for single nucleotide variants (SNV), that are the vast
majority of variants called (Table 1). Inspecting the location of the
called mutations, an abnormal accumulation of mutations was found in
some parts of chromosome 7 and 19, where mutation rates were much higher
than in any other chromosome and what was causing the differences in
mutations called in the whole genome (a 1.75 fold change is shown in
table 1, with an increased mutational burden in the untreated samples).
That

is why, for investigating the tumour burden of the samples, all
chromosomes but those were used.

<img src="media/image11.png" style="width:3.15167in;height:3.41899in" />

**Figure 8.** Mutational burden of mouse embryonic fibroblasts, DMBA
treated (n=3) and control (n=2). Number of Somatic Variants per megabase
called from whole genome sequencing data against a normal sample of
cultured fibroblasts.

<img src="media/image12.png" style="width:5.91181in;height:4.18958in" />
The number of variants called against a sample of untreated MEF in all
chromosomes but chromosomes 7 and 19 is shown in figure 8. Differences
are observed between control and DMBA-treated MEF. This set of variants
was used to determine the mutational signature left by DMBA on MEF.
These signatures, extracted using nonnegative matrix factorization
(NMF), are shown in figure 9.

### 

**Figure 9.** Mutational signatures extracted from Mouse embryonic
fibroblasts treated with DMBA and control using nonnegative matrix
factorization (NMF). Values indicate percentage of mutations for each
SBS and context.

Using matrix factorization algorithms lead to the extraction of 3
mutational signatures from the mutational spectrum of all five mouse
embryonic fibroblasts samples, three treated with DMBA and two
untreated. One was labelled as “DMBA-like signature” for being the one
that correlated the most with signatures found in the literature for
this type of treatment.

<img src="media/image13.png" style="width:3.84344in;height:2.04052in" />

**Figure 10.** Signature contribution of the three extracted signatures
(DMBA-like, Background 1 and Background 2) to the mutational burden in
each MEF sample.

This labelling of the signature as DMBA-like signature is confirmed when
looking at the contribution of each signature to the final number of
somatic mutations in all five samples (figure 10), because it is seen
that the DMBA-like signature is predominantly found in DMBA-treated
samples and not so often found in untreated samples. We can consider,
therefore, that the effect of DMBA leaves this specific signature in the
treated fibroblasts.

**4.4 Fitting the signature found in MEF to the skin samples**
--------------------------------------------------------------

After extracting the signatures present in MEF, they can be fitted,
using regression models, to investigate to which extent are they present
in the skin of both naked mole rat and mouse. Also, other signatures
extracted in the COSMIC project were considered, SBS7 and SBS38, as they
are linked to somatic signatures found in skin due to UV light exposure
and in cSCC.

<img src="media/image14.png" style="width:4.70139in;height:2.00972in" />

**Figure 11.** Signature fitting of the signatures extracted from MEF
and the COSMIC signatures SBS38 and SBS7. Percentages indicate relative
contribution to the total number of somatic SBS called, colours indicate
the different signatures, the three signatures extracted from MEF
(DMBA-like, background 1 and background 2) and two COSMIC signatures
(SBS7 and SBS38).

Inspecting the change in proportions found in the two conditions, an
increase in the share of mutations caused by the DMBA-like signature is
observed in the treated with respect to the control samples in both
species. Other differences are observed, in the background signatures of
MEF and in the COSMIC SBS signatures, in which there is an increase in
treated mouse skin in comparison to control mouse skin and a decrease in
treated NMR skin in comparison to control NMR skin.

 **4.5 Mutational signature in variant allelic frequency clusters**
------------------------------------------------------------------

<img src="media/image15.png" style="width:4.67708in;height:3.66875in" />
Variant allelic frequencies (VAF) of somatic variants found in the NMR
and mouse treated skin range between the very low values, present in 10
or 20% of the sequenced cells, and the higher ranges detected, present
in 40 or 50% of the sequenced cells (Figure 2). Considering the DMBA
treatment is applied and then a proliferation factor (TPA), it is safe
to assume DMBA-affected cells grow at a similar rate and mutations
appear at the same moment in cell fate, and therefore a specific range
of allelic frequency should contain most of the DMBA-caused variants.

**Figure 12.** Cosine similarity of subtracted DMBA control mutational
signatures and the DMBA signature in one percentile variant allelic
frequencies. Density plot of the distribution of VAF in treated skin of
NMR and mouse and heatmap of similarity to the DMBA signature in
McCreery et al., 2015 in each VAF cluster. Each interval of the heatmap
represents a VAF cluster, containing 1% of the somatic mutations called.

Figure 12 shows a difference in the similarity of the mutational
spectrums to the DMBA signature. A lower similarity is seen in both
species in the lowest VAF clusters, below 0.1 VAF. And a higher
similarity is seen after the VAF distribution peaks in both species,
with more VAF clusters with a higher similarity in mouse than in NMR.
The highest similarity is observed in NMR in the VAF cluster that
contains the mutations with the highest VAF, from 0.2 to 0.25.

**4.6 Nonsynonymous mutations in cancer driver genes**
------------------------------------------------------

Not all somatic SBS in the treated skin of both species influence
protein synthesis. Of those that affect the coding sequence of a gene,
some have a synonymous effect in protein sequence, while others can
disrupt the synthesis of the peptide, affect the aminoacidic sequence or
alter intronic sequences, among other effects. Considering only those
SBS that affect the coding sequence of a gene and do not have a
synonymous effect on protein synthesis, the identified mutations are
displayed on figure 13.

<img src="media/image16.png" style="width:7.66458in;height:3.73403in" />

**Figure 13.** Human cancer driver genes with nonsynonymous mutations,
and their estimated presence in the sample from allelic frequencies.
Each circle represents a single base substitution (SBS) with a
nonsynonymous effect, each color represents the gene name of the mutated
gene, size of the circle is the allelic frequency of the SBS and labels
are displayed for mutations with an allelic frequency above 0.1,
displaying the gene name of the mutated gene. Positions in the plot of
the circles are randomized for visualization purposes and, therefore,
coincidence of two mutations in the same coordinates is only an effect
of this randomization. Different squares represent different 1
mm<sup>3</sup> biopsy. An asterisk indicates the skin mouse biopsy
collected from a papilloma, which is the biopsy containing the
characteristic SBS found in the literature in a DMBA treatment: Hras
Q61L.

Figure 13 includes a visualization of nonsynonymous SBS that affect
described cancer related genes and their presence in each 1
mm<sup>3</sup> biopsy. A handful of this type of mutations are found in
each NMR sample displayed in the figure. The same is true for each mouse
biopsy, with the remarking difference that a higher number of these type
of mutations is found in the papilloma sample, a sample that also
includes a very well described change due to the DMBA treatment, Hras
Q61L.

<img src="media/image17.png" style="width:8.15972in;height:3.15208in" />

**Figure 14.** Nonsynonymous SBS in cancer related genes per biopsy.
Each coloured square represents the presence of a nonsynonymous SBS.
Each row represents a different biopsy named by its condition (treated
or control) and the sample id. The colour of the square represents the
species. Columns are ordered based on the number of biopsies containing
a nonsynonymous mutation in that gene; at the left, mutations present in
a higher number of biopsies, at the right, mutations present in only one
biopsy.

It is relevant to point out which cancer related genes are more commonly
mutated in each group (species and condition). This, and in which
biopsies is present is shown in figure 13. In mouse, 19 cancer related
genes are altered in more than one biopsy: *EPS15, PWWP2A, FGFR1OP,
MLLT10, CAMTA1, CHD4, DROSHA, FGFR2, FLT3, HMGA, KMT2D, MECOM, MITF,
NOTCH2, PRDM1, PTPRB, RANBP17, TERT, TRAP;* ordered from the most
recurrent (in 4 out of 5 samples), to the less recurrent (in 2 out of 5
samples). In treated NMR skin, 7 genes are altered in more than one
biopsy: *APC* (3/6 biopsies), *GMPS, LRP1B, ATM, MAP3K13, MTOR, TPM3*
(2/6 biopsies). In control NMR skin, 3 genes are altered in more than 1
biopsy: *FANCG, GPC3, PCSK7;* all mutated in 2 out of 3 biopsies.

**5. Discussion** 
=================

Cancer is one of the deadliest diseases in the world. Driven by
disruptions in cell functioning, the previous steps to it are the
accumulation of mutations. Studying somatic mutations is, therefore, the
study of cancer. And one way of studying how to treat cancer is studying
organisms that are resistant to the disease. Combining all this, we aim
to study the behaviour of somatic mutations in the cancer-resistant
naked mole rat.

More specifically, a well stablished model for cutaneous squamous cell
carcinoma was used. It consists on the application of a mutagen (DMBA)
and a proinflammatory chemical (TPA) to the skin of the animal. This
leads to appearance of mutations in the skin that can result in mutant
clones and progress into a malignant lesion. Comparing the effect of
this on NMR to the effect of this on mouse by sequencing treated and
control skin from both animals, we aim to investigate on their response
to the treatment.

Looking broadly into somatic mutations in NMR and mouse treated with
DMBA/TPA, a higher mutational burden is found on treated skin biopsies
of both species (figure 1). This burden is calculated in sequencings
subsampled to the same sequencing depth, and the variants have similar
distributions of allelic frequency in all biopsies (figure 2). Meaning
the latter that the variants were present in a similar proportion of
sequenced cells in all samples. Focusing only on single base
substitutions (SBS), a type of mutation called by the software used,
some characteristics of somatic SBS can be analyzed. Firstly, the SBS
type (C&gt;A, C&gt;T, C&gt;G, T&gt;A, T&gt;C, T&gt;G) and the immediate
context (the 5´and 3´nucleotides to the SBS) are used to classify SBSs
into 96 different types. The proportion of mutations that is included
into each type is known as the mutational spectrum of the sample. And,
because mutational processes tend to have biases for certain SBS types
and contexts, they can be used to infer the mutational processes active
using the mathematical framework of mutational signatures (Alexandrov et
al., 2012). The treatment used for this project has a bias for T&gt;A
transversions in C.G, G.G and T.G contexts (McCreery et al., 2015;
Nassar et al., 2015). Therefore, looking at number of T&gt;A
transversions in those three contexts can help identifying a treatment
effect. In figure 3, T&gt;A transversions increase with the treatment in
the two studied species. And that the first described context (C.G)
increases significantly in NMR, highlighting a possible treatment
effect. Comparing this to existing DMBA/TPA literature, it is consistent
with an increase in T&gt;A mutations, although the reported proportion
of T&gt;A overall is lower in our model (around 15% in NMR and 25% in
mouse) than in the cited article (over 60%), possibly due to the fact
that the samples reported in the article are extracted from papillomas
and the sequenced cells are only epithelial, while our model includes
nonpapilloma samples and all cell types within the volume extracted
(Nassar et al., 2015).

The somatic mutational spectrum of NMR and mouse control and DMBA/TPA
treated is depicted in figures 5 and 6. Where the subtracted mutational
spectrum, calculated by subtracting the proportion of mutations in each
of the 96 mutation types in the control biopsies from the treated
biopsies, represents the SBS types and contexts that are enriched in the
treated skin. This approach is similar to the one followed at
Pleguezuelos-Manzano et al., 2020. This subtracted spectrum highlights
that the SBS types and contexts mainly increased in treated mouse skin
are the same reported in the DMBA literature, T&gt;A in a C.G context
and T&gt;A in a G.G context, with a two percent increase in the treated
with respect to the control skin. This indicates that there is an
increased mutagenesis due to the DMBA treatment in mouse skin. In NMR,
both SBS and contexts are also increased in the treated skin together
with other contexts of T&gt;A. But because this increase is lower than
in mouse and a similar increase is found in other SBS and contexts, the
effect of the treatment is not as clear.

Some studies of somatic mutations using the mutational signature
framework indicate that the use of only the immediate context of SBS may
fall short to illustrate the mutational processes present (Zhang et al.,
2020). And this is also found in other studies where it is proved that
some mutational processes act also leaving a signature in larger
contexts due to their mutational mechanism. For instance, the metabolite
in the cited article forms adenine adducts between the mutated base and
a base three nucleotides away (Pleguezuelos-Manzano et al., 2020).
Because of that, a broader context of T&gt;A transversions was also
considered in this project. Revealing no significance difference in
nucleotides that were not the immediate to the mutation (figure 7).
However, clustering the biopsies based on the cosine similarity of their
immediate context and based on their expanded, 5 nucleotide context, a
better group separation is achieved using the latter, indicating that
differences in the mutational context between the groups lay more than 1
nucleotide away from the mutation.

The next step in the project is characterizing how does DMBA act on
mouse embryonic fibroblasts (MEF), with the goal of trying to find this
same signature on mouse and NMR skin. For this, MEF were cultured and
either treated with DMBA or with DMSO (control), and mutations were
called in the DMBA-treated and control samples using one control sample
as “normal”. Mutational spectrums for all samples were obtained and
three mutational signatures were extracted from those spectrums using
NMF (figure 9). In these three extracted signatures, one had the highest
similarity to the DMBA signature found in the literature and had a
higher presence in DMBA treated MEF (figure 10). Which allowed to label
it as the DMBA signature, and the other signatures extracted as
“Background” signatures. The next step was fitting these signatures
together with other signatures that can be found on the skin, SBS7 and
SBS38. SBS7 and SBS38 are signatures extracted in the project catalogue
of somatic mutations in cancer (COSMIC) from human skin cancers and
whose etiology is linked to UV light exposure (Alexandrov et al., 2020).
In figure 11, the proportion of all fitted mutations attributed to each
signature is shown. And it is seen that there are more mutations
attributed to the DMBA signature in the treated biopsies of both NMR and
mouse, indicating a DMBA treatment effect.

Up until this point, the somatic mutations called in each biopsy were
considered a single unit, as if a population of clonal cells had been
sequenced. But nothing further from reality. The population of sequenced
cells emerges from different cell lineages present in the animal’s skin,
epithelial, mesenchymal, immune cells, all of them are present in the
skin and sequenced. And within each cell type different subclonal
populations can be found. Somatic mutations called can, therefore, be
present in any of these lineages and subclones. Variant allele frequency
(VAF) of the mutations called was used to represent proportion of
sequenced cells containing each mutation. Mutations due to the treatment
(and its mutational signature) should be found in a specific subclonal
population with a specific VAF (Nik-Zainal et al., 2012). In our case,
this is modelled using ranges of VAF. More complex models use a similar
approach to infer the presence of the mutational signature using other
traits of the mutations (Harrigan et al., 2020; Rubanova et al., 2020).
In our case, the bulk of somatic variants called in each species and
treatment were clustered in 100 groups based on their VAF, and for each
group the same approach followed in figures 5 and 6 was used. The
mutational spectrum of the subtraction was obtained. To compare this
spectrum to one of a DMBA treatment, the cosine similarity of this
spectrum and the DMBA spectrum in McCreery et al., 2015 was computed and
plotted in figure 12. The similarity of some VAF clusters to the DMBA
signature was high (up to 0.5) and indicated that the subclones derived
from the treatment were present in those VAF clusters. The signature is
present in mouse in more VAF clusters (observing the number of stripes
with a high similarity in figure 12) than in NMR. And in NMR the
signature is found with the highest similarity in the last VAF cluster
(from 0.2 to 0.25), indicating that the mutations due to the treatment
are present in 40-50% of the sequenced cells in each biopsy. All this
indicates that there is a treatment effect in both species. Moreover,
this effect is only found on certain VAF clusters and the subclones
derived from the treatment represent more than 10% of the sequenced
cells in both species. This last finding can be because, after the DMBA
treatment, TPA was applied, leading to a proliferation of the mutated
cells.

The last trait studied of the somatic mutations called is which genes do
they affect. Especially those that alter the aminoacidic sequence of the
expressed protein. Because, as mentioned at the beginning of the
discussion, the disruption of driver genes can have a determinant effect
on the malignant progression of a tissue. In figures 13 and 14, these
disrupting mutations on cancer driver genes are illustrated. They alter
many genes in each treated biopsy of the two species (figure 13). And,
as expected, more cancer driver genes are hit in the sample obtained
from a papilloma. In which a specific DMBA-related mutation is found,
*Hras* Q61L (Nassar et al., 2015). Other cancer genes carry
nonsynonymous mutations in treated mouse skin such as *NOTCH2* or *KIT.*
And, surprisingly for a cancer-resistant organism, cancer genes are also
altered in NMR such as *MTOR, APC* or *NOTCH2.*

All in all, the expected disruptions are found in the DMBA/TPA treated
mouse skin. An increase in T&gt;A transversions is seen (figure 3), the
DMBA signature is found in the subtracted mutational spectrum (figure
6), the signature found in MEF contributes to more mutations in the
treated skin (figure 11), VAF cluster exhibit a good similarity to the
DMBA signature reported in the literature (figure 12) and cancer related
genes are hit with the treatment, including *HRAS* Q61L (figure 13).
Corresponding to the findings in the literature.

However, there is no existing evidence on the use of a DMBA/TPA
treatment on NMR. And there was not any evidence of papillomas in the
skin of the treated animal. Looking broadly at all somatic mutations, an
increase in T&gt;A transversions was detected on NMR, and a
statistically significant increase of CTG&gt;CAG is found (figure 3).
And, although the DMBA signature is not so clear on the subtracted
mutational spectrum (figure 5), the signature extracted from MEF is
contributing to more mutations in the treated skin (figure 11). To this
it can be added that the signature is found on certain VAF clusters
(figure 12). Concluding that the DMBA is causing mutations on the skin
of the animal. But this effect is shier than in the mouse, as it cannot
be seen looking at the subtracted mutational spectrum. Previous reports
on NMR cancer resistance point to a difference in the DNA repair
machinery of NMR (MacRae et al., 2015). From our results it can be seen
that NMR is not resistant to the effect of chemical mutagens and the
effect of this is similar to the effect found in mouse, meaning that
their repair machineries are acting similarly on the effect of the
mutagen. Other aspects that are attributed to NMR cancer-resistance
cannot be evaluated from our results (Seluanov et al., 2009; Tian et
al., 2013). To sum up, the naked mole-rat was resistant to malignant
lesions caused by the DMBA/TPA treatment, but not to somatic mutations
caused by it. And the effect of the DMBA treatment is observed less
evidently in naked mole-rat compared with mouse.

**Bibliography**
================

Aitken, S. J., Anderson, C. J., Connor, F., Pich, O., Sundaram, V.,
Feig, C., Rayner, T. F., Lukk, M., Aitken, S., Luft, J., Kentepozidou,
E., Arnedo-Pac, C., Beentjes, S. V., Davies, S. E., Drews, R. M., Ewing,
A., Kaiser, V. B., Khamseh, A., López-Arribillaga, E., … Taylor, M. S.
(2020). Pervasive lesion segregation shapes cancer genome evolution.
*Nature*, *583*(7815), 265–270.
https://doi.org/10.1038/s41586-020-2435-1Alexandrov, L. B., Kim, J.,
Haradhvala, N. J., Huang, M. N., Tian Ng, A. W., Wu, Y., Boot, A.,
Covington, K. R., Gordenin, D. A., Bergstrom, E. N., Islam, S. M. A.,
Lopez-Bigas, N., Klimczak, L. J., McPherson, J. R., Morganella, S.,
Sabarinathan, R., Wheeler, D. A., Mustonen, V., Alexandrov, L. B., …
Stratton, M. R. (2020). The repertoire of mutational signatures in human
cancer. *Nature*, *578*(7793), 94–101.
https://doi.org/10.1038/s41586-020-1943-3Alexandrov, L. B., Nik-Zainal,
S., Wedge, D. C., Aparicio, S. A. J. R., Behjati, S., Biankin, A. V.,
Bignell, G. R., Bolli, N., Borg, A., Børresen-Dale, A. L., Boyault, S.,
Burkhardt, B., Butler, A. P., Caldas, C., Davies, H. R., Desmedt, C.,
Eils, R., Eyfjörd, J. E., Foekens, J. A., … Stratton, M. R. (2013).
Signatures of mutational processes in human cancer. *Nature*,
*500*(7463), 415–421. https://doi.org/10.1038/nature12477Alexandrov, L.
B., Nik-zainal, S., Wedge, D. C., Campbell, P. J., & Stratton, M. R.
(2012). Deciphering Signatures of Mutational Processes Operative in
Human Cancer. *CellReports*, *3*(1), 246–259.
https://doi.org/10.1016/j.celrep.2012.12.008Alexandrov, L. B.,
Nik-Zainal, S., Wedge, D. C., Campbell, P. J., & Stratton, M. R. (2013).
Deciphering Signatures of Mutational Processes Operative in Human
Cancer. *Cell Reports*, *3*(1), 246–259.
https://doi.org/10.1016/j.celrep.2012.12.008Amberg, N., Holcmann, M.,
Glitzner, E., Novoszel, P., Stulnig, G., & Sibilia, M. (2015). Mouse
models of nonmelanoma skin cancer. *Methods in Molecular Biology*,
*1267*, 217–250.
https://doi.org/10.1007/978-1-4939-2297-0\_10Baez-Ortega, A., & Gori, K.
(2019). Computational approaches for discovery of mutational signatures
in cancer. *Briefings in Bioinformatics*, *20*(1), 77–88.
https://doi.org/10.1093/bib/bbx082Brandt, M. G., & Moore, C. C. (2019).
Nonmelanoma Skin Cancer. In *Facial Plastic Surgery Clinics of North
America* (Vol. 27, Issue 1, pp. 1–13). W.B. Saunders.
https://doi.org/10.1016/j.fsc.2018.08.001Bray, F., Ferlay, J.,
Soerjomataram, I., Siegel, R. L., Torre, L. A., & Jemal, A. (2018).
Global cancer statistics 2018: GLOBOCAN estimates of incidence and
mortality worldwide for 36 cancers in 185 countries. *CA: A Cancer
Journal for Clinicians*, *68*(6), 394–424.
https://doi.org/10.3322/caac.21492Brunet, J. P., Tamayo, P., Golub, T.
R., & Mesirov, J. P. (2004). Metagenes and molecular pattern discovery
using matrix factorization. *Proceedings of the National Academy of
Sciences of the United States of America*, *101*(12), 4164–4169.
https://doi.org/10.1073/pnas.0308531101Cai, Y., Patel, D. J., Broyde,
S., & Geacintov, N. E. (2010). Base sequence context effects on
nucleotide excision repair. In *Journal of Nucleic Acids* (Vol. 2010).
Hindawi Limited. https://doi.org/10.4061/2010/174252Choi, J. H.,
Besaratinia, A., Lee, D. H., Lee, C. S., & Pfeifer, G. P. (2006). The
role of DNA polymerase ι in UV mutational spectra. *Mutation Research -
Fundamental and Molecular Mechanisms of Mutagenesis*, *599*(1–2), 58–65.
https://doi.org/10.1016/j.mrfmmm.2006.01.003Drost, J., Van Boxtel, R.,
Blokzijl, F., Mizutani, T., Sasaki, N., Sasselli, V., De Ligt, J.,
Behjati, S., Grolleman, J. E., Van Wezel, T., Nik-Zainal, S., Kuiper, R.
P., Cuppen, E., & Clevers, H. (2017). Use of CRISPR-modified human stem
cell organoids to study the origin of mutational signatures in cancer.
*Science*, *358*(6360), 234–238.
https://doi.org/10.1126/science.aao3130Gao, Y. B., Chen, Z. L., Li, J.
G., Hu, X. Da, Shi, X. J., Sun, Z. M., Zhang, F., Zhao, Z. R., Li, Z.
T., Liu, Z. Y., Zhao, Y. Da, Sun, J., Zhou, C. C., Yao, R., Wang, S. Y.,
Wang, P., Sun, N., Zhang, B. H., Dong, J. S., … He, J. (2014). Genetic
landscape of esophageal squamous cell carcinoma. *Nature Genetics*,
*46*(10), 1097–1102. https://doi.org/10.1038/ng.3076Hanahan, D., &
Weinberg, R. A. (2011). Hallmarks of cancer: The next generation. In
*Cell* (Vol. 144, Issue 5, pp. 646–674). Elsevier.
https://doi.org/10.1016/j.cell.2011.02.013Harrigan, C. F., Rubanova, Y.,
Morris, Q., & Selega, A. (2020). TrackSigFreq: Subclonal reconstructions
based on mutation signatures and allele frequencies. *Pacific Symposium
on Biocomputing*, *25*(2020), 238–249.
https://doi.org/10.1142/9789811215636\_0022Hecht, S. S. (2003). Tobacco
carcinogens, their biomarkers and tobacco-induced cancer. In *Nature
Reviews Cancer* (Vol. 3, Issue 10, pp. 733–744). European Association
for Cardio-Thoracic Surgery. https://doi.org/10.1038/nrc1190Hilton, H.
G., Rubinstein, N. D., Janki, P., Ireland, A. T., Bernstein, N., Fong,
N. L., Wright, K. M., Smith, M., Finkle, D., Martin-McNulty, B., Roy,
M., Imai, D. M., Jojic, V., & Buffenstein, R. (2019). Single-cell
transcriptomics of the naked molerat reveals unexpected features of
mammalian immunity. *PLoS Biology*, *17*(11), e3000528.
https://doi.org/10.1371/journal.pbio.3000528Jans, J., Schul, W., Sert,
Y. G., Rijksen, Y., Rebel, H., Eker, A. P. M., Nakajima, S., Van Steeg,
H., De Gruijl, F. R., Yasui, A., Hoeijmakers, J. H. J., & Van Der Horst,
G. T. J. (2005). Powerful skin cancer protection by a CPD-photolyase
transgene. *Current Biology*, *15*(2), 105–115.
https://doi.org/10.1016/j.cub.2005.01.001Lee, B. P., Smith, M.,
Buffenstein, R., & Harries, L. W. (2020). Negligible senescence in naked
mole rats may be a consequence of well-maintained splicing regulation.
*GeroScience*, *42*(2), 633–651.
https://doi.org/10.1007/s11357-019-00150-7Lee, D. D., & Seung, H. S.
(1999). Learning the parts of objects by non-negative matrix
factorization. *Nature*, *401*(6755), 788–791.
https://doi.org/10.1038/44565Li, S., Crawford, F. W., & Gerstein, M. B.
(2020). Using sigLASSO to optimize cancer mutation signatures jointly
with sampling likelihood. *Nature Communications*, *11*(3575).
https://doi.org/10.1038/s41467-020-17388-xLowry, W. E., Flores, A., &
White, A. C. (2016). Exploiting Mouse Models to Study Ras-Induced
Cutaneous Squamous Cell Carcinoma. In *Journal of Investigative
Dermatology* (Vol. 136, Issue 8, pp. 1543–1548). Elsevier B.V.
https://doi.org/10.1016/j.jid.2016.03.017Lyu, X., Garret, J., Rätsch,
G., & Lehmann, K. Van. (2020). Mutational signature learning with
supervised negative binomial non-negative matrix factorization.
*Bioinformatics (Oxford, England)*, *36*(1), i154–i160.
https://doi.org/10.1093/bioinformatics/btaa473MacRae, S. L., Croken, M.
M. K., Calder, R. B., Aliper, A., Milholland, B., White, R. R.,
Zhavoronkov, A., Gladyshev, V. N., Seluanov, A., Gorbunova, V., Zhang,
Z. D., & Vijg, J. (2015). DNA repair in species with extreme lifespan
differences. *Aging*, *7*(12), 1171–1184.
https://doi.org/10.18632/aging.100866Martincorena, I., Roshan, A.,
Gerstung, M., Ellis, P., Van Loo, P., McLaren, S., Wedge, D. C., Fullam,
A., Alexandrov, L. B., Tubio, J. M., Stebbings, L., Menzies, A., Widaa,
S., Stratton, M. R., Jones, P. H., & Campbell, P. J. (2015). High burden
and pervasive positive selection of somatic mutations in normal human
skin. *Science*, *348*(6237), 880–886.
https://doi.org/10.1126/science.aaa6806McCreery, M. Q., Halliwill, K.
D., Chin, D., Delrosario, R., Hirst, G., Vuong, P., Jen, K. Y.,
Hewinson, J., Adams, D. J., & Balmain, A. (2015). Evolution of
metastasis revealed by mutational landscapes of chemically induced skin
cancers. *Nature Medicine*, *21*(12), 1514–1520.
https://doi.org/10.1038/nm.3979Miyata, M., Furukawa, M., Takahashi, K.,
Gonzalez, F. J., & Yamazoe, Y. (2001). Mechanism of
7,12-Dimethylbenz\[a\]anthracene-Induced Immunotoxicity: Role of
Metabolic Activation at the Target Organ. In *Jpn. J. Pharmacol* (Vol.
86).Morris, L. G. T., Kaufman, A. M., Gong, Y., Ramaswami, D., Walsh, L.
A., Turcan, Ş., Eng, S., Kannan, K., Zou, Y., Peng, L., Banuchi, V. E.,
Paty, P., Zeng, Z., Vakiani, E., Solit, D., Singh, B., Ganly, I., Liau,
L., Cloughesy, T. C., … Chan, T. A. (2013). Recurrent somatic mutation
of FAT1 in multiple human cancers leads to aberrant Wnt activation.
*Nature Genetics*, *45*(3), 253–261.
https://doi.org/10.1038/ng.2538Nassar, D., Latil, M., Boeckx, B.,
Lambrechts, D., & Blanpain, C. (2015). Genomic landscape of
carcinogen-induced and genetically induced mouse skin squamous cell
carcinoma. *Nature Medicine*, *21*(8), 946–954.
https://doi.org/10.1038/nm.3878Nik-Zainal, S., Van Loo, P., Wedge, D.
C., Alexandrov, L. B., Greenman, C. D., Lau, K. W., Raine, K., Jones,
D., Marshall, J., Ramakrishna, M., Shlien, A., Cooke, S. L., Hinton, J.,
Menzies, A., Stebbings, L. A., Leroy, C., Jia, M., Rance, R., Mudie, L.
J., … Campbell, P. J. (2012). The life history of 21 breast cancers.
*Cell*, *149*(5), 994–1007.
https://doi.org/10.1016/j.cell.2012.04.023Omichessan, H., Severi, G., &
Perduca, V. (2019). Computational tools to detect signatures of
mutational processes in DNA from tumours: A review and empirical
comparison of performance. *PLoS ONE*, *14*(9), 1–28.
https://doi.org/10.1371/journal.pone.0221235Petruseva, I. O., Evdokimov,
A. N., & Lavrik, O. I. (2017). Genome Stability Maintenance in Naked
Mole-Rat. *Acta Naturae*, *9*(4), 31–41.
http://www.ncbi.nlm.nih.gov/pubmed/29340215Pitolli, C., Wang, Y., Candi,
E., Shi, Y., Melino, G., & Amelio, I. (2019). P53-mediated tumor
suppression: DNA-damage response and alternative mechanisms. In
*Cancers* (Vol. 11, Issue 12). MDPI AG.
https://doi.org/10.3390/cancers11121983Pleguezuelos-Manzano, C.,
Puschhof, J., Rosendahl Huber, A., van Hoeck, A., Wood, H. M., Nomburg,
J., Gurjao, C., Manders, F., Dalmasso, G., Stege, P. B., Paganelli, F.
L., Geurts, M. H., Beumer, J., Mizutani, T., Miao, Y., van der Linden,
R., van der Elst, S., Ambrose, J. C., Arumugam, P., … Clevers, H.
(2020). Mutational signature in colorectal cancer caused by genotoxic
pks + E. coli. *Nature*, *580*(7802), 269–273.
https://doi.org/10.1038/s41586-020-2080-8Poplin, R., Ruano-Rubio, V.,
DePristo, M., Fennell, T., Carneiro, M., Van der Auwera, G., Kling, D.,
Gauthier, L., Levy-Moonshine, A., Roazen, D., Shakir, K., Thibault, J.,
Chandran, S., Whelan, C., Lek, M., Gabriel, S., Daly, M., Neale, B.,
MacArthur, D., & Banks, E. (2017). Scaling accurate genetic variant
discovery to tens of thousands of samples. *BioRxiv*, 201178.
https://doi.org/10.1101/201178Rajapaksa, K. S. (2010). Ovarian
Metabolism of Xenobiotics. In *Comprehensive Toxicology, Second Edition*
(Vol. 11, pp. 457–467). Elsevier Inc.
https://doi.org/10.1016/B978-0-08-046884-6.01127-1Robinson, P. S.,
Coorens, T. H. H., Palles, C., Mitchell, E., Abascal, F., Olafsson, S.,
Lee, B., Lawson, A. R. J., Lee-Six, H., Moore, L., Sanders, M. A.,
Hewinson, J., Martin, L., Pinna, C. M. A., Galvotti, S., Campbell, P.
J., Martincorena, I., Tomlinson, I., & Stratton, M. (2020). Elevated
somatic mutation burdens in normal human cells due to defective DNA
polymerases. *BioRxiv*, 2020.06.23.167668.
https://doi.org/10.1101/2020.06.23.167668Ross, J. A., & Nesnow, S.
(1999). Polycyclic aromatic hydrocarbons: Correlations between DNA
adducts and ras oncogene mutations. *Mutation Research - Fundamental and
Molecular Mechanisms of Mutagenesis*, *424*(1–2), 155–166.
https://doi.org/10.1016/S0027-5107(99)00016-0Rubanova, Y., Shi, R.,
Harrigan, C. F., Li, R., Wintersinger, J., Sahin, N., Deshwar, A.,
Dentro, S. C., Leshchiner, I., Gerstung, M., Jolly, C., Haase, K.,
Tarabichi, M., Wintersinger, J., Deshwar, A. G., Yu, K., Gonzalez, S.,
Rubanova, Y., Macintyre, G., … Morris, Q. (2020). Reconstructing
evolutionary trajectories of mutation signature activities in cancer
using TrackSig. *Nature Communications*, *11*(1), 1–12.
https://doi.org/10.1038/s41467-020-14352-7Saunders, C. T., Wong, W. S.
W., Swamy, S., Becq, J., Murray, L. J., & Cheetham, R. K. (2012).
Strelka: Accurate somatic small-variant calling from sequenced
tumor-normal sample pairs. *Bioinformatics*, *28*(14), 1811–1817.
https://doi.org/10.1093/bioinformatics/bts271Seluanov, A., Hine, C.,
Azpurua, J., Feigenson, M., Bozzella, M., Mao, Z., Catania, K. C., &
Gorbunova, V. (2009). Hypersensitivity to contact inhibition provides a
clue to cancer resistance of naked mole-rat. *Proceedings of the
National Academy of Sciences of the United States of America*,
*106*(46), 19352–19357. https://doi.org/10.1073/pnas.0905252106Tian, X.,
Azpurua, J., Hine, C., Vaidya, A., Myakishev-Rempel, M., Ablaeva, J.,
Mao, Z., Nevo, E., Gorbunova, V., & Seluanov, A. (2013).
High-molecular-mass hyaluronan mediates the cancer resistance of the
naked mole rat. *Nature*, *499*(7458), 346–349.
https://doi.org/10.1038/nature12234Wang, N. J., Sanborn, Z., Arnett, K.
L., Bayston, L. J., Liao, W., Proby, C. M., Leigh, I. M., Collisson, E.
A., Gordon, P. B., Jakkula, L., Pennypacker, S., Zou, Y., Sharma, M.,
North, J. P., Vemula, S. S., Mauro, T. M., Neuhaus, I. M., LeBoit, P.
E., Hur, J. S., … Cho, R. J. (2011). Loss-of-function mutations in Notch
receptors in cutaneous and lung squamous cell carcinoma. *Proceedings of
the National Academy of Sciences of the United States of America*,
*108*(43), 17761–17766. https://doi.org/10.1073/pnas.1114669108Xu, C.
(2018). A review of somatic single nucleotide variant calling algorithms
for next-generation sequencing data. *Computational and Structural
Biotechnology Journal*, *16*, 15–24.
https://doi.org/10.1016/j.csbj.2018.01.003Zhang, Y., Xiao, Y., Yang, M.,
& Ma, J. (2020). Cancer mutational signatures representation by
large-scale context embedding. *Bioinformatics (Oxford, England)*,
*36*(1), i309–i316. https://doi.org/10.1093/bioinformatics/btaa433Zhao,
J., Tian, X., Zhu, Y., Zhang, Z., Rydkina, E., Yuan, Y., Zhang, H., Roy,
B., Cornwell, A., Nevo, E., Shang, X., Huang, R., Kristiansen, K.,
Seluanov, A., Fang, X., & Gorbunova, V. (2020). Reply to: Transformation
of naked mole-rat cells. In *Nature* (Vol. 583, Issue 7814, pp. E8–E13).
Nature Research. https://doi.org/10.1038/s41586-020-2411-9
