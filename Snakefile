#A tool that can be useful in snakemake but I have not used here, wildcard_constraint, gives a regular expressionto which the wildcard in the rule has to stick.
#i.e: wildcard_constraint: tumor_sample="AS-\d+"

import yaml
import os
import glob
from itertools import product
import json
import sys
import re

#First, the function read_config(). It will automatically read the file read_config, which is a json file containing info, as a config.

def read_config():
	with open(ref_config) as json_file:
		return json.load(json_file)

#Now, let's choose the species --config reference=hgla or --config reference=mmus.
#Depending on with which species are we working on a different config file, genome and directory will be used for the whole analysis.
#After hours trying, turns out I dont have to add commas here :) (self-reminder)
if config["reference"] == "hgla":
  ref_config="/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Reference/reference-hgla.json"
  genome="/icgc/dkfzlsdf/analysis/B210/references_genome/Heterocephalus_glaber_female.HetGla_female_1.0.dna_rm.toplevel_short_ids.fa"
  full_path="/icgc/dkfzlsdf/analysis/B210/Javi/hgla/"
  dict="/icgc/dkfzlsdf/analysis/B210/references_genome/Heterocephalus_glaber_female.HetGla_female_1.0.dna_rm.toplevel_short_ids.dict"
  chrname2=expand("{chr}", chr=glob_wildcards("/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Reference/chr/split_H_{chr}.bed").chr )
  kb_bed= "/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Reference/NMR_1kb_intervals.bed"
  chrpath="/icgc/dkfzlsdf/analysis/B210/Javi/hgla/Reference/chr/split_H_"
  config.update(read_config())
elif config["reference"] == "mmus":
  ref_config="/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Reference/reference-mmus.json"
  genome="/icgc/dkfzlsdf/analysis/B210/references_genome/Mus_musculus.GRCm38.dna_rm.toplevel_short_IDs.fa"
  full_path="/icgc/dkfzlsdf/analysis/B210/Javi/mmus/"
  dict="/icgc/dkfzlsdf/analysis/B210/references_genome/Mus_musculus.GRCm38.dna_rm.toplevel_short_IDs.dict"
  chrname2=expand("{chr}", chr=glob_wildcards("/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Reference/chr/{chr}.bed").chr )
  kb_bed="/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Reference/mmus_1kb_intervals.bed" #TODO
  chrpath="/icgc/dkfzlsdf/analysis/B210/Javi/mmus/Reference/chr/"
  config.update(read_config())
else:
  quit()


#Now, let's choose the project --config project=17618. This will indicate the folder conatining the samples to analyze
project = config["project"]

#The species path and the project we select will be our path.
path = str(full_path)+str(project)

#The project folder needs to contain the following folders: fastq, alignment, QC, subsampling and mutational_analysis.
#Inside the fastq folder, there are the fastq.gz we want to analyze. From this, we will extract the following wildcards.
#The name of the fastq.gz I am using is the following: AS-475131_TUMOR-LR-lane2_R1.fastq.gz

#The full name of the sample, with the sample id and the lane id
name_full = expand("{ID}",ID=glob_wildcards(str(path)+"/fastq/{id}_R2.fastq.gz").id )

#The sample id (i.e AS-472423) and _ the phenotype of the file, TUMOR is the sample containing the somatic mutations, in our case, a skin sample. NORMAL indicates the sample containing the germline mutations.
name_sample = expand("{ID}",ID=glob_wildcards(str(path)+"/fastq/{id}-LR-{other}_R2.fastq.gz").id )

#For mutect, we select only tumor samples and extract their names as a wildcard.
tumor_sample = expand("{ID}",ID=glob_wildcards(str(path)+"/fastq/{id}_TUMOR-LR-{other}_R2.fastq.gz").id )

#The lane id
lane = expand("{ID}",ID=glob_wildcards(str(path)+"/fastq/{other}-LR-{id}_R2.fastq.gz").id )
#Get only unique values for the lane names
name_lane=list(set(lane))
normal_sample = expand("{ID}",ID=glob_wildcards(str(path)+"/fastq/{id}_NORMAL-LR-"+str(name_lane[1])+"_R2.fastq.gz").id )
normal_sample = str(normal_sample).strip('[]').strip("''")

#Now we will extract information from the config file of the location of the target files for the depth calculation
TARGET_BED = os.path.join(config["referenceInfo"]["path"], config["annotation"]["wesTargetFile"])
TARGET_BED_SORTED=os.path.join(config["referenceInfo"]["path"], config["annotation"]["wesTargetFile"]).replace(".bed", ".sorted.bed")
NONTARGET_BED = os.path.join(config["referenceInfo"]["path"], config["annotation"]["wesTargetFile"]).replace(".bed", ".sorted.complement.bed")
GENOMEFILE =  os.path.join(config["referenceInfo"]["path"], config["referenceInfo"]["genomeFile"])


#As in every Snakemake, we need to set all the output files as inputs in rule all. To avoid the use of wildcards in the final rule of the file.
rule all:
	input:
		expand(path+"/fastq/{name}_R{number}.fastq", name = name_full, number = ["1", "2"]),
		expand(path+"/fastq/{name_full}_R{number}_val_{number}.fq", number=["1","2"], name_full = name_full),
		expand(path+"/alignment/{name_full}.bam", name_full = name_full),
		expand(path+"/alignment/{name_full}_sort.bam", name_full = name_full),
		expand(path+"/alignment/{name_full}_no_dupl_sort.bam", name_full = name_full),
		expand(path+"/alignment/{name_sample}_merged.bam", name_sample = name_sample),
		expand(path+"/alignment/{name_sample}_merged.bam.bai", name_sample = name_sample),
		expand(path+"/QC/{name_sample}.{region}.bed.gz", name_sample=name_sample, region=["target", "nontarget"]),
		expand(path+"/QC/{name_sample}.targets_coverageQC.tsv", name_sample=name_sample),
		expand(path+"/QC/{name_sample}.coverageQC.tsv", name_sample=name_sample),
		path+"/QC/qc-targetcoverage.html",
		expand(path+"/subsampling/{name_sample}_merged_subsampled.bam", name_sample = name_sample),
		expand(path+"/subsampling/{name_sample}_merged_subsampled_RG_LG.bam", name_sample = name_sample),
		expand(path+"/subsampling/{name_sample}_merged_subsampled_RG_LG.bam.bai", name_sample = name_sample),
		expand(path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output.vcf", tumor_sample = tumor_sample, chrname2 = chrname2),
		expand(path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf", tumor_sample = tumor_sample, chrname2 = chrname2),
		expand(path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf.gz", tumor_sample = tumor_sample, chrname2 = chrname2),
		expand(path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered.vcf", tumor_sample = tumor_sample),
		expand(path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered_sort.vcf", tumor_sample = tumor_sample),
		expand(path+"/Mutation_calling/{tumor_sample}/pass.vcf", tumor_sample = tumor_sample),
		expand(path+"/Mutation_calling/1kb_intervals/{tumor_sample}.tsv", tumor_sample = tumor_sample),
		expand(path+"/Mutation_calling/AF/{tumor_sample}.tsv", tumor_sample = tumor_sample) 


#It is also needed to set the order of the rules.
ruleorder: gzip > trim_galore > bwa_mem > samtools_sort > remove_duplicates > samtools_merge > samtools_index > depth_targets > depth_nontargets > calculateTargetSize > targets_coverageQC > coverageQC > rmdreport > samtools_subsample > add_readgroup_subsampled > samtools_index_subsampled > Mutect2 > filter_vcf > index_vcf > merge_vcf > sort_vcf > pass_vcf  > Allelic_F


rule gzip:
	input:
		path+"/fastq/{name_full}_R{number}.fastq.gz"
	output:
		path+"/fastq/{name_full}_R{number}.fastq"
	shell:
		"gzip -dcf {input} > {output}"

rule trim_galore:
	input:
		r1=path+"/fastq/{name_full}_R1.fastq",
		r2=path+"/fastq/{name_full}_R2.fastq",
		path=path
	output:
		path+"/fastq/{name_full}_R1_val_1.fq",
		path+"/fastq/{name_full}_R2_val_2.fq",
	shell:
		"module load trim-galore/0.5.0 ; module load pypy/2.7-6.0.0 ; trim_galore  --output_dir {input.path}/fastq/  --paired {input.r1} {input.r2}  "

rule bwa_mem:
	input:
		R1=path+"/fastq/{name_full}_R1_val_1.fq",
		R2=path+"/fastq/{name_full}_R2_val_2.fq",
		Genome=expand("{gen}",gen=genome)
	output:
		path+"/alignment/{name_full}.bam"
	shell:
		"module load samtools/default ;\
		module load bwa  ; \
		bwa mem -t {threads} {input.Genome} {input.R1} {input.R2} | samtools view -h -b  > {output} "

rule samtools_sort:
	input:
		path+"/alignment/{name_full}.bam"
	output:
		path+"/alignment/{name_full}_sort.bam"
	shell:
		"module load samtools/default ; samtools sort  -O BAM {input} > {output} "


rule remove_duplicates:
	input:
		path+"/alignment/{name_full}_sort.bam"
	output:
		outbam=path+"/alignment/{name_full}_no_dupl_sort.bam",
		metrics=path+"/alignment/{name_full}_dupl_metrics.txt"
	shell:
		"module load gatk/4.0.9.0 ; gatk MarkDuplicates -I {input}  -O {output.outbam} -M {output.metrics}  --REMOVE_DUPLICATES=true --TMP_DIR {path}"


rule samtools_merge:
	input:
		expand(path+"/alignment/{{name_sample}}-LR-{name_lane}_no_dupl_sort.bam", name_lane = name_lane)
	output:
		path+"/alignment/{name_sample}_merged.bam"
	shell:
		"module load samtools/default ; samtools merge  {output} {input} "

rule samtools_index:
	input:
		path+"/alignment/{name_sample}_merged.bam"
	output:
		path+"/alignment/{name_sample}_merged.bam.bai"

	shell:
		"module load samtools; samtools index {input}"

rule targets:
	input:
		target = TARGET_BED
	params:
		stderr = os.path.join(path, "logs", "target_bed.stderr"),
		stdout = os.path.join(path, "logs", "target_bed.stdout")
	output:
		target_sorted = TARGET_BED_SORTED
	shell:
		"module load bedtools && "
		"file {input.target} | grep gzip && TOOL='zcat' || TOOL='cat' ; "
		"${{TOOL}} {input.target} | sort -V -k1,1 -k2,2 | "
		"bedtools merge -i - > {output.target_sorted} "
rule nontargets:
	input:
		targets=TARGET_BED_SORTED,
		genome=GENOMEFILE
	params:
		stderr = os.path.join(path, "logs", "nontarget_bed.stderr"),
		stdout = os.path.join(path, "logs", "nontarget_bed.stdout")

	output:
		nontargets=NONTARGET_BED
	shell:
		"module load bedtools && "
		"bedtools complement -i {input.targets} "
		"-g {input.genome} > {output.nontargets}"

rule depth_targets:
	input:
		bam = path+"/alignment/{name_sample}_merged.bam",
		bai=path+"/alignment/{name_sample}_merged.bam.bai",
		targets = TARGET_BED_SORTED

	params:
		stderr = os.path.join(path, "logs", "{name_sample}.depth_targets.stderr"),
		stdout = os.path.join(path, "logs", "{name_sample}.depth_targets.stdout")
	output:
		path+"/QC/{name_sample}.target.bed.gz"
	shell:
		"module load sambamba ;"
		"sambamba depth region -L  {input.targets} {input.bam} | gzip -c > {output}"

rule depth_nontargets:
	input:
		bam = path+"/alignment/{name_sample}_merged.bam",
		bai=path+"/alignment/{name_sample}_merged.bam.bai",
		targets = NONTARGET_BED

	params:
		stderr = os.path.join(path, "logs", "{name_sample}.depth_nontargets.stderr"),
		stdout = os.path.join(path, "logs", "{name_sample}.depth_nontargets.stdout")
	output:
		path+"/QC/{name_sample}.nontarget.bed.gz"
	shell:
		"module load sambamba ;"
		"sambamba depth region -L  {input.targets} {input.bam} | gzip -c > {output}"

rule calculateTargetSize:
	input: TARGET_BED_SORTED
	run:
		config = read_config()
		config["annotation"]["wesTargetSize"] =  calcBEDlength(TARGET_BED_SORTED)
		with open(config_infile, "w") as fout:
			json.dump(config, fout)

rule targets_coverageQC:
	input:
		bam=path+"/alignment/{name_sample}_merged.bam",
		bai=path+"/alignment/{name_sample}_merged.bam.bai",
		targets = TARGET_BED_SORTED
	params:
		threads = 4,
		COVERAGEQC_BINARY="/tbi/software/x86_64/otp/roddy/plugins/3.5/AlignmentAndQCWorkflows_1.2.73-1/resources/analysisTools/qcPipelineTools/coverageQcD/coverageQc",
		COVERAGEQC_CHROMSIZES= GENOMEFILE,
		stderr = os.path.join(path, "logs", "{name_sample}.otr.stderr"),
		stdout = os.path.join(path, "logs", "{name_sample}.otr.stdout"),
		TARGETSIZE= read_config()["annotation"]["wesTargetSize"]

	output:
		path+"/QC/{name_sample}.targets_coverageQC.tsv"
	shell:
		"module load bedtools/2.24.0 && "
		"module load samtools/1.5 && "
		"bedtools intersect -abam {input.bam} -b {input.targets} | "
		"samtools view -@{params.threads} -b -q 1 -F 1024 | "
		"{params.COVERAGEQC_BINARY} --alignmentFile /dev/stdin --outputFile={output} --processors={params.threads} --basequalCutoff=1 --targetsize={params.TARGETSIZE} --ungappedSizes={params.COVERAGEQC_CHROMSIZES}"

rule coverageQC:
	input:
		path+"/alignment/{name_sample}_merged.bam"
	params:
		threads = 4,
		COVERAGEQC_BINARY="/tbi/software/x86_64/otp/roddy/plugins/3.5/AlignmentAndQCWorkflows_1.2.73-1/resources/analysisTools/qcPipelineTools/coverageQcD/coverageQc",
		COVERAGEQC_CHROMSIZES= GENOMEFILE,
		stderr = os.path.join(path, "logs", "{name_sample}.otr.stderr"),
		stdout = os.path.join(path, "logs", "{name_sample}.otr.stdout"),
	output:
		path+"/QC/{name_sample}.coverageQC.tsv"
	shell:
		"module load bedtools/2.24.0 && "
		"module load samtools/1.5 && "
		"samtools view -@{params.threads} -b -q 1 -F 1024  {input}| "
		"{params.COVERAGEQC_BINARY} --alignmentFile /dev/stdin --outputFile={output} --processors={params.threads} --basequalCutoff=1 --ungappedSizes={params.COVERAGEQC_CHROMSIZES}"

rule rmdreport:
	input:
		expand(path+"/QC/{name_sample}.{region}.bed.gz", region=["target", "nontarget"], name_sample = name_sample),
		expand(path+"/QC/{name_sample}.targets_coverageQC.tsv", name_sample = name_sample),
		expand(path+"/QC/{name_sample}.coverageQC.tsv", name_sample = name_sample)
	params:
		rscript=path+"/QC/qc-targetcoverage.Rmd",
		stderr = os.path.join(path, "logs", "rmdreport.stderr"),
		stdout = os.path.join(path, "logs", "rmdreport.stdout"),
	output:
		html=path+"/QC/qc-targetcoverage.html", # R requires the full path TODO change to defined project dir
		tsv=path+"/QC/view_values.tsv"
	shell:
		 "module load R; "
		  "R -e 'rmarkdown::render(\"{params.rscript}\", output_file=\"{output.html}\")'; "

rule samtools_subsample:
	input:
		bam = path+"/alignment/{name_sample}_merged.bam",
		view_values = path+"/QC/view_values.tsv"
	output:
		path+"/subsampling/{name_sample}_merged_subsampled.bam"
	shell:
		" module load samtools/1.10; i=$( echo {input.bam} | sed \"s/AS/%AS/\" |  cut -f 2 -d \"%\" | cut -f 1 -d \"_\"  ) ; value=$( grep $i {input.view_values} | cut -f 3 | sort | uniq ) ; samtools view -hbs $value {input.bam} > {output} "

rule add_readgroup_subsampled:
	input:
		bam = path+"/subsampling/{name_sample}_merged_subsampled.bam",
		snakedata = path+"/fastq/snakedata.tsv" #You can find it at /icgc/dkfzlsdf/midterm/THEPROJECTNUMBER/data/RUNNUMBER (/icgc/dkfzlsdf/midterm/018431/data/200520_A0032..)
	output:
		path+"/subsampling/{name_sample}_merged_subsampled_RG_LG.bam"
	shell:
		'module load gatk/4.0.9.0 ; i=$( echo {input.bam} | sed \"s/AS/%AS/\" |  cut -f 2 -d \"%\" | cut -f 2 -d \"_\"  )  ; \
		gatk AddOrReplaceReadGroups  -ID ID  -PL Illumina -LB  ID   -PU ID  -SM $i  -I  {input.bam} -O  {output} '

rule samtools_index_subsampled:
	input:
		path+"/subsampling/{name_sample}_merged_subsampled_RG_LG.bam"
	output:
		path+"/subsampling/{name_sample}_merged_subsampled_RG_LG.bam.bai"

	shell:
		"module load samtools; samtools index {input}"

rule Mutect2:
	input:
		tumor=path+"/subsampling/{tumor_sample}_TUMOR_merged_subsampled_RG_LG.bam",
		normal=path+"/subsampling/"+normal_sample+"_NORMAL_merged_subsampled_RG_LG.bam",
		ind_tumor=path+"/subsampling/{tumor_sample}_TUMOR_merged_subsampled_RG_LG.bam.bai",
		ind_normal=path+"/subsampling/"+normal_sample+"_NORMAL_merged_subsampled_RG_LG.bam.bai",
		ref=genome,
		chr=chrpath+"{chrname2}.bed",
		path=path+"/Mutation_calling"
	output:
		path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output.vcf"
	shell:
		" mkdir -p {input.path}/{wildcards.tumor_sample}/Mutect2; module load gatk/4.0.9.0 ; module load samtools ;  gatk Mutect2 -L {input.chr}  -R  {input.ref}  -I {input.tumor}  -I {input.normal} -tumor TUMOR -normal NORMAL  -O {output} "

rule filter_vcf:
	input:
		path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output.vcf"
	output:
		path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf"
	shell:
		" module load gatk/4.0.9.0 ; gatk  FilterMutectCalls -V {input} -O {output} "

rule index_vcf:
	input:
		path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf"
	output:
		gz=path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf.gz",
		tbi=path+"/Mutation_calling/{tumor_sample}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf.gz.tbi"
	shell:
		"export PERL5LIB=/tbi/software/x86_64/vcftools/vcftools-0.1.12b/el7/lib/perl5/site_perl/   ;\
		export PATH=${{PATH}}:/home/f528r/tabix-0.2.6/ ;\
		bgzip -c  {input}  > {output.gz} ;\
		tabix -p vcf {output.gz}  "

rule merge_vcf:
	input:
		vcf=expand(path+"/Mutation_calling/{{tumor_sample}}/Mutect2/{chrname2}_Mutect2_output_filtered.vcf.gz", chrname2=chrname2),
		path=path+"/Mutation_calling/"
	output:
		path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered.vcf"
	shell:
		" module load vcftools/default ;\
		export PERL5LIB=/tbi/software/x86_64/vcftools/vcftools-0.1.12b/el7/lib/perl5/site_perl/   ;\
		export PATH=${{PATH}}:/home/f528r/tabix-0.2.6/ ;\
		  vcf-merge {input.vcf} > {input.path}/{wildcards.tumor_sample}/Mutect2/merged.vcf ;\
		vcf-sort -c {input.path}/{wildcards.tumor_sample}/Mutect2/merged.vcf  > {output} ;\
		rm {input.path}/{wildcards.tumor_sample}/Mutect2/merged.vcf "

rule sort_vcf:
	input:
		vcf=path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered.vcf",
		dict= dict
	output:
		path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered_sort.vcf"
	shell:
		" module load gatk/4.0.9.0 ;  gatk SortVcf --SEQUENCE_DICTIONARY={input.dict}  --INPUT={input.vcf} --OUTPUT={output} "

rule pass_vcf:
	input:
		path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered_sort.vcf"
	output:
		path+"/Mutation_calling/{tumor_sample}/pass.vcf"
	shell:
		'egrep \"PASS|#CHROM|##\" {input} > {output} '
		
rule kb_interval:
  input:
    vcf=path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered_sort.vcf",
    path=path+"/Mutation_calling",
    bed=kb_bed #I need the 1kb bed of mmus
  output:
    tsv=path+"/Mutation_calling/1kb_intervals/{tumor_sample}.tsv",
    bed=path+"/Mutation_calling/{tumor_sample}/1kb_intervals/Mutect2_PASS.bed"
  shell:
    "module load bedtools; mkdir -p {input.path}/{wildcards.tumor_sample}/1kb_intervals ; grep PASS {input.vcf} | awk \'{{if (length($4)==1 && length($5)==1) print}}\' | awk \'{{ print $1\"\\t\"$2\"\\t\"$2\"\\t\"$4\">\"$5\"\\t.\\t+\"}}\' > {output.bed} ; bedtools intersect -a {input.bed} -b {output.bed} -wa -c -loj > {input.path}/1kb_intervals/1kb_intervals_depth.txt ; cp {input.path}/1kb_intervals/1kb_intervals_depth.txt  {output.tsv}"
    
rule Allelic_F:
  input:
    vcf=path+"/Mutation_calling/{tumor_sample}/Mutect2_output_filtered_sort.vcf",
    path=path+"/Mutation_calling"
  output:
    path+"/Mutation_calling/AF/{tumor_sample}.tsv"
  shell:    
    " mkdir -p {input.path}/{wildcards.tumor_sample}/AF ; grep PASS {input.vcf} | awk \'{{if (length($4)==1 && length($5)==1) print}}\' | sed \'s/\t\.//g\' | cut -f 1,2,3,4,8 | sed \'s/\:/\t/g\' | awk \'{{ print $1\"\t\"$2\"\t\"$3\">\"$4\"\t\"$7 }}\' > {output} "
    
    
