############################################
## -------------------------------------- ##
## ------------- Description ------------ ##
## -------------------------------------- ##
############################################

#Input Notes:
#This is the second step of this pipeline and is starting using fastq data files that were organized into folders using the 1_Obtain_Virome_Data.sh
#All fastq data files that correspond to samples with index ending in 23 and 25 as indicated in metadata file are control reads and should be put in a directory ./Raw_Control_Reads
 #Within the Raw_Control_Reads directory each fastq forward and corresponding reverse read was placed in a folder labed with its SampleID
 #Sample file names were changed to be "SampleID"_forward.fq (e.g. for sample HV_001_25 reads it would be ./Raw_Control_Reads/HV_001_25_forward.fq and ./Raw_Control_Reads/HV_001_25_reverse.fq
#All other fastq data files that are not control reads were put in a different directory ./Raw_Reads
 #Within the Raw_Reads directory each fastq forward and corresponding reverse read was placed in a folder labed with its SampleID
 #Sample file names were changed to be "SampleID"_forward.fq (e.g. for sample HV_001_01 reads it would be ./Raw_Reads/HV_001_01_forward.fq and ./Raw_Reads/HV_001_01_reverse.fq

#Output Notes:
#This pipeline will generate a fasta file containing virome contigs (unmapped_1000_contigs.fa)
#The next step in this pipeline is The next step in pipeline is 3. Contig Mapping (3_Virome_Contig_Mapping.sh)

#General Notes:
#This pipeline is designed to be run using the Holland Computing Center at the University of Nebraska. Some tool comands may differ depending on installation of the tool. Please refer to the listed githubs for each tool used as mentioned in script for further information if issues arise 
#Some file locations may differ from yours so this needs to be changed accordingly. This script is designed to be run all in one folder
#This script is designed to be used for viral metagenome data but can be altered to be used for any type of metagenome data by changing host and contamination removal steps

#######################################################################
## ----------------------------------------------------------------- ##
## -- Quality Assessment of raw sequenceing reads in fastq format -- ##
## ----------------------------------------------------------------- ##
#######################################################################

#This pipeline should be preformed on all samples in Raw_Reads directory generated using the 1_Obtain_Virome_Data.sh script. Here we will use one sample (HV_001_01) as an example. Please repeat each step for all samples.  

## --------------------------------------
## --------- Quality Evaluation ---------
## --------------------------------------

#Fastq data files were initially quality checked using FastQC 
#FastQC is available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

fastqc Raw_Reads/HV_001_01_forward.fq
fastqc Raw_Reads/HV_001_01_reverse.fq

## ---------------------------------------
## ------------ Phi X removal ------------
## ---------------------------------------

#Phi X contaminating reads were removed from fastq data files using BBduk from the BBmap package
#BBmap can be found at: https://github.com/BioInfoTools/BBMap
#BBmap guide is avalible at https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/ 

mkdir Phi_Remove_Trimmed/

bbduk.sh in=Raw_Reads/HV_001_01_forward.fq out=Phi_Remove/HV_001_01_filt_forward.fq k=31 ref=artifacts,phix ordered cardinality
bbduk.sh in=Raw_Reads/HV_001_01_reverse.fq out=Phi_Remove/HV_001_01_filt_reverse.fq k=31 ref=artifacts,phix ordered cardinality

## --------------------------------------
## ---------- Quality trimming ----------
## --------------------------------------

#Phix removed data files were trimmed to remove low quality reads and short reads and return matching paired end reads using the Sickle
#Sickle is avalible at https://github.com/najoshi/sickle

mkdir Phi_Remove/

sickle pe -t sanger -f Phi_Remove/HV_001_01_filt_forward.fq -r Phi_Remove/HV_001_01_filt_reverse.fq -o Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R1.fastq -p Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R2.fastq -s Phi_Remove_Trimmed/HV_001_01_S4_trimmed_single.fastq -q 30 -l 75

## -----------------------------------------------
## - Bacterial 16S rDNA contamination evaluation -
## -----------------------------------------------

#16S reads were obtained from the SILVA database v.138.1 to be used as a reference

mkdir 16S_Ref/
mkdir 16S_Mapped/

wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz ./16S/

#BBmap was used to make a reference using the 16S rDNA reads obtained 

bbmap.sh ref=./16S_Ref/SILVA_138.1_SSURef_tax_silva.fasta.gz -Xmx23g

#Trimmed Data files were mapped to the 16S rDNA reference using BBMap and flags and parameters described for high precision mapping of contamination detection as suggested in the BBMap Guide; percentages of mapped reads were determined

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=./16S_Ref/ qtrim=rl trimq=10 untrim -Xmx23g in=Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R1.fastq outm=16S_Mapped/HV_001_01_16S_R1.fq
bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=./16S_Ref/ qtrim=rl trimq=10 untrim -Xmx23g in=Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R2.fastq outm=16S_Mapped/HV_001_01_16S_R2.fq


## ----------------------------------------
## - Removal of human contaminating reads -
## ----------------------------------------

#The hg19 human genome was retreived and indexed using BBmap

mkdir Human_Ref/
mkdir Human_Remove/

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz ./Human_Ref/
tar zvfx ./Human_Ref/chromFa.tar.gz
cat ./Human_Ref/*.fa > ./Human_Ref/hg19.fa
rm ./Human_Ref/chr*.fa

bbmap.sh ref=./Human_Ref/hg19.fa -Xmx23g

#Trimmed Data files were mapped to the hg19 human indexed genome using with standard operational flags for high precision mapping with low sensitivity in order to lower the risk of false positive mapping as suggested in the BBMap Guide

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=./Human_Ref qtrim=rl trimq=10 untrim -Xmx23g in=Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R1.fastq outu=Human_Remove/HV_001_01_clean_R1.fq outm=Human_Remove/HV_001_01_human_R1.fq
bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=./Human_Ref qtrim=rl trimq=10 untrim -Xmx23g in=Phi_Remove_Trimmed/HV_001_01_S4_trimmed_pair_R2.fastq outu=Human_Remove/HV_001_01_clean_R2.fq outm=Human_Remove/HV_001_01_human_R2.fq

## --------------------------------------
## --------- Quality Trimming 2 ---------
## --------------------------------------

# Trimming was prefored on human contamination removed reads to make paired end reads the same length

mkdir Human_Remove_Trimmed/ 

sickle pe -t sanger -f Human_Remove/HV_001_01_clean.R1.fq -r Human_Remove/HV_001_01_clean.R2.fq -o Human_Remove_Trimmed/HV_001_01_final_R1.fastq -p Human_Remove_Trimmed/HV_001_01_final_R2.fastq -s Human_Remove_Trimmed/HV_001_01_final_single.fastq

###########################################
## ------------------------------------- ##
## ---- Denovo Virome Meta-Assembly ---- ##
## ------------------------------------- ##
###########################################

## --------------------------------------
## -------------- Assembly --------------
## --------------------------------------

#A data set was formed by combining all sample (not negative control) raw filtered reads and running an assembly on those reads to make one meta-assembly:

#combine all forward reads and reverse reads for an overall virome assembly (for all non-control samples)

cat Human_Remove_Trimmed/*_final_R1.fastq > forward_overall.fastq
cat Human_Remove_Trimmed/*_final_R2.fastq > reverse_overall.fastq

#Megahit v.1.2.8 was used for denovo assembly and can be found at https://github.com/voutcn/megahit

megahit -1 forward_overall.fastq -2 reverse_overall.fastq -o MEGAHIT_OUTPUT_ALL -t 36 --presets meta-large

## --------------------------------------
## --- Quality evaluation of assembly ---
## --------------------------------------

#QUAST v.5.0.2 was run on the resulting assembly output to assess the meta assembly quality
#QUAST can be found at: https://github.com/ablab/quast

quast.py ./MEGAHIT_OUTPUT_ALL/final_contigs.fa -o QUAST_OUTPUT

##################################################
## -------------------------------------------- ##
## - Generate Greater Than >1kb Meta-Assembly - ##
## -------------------------------------------- ##
##################################################


#All contigs smaller than 1000 bp in legnth were removed from the meta-assembly
#removesmalls.pl perl script can be found at: https://github.com/drtamermansour/p_asteroides/blob/master/scripts/removesmalls.pl

cp MEGAHIT_OUTPUT_ALL/final_contigs.fa > MEGAHIT_OUTPUT_ALL/ALL_CONTIGS.fa 
perl removesmalls.pl 1000 MEGAHIT_OUTPUT_ALL/ALL_CONTIGS.fa  > ./ALL_1000_CONTIGS.fa

#######################################################
## ------------------------------------------------- ##
## - Denovo Negative Control Samples Meta-Assembly - ##
## ------------------------------------------------- ##
#######################################################

## --------------------------------------
## -------------- Assembly --------------
## --------------------------------------

#An additional data set was made that contained only negative control reads that can be used to remove any contaminating reads from the overall assembly

#cp all control raw reads (all samples with index ending in 23 and 25 as indicated in metadata file) into a folder and combine all forward and reverse reads

mkdir Controls/

cat Raw_Control_Reads/*forward.fq > ./Controls/control_forward.fq
cat Raw_Control_Reads/*reverse.fq > ./Controls/control_reverse.fq

megahit -1 ./Controls/control_forward.fq -2 ./Controls/control_reverse.fq.gz  -o CON_MEGAHIT_OUTPUT -t 36 --presets meta-large

## --------------------------------------
## --- Quality evaluation of assembly ---
## --------------------------------------

cp CON_MEGAHIT_OUTPUT/final_contigs.fa > CON_MEGAHIT_OUTPUT/CON_CONTIGS.fa 
perl removesmalls.pl 1000 CON_MEGAHIT_OUTPUT/CON_CONTIGS.fa > ./CON_1000_CONTIGS.fa
 
#############################################
## --------------------------------------- ##
## - Removal of Negative Control Contigs - ##
## --------------------------------------- ##
#############################################

#Negative control contigs were used as a refernce and any contigs from the overall assembly that mapped to control reads were removed 

#BWA v.0.7 was used for mapping of overall contigs to the negative contigs
#BWA can be found at: https://github.com/lh3/bwa

bwa index CON_1000_CONTIGS.fa
bwa mem CON_1000_CONTIGS.fa ALL_1000_CONTIGS.fa > ALL_1000_CONTIGS.sam

#Samtools was used to manipulate the sam file and pull out all unmapped contigs to be used for downstream processing
#Samtools can be found at: https://github.com/samtools/samtools

samtools faidx CON_1000_CONTIGS.fa
samtools import CON_1000_CONTIGS.fa.fai ALL_1000_CONTIGS.sam ALL_1000_CONTIGS.bam
samtools sort -O BAM -o ALL_1000_CONTIGS.sorted.bam ALL_1000_CONTIGS.bam
samtools index ALL_1000_CONTIGS.sorted.bam
samtools idxstats ALL_1000_CONTIGS.sorted.bam > ALL_1000_CONTIGS.sorted.idxstats.txt

#Pull out mapped reads

samtools view -b -F 4 -b ALL_1000_CONTIGS.sorted.bam > mapped_1000_contigs.bam
samtools fasta mapped_1000_contigs.bam > mapped_1000_contigs.fa
samtools view -c -F 4 ALL_1000_CONTIGS.sorted.bam >> num_mapped_reads_all_1000.txt

#Pull out ummapped reads

samtools view -b -f 4 -b ALL_1000_CONTIGS.sorted.bam > unmapped_1000_contigs.bam
samtools fasta unmapped_1000_contigs.bam > unmapped_1000_contigs.fa
samtools view -c -f 4 ALL_1000_CONTIGS.sorted.bam >> num_unmapped_reads_all_1000.txt

#### The output of unmapped_1000_contigs.fa will be used for downstream processing ####
#### The next step in pipeline is 3. Contig Mapping (3_Virome_Contig_Mapping.sh) ####

