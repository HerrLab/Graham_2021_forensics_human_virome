############################################
## -------------------------------------- ##
## ------------- Description ------------ ##
## -------------------------------------- ##
############################################

#Input Notes:
#This is the third step of this pipeline is designed to be used after completing 2_Virome_Assembly.sh
#This pipeline should begin using the fasta file unmapped_1000_contigs.fa generated using 2_Virome_Assembly.sh

#Output Notes:
#This pipeline will generate a read abundance table that can be used for input into Phyloseq in R (counts2_bowtie.txt)
#The next step in this pipeline is the next step in pipeline is 4. Contig Viral Annotation (4_Contig_Annotation.sh)

#General Notes:
#This pipeline is designed to be run using the Holland Computing Center at the University of Nebraska. Some tool commands may differ depending on installation of the tool. Please refer to the listed Githubs for each tool used as mentioned in script for further information if issues arise 
#Some file locations may differ from yours so this needs to be changed accordingly. This script is designed to be run all in one folder
#This script is designed to be used post the first and second steps of the pipeline, however any contig assembly fasta file can be used with this script

############################################
## -------------------------------------- ##
## --------- Build Contig Index --------- ##
## -------------------------------------- ##
############################################

#The global alignment tool Bowtie2 v.2.3.5 was used for indexing and mapping
#Bowtie2 can be found at: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

bowtie2-build -f unmapped_1000_contigs.fa VContigs_index

############################################
## -------------------------------------- ##
## ----------- Contig Mapping ----------- ##
## -------------------------------------- ##
############################################

#Mapping was performed for each sample (not for negative controls) so repeat this step for each individual sample. Here we use the example of using sample HV_001_01

bowtie2 --end-to-end -x VContigs_index -1 Raw_Reads/HV_001_001_forward.fq.gz -2 Raw_Reads/HV_001_01_reverse.fq.gz | samtools view -F 4 -o HV_001_01.bam

#Samtools was used for subsiquent manipluation of the sam file output from bowtie2
#Samtools can be found at: https://github.com/samtools/samtools

samtools sort -O BAM -o HV_001_01.sorted.bam HV_001_01.bam 
samtools index HV_001_01.sorted.bam
samtools idxstats HV_001_01.sorted.bam > HV_001_01.sorted.bam.idxstats.txt

#################################################
## ------------------------------------------- ##
## - Contig Mapping to Abundance Table for R - ##
## ------------------------------------------- ##
#################################################

#Once all samples have been run through this pipeline you can convert the mapped abundance statistics into a read abundance table that can be used as an OTU table in Phyloseq in R. 
#get_count_table.py python script avalible at: https://github.com/edamame-course/Metagenome

python2 get_count_table.py *.idxstats.txt > counts_bowtie.txt
sed 's/.idxstats.txt//g' counts_bowtie.txt > counts2_bowtie.txt
