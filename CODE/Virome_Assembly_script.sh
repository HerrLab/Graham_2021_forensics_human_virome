
#################################################################
##                                                             ##
## Quality Assessment of raw sequenceing reads in fastq format ##
##                                                             ##
#################################################################

## Quality Evaluation ##

#Fastq data files were initially quality checked using FastQC 
#FastQC is available at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

fastqc forward.fq
fastqc reverse.fq

## Phi X removal ##

#Phi X contaminating reads were removed from fastq data files using BBduk from the BBmap package
#BBmap guide is avalible at https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/ 

bbduk.sh in=forward.fq out=filt_forward.fq k=31 ref=artifacts,phix ordered cardinality
bbduk.sh in=reverse.fq out=filt_reverse.fq k=31 ref=artifacts,phix ordered cardinality

## Quality trimming ##

#Phix removed data files were trimmed to remove low quality reads and short reads and return matching paired end reads using the Sickle
#Sickle is avalible at https://github.com/najoshi/sickle

sickle pe -t sanger -f filt_forward.fq -r filt_reverse.fq -o S4_trimmed_pair_R1.fastq -p S4_trimmed_pair_R2.fastq -s S4_trimmed_single.fastq -q 30 -l 75

## Bacterial 16S rDNA contamination evaluation ##

#16S reads were obtained from the SILVA database v.138.1 to be used as a reference

wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz $WORK/16S/

#BBmap was used to make a reference using the 16S rDNA reads obtained 

bbmap.sh ref=$WORK/16S/SILVA_138.1_SSURef_tax_silva.fasta.gz -Xmx23g

#Trimmed Data files were mapped to the 16S rDNA reference using BBMap and flags and parameters described for high precision mapping of contamination detection as suggested in the BBMap Guide; percentages of mapped reads were determined

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=$WORK/16S/ qtrim=rl trimq=10 untrim -Xmx23g in=S4_trimmed_pair_R1.fastq outm=16S.R1.fq
bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=$WORK/16S/ qtrim=rl trimq=10 untrim -Xmx23g in=S4_trimmed_pair_R2.fastq outm=16S.R2.fq

## Removal of human contaminating reads ##

#The hg19 human genome was retreived and indexed using BBmap

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz $WORK/Human/
tar zvfx $WORK/Human/chromFa.tar.gz
cat $WORK/Human/*.fa > $WORK/Human/hg19.fa
rm $WORK/Human/chr*.fa
bbmap.sh ref=$WORK/Human/hg19.fa -Xmx23g

#Trimmed Data files were mapped to the hg19 human indexed genome using with standard operational flags for high precision mapping with low sensitivity in order to lower the risk of false positive mapping as suggested in the BBMap Guide

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$WORK/Human qtrim=rl trimq=10 untrim -Xmx23g in=S4_trimmed_pair_R1.fastq outu=clean.R1.fq outm=human.R1.fq
bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$WORK/Human qtrim=rl trimq=10 untrim -Xmx23g in=S4_trimmed_pair_R2.fastq outu=clean.R2.fq outm=human.R2.fq

## Trimming to make paired end reads the same length ##

sickle pe -t sanger -f clean.R1.fq -r clean.R2.fq -o final_R1.fastq -p final_R2.fastq -s final_single.fastq

##############
##          ##
## Assembly ##
##          ##
##############
