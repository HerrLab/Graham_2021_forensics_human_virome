############################################
## -------------------------------------- ##
## ------------- Description ------------ ##
## -------------------------------------- ##
############################################

#Input Notes:
#This is the fourth step of this pipeline is designed to be used after completing 2_Virome_Assembly.sh
  # You can perform this step either before or after the third step:  3. Virome Contig Mapping (3_Virome_Contig_Mapping.sh)
#This pipeline should begin using the fasta file unmapped_1000_contigs.fa generated using 2_Virome_Assembly.sh

#Output Notes:
#This pipeline will result in multiple taxonomic classification tables for each of the 1298 contigs. These were then put together into one taxonomy file named vc_tax.csv
#The next step in this pipeline is the next step in pipeline is 5. Phyloseq Overall Taxonomy (5_Phyloseq_Overall_Taxonomy.R)

#General Notes:
#This pipeline is designed to be run using the Holland Computing Center at the University of Nebraska. Some tool commands may differ depending on installation of the tool. Please refer to the listed Githubs for each tool used as mentioned in script for further information if issues arise 
#Some file locations may differ from yours so this needs to be changed accordingly. This script is designed to be run all in one folder
#This script is designed to be used post the first and second step of the pipeline, however any viral contig assembly fasta file can be used with this script

############################################
## -------------------------------------- ##
## ---------------- CheckV -------------- ##
## -------------------------------------- ##
############################################

#CheckV v.0.7 is used to identify any contigs that contain at least one viral gene using current viral reference databases
#CheckV v.0.7 can be found at: https://pypi.org/project/checkv/

checkv end_to_end unmapped_1000_contigs.fa ./CHECKV_OUTPUT -t 8

#From the output file of CHECKV_OUTPUT/quality_summary.tsv all contigs that had "no viral genes detected" were removed and all others were retained. 
#This resulted in 1400 contigs containing at least one viral gene. The names of the contigs were put into a file checkv_contig_list.txt

#A fasta file was then made containing only the viral gene containing contigs (as identifed in CheckV) and was named as vir_fin.fa
grep -A 1 -w -f checkv_contig_list.txt unmapped_1000_contigs.fa >> vir_fin.fa

############################################
## -------------------------------------- ##
## --------------- Kraken2 -------------- ##
## -------------------------------------- ##
############################################

#Kraken2 v.2.0.8-beta was used to annotate the checkv identified viral contigs in the vir_fin.fa file
#Kraken2 can be found at: http://ccb.jhu.edu/software/kraken2/index.shtml

mkdir KRAKEN/
kraken2-build --download-taxonomy --db KRAKEN/
kraken2-build --download-library viral --db KRAKEN/
kraken2-build --build --db KRAKEN/
kraken2 --use-names  --db vir_fin.fa > kraken_output.txt

#Will result in a kraken_output.txt file that will be used in conjunction  with the output files of kaiju, demovir, and independent blast searches

############################################
## -------------------------------------- ##
## ---------------- Kaiju --------------- ##
## -------------------------------------- ##
############################################

#Kaiju v.1.7 was used to annotate the Checkv identified viral contigs in the vir_fin.fa file
#Kaiju can be found at: https://kaiju.binf.ku.dk

kaiju -t $NODES -f $KAIJU_DB_VIRUSES -i vir_fin.fa -x -m 11 -a greedy -e 5 -E 0.001 -v -o kaiju_virus.out

kaiju-addTaxonNames -t $NODES -n $NAMES -v  -r superkingdom,phylum,class,order,family,genus,species -i kaiju_virus.out -o kaiju_virus_names.tsv

kaiju -t $NODES -f $KAIJU_DB -i vir_fin.fa -x -m 11 -a greedy -e 5  -E 0.001 -v -o kaiju_DB.out

kaiju-addTaxonNames -t $NODES -n $NAMES -v -r superkingdom,phylum,class,order,family,genus,species -i kaiju_DB.out -o kaiju_DB_names.tsv

#Will result in both kaiju_virus_names.tsv and kaiju_DB_names.tsv files that will be used in conjunction  with the output files of demovir, kraken2, and independent blast searches


############################################
## -------------------------------------- ##
## --------------- Demovir -------------- ##
## -------------------------------------- ##
############################################

#Demovir was used to annotate the Checkv identified viral contigs in the vir_fin.fa file
#Demovir can be found at: https://github.com/feargalr/Demovir

#For those using the HCC at UNL make sure to module load the following to use Demovir:
  #module load usearch/11.0
  #module load prodigal/2.6
  #module load ​R/3.6

mkdir Demovir
#git clone https://github.com/feargalr/Demovir.git
cd Demovir/
chmod +x *.sh
#download database from: https://figshare.com/articles/NR_Viral_TrEMBL/5822166
./format_db.sh

prodigal -a AA.fasta -i vir_fin.fa -p meta &> /dev/null
usearch -ublast AA.fasta -db uniprot_trembl.viral.udb -evalue 1e-5 -blast6out trembl_ublast.viral.txt -trunclabels &> /dev/null
sort -u -k1,1 trembl_ublast.viral.txt > trembl_ublast.viral.u.txt
rm AA.fasta
cut -f 1,2 trembl_ublast.viral.u.txt | sed 's/_[0-9]\+\t/\t/' | cut -f 1 | paste trembl_ublast.viral.u.txt - > trembl_ublast.viral.u.contigID.txt
rm trembl_ublast.viral.u.txt trembl_ublast.viral.txt
Rscript demovir.R
rm trembl_ublast.viral.u.contigID.txt

#Will result in a demovir_assignments.txt file that will be used in conjunction with the output files of kaiju, kraken2, and independent blast searches

############################################
## -------------------------------------- ##
## --------------- BlastN --------------- ##
## -------------------------------------- ##
############################################

#Each contig was run through Blast using the BlastN algorithm. Only hits with >10% coverage were considered true results. The top his was recorded in an excel sheet manually

#All annotation outputs were recorded in an excel file and the best annotation for each contig was determined as explained in the materials in methods. This resulted in the vc_tax.csv file.
#Contigs that had 100% syntony to the human or animal genomes and any proviral genomes were removed from the 1400 set of contigs thus resulting in a final set of 1298 contigs. 
