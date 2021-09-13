# Graham_2021_forensics_human_virome

Github Repository for the paper on the use of the human virome as a forensic marker

This work is [deposited to bioRxiv]().

To cite this work or code:
Graham E, et al. To be updated upon acceptance for publication.

_ASSEMBLY_
This directory contains the meta-assembly generated using script 02_Virome_Assembly.sh. It was split into two parts due to its large size in order to upload onto Github. To put the assembly back together into one file first unzip both files and then run the following script:
cat unmapped_1000_contigs1.fa unmapped_1000_contigs2.fa >> unmapped_1000_contigs.fa
