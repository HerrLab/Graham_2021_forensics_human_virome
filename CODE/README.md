# Graham_2021_forensics_human_virome

Github Repository for the paper on the use of the human virome as a forensic marker

This work is [deposited to bioRxiv]().

To cite this work or code:
Graham E, et al. To be updated upon acceptance for publication.

__ABSTRACT__
The use of skin virome for human identification purposes offers a unique approach to instances where a viable and statistically relevant human DNA profile is unavailable. The skin virome may act as an alternative DNA profile and/or an additional form of probative genetic material. To date, no study has attempted to investigate the human virome over a time series across various physical locations of the body to identify its potential as a tool for human identification. For this study, we set out to evaluate the stability, diversity, and individualization of the human skin virome. An additional goal was to identify viral signatures that can be used in conjunction with traditional forensic STR loci. In order to accomplish this, human virome metagenomes were collected and sequenced from 42 individuals at three anatomical locations (left hand, right hand, and scalp) across multiple collection periods over a 6-month window of time. Assembly dependent and independent bioinformatic approaches were employed, along with a database-based assessment, which resulted in three sets of stable putative viral markers. In total, with the three sets combined, 59 viral species and uncharacterized viral genome assemblies were identified as being significantly stable (P=5.3x10^-15). Viral diversity, based on presence or absence, is significantly different across subjects (P<0.001). Here we demonstrate that not only is the human virome applicable to be used for human identification, but we have identified many viral signatures that can be used for forensic applications, thus providing a foundation to the novel field of forensic virology. 

__FUNDING__
This work was completed using the [Holland Computing Center](https://hcc.unl.edu/) of the University of Nebraska, which receives support from the Nebraska Research Initiative. This work was supported by the Department of Justice [grant numbers 2017-IJ-CX-0025 and 2019-R2-CX-0048]. All of the funding agencies had no role in study design, data collection and interpretation, or the decision to submit the work for publication.

More info:
[Herr Lab Website](http://herrlab.com/);
[Fernando Lab Website](https://fernandolab.unl.edu/)

 ## CODE DIRECTORY ##

This directory contains all scripts used in this study. At the beginning of each script is a description of what files are needed for input to execute the script, a brief general description of the script, and a description of all expected outputs. For any issues or questions or issues experienced in scripts please contact Ema Graham (ema.graham@huskers.unl.edu). All R-markdown scripts were designed to be run using R v.3.6.3. 

__ORDER OF SCRIPTS__

1. Obtain Virome Data (01_Obtain_Virome_Data.sh)
2. Virome Assembly (02_Virome_Assembly.sh)
3. Virome Contig Mapping (03_Virome_Contig_Mapping.sh)
4. Contig Annotation (04_Contig_Annotation.sh)
5. Phyloseq Overall Taxonomy (05_Phyloseq_Overall_Taxonomy_Fig1.Rmd)
6. Relative Abundance by Location within Subject (06_Relative_Abundance_Fig2.Rmd)
7. Stable Viral Families (07_Stable_Families_Fig3.Rmd)
8. Generation of Set A Stable Markers (08_Set_A.Rmd)
9. Mapping of Reads to NCBI Reference Genomes for Set B (09_NCBI_Genomes_Mapping_SetB.sh)
10. Generation of Set B Stable Markers (10_Set_B.Rmd)
11. Generation of Set C Stable Markers (11_Set_C.Rmd)
12. Number of Markers Identified as being Stable by Location and Set (12_Upset_Plots_Fig4.Rmd)
13. Heatmap Profile using all Identified Stable Markers (13_Heatmap_Profile_Fig5.Rmd)
14. Stability Statistics of Marker Sets (14_Within_vs_Between_Dissimilarity_Fig6.Rmd)
15. Alpha and Beta Diversity of the Overall Marker Set 15_Alpha_Beta_Diversity_Fig7.Rmd

__SCRIPTS USED TO MAKE MANUSCRIPT FIGURES__

- Figure 1: 05_Phyloseq_Overall_Taxonomy_Fig1.Rmd
- Figure 2: 06_Relative_Abundance_Fig2.Rmd
- Figure 3: 07_Stable_Families_Fig3.Rmd
- Figure 4: 12_Upset_Plots_Fig4.Rmd
- Figure 5: 13_Heatmap_Profile_Fig5.Rmd
- Figure 6: 14_Within_vs_Between_Dissimilarity_Fig6.Rmd
- Figure 7: 15_Alpha_Beta_Diversity_Fig7.Rmd
- Supplementary Figure 2: 05_Phyloseq_Overall_Taxonomy_Fig1.Rmd
- Supplementary Figure 3: 05_Phyloseq_Overall_Taxonomy_Fig1.Rmd
