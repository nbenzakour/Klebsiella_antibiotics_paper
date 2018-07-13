# Klebsiella_antibiotics_paper

This repository contains the R scripts used to generate some of the supporting information for ourpaper looking at porine defects in *Klebsiella pneumoniae*, as well some additional figures that were not used but which can provide useful additional insights (manuscript under review, soon available on BioRXiv) .

Briefly, we collected a total of 1,557 draft and complete *K. pneumoniae* genomes publicly available in Genbank (Feb 2017) to investigate the significance of mutations affecting major porines OmpK35 and OmpK36 in the context of a wider population. 

data/ contains screening data stored in KP_screening_data.csv, with the following descriptors:

| Identifier        | Description                                                                               |  
|-------------------|-------------------------------------------------------------------------------------------| 
| Id                | Isolate assembly identifier                                                               | 
| genome_size       | Genome size in bp                                                                         | 
| contigs           | Number of contigs                                                                         | 
| n50               | N50 in bp                                                                                 | 
| ST                | Multi-locus sequence type                                                                 | 
| Insertion         | Di-nucleotide insertion of GD or TD in OmpK36 loop L3                                     | 
| ompK36            | Presence or absence (1 or 0)                                                              | 
| ompK35            | Presence or absence (1 or 0)                                                              | 
| Source            | Source of isolates as described in original record                                        | 
| Source_edited     | Source of isolates as described in original record - cleaned up (typo, font case, etc.)   | 
| Source_summarised | Source of isolates summarised in 3 major categories (blood, urine and other)              | 
| Host              | Host origin (about half of the isolates don't have this information)                      | 
| Country           | Country of Origin                                                                         | 
| Year              | Date of isolation simplified as year only                                                 | 
| resistance_genes  | Resistance profile as predicted by ResFinder                                              | 
| plasmids          | Plasmids replicons identified by PlasmidFinder (not used in present analysis)             | 
| Organism          | Organism name                                                                             | 
| Strain            | Strain name                                                                               | 

scripts/ contains the following scripts:

* KP_associations_EDA.md
* KP_population_EDA.md
* KP_resistance_determinants_EDA.md

images/ contains images generated:

* Fig S7 KP_unnamed-chunk-2-1.png
* Fig S8 KPassoc_unnamed-chunk-8-1.png
* Fig S9 KPassoc_unnamed-chunk-9-1.png & KPassoc_unnamed-chunk-9-2.png
* Fig S10 KPassoc_unnamed-chunk-10-1.png & KPassoc_unnamed-chunk-10-2.png
* Fig S11 KPres_unnamed-chunk-7-1.png