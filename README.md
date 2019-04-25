# HEh2

This directory includes scripts that were used to produce results for the following manuscript:

Hernandez RD, Uricchio LH, Hartman K, Ye C, Dahl A, Zaitlen N. Ultra-rare variants drive substantial cis-heritability of human gene expression (https://www.biorxiv.org/content/early/2017/12/15/219238).

Unfortunately, much of this code is not plug and play. Several scripts are geared toward a particular data structure with hard coded directories. However, the code is being made available to make it clear exactly what was done to produce the results. Any questions should be addressed to Ryan Hernandez <ryan.hernandez@me.com>. 

Note that we have compared our implementation of Haseman-Elston regression with the implementation in GCTA, and found them to report estimates that are identical within 8 decimal points (likely owing to differences in single vs double point precision). Since GCTA is widely used, we recommend using that package, but make our code available to ensure openness of our research methods.

Data:

WGS Genotype data in VCF format are available from 1000 Genomes Project: http://www.internationalgenome.org/data#download

RNA seq data are available in two ways:
1) Summary FPKM across genes and individuals available in this repo: expressions.genes.eur.matrix.eqtl.txt.gz
2) Directly from GEUVADIS: https://www.ebi.ac.uk/Tools/geuvadis-das/

Subject IDs:
Unrelated individuals and population labels and sex listed in: 20140610_samples_to_analyse_v10.txt.gz


Scripts:

1) Process VCF to 012 format and produce variant position list for each gene can be run in Perl using: getGene012.pl

2) Run H-E on 012/pos/expression data/Covariats: HEh2.R

   To run this code, use command line R:
   Rscript HEh2.R help

3) Functions for running AI implementation of REML: reml_func.R

4) Plot data in R: finalFigures_NatGen4.R