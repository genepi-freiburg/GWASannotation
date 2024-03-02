#31Jan2024
GOAL: create an annotation pipeline for GWAS loci, by editing/ upgrading the ProGeM pipline.
- botom up: add coloc with eQTL and pQTLs
- top dow: add PoPS 
optparse was added to scripts/ProGeM_settings.r to facilitate using the pipeline in parallel for different tratis.

<<<<<<< HEAD
INPUT:
- sumstats RDS file in hg38 build with the following columns 

| Column | Description |
| --- | --- |
| Name | CHR:POS:A2:A1 (example: *chr1:10177:A:AC*) |
| rsID |  Reference SNP cluster ID (example *rs367896724*) |
| CHR | chromosome number (example: *1*) |
| POS | SNP position (example: *10177*) |
| A1 | effect allele (example: *AC*) |
| A2 | other allele (example: *A*) |
| BETA | beta |
| SE | se|
| nlog10P | -log10(P) |
| AF | frequency of A1 |
| N | sample size |
*script 00_process_sumstats.R can be used to create an input file from a gz and tabix regenie output
* analysis will be conducted in hg38!!! 

- output path: full output path (including the basename of the output files)
