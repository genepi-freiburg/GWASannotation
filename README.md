### 31Jan2024
## GOAL: create an annotation pipeline for GWAS loci, by editing/ upgrading the ProGeM pipline.
- botom up: add coloc with eQTL and pQTLs
- top dow: add PoPS 


**GWASAnno_main.R parameters:**  
| Column | Description |
| --- | --- |
| --GWAS_RDS | GWAS summary stats .RSD [required] |
| --output_path | Output folder path [required] |
| --bfile | Bfile to use for MAGMA (if not provided UKB_14K_hg38 will be used) |
| --magma_annotated_genes | Magma_annotated_genes to use for MAGMA (if not provided it will be created)
| --eQTL_datasets_coloc | Comma-separated eQTL datasets to use for coloc (default=GTEXv8)
| --eQTL_tissues_interest_coloc | Comma-separated tissues of interest to be selected from the eQTL datasets -check 04_utils/tissues_eQTL.txt
| --pQTL_datasets_coloc | Comma-separated pQTL datasets to use for coloc (default=c('Icelanders_pGWAS','UKB_PPP_EUR')
| --coloc_input_path | Path with coloc results for the selected databases (default=NA and runs coloc analysis)
| --r2_thresh | Threshold for r2 values (default=0.8)
| --interval_kb | Genomic interval (in kb) either side of the sentinel SNP (default=500)
| --LD_region_overhang_kb | Genomic interval (in kb) of the overhang of the left-most and right-most proxy/sentinel variant at each locus (default=5)
| --sumstats_2_max_nlog10P_thresh | Sumstats_2_max_nlog10P will be filtered for > this threshold (default=-log10(5e-8))
| --PP.H4.abf_thresh | PP.H4.abf will be filtered for > this threshold (default=0.8)


**--GWAS_RDS format:** 
- Sumstats RDS file built using genome build hg38, with the following columns:

| Column | Description |
| -- | --- |
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

* script 00_process_sumstats.R can be used to create an input file from a gz and tabix regenie output
* analysis will be conducted in hg38!!! 

