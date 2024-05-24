## GOAL: create an annotation pipeline for GWAS loci, by editing/ upgrading the ProGeM pipline.
- botom up: add coloc with eQTL and pQTLs
- top dow: add PoPS 


**GWASAnno_main.R parameters:**  
| Column | Description |
| --- | --- |
| --GWAS_RDS | GWAS summary stats .RSD [required] |
| --output_path | Output folder path [required] |
| --bfile | bfile to use for selecting proxies and for MAGMA (if not provided UKB_14K_hg38 will be used) |
| --magma_annotated_genes | Magma_annotated_genes to use for MAGMA (if not provided it will be created)
| --eQTL_datasets_coloc | Comma-separated eQTL datasets to use for coloc (default=GTEXv8)
| --eQTL_tissues_interest_coloc | Comma-separated tissues of interest to be selected from the eQTL datasets -check 04_utils/tissues_eQTL.txt
| --pQTL_datasets_coloc | Comma-separated pQTL datasets to use for coloc (default=c('Icelanders_pGWAS','UKB_PPP_EUR')
| --coloc_input_path | Path with coloc results for the selected databases (default=NA and runs coloc analysis)
| --r2_thresh | Threshold for r2 values (default=0.8)
| --interval_kb | Genomic interval (in kb) either side of the sentinel SNP (default=500), used to define genes that nearby the index SNP 
| --LD_region_overhang_kb | Genomic interval (in kb) of the overhang of the left-most and right-most proxy/index SNP at each locus (default=5), used to identify genes with variants in LD with index or proxy SNPs
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
| A1 | alternative/ effect allele (example: *AC*) |
| A2 | reference/ other allele (example: *A*) |
| BETA | beta |
| SE | se|
| nlog10P | -log10(P) |
| AF | frequency of A1 |
| N | sample size |

* scripts 04_utils/00_process_sumstats_REGENIE.R and 04_utils/00_process_sumstats_metal.R can be used to create an input file from a gz and tabix regenie output
* analysis will be conducted in hg38!!! 

**--Post processing - scoring system:**  
*script 00_scripts/08_postprocessing_ProGeM.R can be used with the following options in case an alternative scoring is desired  

| Column | Description |
| --- | --- |
| --output_path | output_path for GWASAnno full path [required] |
| --output_file_name |  output file name [default='GWASAnno_summary.txt]' |
| --eQTL_tissues_interest_coloc | Comma-separated tissues of interest previously selected from the eQTL datasets |
| --nearest | nearest score [default=1] |
| --second_nearest | second_nearest score [default=0.99] |
| --third_nearest | third_nearest score [default=0.98] |
| --LD_overlapping | LD_overlapping score [default=1] |
| --lead_IMPACT | VEP lead IMPACT score [default=1.02] |
| --coloc_eQTL_tissues_interest | coloc_eQTL_tissues_interest score (if --eQTL_tissues_interest_coloc are not defined, --coloc_eQTL_tissues_interest is not used) [default=1.04]|
| --coloc_eQTL | coloc_eQTL score (if --eQTL_tissues_interest_coloc is defined, this score is only applied if NO coloc_eQTL_tissues_interest is identified) [default=1.03]|
| --lead_eQTL | lead_eQTL or proxie_eQTL score (this score is only applied if NO coloc_eQTL is identified) [default=0.75]|
| --coloc_pQTL | coloc_pQTL score [default=1.04] |
| --PoPS_top1 | PoPS_top1 score [default=0] |
| --PoPS_top2 | PoPS_top2 score [default=0] |
| --PoPS_top3 | PoPS_top3 score [default=0] |
