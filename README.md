## GOAL: prioritize the most probable causal genes for each independent GWAS association signal
This annotation pipeline integrates multiple lines of evidence, building upon and enhancing the ProGeM framework (Stacey et al. 2019). Our approach incorporates gene proximity, functional variant annotations using Ensembl 101 Variant Effect Predictor (VEP) (McLaren et al. 2016), genetic colocalization with expression and protein quantitative trait loci (eQTL and pQTL), and linkage disequilibrium (LD)-based gene overlap, to generate a comprehensive scoring system for gene prioritization.  

**Further information can be found in this publication:**  
Monteiro-Martins S, et al. Genetic screens of imaging-derived kidney volumes identify genes linked to kidney function. Kidney Int. 2025 Oct 10;S0085-2538(25)00781-1.  
[PMID: 41077127](https://www.kidney-international.org/article/S0085-2538(25)00781-1/fulltext)  


**R packages used**
- GenomicAlignments
- GenomicFeatures
- GenomicRanges
- parallel
- coloc
- data.table
- readxl
- writexl
- dplyr
- optparse
- genepicoloc (currently only available internally).  
_*R version 4.6.0*_

**VEP version 101**.  
_*INTERNAL: VEP requires DBI! Please confirm you have DBI installed (you can check by typing perl -e 'use DBI')._

**Analysis will be conducted in hg38!** 

### GWASAnno_main.R parameters: 
| Column | Description |
| --- | --- |
| --GWAS_RDS | GWAS summary stats (.RSD) [required] |
| --output_path | Complete directory path, including output file prefix [required] |
| --sumstats_type | GWAS summary type: quant or cc (default=quant) |
| --GWAS_max_nlog10P_thresh |GWAS summary statistics will be filtered for > this threshold  (default=-log10(5e-8)) |
| --bfile | Plink bfile (bed, bim, fam) to use as LDref for selecting proxies (if not provided UKB_14K_hg38 will be used) |
| --r2_thresh | Threshold for r2 values (default=0.8)
| --eQTL_datasets_coloc | Comma-separated eQTL datasets to use for coloc. Possible options: GTEXv8, eQTLGen, Kidney_eQTL (default=GTEXv8)
| --eQTL_tissues_interest_coloc | Comma-separated tissues of interest to be selected from the eQTL datasets -check 04_utils/tissues_eQTL.txt (default=NA)
| --pQTL_datasets_coloc | Comma-separated pQTL datasets to use for coloc. Possible datasets: Icelanders_pGWAS,UKB_PPP_EUR (default=c('Icelanders_pGWAS','UKB_PPP_EUR')
| --coloc_input_path | Path with coloc results for the selected datasets (default=NA and runs coloc analysis)
| --QTL_coloc_max_nlog10P_thresh | eQTL and pQTL datasets used for colocalization will be filtered for > this threshold  (default=-log10(5e-8))
| --PP.H4.abf_thresh | PP.H4.abf will be filtered for > this threshold (default=0.8)
| --interval_kb | Genomic interval (in kb) either side of the sentinel SNP, used to define genes nearby the index SNP (nearest genes) (default=500)
| --LD_region_overhang_kb | Genomic interval (in kb) of the overhang of the left-most and right-most proxy/index SNP at each locus, used to identify genes with variants in LD with index or proxy SNPs (default=5)

  
  
### --GWAS_RDS format:
- Sumstats RDS file built using genome build hg38, with the following columns:

| Column | Description |
| -- | --- |
| Name | CHR:POS:REF:ALT (fixed format to match the variant identifier in coloc datasets; example: *chr1:10177:A:AC*) |
| rsID |  Variant identifier that will be shown on the output (need to be in concordance with the LDref Plink bfile; example *rs367896724*) |
| CHR | chromosome number (example: *1*) |
| POS | SNP position (example: *10177*) |
| A1 | alternative/ effect allele (example: *AC*) |
| A2 | reference/ other allele (example: *A*) |
| BETA | beta |
| SE | se|
| nlog10P | -log10(P) |
| AF | frequency of A1 |
| N | sample size |

_*INTERNAL: some scripts can be found on 05_internal/ that can be used to help create the input file for GWASannotation pipeline. However, these step is not part of the pipeline! Make sure your input file for --GWAS_RDS fits the description_  
_**Before running the pipeline check if your GWAS$Name overlaps the Names from the colocalization datasets**_
<pre> 
  c <- data.table::fread("/data/public_resources/CKDGen/preprocessing/Wuttke2019/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt_hg38.gz")
  sumstats <- readRDS("/path/to/sumstats_hg38_dedup.RDS")
  length(unique(sumstats$Name))
  #Intersection between sumstats and coloc
  int <- sumstats[which(sumstats$Name %in% c$Name_hg38),]
  length(unique(int$Name))

  #check with inverted allele order
  split_data <- do.call(rbind, strsplit(sumstats$rsID, ":"))
  sumstats$A1 <- split_data[, 3]
  sumstats$A2 <- split_data[, 4]
  sumstats$Name_inv <- paste0("chr",sumstats$CHR,":", sumstats$POS,":", sumstats$A2,":", sumstats$A1)
  int_inv <- sumstats[which(sumstats$Name_inv %in% c$Name_hg38),]
  length(unique(int_inv$Name))
</pre>

### Post processing - scoring system:
*The script 00_scripts/06_postprocessing.R can be used with the following options in case an alternative scoring is desired  

| Column | Description |
| --- | --- |
| --output_path | output_path for GWASAnnotation full path (same as previously) [required] |
| --output_file_name |  output file name [default='GWASAnno_summary.txt]' |
| --eQTL_tissues_interest_coloc | Comma-separated tissues of interest selected from the eQTL datasets (this selection can be altered as long as the tissues selected are part of the eQTL datasets used on GWASAnno_main.R)
| --nearest | nearest score [default=1] |
| --second_nearest | second_nearest score [default=0.99] |
| --third_nearest | third_nearest score [default=0.98] |
| --LD_overlapping | LD_overlapping score [default=1] |
| --lead_IMPACT | VEP lead IMPACT score [default=1.02] |
| --coloc_eQTL_tissues_interest | coloc_eQTL_tissues_interest score (if --eQTL_tissues_interest_coloc are not defined, --coloc_eQTL_tissues_interest is not used) [default=1.04]|
| --coloc_eQTL | coloc_eQTL score (if --eQTL_tissues_interest_coloc is defined, this score is only applied if NO coloc_eQTL_tissues_interest is identified) [default=1.03]|
| --lead_eQTL | lead_eQTL or proxie_eQTL score (this score is only applied if NO coloc_eQTL is identified) [default=0.75]|
| --coloc_pQTL | coloc_pQTL score [default=1.04] |

### Examples  
* Run GWASAnnotation without --eQTL_tissues_interest_coloc
<pre> 
output_path="/path/to/output/test"
input_file="/path/to/sumstats_hg38_dedup.RDS"
bfile="/path/to/plink_b_file"
sbatch -p TempCompute --nodelist="imbi12" --output /path/to/logs/GWASAnno_test.txt --job-name=GWASAnno --mail-type=END --mail-user=email@domain.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R  \
--GWAS_RDS $input_file \
--bfile $bfile \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--output_path $output_path "
</pre>

* Run 06_postprocessing selecting --eQTL_tissues_interest_coloc
<pre> 
output_path="/path/to/output/test"
sbatch -p TempCompute --nodelist="imbi12" --output /path/to/logs/GWASAnno_test.txt  --job-name=Anno --mail-type=END --mail-user=email@domain.de --wrap="Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/06_postprocessing.R  \
    --output_path $output_path \
    --output_file_name "GWASAnno_summary_v2.txt" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Liver,Kidney_eQTL.TubsigeQTLs" "
</pre>

* Run GWASAnnotation using previously generated colocalization results.
<pre> 
output_path="/path/to/output/test"
input_file="/path/to/sumstats_hg38_dedup.RDS"
bfile="/path/to/plink_b_file"
sbatch -p TempCompute --nodelist="imbi12" --output /path/to/logs/GWASAnno_test.txt --job-name=GWASAnno --mail-type=END --mail-user=email@domain.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R  \
--GWAS_RDS $input_file \
--bfile $bfile \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--coloc_input_path "/path/to/coloc/output/" \
--output_path $output_path "
</pre>

