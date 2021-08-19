# SuSiEx

**SuSiEx** is a Python based command line tool that performs cross-ethnic fine-mapping using GWAS summary statistics and LD reference panels. The method is built on the Sum of Single Effects (SuSiE) model:

Gao Wang, Abhishek Sarkar, Peter Carbonetto, Matthew Stephens. A simple new approach to variable selection in regression, with application to genetic fine-mapping. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, https://doi.org/10.1111/rssb.12388, 2020.


## Getting Started

- Clone this repository using the following git command:
   
    `git clone https://github.com/getian107/SuSiEx.git`

    Alternatively, download the source files from the github website (https://github.com/getian107/SuSiEx).
    
- **SuSiEx** requires Python package **scipy** (https://www.scipy.org/) installed.

- Once Python and its dependencies have been installed, running

    `./SuSiEx.py --help` or `./SuSiEx.py -h`

    will print a list of command-line options.


## Using SuSiEx

`
python SuSiEx.py --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --ref_file=REF_FILE --ld_file=LD_MATRIX_FILE --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILENAME --chr=CHR --bp=BP --chr_col=CHR_COL --snp_col=SNP_COL --bp_col=BP_COL --a1_col=A1_COL --a2_col=A2_COL --eff_col=EFF_COL --se_col=SE_COL --pval_col=PVAL_COL --plink=PLINK [--keep-ambig=KEEP_AMBIGUOUS_SNPS --maf=MAF_THRESHOLD --n_sig=NUMBER_OF_SIGNALS --level=LEVEL --min_purity=MINIMUM_PURITY --mult-step=MULT_STEP_FITTING --pval_thresh=MARGINAL_PVAL_THRESHOLD --max_iter=MAXIMUM_ITERATIONS --tol=TOLERANCE]
`

- SUM_STATS_FILE (required): Full path and filename of the GWAS summary statistics for each population, separated by comma. Each file must contain a header line. The column corresponding to the effect size estimates must have a header of BETA or OR, indicating whether the effect estimates are regression coefficients or odds ratios.

- GWAS_SAMPLE_SIZE (required): Sample size of the GWAS for each population, in the order corresponding to the GWAS summary statistics files, separated by comma.

- REF_FILE (required): Full path and filename prefix of the reference panel for each population in PLINK binary format (.bed/.bim/.fam), in the order corresponding to the GWAS summary statistics files, separated by comma. The format for SNP IDs (e.g., rs IDs or chr:bp:a1:a2) and version of genome build must be consistent across GWAS summary statistics files and reference panels. A curated 1000 Genomes phase 3 reference panle in PLINK format can be download from the MAGMA website: https://ctg.cncr.nl/software/magma.

- LD_MATRIX_FILE (required): Full path and filename prefix of the LD matrix to be computed from the reference panels for each population, in the order corresponding to the GWAS summary statistics files, separated by comma.

- OUTPUT_DIR (required): Full path of the output directory.

- OUTPUT_FILENAME (required): Prefix of the output files.

- CHR (required): Chromosome of the fine-mapping region.

- BP (required): The start and end base pair coordinate of the fine-mapping region, separated by comma.

- CHR_COL (required): The column number of the chromosome code in each GWAS summary statistics file, separated by comma.

- SNP_COL (required): The column number of the SNP IDs in each GWAS summary statistics file, separated by comma.

- BP_COL (required): The column number of the base pair coordinate in each GWAS summary statistics file, separated by comma.

- A1_COL (required): The column number of the A1 allele (effective allele) in each GWAS summary statistics file, separated by comma.

- A2_COL (required): The column number of the A2 allele (alternative allele) in each GWAS summary statistics file, separated by comma.

- EFF_COL (required): The column number of the effect size estimate (beta or odds ratio) in each GWAS summary statistics file, separated by comma.

- SE_COL (required): The column number of the standard error of the effect size estimate in each GWAS summary statistics file, separated by comma.

- PVAL_COL (required): The column number of the p-value in each GWAS summary statistics file, separated by comma.

- PLINK (required): The full path and filename of PLINK.

- KEEP_AMBIGUOUS_SNPS (optional): If False, all ambiguous SNPs (A/T and G/C SNPs) will be removed; if True, ambiguous SNPs with A1 and A2 matching exactly with the A1 and A2 in the reference panel will be retained. Default is False.

- MAF_THRESHOLD (optional): Minor allele frequency (MAF) threshold applied to the reference panel. Default is 0.005.

- NUMBER_OF_SIGNALS (optional): Maximum number of causal signals in the SuSiEx model. Default is 10.

- LEVEL (optional): Coverage level of the credible set. Default is 95%.

- MINIMUM_PURITY (optional): Minimum purity of the credible set. Credible sets with purity below this specified value will be filtered out. Default is 0.5.

- MULT_STEP_FITTING (optional): Use the multi-step modeling fitting approach. The model is first fitted using 5 signals. If the algorithm doesn't converge, the maximum number of signals is progressively reduced from 5 to 1 until convergence. If 5 credible sets are identified, the maximum number of signals is increased to 10. If the algorithm doesn't converge with 10 signals, the maximum number of signals is then progressively reduced from 10 to 5 until convergence.

- MARGINAL_PVAL_THRESHOLD (optional): Filtering threshold for the marginal p-value. Credible sets containing no marginal p-value below this specified value will be filtered out. Default is 1e-05.

- MAXIMUM_ITERATIONS (optional): Maximum number of iterations allowed for the model fitting algorithm. Default is 100.

- TOLERANCE (optional): Tolerance for the convergence of the variational algorithm. Default is 1e-04.


## Output

A `.summary` file, a `.cs` file and a `.snp` file will be written to the specified output directory.
The `.snp` file contains a list of SNPs that are used in the fine-mapping algorithm.
These are the SNPs that are located in the specified fine-mapping region, available in the GWAS summary statistics and the reference panel for at least one population, and survived the minor allele frequency filtering.

If the varitional algorithm did not converge, "FAIL" will be written to both the `.summary` file and the `.cs` file.
If no credible set was identified at the specified coverage level after purity and marginal p-value filtering, "NULL" will be written to both the `.summary` file and the `.cs` file.

Otherwise, the `.summary` file contains credible set level information, which has a header line of the fine-mapping region and the following columns:

- CS_ID: ID of the credible set.

- CS_LENGTH: Size (number of SNPs) of the credible set.

- CS_PURITY: Purity of the credible set.

- MAX_PIP_SNP: SNP in the credible set that had the largest posterior inclusion probability (PIP).

- BP: The base pair coordinate of the MAX_PIP_SNP.

- REF_ALLELE: The reference allele of the MAX_PIP_SNP in each population, separated by comma.

- ALT_ALLELE: The alternative allele of the MAX_PIP_SNP in each population, separated by comma.

- REF_FRQ: Frequency of the reference allele in each reference panel, separated by comma.

- BETA: Marginal per-allele effect size of the MAX_PIP_SNP with respect to the reference allele in each population, separated by comma.

- SE: The standard error of the marginal per-allele effect size of the MAX_PIP_SNP in each population, separated by comma.

- -LOG10P: -log10 of the marginal p-value of the MAX_PIP_SNP in each population, separated by comma.

- MAX_PIP: Maximum posterior inclusion probability (PIP) in the credible set.


The `.cs` file contains information for all the SNPs included in credible sets and has the following columns:

- CS_ID: ID of the credible set.

- SNP: SNP identifier.

- BP: The base pair coordinate of the SNP.

- REF_ALLELE: The reference allele of the SNP in each population, separated by comma.

- ALT_ALLELE: The alternative allele of the SNP in each population, separated by comma.

- REF_FRQ: Frequency of the reference allele in each reference panel, separated by comma.

- BETA: Marginal per-allele effect size of the SNP with respect to the reference allele in each population, separated by comma.

- SE: The standard error of the marginal per-allele effect size of the SNP in each population, separated by comma.

- -LOG10P: -log10 of the marginal p-value of the SNP in each population, separated by comma.

- CS_PIP: Posterior inclusion probability (PIP) of the SNP.

- OVRL_PIP: Posterior inclusion probability (PIP) of the SNP in any of the credible set.


## Example

```
python SuSiEx.py \
    --sst_file=${sst_eur},${sst_afr} \
    --n_gwas=${n_eur},${n_afr} \
    --ref_file=${ref_eur},${ref_afr} \
    --ld_file=${ld_eur},${ld_afr} \
    --out_dir=${out_dir} \
    --out_name=${out_name} \
    --chr=1 \
    --bp=94813205,95812998 \
    --chr_col=1,1 \
    --snp_col=2,2 \
    --bp_col=3,3 \
    --a1_col=4,4 \
    --a2_col=5,5 \
    --eff_col=8,8 \
    --se_col=9,9 \
    --pval_col=10,10 \
    --keep-ambig=False \
    --mult-step=True \
    --maf=0.01 \
    --plink=${plink_dir}/plink_v1.90
```


## Support

Please direct questions or bug reports to Tian Ge (tge1@mgh.harvard.edu).

