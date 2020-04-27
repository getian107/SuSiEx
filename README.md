# SuSiEx

**SuSiEx** is a Python based command line tool that performs cross-ethnic fine-mapping using GWAS summary statistics and external LD reference panels. The method is built on the sum of single effects (SuSiE) model:

Gao Wang, Abhishek Sarkar, Peter Carbonetto, Matthew Stephens. A simple new approach to variable selection in regression, with application to genetic fine-mapping. *bioRxiv*, https://www.biorxiv.org/content/10.1101/501114v2, 2019.


## Getting Started

- Clone this repository using the following git command:
   
    `git clone https://github.com/getian107/SuSiEx.git`

    Alternatively, download the source files from the github website (https://github.com/getian107/SuSiEx).
    
- **SuSiEx** requires Python packages **scipy** (https://www.scipy.org/) installed.

- Once Python and its dependencies have been installed, running

    `./SuSiEx.py --help` or `./SuSiEx.py -h`

    will print a list of command-line options.


## Using SuSiEx

`
python SuSiEx.py --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --ld_file=LD_MATRIX --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILENAME [--n_sig=NUMBER_OF_SIGNALS --level=LEVEL --min_purity=MINIMUM_PURITY --pval_thresh=MARGINAL_PVAL_THRESHOLD --max_iter=MAXIMUM_ITERATIONS --tol=TOLERANCE]
`
- SUM_STATS_FILE (required): Full paths and file names of the GWAS summary statistics, separated by comma.

- GWAS_SAMPLE_SIZE (required): Sample sizes of the GWAS, in the order corresponding to the GWAS summary statistics files, separated by comma.

- LD_MATRIX (required): Full paths and file names of the LD matrices computed from reference panels, in the order corresponding to the GWAS summary statistics files, separated by comma.

- OUTPUT_DIR (required): Output directory.

- OUTPUT_FILENAME (required): Prefix of the output files.

- NUMBER_OF_SIGNALS (optional): Maximum number of causal signals in the SuSiEx model. Default is 10.

- LEVEL (optional): Coverage level of the credible set. Default is 95%.

- MINIMUM_PURITY (optional): Minimum purity of the credible set. Credible sets with purity below this specified value will be filtered out. Default is 0.5.

- MARGINAL_PVAL_THRESHOLD (optional): Filtering threshold for the marginal p-value. Credible sets containing no marginal p-value below this specified value will be filtered out. Default is 1e-05.

- MAXIMUM_ITERATIONS (optional): Maximum number of iterations allowed for the model fitting algorithm. Default is 100.

- TOLERANCE (optional): Tolerance for the convergence of the variational algorithm. Default is 1e-04.


## Output

A `.summary` file and a `.cs` file will be written to the specified output directory.
If the varitional algorithm did not converge, "FAIL" will be written to both files.
If no credible set was identified at the specified coverage level, "NULL" will be written to both files.
Otherwise, the `.summary` file contains the following columns:

- CS_ID: ID of the credible set.

- CS_LENGTH: size (number of SNPs) of the credible set.

- CS_PURITY: purity of the credible set.

- MAX_PIP: maximum posterior inclusion probability (PIP) in the credible set.

The `.cs` file contains the following columns:

- CS_ID: ID of the credible set.

- SNP: SNP identifier.

- CS_PIP: posterior inclusion probability (PIP) of the SNP.

- OVRL_PIP: posterior inclusion probability (PIP) of the SNP in any of the credible set.


## Example

```
python SuSiEx.py \
    --sst_file=${sst_eur},${sst_afr} \   
    --n_gwas=${n_eur},${n_afr} \    
    --ld_file=${ld_eur},${ld_afr} \   
    --out_dir=${out_dir} \   
    --out_name=${out_name}
```


## Support

Please direct questions or bug reports to Tian Ge (tge1@mgh.harvard.edu).

