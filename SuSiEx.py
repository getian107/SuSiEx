#!/usr/bin/env python

"""
SuSiEx: trans-ethnic fine-mapping

Usage:

Simulation mode:

python SuSiEx.py --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --ld_file=LD_MATRIX_FILE --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILENAME --sim=True
		 [--n_sig=NUMBER_OF_SIGNALS --level=LEVEL --min_purity=MINIMUM_PURITY --pval_thresh=MARGINAL_PVAL_THRESHOLD
		  --max_iter=MAXIMUM_ITERATIONS --tol=TOLERANCE]

Non-simulation mode:

python SuSiEx.py --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --ref_file=REF_FILE --ld_file=LD_MATRIX_FILE --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILENAME --sim=False
                 --chr=CHR --bp=BP --chr_col=CHR_COL --snp_col=SNP_COL --bp_col=BP_COL --a1_col=A1_COL --a2_col=A2_COL
                 --eff_col=EFF_COL --pval_col=PVAL_COL --plink=PLINK
		 [--keep-ambig=KEEP_AMBIGUOUS_SNPS --maf=MAF_THRESHOLD --n_sig=NUMBER_OF_SIGNALS --level=LEVEL --min_purity=MINIMUM_PURITY
		  --pval_thresh=MARGINAL_PVAL_THRESHOLD --max_iter=MAXIMUM_ITERATIONS --tol=TOLERANCE]

"""


import sys
import getopt

import parse_genet
import parse_genet_sim
import SuSiE
import parse_pip
import write_cs
import scipy as sp


def parse_param():
    long_opts_list = ['sst_file=', 'n_gwas=', 'ref_file=', 'ld_file=', 'out_dir=', 'out_name=', 'sim=',
                      'chr=', 'bp=', 'chr_col=', 'snp_col=', 'bp_col=', 'a1_col=', 'a2_col=',
                      'eff_col=', 'pval_col=', 'plink=', 'keep-ambig=', 'maf=',
		      'n_sig=', 'level=', 'min_purity=', 'pval_thresh=', 'max_iter=', 'tol=', 'help']

    param_dict = {'sst_file': None, 'n_gwas': None, 'ref_file': None, 'ld_file': None, 'out_dir': None, 'out_name': None, 'sim': 'True',
                  'chr': None, 'bp': None, 'chr_col': None, 'snp_col': None, 'bp_col': None, 'a1_col': None, 'a2_col': None,
                  'eff_col': None, 'pval_col': None, 'plink': None, 'keep-ambig': 'False', 'maf': 0.005,
		  'n_sig': 10, 'level': 0.95, 'min_purity': 0.5, 'pval_thresh': 1e-5, 'max_iter': 100, 'tol': 1e-4}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--sst_file": param_dict['sst_file'] = arg.split(',')
            elif opt == "--n_gwas": param_dict['n_gwas'] = map(int, arg.split(','))
            elif opt == "--ref_file": param_dict['ref_file'] = arg.split(',')
            elif opt == "--ld_file": param_dict['ld_file'] = arg.split(',')
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--sim": param_dict['sim'] = arg
            elif opt == "--chr": param_dict['chr'] = int(arg)
            elif opt == "--bp": param_dict['bp'] = map(int, arg.split(','))
            elif opt == "--chr_col": param_dict['chr_col'] = map(int, arg.split(','))
            elif opt == "--snp_col": param_dict['snp_col'] = map(int, arg.split(','))
            elif opt == "--bp_col": param_dict['bp_col'] = map(int, arg.split(','))
            elif opt == "--a1_col": param_dict['a1_col'] = map(int, arg.split(','))
            elif opt == "--a2_col": param_dict['a2_col'] = map(int, arg.split(','))
            elif opt == "--eff_col": param_dict['eff_col'] = map(int, arg.split(','))
            elif opt == "--pval_col": param_dict['pval_col'] = map(int, arg.split(','))
            elif opt == "--plink": param_dict['plink'] = arg
            elif opt == "--keep-ambig": param_dict['keep-ambig'] = arg
            elif opt == "--maf": param_dict['maf'] = float(arg)
            elif opt == "--n_sig": param_dict['n_sig'] = int(arg)
            elif opt == "--level": param_dict['level'] = float(arg)
            elif opt == "--min_purity": param_dict['min_purity'] = float(arg)
            elif opt == "--pval_thresh": param_dict['pval_thresh'] = float(arg)
            elif opt == "--max_iter": param_dict['max_iter'] = int(arg)
            elif opt == "--tol": param_dict['tol'] = float(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['sst_file'] == None:
        print('* Please provide GWAS summary statistics using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Please provide the sample size of each GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['ld_file'] == None:
        print('* Please provide the directory and filename prefix of the LD matrix for each GWAS summary statistics file using --ld_file\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Please specify the prefix of the output filenames using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sst_file']) != len(param_dict['n_gwas']) or
          len(param_dict['sst_file']) != len(param_dict['ld_file'])):
        print('* Length of sst_file, n_gwas and ld_file does not match\n')
        sys.exit(2)

    n_pop = len(param_dict['sst_file'])

    if param_dict['sim'] == 'False':
        if param_dict['ref_file'] == None or len(param_dict['ref_file']) != n_pop:
            print('* Please provide a reference panel for each GWAS using --ref-file\n')
            sys.exit(2)
        elif param_dict['chr'] == None:
            print('* Please provide the chromosome code of the fine-mapping region using --chr\n')
            sys.exit(2)
        elif param_dict['bp'] == None or len(param_dict['bp']) != 2:
            print('* Please provide the start and end base pair coordinate of the fine-mapping region using --bp\n')
            sys.exit(2)
        elif param_dict['chr_col'] == None or len(param_dict['chr_col']) != n_pop:
            print('* Please provide the column number of the chromosome code in each GWAS summary statistics file using --chr_col\n')
            sys.exit(2)
        elif param_dict['snp_col'] == None or len(param_dict['snp_col']) != n_pop:
            print('* Please provide the column number of the SNP IDs in each GWAS summary statistics file using --snp_col\n')
            sys.exit(2)
        elif param_dict['bp_col'] == None or len(param_dict['bp_col']) != n_pop:
            print('* Please provide the column number of the base pair coordinate in each GWAS summary statistics file using --bp_col\n')
            sys.exit(2)
        elif param_dict['a1_col'] == None or len(param_dict['a1_col']) != n_pop:
            print('* Please provide the column number of the A1 allele in each GWAS summary statistics file using --a1_col\n')
            sys.exit(2)
        elif param_dict['a2_col'] == None or len(param_dict['a2_col']) != n_pop:
            print('* Please provide the column number of the A2 allele in each GWAS summary statistics file using --a2_col\n')
            sys.exit(2)
        elif param_dict['eff_col'] == None or len(param_dict['eff_col']) != n_pop:
            print('* Please provide the column number of the effect size estimate (beta or odds ratio) in each GWAS summary statistics file using --eff_col\n')
            sys.exit(2)
        elif param_dict['pval_col'] == None or len(param_dict['pval_col']) != n_pop:
            print('* Please provide the column number of the p-value in each GWAS summary statistics file using --pval_col\n')
            sys.exit(2)
        elif param_dict['plink'] == None:
            print('* Please provide the directory and filename of PLINK using --plink\n')
            sys.exit(2)
        elif param_dict['keep-ambig'] == 'False':
            print('* All ambiguous SNPs will be removed\n')
        elif param_dict['keep-ambig'] == 'True':
            print('* Ambiguous SNPs will be retained if A1/A2 in the summary statistics match A1/A2 in the reference panel; ' +
                     'Use --keep-ambig=False to remove all ambiguous SNPs if not certain about allele matching between summary statistics and the reference\n')


    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    n_pop = len(param_dict['sst_file'])

    if param_dict['sim'] == 'False':
        for pp in range(n_pop):
            parse_genet.calc_ld(param_dict['ref_file'][pp], param_dict['ld_file'][pp], param_dict['plink'], param_dict['chr'], param_dict['bp'], param_dict['maf'])

        ref_dict = {}
        for pp in range(n_pop):
            ref_dict[pp] = parse_genet.parse_ref(param_dict['ref_file'][pp], param_dict['ld_file'][pp]+'_ref', param_dict['ld_file'][pp]+'_frq')

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet.parse_sumstats(param_dict['sst_file'][pp], ref_dict[pp], param_dict['chr'], param_dict['bp'],
                           param_dict['chr_col'][pp], param_dict['snp_col'][pp], param_dict['bp_col'][pp], param_dict['a1_col'][pp], param_dict['a2_col'][pp],
                           param_dict['eff_col'][pp], param_dict['pval_col'][pp], param_dict['n_gwas'][pp], param_dict['keep-ambig'])

        ld_dict = {}
        for pp in range(n_pop):
            ld_dict[pp] = parse_genet.parse_ld(param_dict['ld_file'][pp]+'.ld', ref_dict[pp], sst_dict[pp])

        snp_dict, beta, tau_sq, pval_min, ind, ld = parse_genet.align_sumstats(sst_dict, ld_dict, n_pop)

        for pp in range(n_pop):
            parse_genet.clean_files(param_dict['ld_file'][pp])

    else:
        print('### Simulation Mode ###')

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet_sim.parse_sumstats(param_dict['sst_file'][pp], param_dict['n_gwas'][pp])

        ld_dict = {}
        for pp in range(n_pop):
            ld_dict[pp] = parse_genet_sim.parse_ld(param_dict['ld_file'][pp])

        snp_dict, beta, tau_sq, pval_min, ind, ld = parse_genet_sim.align_sumstats(sst_dict, ld_dict, n_pop)



    alpha, b, b_sq, sigma_sq, elbo_new, n_iter, flag = \
	SuSiE.SUSIE_sst_xethn(beta, ind, param_dict['n_gwas'], ld, tau_sq, param_dict['n_sig'], param_dict['max_iter'], param_dict['tol'])

    n_cs, alpha, cs_bin, cs_purity, pip = \
	parse_pip.pip(flag, alpha, ld, pval_min, param_dict['level'], param_dict['min_purity'], param_dict['pval_thresh'])

    write_cs.write_cs(param_dict['chr'], param_dict['bp'], n_pop, n_cs, snp_dict, beta, ind, param_dict['n_gwas'], 
        alpha, cs_bin, cs_purity, pip, param_dict['out_dir'], param_dict['out_name'])

    print('... Done ...\n')


if __name__ == '__main__':
    main()


