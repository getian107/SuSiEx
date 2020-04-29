#!/usr/bin/env python

"""
SuSiEx: trans-ethnic fine-mapping

Usage:
python SuSiEx.py --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --ld_file=LD_MATRIX --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILENAME
		 [--n_sig=NUMBER_OF_SIGNALS --level=LEVEL --min_purity=MINIMUM_PURITY --pval_thresh=MARGINAL_PVAL_THRESHOLD
		  --max_iter=MAXIMUM_ITERATIONS --tol=TOLERANCE]

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
    long_opts_list = ['sst_file=', 'n_gwas=', 'ld_file=', 'out_dir=', 'out_name=',
		      'n_sig=', 'level=', 'min_purity=', 'pval_thresh=', 'max_iter=', 'tol=',
		      'sim=', 'ref_file=', 'help']

    param_dict = {'sst_file': None, 'n_gwas': None, 'ld_file': None, 'out_dir': None, 'out_name': None,
		  'n_sig': 10, 'level': 0.95, 'min_purity': 0.5, 'pval_thresh': 1e-5, 'max_iter': 100, 'tol': 1e-4,
		  'sim': 'True', 'ref_file': None}

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
            elif opt == "--ld_file": param_dict['ld_file'] = arg.split(',')
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--n_sig": param_dict['n_sig'] = int(arg)
            elif opt == "--level": param_dict['level'] = float(arg)
            elif opt == "--min_purity": param_dict['min_purity'] = float(arg)
            elif opt == "--pval_thresh": param_dict['pval_thresh'] = float(arg)
            elif opt == "--max_iter": param_dict['max_iter'] = int(arg)
            elif opt == "--tol": param_dict['tol'] = float(arg)
            elif opt == "--sim": param_dict['sim'] = arg
	    elif opt == "--ref_file": param_dict['ref_file'] = arg.split(',')
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['sst_file'] == None:
        print('* Please provide GWAS summary statistics file(s) using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Please provide the sample size(s) of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['ld_file'] == None:
        print('* Please provide LD matrix file(s) using --sst_file\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Please specify the output filename using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sst_file']) != len(param_dict['n_gwas']) or 
          len(param_dict['sst_file']) != len(param_dict['ld_file'])):
        print('* Length of sst_file, n_gwas and ld_file does not match\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    n_pop = len(param_dict['sst_file'])

    if param_dict['sim'] == 'False':
        ref_dict = {}
        for pp in range(n_pop):
            ref_dict[pp] = parse_genet.parse_ref(param_dict['ref_file'][pp])

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet.parse_sumstats(param_dict['sst_file'][pp], ref_dict[pp], param_dict['n_gwas'][pp])

        ld_dict = {}
        for pp in range(n_pop):
            ld_dict[pp] = parse_genet.parse_ld(param_dict['ld_file'][pp], ref_dict[pp], sst_dict[pp])

        snp_dict, beta, tau_sq, pval_min, ind, ld = parse_genet.align_sumstats(sst_dict, ld_dict, n_pop)

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

    write_cs.write_cs(n_cs, snp_dict, alpha, cs_bin, cs_purity, pip, param_dict['out_dir'], param_dict['out_name'])

    print('... Done ...\n')


if __name__ == '__main__':
    main()


