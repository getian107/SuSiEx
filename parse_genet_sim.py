#!/usr/bin/env python

"""
Parse GWAS summary statistics and LD matrices

"""


import scipy as sp
from scipy.stats import norm


def parse_sumstats(sst_file, n_subj):
    print('... parse sumstats file: %s ...' % sst_file)

    sst_dict = {'SNP':[], 'BETA':[], 'P': [], 'MISS':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            sst_dict['SNP'].append(ll[1])
            beta = ll[10]
            pval = ll[11]

            if beta == 'NA' or pval == 'NA':
                sst_dict['BETA'].append(0.0)
                sst_dict['P'].append(1.0)
                sst_dict['MISS'].append(True)
            else:
                sst_dict['BETA'].append(float(beta)/sp.sqrt(n_subj))
                sst_dict['P'].append(float(pval))
                sst_dict['MISS'].append(False)

    sst_dict['BETA'] = sp.array(sst_dict['BETA'], ndmin=2).T
    sst_dict['P'] = sp.array(sst_dict['P'], ndmin=2).T
    sst_dict['MISS'] = sp.array(sst_dict['MISS'], ndmin=2).T

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))

    return sst_dict


def parse_ld(ld_file):
    print('... parse LD file: %s ...' % ld_file)

    with open(ld_file) as ff:
        ld = [[float(val) for val in (line.strip()).split()] for line in ff]

    ld = sp.array(ld)
    ld[sp.isnan(ld)] = 0.0
    n_snp = len(ld)
    for ii in range(n_snp):
        ld[ii,ii] = 1.0

    print('... %s LD matrix read from %s ...' % (sp.shape(ld), ld_file))

    return ld


def align_sumstats(sst_dict, ld_dict, n_pop):
    print('... align summary statistics and LD references across populations ...')

    snp_dict = {}
    snp_dict['SNP'] = sst_dict[0]['SNP']
    n_snp = len(snp_dict['SNP'])

    beta = sp.zeros((n_snp,n_pop))
    pval = sp.zeros((n_snp,n_pop))
    ind = sp.zeros((n_snp,n_pop), dtype=bool)
    ld = sp.zeros((n_snp,n_snp,n_pop))

    for pp in range(n_pop):
        beta[:,pp,None] = sst_dict[pp]['BETA']
        pval[:,pp,None] = sst_dict[pp]['P']
        ind[:,pp,None] = sst_dict[pp]['MISS']
        ld[:,:,pp] = ld_dict[pp]

    tau_sq = sp.amax(beta**2,axis=0)
    pval_min = sp.amin(pval,axis=1).reshape(len(pval),1)
    
    return snp_dict, beta, tau_sq, pval_min, ind, ld


