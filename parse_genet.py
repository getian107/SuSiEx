#!/usr/bin/env python

"""
Parse GWAS summary statistics and LD matrices

"""


import scipy as sp
from scipy.stats import norm


def parse_ref(ref_file, chrom):
    print('... parse reference file: %s ...' % (ref_file + '.bim'))

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[]}
    with open(ref_file + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom:
                ref_dict['CHR'].append(chrom)
                ref_dict['SNP'].append(ll[1])
                ref_dict['BP'].append(int(ll[3]))
                ref_dict['A1'].append(ll[4])
                ref_dict['A2'].append(ll[5])

    print('... %d SNPs on chromosome %d read from %s ...' % (len(ref_dict['SNP']), chrom, ref_file))
    return ref_dict


def parse_sumstats(ref_dict, sst_file, n_subj):
    print('... parse sumstats file: %s ...' % sst_file)

    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            sst_dict['SNP'].append(ll[0])
            sst_dict['A1'].append(ll[1])
            sst_dict['A2'].append(ll[2])

    print('... %d SNPs read from %s ...' % (len(sst_dict['SNP']), sst_file))


    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2']))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1']], [mapping[aa] for aa in sst_dict['A2']])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2']], [mapping[aa] for aa in sst_dict['A1']]))

    comm_snp = ref_snp & sst_snp

    print('... %d common SNPs in the reference and sumstats ...' % len(comm_snp))


    n_sqrt = sp.sqrt(n_subj)
    sst_eff = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = sp.log(float(ll[3]))

                p = max(float(ll[4]), 1e-323)

                if sp.isnan(beta) or sp.isnan(p):
                    beta_std = 0.0
                else:
                    beta_std = sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt

                sst_eff.update({snp: beta_std})

            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = sp.log(float(ll[3]))
 
                p = max(float(ll[4]), 1e-323)

                if sp.isnan(beta) or sp.isnan(p):
                    beta_std = 0.0
                else:
                    beta_std = -1*sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt

                sst_eff.update({snp: beta_std})


    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'BETA':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            sst_dict['SNP'].append(snp)
            sst_dict['CHR'].append(ref_dict['CHR'][ii])
            sst_dict['BP'].append(ref_dict['BP'][ii])
            sst_dict['A1'].append(ref_dict['A1'][ii])
            sst_dict['A2'].append(ref_dict['A2'][ii])
            sst_dict['BETA'].append(sst_eff[snp])

    return sst_dict


def parse_ld(ld_file):

    with open(ld_file) as ff:
        ld = [[float(val) for val in (line.strip()).split()] for line in ff]

    return ld





