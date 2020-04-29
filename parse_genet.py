#!/usr/bin/env python

"""
Parse GWAS summary statistics and LD matrices

"""


import scipy as sp
from scipy.stats import norm


def parse_ref(ref_file):
    print('... parse reference file: %s ...' % (ref_file + '.bim'))

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[]}
    with open(ref_file + '.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            ref_dict['CHR'].append(ll[0])
            ref_dict['SNP'].append(ll[1])
            ref_dict['BP'].append(int(ll[3]))
            ref_dict['A1'].append(ll[4])
            ref_dict['A2'].append(ll[5])

    print('... %d SNPs read from %s ...' % (len(ref_dict['SNP']), ref_file))

    return ref_dict


def parse_sumstats(sst_file, ref_dict, n_subj):
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


    mapping = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    ATGC = ['A', 'T', 'G', 'C']

    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2']))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1'] if aa in ATGC], [mapping[aa] for aa in sst_dict['A2'] if aa in ATGC])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2'] if aa in ATGC], [mapping[aa] for aa in sst_dict['A1'] if aa in ATGC]))

    comm_snp = ref_snp & sst_snp

    print('... %d common SNPs in the reference and sumstats ...' % len(comm_snp))


    n_sqrt = sp.sqrt(n_subj)
    sst_eff = {}
    sst_pval = {}
    sst_miss = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]
            if a1 not in ATGC or a2 not in ATGC:
                continue

            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if ll[3] == 'NA' or ll[4] == 'NA':
                    sst_eff.update({snp: 0.0})
                    sst_pval.update({snp: 1.0})
                    sst_miss.update({snp: True})
                else:
                    if 'BETA' in header:
                        beta = float(ll[3])
                    elif 'OR' in header:
                        beta = sp.log(float(ll[3]))

                    pval = max(float(ll[4]), 1e-323)

                    sst_eff.update({snp: sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                    sst_pval.update({snp: pval})
                    sst_miss.update({snp: False})
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if ll[3] == 'NA' or ll[4] == 'NA':
                    sst_eff.update({snp: 0.0})
                    sst_pval.update({snp: 1.0})
                    sst_miss.update({snp: True})
                else:
                    if 'BETA' in header:
                        beta = float(ll[3])
                    elif 'OR' in header:
                        beta = sp.log(float(ll[3]))

                    pval = max(float(ll[4]), 1e-323)

                    sst_eff.update({snp: -1*sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                    sst_pval.update({snp: pval})
                    sst_miss.update({snp: False})                    


    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'BETA':[], 'P':[], 'MISS':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            sst_dict['CHR'].append(ref_dict['CHR'][ii])
            sst_dict['SNP'].append(ref_dict['SNP'][ii])
            sst_dict['BP'].append(ref_dict['BP'][ii])
            sst_dict['A1'].append(ref_dict['A1'][ii])
            sst_dict['A2'].append(ref_dict['A2'][ii])
            sst_dict['BETA'].append(sst_eff[snp])
            sst_dict['P'].append(sst_pval[snp])
            sst_dict['MISS'].append(sst_miss[snp])

    sst_dict['BETA'] = sp.array(sst_dict['BETA'], ndmin=2).T
    sst_dict['P'] = sp.array(sst_dict['P'], ndmin=2).T
    sst_dict['MISS'] = sp.array(sst_dict['MISS'], ndmin=2).T

    return sst_dict


def parse_ld(ld_file, ref_dict, sst_dict):
    print('... parse LD file: %s ...' % ld_file)

    with open(ld_file) as ff:
        ld = [[float(val) for val in (line.strip()).split()] for line in ff]

    ld = sp.array(ld)
    ld[sp.isnan(ld)] = 0.0
    n_snp = len(ld)
    for ii in range(n_snp):
        ld[ii,ii] = 1.0

    idx = [ii for (ii, snp) in enumerate(ref_dict['SNP']) if snp in sst_dict['SNP']]
    if len(idx) != 0:
        ld = ld[sp.ix_(idx,idx)]
    else:
        ld = sp.array([])

    print('... %s LD matrix read from %s ...' % (sp.shape(ld), ld_file))

    return ld


def align_sumstats(sst_dict, ld_dict, n_pop):
    print('... align summary statistics and LD references across populations ...')

    snp = []
    for pp in range(1,n_pop):
        snp = list(set(snp+sst_dict[pp]['SNP']))

    n_snp = len(snp)
    bp = []
    for jj in range(n_snp):
        for pp in range(n_pop):
            if snp[jj] in sst_dict[pp]['SNP']:
                idx = sst_dict[pp]['SNP'].index(snp[jj])
                bp.append(sst_dict[pp]['BP'][idx])
                break

    idx = sp.argsort(bp)
    snp_dict = {}
    snp_dict['SNP'] = [snp[jj] for jj in idx]


    beta = sp.zeros((n_snp,n_pop))
    pval = sp.zeros((n_snp,n_pop))
    ind = sp.zeros((n_snp,n_pop), dtype=bool)
    ld = sp.zeros((n_snp,n_snp,n_pop))

    for pp in range(n_pop):
        idx = [jj for jj in range(n_snp) if snp_dict['SNP'][jj] in sst_dict[pp]['SNP']]
        idx_pp = [sst_dict[pp]['SNP'].index(snp_dict['SNP'][jj]) for jj in idx] 

        beta[idx,pp,None] = sst_dict[pp]['BETA'][idx_pp]
        pval[idx,pp,None] = sst_dict[pp]['P'][idx_pp]
        ind[idx,pp,None] = sst_dict[pp]['MISS'][idx_pp]
        ld[sp.ix_(idx,idx,[pp])] = ld_dict[pp][sp.ix_(idx_pp,idx_pp)].reshape(len(idx),len(idx),1)

    tau_sq = sp.amax(beta**2,axis=0)
    pval_min = sp.amin(pval,axis=1).reshape(len(pval),1)

    return snp_dict, beta, tau_sq, pval_min, ind, ld


