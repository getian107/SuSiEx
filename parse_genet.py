#!/usr/bin/env python

"""
Parse GWAS summary statistics and LD matrices

"""


import scipy as sp
from scipy.stats import norm
import gzip
import subprocess


def calc_ld(ref_file, ld_file, plink, chrom, bp, maf):
    print('... calculate LD matrix: %s ...' % (ld_file+'.ld'))

    SNP = []
    with open(ref_file+'.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == chrom and int(ll[3]) >= bp[0] and int(ll[3]) <= bp[1]:
                SNP.append(ll[1])

    with open(ld_file+'.snp', 'w') as ff:
        for snp in SNP:
            ff.write('%s\n' % snp)


    cmd = '%s --bfile %s --keep-allele-order --chr %d --extract %s --maf %f --make-bed --out %s' \
        % (plink, ref_file, chrom, ld_file+'.snp', maf, ld_file+'_ref')
    subprocess.check_output(cmd, shell=True)

    cmd = '%s --bfile %s --keep-allele-order --r square --out %s' %(plink, ld_file+'_ref', ld_file)
    subprocess.check_output(cmd, shell=True)

    cmd = '%s --bfile %s --keep-allele-order --freq --out %s' %(plink, ld_file+'_ref', ld_file+'_frq')
    subprocess.check_output(cmd, shell=True)


def parse_ref(ld_file_ref, ld_file_frq):
    print('... parse reference file: %s ...' % (ld_file_ref+'.bim'))

    ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'FRQ':[]}
    with open(ld_file_ref+'.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            ref_dict['CHR'].append(int(ll[0]))
            ref_dict['SNP'].append(ll[1])
            ref_dict['BP'].append(int(ll[3]))
            ref_dict['A1'].append(ll[4])
            ref_dict['A2'].append(ll[5])

    with open(ld_file_frq+'.frq') as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            ref_dict['FRQ'].append(ll[4])

    print('... %d SNPs in the fine-mapping region read from %s ...' % (len(ref_dict['SNP']), ld_file_ref+'.bim'))

    return ref_dict


def parse_sumstats(sst_file, ref_dict, chrom, bp, chr_col, snp_col, bp_col, a1_col, a2_col, eff_col, se_col, pval_col, n_subj, ambig):
    print('... parse sumstats file: %s ...' % sst_file)

    ATGC = ['A', 'T', 'G', 'C']
    sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[chr_col-1]) == chrom and int(ll[bp_col-1]) >= bp[0] and int(ll[bp_col-1]) <= bp[1] and \
                ll[a1_col-1] in ATGC and ll[a2_col-1] in ATGC:
                sst_dict['SNP'].append(ll[snp_col-1])
                sst_dict['A1'].append(ll[a1_col-1])
                sst_dict['A2'].append(ll[a2_col-1])

    print('... %d SNPs in the fine-mapping region read from %s ...' % (len(sst_dict['SNP']), sst_file))


    mapping = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2']))

    sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1']], [mapping[aa] for aa in sst_dict['A2']])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2']], [mapping[aa] for aa in sst_dict['A1']]))

    comm_snp = ref_snp & sst_snp


    n_sqrt = sp.sqrt(n_subj)
    sst_eff = {}
    sst_pval = {}
    sst_logp = {}
    sst_miss = {}
    with open(sst_file) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            chr = int(ll[chr_col-1]); snp = ll[snp_col-1]; a1 = ll[a1_col-1]; a2 = ll[a2_col-1]
            if chr != chrom or a1 not in ATGC or a2 not in ATGC:
                continue

            if (a1=='A' and a2=='T') or (a1=='T' and a2=='A') or (a1=='G' and a2=='C') or (a1=='C' and a2=='G'):
                if ambig == 'TRUE' and (snp, a1, a2) in comm_snp:
                    if ll[eff_col-1] == 'NA' or ll[pval_col-1] == 'NA':
                        sst_eff.update({snp: 0.0})
                        sst_pval.update({snp: 1.0})
                        sst_logp.update({snp: 0.0})
                        sst_miss.update({snp: True})
                    else:
                        if header[eff_col-1] == 'BETA':
                            beta = float(ll[eff_col-1])
                        elif header[eff_col-1] == 'OR':
                            beta = sp.log(float(ll[eff_col-1]))

                        se = float(ll[se_col-1])
                        pval = max(float(ll[pval_col-1]), 1e-323)

                        # sst_eff.update({snp: sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                        sst_eff.update({snp: beta/se/n_sqrt})
                        sst_pval.update({snp: pval})
                        sst_logp.update({snp: -1*sp.log10(2)-norm.logcdf(-1*abs(beta/se))/sp.log(10)})
                        sst_miss.update({snp: False})

                elif ambig == 'TRUE' and (snp, a2, a1) in comm_snp:
                    if ll[eff_col-1] == 'NA' or ll[pval_col-1] == 'NA':
                        sst_eff.update({snp: 0.0})
                        sst_pval.update({snp: 1.0})
                        sst_logp.update({snp: 0.0})
                        sst_miss.update({snp: True})
                    else:
                        if header[eff_col-1] == 'BETA':
                            beta = float(ll[eff_col-1])
                        elif header[eff_col-1] == 'OR':
                            beta = sp.log(float(ll[eff_col-1]))

                        se = float(ll[se_col-1])
                        pval = max(float(ll[pval_col-1]), 1e-323)

                        # sst_eff.update({snp: -1*sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                        sst_eff.update({snp: -1*beta/se/n_sqrt})
                        sst_pval.update({snp: pval})
                        sst_logp.update({snp: -1*sp.log10(2)-norm.logcdf(-1*abs(beta/se))/sp.log(10)})
                        sst_miss.update({snp: False})

            elif (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if ll[eff_col-1] == 'NA' or ll[pval_col-1] == 'NA':
                    sst_eff.update({snp: 0.0})
                    sst_pval.update({snp: 1.0})
                    sst_logp.update({snp: 0.0})
                    sst_miss.update({snp: True})
                else:
                    if header[eff_col-1] == 'BETA':
                        beta = float(ll[eff_col-1])
                    elif header[eff_col-1] == 'OR':
                        beta = sp.log(float(ll[eff_col-1]))

                    se = float(ll[se_col-1])
                    pval = max(float(ll[pval_col-1]), 1e-323)

                    # sst_eff.update({snp: sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                    sst_eff.update({snp: beta/se/n_sqrt})
                    sst_pval.update({snp: pval})
                    sst_logp.update({snp: -1*sp.log10(2)-norm.logcdf(-1*abs(beta/se))/sp.log(10)})
                    sst_miss.update({snp: False})

            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if ll[eff_col-1] == 'NA' or ll[pval_col-1] == 'NA':
                    sst_eff.update({snp: 0.0})
                    sst_pval.update({snp: 1.0})
                    sst_logp.update({snp: 0.0})
                    sst_miss.update({snp: True})
                else:
                    if header[eff_col-1] == 'BETA':
                        beta = float(ll[eff_col-1])
                    elif header[eff_col-1] == 'OR':
                        beta = sp.log(float(ll[eff_col-1]))

                    se = float(ll[se_col-1])
                    pval = max(float(ll[pval_col-1]), 1e-323)

                    # sst_eff.update({snp: -1*sp.sign(beta)*abs(norm.ppf(pval/2.0))/n_sqrt})
                    sst_eff.update({snp: -1*beta/se/n_sqrt})
                    sst_pval.update({snp: pval})
                    sst_logp.update({snp: -1*sp.log10(2)-norm.logcdf(-1*abs(beta/se))/sp.log(10)})
                    sst_miss.update({snp: False})


    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'FRQ':[], 'BETA':[], 'P':[], 'LOGP':[], 'MISS':[]}
    for (ii, snp) in enumerate(ref_dict['SNP']):
        if snp in sst_eff:
            sst_dict['CHR'].append(ref_dict['CHR'][ii])
            sst_dict['SNP'].append(ref_dict['SNP'][ii])
            sst_dict['BP'].append(ref_dict['BP'][ii])
            sst_dict['A1'].append(ref_dict['A1'][ii])
            sst_dict['A2'].append(ref_dict['A2'][ii])
            sst_dict['FRQ'].append(ref_dict['FRQ'][ii])
            sst_dict['BETA'].append(sst_eff[snp])
            sst_dict['P'].append(sst_pval[snp])
            sst_dict['LOGP'].append(sst_logp[snp])
            sst_dict['MISS'].append(sst_miss[snp])

    sst_dict['FRQ'] = sp.array(sst_dict['FRQ'], ndmin=2).T
    sst_dict['BETA'] = sp.array(sst_dict['BETA'], ndmin=2).T
    sst_dict['P'] = sp.array(sst_dict['P'], ndmin=2).T
    sst_dict['LOGP'] = sp.array(sst_dict['LOGP'], ndmin=2).T
    sst_dict['MISS'] = sp.array(sst_dict['MISS'], ndmin=2).T

    print('... %d common SNPs in the reference and sumstats ...' % len(sst_dict['SNP']))

    return sst_dict


def parse_ld(ld_file, ref_dict, sst_dict):
    print('... parse LD file: %s ...' % ld_file)

    with gzip.open(ld_file) as ff:
        ld = [[float(val) for val in (line.strip()).split()] for line in ff]

    ld = sp.array(ld)
    ld[sp.isnan(ld)] = 0.0
    n_snp = len(ld)
    for ii in range(n_snp):
        ld[ii,ii] = 1.0

    idx = [ii for (ii,snp) in enumerate(ref_dict['SNP']) if snp in sst_dict['SNP']]
    if len(idx) != 0:
        ld = ld[sp.ix_(idx,idx)]
    else:
        ld = sp.array([])

    print('... %s LD matrix read from %s ...' % (sp.shape(ld), ld_file))

    return ld


def align_sumstats(sst_dict, ld_dict, n_pop):
    print('... align summary statistics and LD references across populations ...')

    snp = []
    for pp in range(n_pop):
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
    snp_dict['BP'] = [bp[jj] for jj in idx]

    print('... a total of %d SNPs across populations ...' % len(snp_dict['SNP']))


    snp_dict['A1'] = sp.array([['NA']*n_pop]*n_snp)
    snp_dict['A2'] = sp.array([['NA']*n_pop]*n_snp)
    snp_dict['FRQ'] = sp.zeros((n_snp,n_pop))
    beta = sp.zeros((n_snp,n_pop))
    pval = sp.ones((n_snp,n_pop))
    logp = sp.zeros((n_snp,n_pop))
    ind = sp.ones((n_snp,n_pop), dtype=bool)
    ld = sp.zeros((n_snp,n_snp,n_pop))
    for pp in range(n_pop):
        for ii in range(n_snp):
            ld[ii,ii,pp] = 1.0

    for pp in range(n_pop):
        idx = [jj for jj in range(n_snp) if snp_dict['SNP'][jj] in sst_dict[pp]['SNP']]
        idx_pp = [sst_dict[pp]['SNP'].index(snp_dict['SNP'][jj]) for jj in idx] 

        snp_dict['A1'][idx,pp] = [sst_dict[pp]['A1'][jj] for jj in idx_pp]
        snp_dict['A2'][idx,pp] = [sst_dict[pp]['A2'][jj] for jj in idx_pp]
        snp_dict['FRQ'][idx,pp,None] = sst_dict[pp]['FRQ'][idx_pp]
        beta[idx,pp,None] = sst_dict[pp]['BETA'][idx_pp]
        pval[idx,pp,None] = sst_dict[pp]['P'][idx_pp]
        logp[idx,pp,None] = sst_dict[pp]['LOGP'][idx_pp]
        ind[idx,pp,None] = sst_dict[pp]['MISS'][idx_pp]
        ld[sp.ix_(idx,idx,[pp])] = ld_dict[pp][sp.ix_(idx_pp,idx_pp)].reshape(len(idx),len(idx),1)

    tau_sq = sp.amax(beta**2,axis=0)
    pval_min = sp.amin(pval,axis=1).reshape(len(pval),1)

    return snp_dict, beta, tau_sq, pval, logp, pval_min, ind, ld


def clean_files(ld_file):
    subprocess.check_output('rm '+ld_file+'.snp', shell=True)
    subprocess.check_output('rm '+ld_file+'_ref.bed', shell=True)
    subprocess.check_output('rm '+ld_file+'_ref.bim', shell=True)
    subprocess.check_output('rm '+ld_file+'_ref.fam', shell=True)
    subprocess.check_output('rm '+ld_file+'_ref.log', shell=True)
    subprocess.check_output('rm '+ld_file+'_frq.frq', shell=True)
    subprocess.check_output('rm '+ld_file+'_frq.log', shell=True)
    subprocess.check_output('rm '+ld_file+'.ld', shell=True)
    subprocess.check_output('rm '+ld_file+'.log', shell=True)


