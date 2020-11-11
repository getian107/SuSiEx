#!/usr/bin/env python

"""
Calculate credible sets and PIPs

"""


import scipy as sp


def pip(flag, alpha, ld, pval_min, level, min_purity, pval_thresh):
    print('... calculate credible sets and PIPs ...')

    if flag == False or sp.sum(sp.isnan(alpha)) > 0:
        n_cs = sp.nan
        alpha = sp.nan
        cs_bin = sp.nan
        cs_purity = sp.nan
        pip = sp.nan

        return n_cs, alpha, cs_bin, cs_purity, pip

    (n_snp,n_snp,n_pop) = sp.shape(ld)
    (n_snp,n_cs) = sp.shape(alpha)
    cs_bin = sp.zeros((n_snp,n_cs), dtype=bool)
    cs_purity = sp.zeros((n_cs,1))

    for ll in range(n_cs):
        idx = sp.argsort(-alpha[:,ll])
        pip_sort = alpha[idx,ll]
        sum_pip = sp.cumsum(pip_sort)
        cs_size = sp.argmax(sum_pip>=level)
        cs_idx = idx[0:cs_size+1]
        cs_bin[cs_idx,ll] = True
        abs_corr = sp.zeros((n_pop,1))
        for pp in range(n_pop):
            abs_corr_mat = abs(ld[sp.ix_(cs_idx,cs_idx,[pp])])
            abs_corr[pp] = sp.amin(abs_corr_mat)

        cs_purity[ll] = sp.amax(abs_corr)

    flt_in = cs_purity>=min_purity
    for ll in range(n_cs):
        if sp.amin(pval_min[cs_bin[:,ll]]) > pval_thresh:
            flt_in[ll] = False

    flt_in = sp.where(flt_in)[0]
    if len(flt_in) == 0:
        n_cs = 0
        alpha = sp.nan
        cs_bin = sp.nan
        cs_purity = sp.nan
        pip = sp.nan

        return n_cs, alpha, cs_bin, cs_purity, pip
    else:
        alpha = alpha[:,flt_in]
        cs_bin = cs_bin[:,flt_in]
        cs_purity = cs_purity[flt_in]

    (cs_bin,idx) = sp.unique(cs_bin,axis=1,return_index=True)
    n_cs = len(idx)
    alpha = alpha[:,idx]
    cs_purity = cs_purity[idx]
    pip = 1-sp.prod(1-alpha,axis=1).reshape(n_snp,1)

    return n_cs, alpha, cs_bin, cs_purity, pip


