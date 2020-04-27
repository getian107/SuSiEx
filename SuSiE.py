#!/usr/bin/env python

"""
Posterior inference for the cross-ethnic sum of single effects (SuSiE) regression model
using GWAS summary statistics

"""


import scipy as sp
import SER


def SUSIE_sst_xethn(beta, ind, n, D, tau_sq, l, max_iter, tol):
    print('... fit SuSiEx model ...')

    (p,s) = sp.shape(beta)

    alpha = sp.zeros((p,l))
    b = sp.zeros((p,l,s))
    b_sq = sp.zeros((p,l,s))
    loglik = sp.zeros((l,1))
    pst_loglik = sp.zeros((l,1))

    sigma_sq = sp.ones((s,1))

    bhat = sp.zeros((p,s))
    for ss in range(s):
        bhat[:,ss] = sp.sum(sp.dot(D[:,:,ss],b[:,:,ss]),axis=1)

    beta_ll = sp.zeros((p,s))
    elbo_old = -sp.inf
    flag = False

    for n_iter in range(1,max_iter+1):
        print('--- iter-%i ---' % n_iter)

        for ll in range(l):
            for ss in range(s):
                bhat[:,ss] = bhat[:,ss]-sp.dot(D[:,:,ss],b[:,ll,ss])
                beta_ll[:,ss] = beta[:,ss]-bhat[:,ss]

            alpha[:,ll,None], b[:,ll,:], b_sq[:,ll,:], loglik[ll], pst_loglik[ll] = SER.SER_sst_xethn(beta_ll, ind, n, sigma_sq, tau_sq)

            for ss in range(s):
                bhat[:,ss] = bhat[:,ss]+sp.dot(D[:,:,ss],b[:,ll,ss])

        erss = sp.zeros((s,1))
        for ss in range(s):
            erss[ss] = n[ss]-2.0*n[ss]*sp.sum(sp.dot(b[:,:,ss].T,beta[:,ss]))+n[ss]*sp.sum(sp.dot(sp.dot(b[:,:,ss].T,D[:,:,ss]),b[:,:,ss])) \
                       -n[ss]*sp.sum(sp.diag(sp.dot(sp.dot(b[:,:,ss].T,D[:,:,ss]),b[:,:,ss])))+n[ss]*sp.sum(b_sq[:,:,ss])
            sigma_sq[ss] = 1.0/n[ss]*erss[ss]

        elbo_new = sp.sum(loglik)-sp.sum(pst_loglik)
        for ss in range(s):
            elbo_new = elbo_new-0.5*n[ss]*sp.log(2*sp.pi*sigma_sq[ss])-1.0/(2.0*sigma_sq[ss])*erss[ss]

        if abs(elbo_new-elbo_old) < tol:
            flag = True
            break
        else:
            elbo_old = elbo_new

    return alpha, b, b_sq, sigma_sq, elbo_new, n_iter, flag


