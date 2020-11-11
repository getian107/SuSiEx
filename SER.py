#!/usr/bin/env python

"""
Posterior inference for the cross-ethnic single effect regression (SER) model
using GWAS summary statistics

"""


import scipy as sp


def SER_sst_xethn(beta, ind, n, sigma_sq, tau_sq):

    (p,s) = sp.shape(beta)

    v_sq = sp.zeros((p,s))
    z = sp.zeros((p,s))
    logBF = sp.zeros((p,s))

    mu = sp.zeros((p,s))
    mu_sq = sp.zeros((p,s))
    phi_sq = sp.zeros((p,s))

    for ss in range(s):
        for jj in range(p):
            if ind[jj,ss] == True:
                v_sq[jj,ss] = 1e6
            else:
                v_sq[jj,ss] = sigma_sq[ss]/n[ss]

            z[jj,ss] = beta[jj,ss]/sp.sqrt(v_sq[jj,ss])
            logBF[jj,ss] = 0.5*sp.log(v_sq[jj,ss]/(v_sq[jj,ss]+tau_sq[ss]))+0.5*z[jj,ss]**2*tau_sq[ss]/(v_sq[jj,ss]+tau_sq[ss])

            phi_sq[jj,ss] = 1.0/(1.0/v_sq[jj,ss]+1.0/tau_sq[ss])
            mu[jj,ss] = phi_sq[jj,ss]/v_sq[jj,ss]*beta[jj,ss]
            mu_sq[jj,ss] = phi_sq[jj,ss]+mu[jj,ss]**2

    w = 1.0/p*sp.exp(sp.sum(logBF,axis=1)-sp.sum(sp.amax(logBF,axis=0)))
    alpha = w.reshape(p,1)/sp.sum(w)

    b = sp.repeat(alpha,s,axis=1)*mu
    b_sq = sp.repeat(alpha,s,axis=1)*mu_sq

    loglik = sp.sum(sp.amax(logBF,axis=0))+sp.log(sp.sum(w))
    for ss in range(s):
        loglik = loglik-0.5*n[ss]*sp.log(2*sp.pi*sigma_sq[ss])-1.0/(2.0*sigma_sq[ss])*n[ss]

    pst_loglik = 0
    for ss in range(s):
        pst_loglik = pst_loglik-0.5*n[ss]*sp.log(2*sp.pi*sigma_sq[ss]) \
                     -1.0/(2.0*sigma_sq[ss])*(n[ss]-2.0*n[ss]*sp.inner(b[:,ss],beta[:,ss])+n[ss]*sp.sum(b_sq[:,ss]))

    return alpha, b, b_sq, loglik, pst_loglik


