#!/usr/bin/env python

"""
Write credible sets

"""


import scipy as sp


def write_cs(n_cs, snp_dict, alpha, cs_bin, cs_purity, pip, out_dir, out_name):

    print('... write credible sets to file: %s ...' % (out_dir + '/' + out_name + '.cs'))
    with open(out_dir + '/' + out_name + '.cs', 'w') as ff:
        if sp.isnan(n_cs):
            ff.write('%s' % 'FAIL')
        elif n_cs == 0:
            ff.write('%s' % 'NULL')
        else:
            cs_len = sp.sum(cs_bin,axis=0)
            max_pip = sp.zeros((n_cs,1))
            ff.write('%s\t%s\t%s\t%s\n' % ('CS_ID', 'SNP', 'CS_PIP', 'OVRL_PIP'))
            for ll in range(n_cs):
                idx = [jj for jj in range(len(pip)) if cs_bin[jj,ll]]
                snp = [snp_dict['SNP'][jj] for jj in idx]
                cs_pip = [alpha[jj,ll] for jj in idx]
                max_pip[ll] = sp.amax(cs_pip)
                ovrl_pip = [pip[jj] for jj in idx]

                for csid, snp_ll, cs_pip_ll, ovrl_pip_ll in zip([ll+1]*cs_len[ll], snp, cs_pip, ovrl_pip):
                    ff.write('%i\t%s\t%.6f\t%.6f\n' % (csid, snp_ll, cs_pip_ll, ovrl_pip_ll))


    print('... write summary info to file: %s ...' % (out_dir + '/' + out_name + '.summary'))
    with open(out_dir + '/' + out_name + '.summary', 'w') as ff:
        if sp.isnan(n_cs):
            ff.write('%s' % 'FAIL')
        elif n_cs == 0:
            ff.write('%s' % 'NULL')
        else:
            ff.write('%s\t%s\t%s\t%s\n' % ('CS_ID', 'CS_LENGTH', 'CS_PURITY', 'MAX_PIP'))
            for ll in range(n_cs):
                ff.write('%i\t%i\t%.3f\t%.6f\n' % (ll+1, cs_len[ll], cs_purity[ll], max_pip[ll]))


