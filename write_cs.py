#!/usr/bin/env python

"""
Write credible sets

"""


import scipy as sp


def arr2e(xx):
    return ','.join(tuple(map("{:.3e}".format, xx)))

def arr2f(xx):
    return ','.join(tuple(map("{:.3f}".format, xx)))

def list2txt(xx):
    return ','.join(xx)


def write_cs(reg_chr, reg_bp, n_pop, n_cs, snp_dict, beta, ind, pval, logp, n, alpha, cs_bin, cs_purity, pip, out_dir, out_name):

    print('... write credible sets to file: %s ...' % (out_dir + '/' + out_name + '.summary/.cs/.snp'))

    with open(out_dir + '/' + out_name + '.snp', 'w') as ff:
        for snp in snp_dict['SNP']:
            ff.write('%s\n' % snp)

    with open(out_dir + '/' + out_name + '.cs', 'w') as ff, open(out_dir + '/' + out_name + '.summary', 'w') as ss:
        ss.write('# chr%i:%s-%s\n' % (reg_chr, reg_bp[0], reg_bp[1]))

        if sp.isnan(n_cs):
            ff.write('%s' % 'FAIL')
            ss.write('%s' % 'FAIL')
        elif n_cs == 0:
            ff.write('%s' % 'NULL')
            ss.write('%s' % 'NULL')
        else:
            cs_len = sp.sum(cs_bin,axis=0)

            ff.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                ('CS_ID', 'SNP', 'BP', 'REF_ALLELE', 'ALT_ALLELE', 'REF_FRQ', 'BETA', 'SE', '-LOG10P', 'CS_PIP', 'OVRL_PIP'))
            ss.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
                ('CS_ID', 'CS_LENGTH', 'CS_PURITY', 'MAX_PIP_SNP', 'BP', 'REF_ALLELE', 'ALT_ALLELE', 'REF_FRQ', 'BETA', 'SE', '-LOG10P', 'MAX_PIP'))

            for ll in range(n_cs):
                idx = [jj for jj in range(len(pip)) if cs_bin[jj,ll]]
                ind_ll = ind[idx,:]
                snp = [snp_dict['SNP'][jj] for jj in idx]
                bp = [snp_dict['BP'][jj] for jj in idx]
                a1 = snp_dict['A1'][idx,:]
                a2 = snp_dict['A2'][idx,:]
                frq = snp_dict['FRQ'][idx,:]; frq[ind_ll] = 0.5
                eff = beta[idx,:]/sp.sqrt(2*frq*(1-frq)); eff[ind_ll] = sp.nan
                se = 1/sp.sqrt(sp.array(n))/sp.sqrt(2*frq*(1-frq)); se[ind_ll] = sp.nan
                frq[ind_ll] = sp.nan
                p = logp[idx,:]; p[ind_ll] = sp.nan
                cs_pip = alpha[idx,ll]
                ovrl_pip = pip[idx]

                for csid, snp_ll, bp_ll, a1_ll, a2_ll, frq_ll, eff_ll, se_ll, p_ll, cs_pip_ll, ovrl_pip_ll in \
                    zip([ll+1]*cs_len[ll], snp, bp, a1, a2, frq, eff, se, p, cs_pip, ovrl_pip):
                    ff.write('%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n' % \
                    (csid, snp_ll, bp_ll, list2txt(a1_ll), list2txt(a2_ll), arr2f(frq_ll), arr2e(eff_ll), arr2e(se_ll), arr2f(p_ll), cs_pip_ll, ovrl_pip_ll))

                idx_max = sp.argmax(cs_pip)
                ss.write('%i\t%i\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\n' % \
                (ll+1, cs_len[ll], cs_purity[ll], snp[idx_max], bp[idx_max], list2txt(a1[idx_max,:]), list2txt(a2[idx_max,:]), 
                arr2f(frq[idx_max,:]), arr2e(eff[idx_max,:]), arr2e(se[idx_max,:]), arr2f(p[idx_max,:]), cs_pip[idx_max]))


