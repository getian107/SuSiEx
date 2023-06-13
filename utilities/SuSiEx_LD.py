#!/usr/bin/env python

"""
SuSiEx_LD: LD calculation for SuSiEx

Usage:

python SuSiEx.py --ref_file=REF_FILE --ld_file=LD_MATRIX_FILE --chr=CHR --bp=BP --plink=PLINK --maf=MAF_THRESHOLD

"""


import sys
import getopt
import subprocess


def parse_param():
    long_opts_list = ['ref_file=', 'ld_file=', 'chr=', 'bp=', 'plink=', 'maf=']

    param_dict = {'ref_file': None, 'ld_file': None, 'chr': None, 'bp': None, 'plink': None, 'maf': 0.005}

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
            elif opt == "--ref_file": param_dict['ref_file'] = arg
            elif opt == "--ld_file": param_dict['ld_file'] = arg
            elif opt == "--chr": param_dict['chr'] = int(arg)
            elif opt == "--bp": param_dict['bp'] = map(int, arg.split(','))
            elif opt == "--plink": param_dict['plink'] = arg
            elif opt == "--maf": param_dict['maf'] = float(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_file'] == None:
        print('* Please provide the directory and filename prefix of the reference panel using --ref_file\n')
        sys.exit(2)
    elif param_dict['ld_file'] == None:
        print('* Please provide the directory and filename prefix of the LD matrix to be computed using --ld_file\n')
        sys.exit(2)
    elif param_dict['chr'] == None:
        print('* Please provide the chromosome code of the fine-mapping region using --chr\n')
        sys.exit(2)
    elif param_dict['bp'] == None or len(param_dict['bp']) != 2:
        print('* Please provide the start and end base pair coordinate of the fine-mapping region using --bp\n')
        sys.exit(2)
    elif param_dict['plink'] == None:
        print('* Please provide the directory and filename of PLINK using --plink\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    print('... calculate LD matrix ...')

    SNP = []
    with open(param_dict['ref_file']+'.bim') as ff:
        for line in ff:
            ll = (line.strip()).split()
            if int(ll[0]) == param_dict['chr'] and int(ll[3]) >= param_dict['bp'][0] and int(ll[3]) <= param_dict['bp'][1]:
                SNP.append(ll[1])

    with open(param_dict['ld_file']+'.snp', 'w') as ff:
        for snp in SNP:
            ff.write('%s\n' % snp)


    cmd = '%s --bfile %s --keep-allele-order --chr %d --extract %s --maf %f --make-bed --out %s' \
        % (param_dict['plink'], param_dict['ref_file'], param_dict['chr'], param_dict['ld_file']+'.snp', param_dict['maf'], param_dict['ld_file']+'_ref')
    subprocess.check_output(cmd, shell=True)

    cmd = '%s --bfile %s --keep-allele-order --r square bin4 --out %s' %(param_dict['plink'], param_dict['ld_file']+'_ref', param_dict['ld_file'])
    subprocess.check_output(cmd, shell=True)

    cmd = '%s --bfile %s --keep-allele-order --freq --out %s' %(param_dict['plink'], param_dict['ld_file']+'_ref', param_dict['ld_file']+'_frq')
    subprocess.check_output(cmd, shell=True)

    subprocess.check_output('rm '+param_dict['ld_file']+'.snp', shell=True)
    subprocess.check_output('rm -rf '+param_dict['ld_file']+'.nosex', shell=True)
    subprocess.check_output('rm -rf '+param_dict['ld_file']+'_*nosex', shell=True)
    subprocess.check_output('rm '+param_dict['ld_file']+'_ref.bed', shell=True)
    subprocess.check_output('rm '+param_dict['ld_file']+'_ref.fam', shell=True)
    subprocess.check_output('rm '+param_dict['ld_file']+'_ref.log', shell=True)
    subprocess.check_output('rm '+param_dict['ld_file']+'_frq.log', shell=True)
    subprocess.check_output('rm '+param_dict['ld_file']+'.log', shell=True)

    print('... Done ...\n')


if __name__ == '__main__':
    main()



