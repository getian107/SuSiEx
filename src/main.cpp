/*
 * main.cpp
 *
 *  Created on: Dec 31, 2021
 *      Author: kyuan
 */


# include "data.hpp"
# include "model.hpp"

# include <unistd.h>
# include <getopt.h>
# include <fstream>
# include <cstdlib>
# include <map>
# include "omp.h"

static struct option susiexOption[] =
{
	{"sst_file",		required_argument,	NULL,	's'},
	{"n_gwas",			required_argument,	NULL,	'n'},
	{"ref_file",		required_argument,	NULL,	'r'},
	{"ld_file",			required_argument,	NULL,	'l'},
	{"out_dir",			required_argument,	NULL,	'd'},
	{"out_name",		required_argument,	NULL,	'o'},
	{"chr",				required_argument,	NULL,	'c'},
	{"bp",				required_argument,	NULL,	'b'},
	{"chr_col",			required_argument,	NULL,	'C'},
	{"snp_col",			required_argument,	NULL,	'S'},
	{"bp_col",			required_argument,	NULL,	'P'},
	{"a1_col",			required_argument,	NULL,	'A'},
	{"a2_col",			required_argument,	NULL,	'B'},
	{"eff_col",			required_argument,	NULL,	'E'},
	{"se_col",			required_argument,	NULL,	'T'},
	{"pval_col",		required_argument,	NULL,	'V'},
	{"plink",			required_argument,	NULL,	'p'},
	{"keep-ambig",		required_argument,	NULL,	'k'},
	{"maf",				required_argument,	NULL,	'm'},
	{"n_sig",			required_argument,	NULL,	'i'},
	{"level",			required_argument,	NULL,	'L'},
	{"min_purity",		required_argument,	NULL,	'u'},
	{"mult-step",		required_argument,	NULL,	'M'},
	{"pval_thresh",		required_argument,	NULL,	'a'},
	{"max_iter",		required_argument,	NULL,	'I'},
	{"tol",				required_argument,	NULL,	't'},
	{"threads",			required_argument,	NULL,	'X'},
	{"plink_mem",		required_argument,	NULL,	'K'},
	{"help",			no_argument,		NULL,	'h'},
};

void susiexHelp()
{
	std::cout << "SuSiEx" << std::endl;
	std::cout << "Usage: SuSiEx <[options]>" << std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "\t--sst_file      /  -s       SUM_STATS_FILE (required): Full path and filename of the GWAS summary statistics for each population, separated by comma. Each file must contain a header line. The column corresponding to the effect size estimates must have a header of BETA or OR, indicating whether the effect estimates are regression coefficients or odds ratios." << std::endl;
	std::cout << "\t--n_gwas        /  -n       GWAS_SAMPLE_SIZE (required): Sample size of the GWAS for each population, in the order corresponding to the GWAS summary statistics files, separated by comma." << std::endl;
	std::cout << "\t--ref_file      /  -r       REF_FILE (required): Full path and filename prefix of the reference panel for each population in PLINK binary format (.bed/.bim/.fam), in the order corresponding to the GWAS summary statistics files, separated by comma. The format for SNP IDs (e.g., rs IDs or chr:bp:a1:a2) and version of genome build must be consistent across GWAS summary statistics files and reference panels. A curated 1000 Genomes phase 3 reference panle in PLINK format can be download from the MAGMA website: https://ctg.cncr.nl/software/magma." << std::endl;
	std::cout << "\t--ld_file       /  -l       LD_MATRIX_FILE (required): Full path and filename prefix of the LD matrix to be computed from the reference panels for each population, in the order corresponding to the GWAS summary statistics files, separated by comma." << std::endl;
	std::cout << "\t--out_dir       /  -d       OUTPUT_DIR (required): Full path of the output directory." << std::endl;
	std::cout << "\t--out_name      /  -o       OUTPUT_FILENAME (required): Prefix of the output files." << std::endl;
	std::cout << "\t--chr           /  -c       CHR (required): Chromosome of the fine-mapping region." << std::endl;
	std::cout << "\t--bp            /  -b       BP (required): The start and end base pair coordinate of the fine-mapping region, separated by comma." << std::endl;
	std::cout << "\t--chr_col       /  -C       CHR_COL (required): The column number of the chromosome code in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--snp_col       /  -S       SNP_COL (required): The column number of the SNP IDs in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--bp_col        /  -P       BP_COL (required): The column number of the base pair coordinate in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--a1_col        /  -A       A1_COL (required): The column number of the A1 allele (effective allele) in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--a2_col        /  -B       A2_COL (required): The column number of the A2 allele (alternative allele) in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--eff_col       /  -E       EFF_COL (required): The column number of the effect size estimate (beta or odds ratio) in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--se_col        /  -T       SE_COL (required): The column number of the standard error of the effect size estimate in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--pval_col      /  -V       PVAL_COL (required): The column number of the p-value in each GWAS summary statistics file, separated by comma." << std::endl;
	std::cout << "\t--plink         /  -p       PLINK (required): The full path and filename of PLINK." << std::endl;
	std::cout << "\t--keep-ambig    /  -k       KEEP_AMBIGUOUS_SNPS (optional): If False, all ambiguous SNPs (A/T and G/C SNPs) will be removed; if True, ambiguous SNPs whose A1 and A2 can be matched with the A1 and A2 in the reference panel via allele flip will be retained. Default is False." << std::endl;
	std::cout << "\t--maf           /  -m       MAF_THRESHOLD (optional): Minor allele frequency (MAF) threshold applied to the reference panel. Default is 0.005." << std::endl;
	std::cout << "\t--n_sig         /  -i       NUMBER_OF_SIGNALS (optional): Maximum number of causal signals in the SuSiEx model. Default is 5." << std::endl;
	std::cout << "\t--level         /  -L       LEVEL (optional): Coverage level of the credible set. Default is 95%." << std::endl;
	std::cout << "\t--min_purity    /  -u       MINIMUM_PURITY (optional): Minimum purity of the credible set. Credible sets with purity below this specified value will be filtered out. Default is 0.5." << std::endl;
	std::cout << "\t--mult-step     /  -M       MULT_STEP_FITTING (optional): Use the multi-step modeling fitting approach. The model is first fitted using 5 signals. If the algorithm doesn't converge, the maximum number of signals is progressively reduced from 5 to 1 until convergence. If 5 credible sets are identified, the maximum number of signals is increased to 10. If the algorithm doesn't converge with 10 signals, the maximum number of signals is then progressively reduced from 10 to 5 until convergence." << std::endl;
	std::cout << "\t--pval_thresh   /  -a       MARGINAL_PVAL_THRESHOLD (optional): Filtering threshold for the marginal p-value. Credible sets containing no marginal p-value below this specified value will be filtered out. Default is 1e-05." << std::endl;
	std::cout << "\t--max_iter      /  -I       MAXIMUM_ITERATIONS (optional): Maximum number of iterations allowed for the model fitting algorithm. Default is 100." << std::endl;
	std::cout << "\t--tol           /  -t       TOLERANCE (optional): Tolerance for the convergence of the variational algorithm. Default is 1e-04." << std::endl;
	std::cout << "\t--threads       /  -X       N_THREADS (optional): Number of threads for computation. Default is 1." << std::endl;
	std::cout << "\t--plink_mem     /  -K       Maximum memory in MB used by PLINK. Default is 5000" << std::endl;
	std::cout << "\t--help          /  -h       Print this help." << std::endl;
	exit(0);
}

int main(int argc, char **argv)
{
	softpar par;
	if(argc <= 1)
		susiexHelp();
	char optstring[] = "s:n:r:l:d:o:c:b:C:S:P:A:B:E:T:V:p:k:m:i:L:u:M:a:I:t:X:K:h";
	char opt;
	std::map<std::string, bool> str2bool;
	str2bool["True"] = true;
	str2bool["TRUE"] = true;
	str2bool["T"] = true;
	str2bool["False"] = false;
	str2bool["FALSE"] = false;
	str2bool["F"] = false;
	std::map<char, int> parCheck;
	std::map<std::string, std::string> parAll;
	int paridx(-1);
	while((opt = getopt_long_only(argc, argv, optstring, susiexOption, &paridx)) != -1)
	{
		switch (opt)
		{
			case 's': par.split_str(optarg, "sst"); break;
			case 'n': par.split_int(optarg, "n"); break;
			case 'r': par.split_str(optarg, "ref"); break;
			case 'l': par.split_str(optarg, "ld"); break;
			case 'd': par.out_dir = optarg; break;
			case 'o': par.out_name = optarg; break;
			case 'c': par.chr = optarg; break;
			case 'b': par.split_int(optarg, "stend"); break;
			case 'C': par.split_int(optarg, "chr"); break;
			case 'S': par.split_int(optarg, "snp"); break;
			case 'P': par.split_int(optarg, "bp"); break;
			case 'A': par.split_int(optarg, "a1"); break;
			case 'B': par.split_int(optarg, "a2"); break;
			case 'E': par.split_int(optarg, "eff"); break;
			case 'T': par.split_int(optarg, "se"); break;
			case 'V': par.split_int(optarg, "pval"); break;
			case 'p': par.plink = optarg; break;
			case 'k':
				if(str2bool.count(optarg))
					par.keep_ambig = str2bool[optarg];
				else
				{
					std::cerr << "* Cannot identify bool parameter --keep-ambig " << optarg << std::endl;
					std::cerr << "  Please use T, True, TRUE, F, False, or FALSE" << std::endl;
					exit(2);
				}
				break;
			case 'm': par.maf = atof(optarg); break;
			case 'i': par.n_sig = atol(optarg); break;
			case 'L': par.level = atof(optarg); break;
			case 'u': par.min_purity = atof(optarg); break;
			case 'M':
				if(str2bool.count(optarg))
					par.mult_step = str2bool[optarg];
				else
				{
					std::cerr << "* Cannot identify bool parameter --mult-step " << optarg << std::endl;
					std::cerr << "  Please use T, True, TRUE, F, False, or FALSE" << std::endl;
					exit(2);
				}
				break;
			case 'a': par.pth = atof(optarg); break;
			case 'I': par.max_iter = atol(optarg); break;
			case 't': par.tol = atof(optarg); break;
			case 'X': par.nthreads = atol(optarg); break;
			case 'K': par.plink_mem = atol(optarg); break;
			case 'h': susiexHelp(); break;
			case '?': std::cout << "Cannot identify parameter " << opt << std::endl;
						susiexHelp();
		}
		++parCheck[opt];
		if(parCheck[opt] > 1)
		{
			std::cerr << "* Please provide only one argument of option ";
			if(paridx < 0)
				std::cerr << "-" << opt << std::endl;
			else
				std::cerr << "--" << susiexOption[paridx].name << std::endl;
			exit(2);
		}
		if(paridx < 0)
			parAll[(std::string)"-" + opt] = optarg;
		else
			parAll["--" + (std::string)susiexOption[paridx].name] = optarg;
	}
	if(parCheck.count('s') == 0)
	{
		std::cerr << "* Please provide GWAS summary statistics using --sst_file" << std::endl;
        exit(2);
	}
	if(parCheck.count('n') == 0)
	{
		std::cerr << "* Please provide the sample size of each GWAS using --n_gwas" << std::endl;
        exit(2);
	}
	if(parCheck.count('l') == 0)
	{
		std::cerr << "* Please provide the directory and filename prefix of the LD matrix for each GWAS summary statistics file using --ld_file" << std::endl;
        exit(2);
	}
	if(parCheck.count('d') == 0)
	{
		std::cerr << "* Please specify the output directory using --out_dir" << std::endl;
        exit(2);
	}
	if(parCheck.count('o') == 0)
	{
		std::cerr << "* Please specify the prefix of the output filenames using --out_name" << std::endl;
        exit(2);
	}
	if(par.sst_file.size() != par.n_gwas.size() || par.sst_file.size() != par.ld_file.size())
	{
		std::cerr << "* Length of sst_file, n_gwas and ld_file does not match" << std::endl;
        exit(2);
	}
	npop = par.sst_file.size();
	par.precmp = true;
	for(int i = 0 ; i < npop ; ++i)
	{
		std::string path;
		path = par.ld_file[i] + ".ld.bin";
		std::ifstream fpld(path.c_str());
		path = par.ld_file[i] + "_frq.frq";
		std::ifstream fpfrq(path.c_str());
		path = par.ld_file[i] + "_ref.bim";
		std::ifstream fpref(path.c_str());
		if(!fpld || !fpfrq || !fpref)
		{
			par.precmp = false;
			break;
		}
	}
	if(!par.precmp && (parCheck.count('r') == 0 || par.ref_file.size() != npop))
	{
		std::cerr << "* Please provide precomputed LD and frequency files or a reference panel for each GWAS using --ref-file" << std::endl;
        exit(2);
	}
	if(parCheck.count('c') == 0)
	{
		std::cerr << "* Please provide the chromosome code of the fine-mapping region using --chr" << std::endl;
        exit(2);
	}
	if(parCheck.count('b') == 0)
	{
		std::cerr << "* Please provide the start and end base pair coordinate of the fine-mapping region using --bp" << std::endl;
        exit(2);
	}
	if(parCheck.count('C') == 0 || par.chr_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the chromosome code in each GWAS summary statistics file using --chr_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('S') == 0 || par.snp_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the SNP IDs in each GWAS summary statistics file using --snp_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('P') == 0 || par.bp_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the base pair coordinate in each GWAS summary statistics file using --bp_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('A') == 0 || par.a1_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the A1 allele in each GWAS summary statistics file using --a1_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('B') == 0 || par.a2_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the A2 allele in each GWAS summary statistics file using --a2_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('E') == 0 || par.eff_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the effect size estimate (beta or odds ratio) in each GWAS summary statistics file using --eff_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('T') == 0 || par.se_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the standard error of the effect size estimate (beta or log odds ratio) in each GWAS summary statistics file using --col_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('V') == 0 || par.pval_col.size() != npop)
	{
		std::cerr << "* Please provide the column number of the p-value in each GWAS summary statistics file using --pval_col" << std::endl;
        exit(2);
	}
	if(parCheck.count('p') == 0)
	{
		std::cerr << "* Please provide the directory and filename of PLINK using --plink" << std::endl;
        exit(2);
	}

	if(par.keep_ambig == false)
		std::cout << "* All ambiguous SNPs will be removed" << std::endl;
	std::cout << std::endl;

	omp_set_num_threads(par.nthreads);

	par.show();

	dataset dat;
	dat.load(par);

	susiex model(npop, dat.nsnp, dat, par);

	model.susie_sst_xethn();
}

