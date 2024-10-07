/*
 * data.cpp
 *
 *  Created on: Dec 27, 2021
 *      Author: kyuan
 */

# include "data.hpp"

# include <fstream>
# include <sstream>
# include <cstdlib>
# include <map>
# include <algorithm>
# include <cmath>
# include <sys/types.h>
# include <sys/stat.h>
# include <unistd.h>
# include <fcntl.h>
# include <cstdio>

int npop;

softpar snp::par;

long FdGetFileSize(int fd)
{
    struct stat stat_buf;
    int rc = fstat(fd, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

void dataset::load(const softpar & par)
{
	//Configure the comparison functions
	if(par.key_by == 0)
	/*
	 * key by rs ID
	 * suitable for array data
	 */
	{
		var_cmp_func = cmp_rs;
		var_eql_func = eql_rs;
	}
	else if(par.key_by == 1)
	/*
	 * key by pos & ref
	 * suitable for sequencing data or imputed data with same reference panel
	 */
	{
		var_cmp_func = cmp_pos_ref;
		var_eql_func = eql_pos_ref;
	}
	else if(par.key_by == 2)
	/*
	 * key by pos & rs
	 * more complex scenario
	 */
	{
		var_cmp_func = cmp_pos_rs;
		var_eql_func = eql_pos_rs;
	}

	//Compute LD if it is not available
	if(par.precmp == false)
		for(int i = 0 ; i < npop ; ++i)
		{
			std::cout << "... calculate LD matrix: " << par.ld_file[i] << \
					".ld ..." << std::endl;

			std::string path;
			path = par.ref_file[i] + ".bim";
			std::ifstream fpibim(path.c_str());
			path = par.ld_file[i] + ".snp";
			std::ofstream fposnp(path.c_str());
			std::string chr, var;
			long bp;
			std::string line;
			while(getline(fpibim, line))
			{
				std::istringstream cl(line);
				cl >> chr >> var >> line >> bp;
				if(chr == par.chr && bp >= par.start && bp <= par.end)
					fposnp << var << std::endl;
			}

			std::ostringstream cmd;
			cmd << par.plink << " --bfile " << par.ref_file[i] << \
                    " --threads " << par.nthreads << \
                    " --memory " << par.plink_mem << \
					" --keep-allele-order --chr " << par.chr << \
					" --extract " << par.ld_file[i] << ".snp --maf " << \
					par.maf << " --make-bed --out " << par.ld_file[i] << "_ref";
			system(cmd.str().c_str());

			cmd.str("");
			cmd << par.plink << " --bfile " << par.ld_file[i] << \
					"_ref --keep-allele-order --r square bin4 --out " << par.ld_file[i] << \
                    " --threads " << par.nthreads << \
                    " --memory " << par.plink_mem;
			system(cmd.str().c_str());

			cmd.str("");
			cmd << par.plink << " --bfile " << par.ld_file[i] << \
					"_ref --keep-allele-order --freq --out " << par.ld_file[i] + "_frq" << \
                    " --threads " << par.nthreads << \
                    " --memory " << par.plink_mem;
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << ".snp";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << "_ref.bed";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << "_ref.fam";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << "_ref.log";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << "_frq.log";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm " << par.ld_file[i] << ".log";
			system(cmd.str().c_str());
			cmd.str("");
			cmd << "rm -rf " << par.ld_file[i] << "*nosex";
			system(cmd.str().c_str());
		}

	std::vector<std::vector<ldref> > vars;
	std::vector<std::vector<sumstats> > sum;
	vars.resize(npop);
	sum.resize(npop);
	memset(tau_sq, 0, sizeof(double) * npop);
	//Load LD variant information from reference bim and frequency files
	for(int i = 0 ; i < npop; ++i)
	{
		std::string path;
		path = par.ld_file[i] + "_ref.bim";
		std::ifstream fpref(path.c_str());
		path = par.ld_file[i] + "_frq.frq";
		std::ifstream fpfrq(path.c_str());

		std::cout << "... parse reference file: " << par.ld_file[i] << \
				"_ref.bim ..." << std::endl;

		std::vector<ldref> & curVars(vars[i]);
		std::string line, line2, tmp;
		ldref t;
		varCoord frqCoord;
		t.idx = 0;
		getline(fpfrq, line2);
		while(getline(fpref, line))
		{
			std::istringstream cl(line);
			std::string curChr;
			cl >> curChr >> t.coord.id >> tmp >> t.coord.pos >> t.coord.a2 >> t.coord.a1;
			if(curChr != par.chr)
			{
				std::cerr << "Error: Line " << t.idx + 1 << " in reference file: " << \
						par.ld_file[i] << "_ref.bim contain variants not in \"" << \
						par.chr << "\"" << std::endl;
				std::cerr << "\t" << line << std::endl;
				exit(1);
			}
			getline(fpfrq, line2);
			std::istringstream cl2(line2);
			cl2 >> curChr >> frqCoord.id >> frqCoord.a2 >> frqCoord.a1 >> t.frq;
			if(curChr != par.chr)
			{
				std::cerr << "Error: Line " << t.idx + 1 << " in frequency file: " << \
						par.ld_file[i] << "_frq.frq contain variants not in \"" << \
						par.chr << "\"" << std::endl;
				std::cerr << "\t" << line2 << std::endl;
				exit(1);
			}
			if(t.coord.id != frqCoord.id || t.coord.a1 != frqCoord.a1 || \
					t.coord.a2 != frqCoord.a2)
			{
				std::cerr << "Error: Inconsistent variant information in line " << \
						t.idx + 1 << " between \n" << \
						"reference file: " << par.ld_file[i] << "_ref.bim\n" << \
						'\t' << line << \
						"frequency file: " << par.ld_file[i] << "_frq.frq\n" << \
						'\t' << line2 << std::endl;
				exit(1);
			}
			curVars.push_back(t);
			++t.idx;
		}
		std::cout << "... " << curVars.size() << " SNPs in the fine-mapping region read from " << \
				par.ld_file[i] << "_ref.bim ..." << std::endl;

		//Sort and check redundancy
		sort(curVars.begin(), curVars.end());
		int nvar(curVars.size());
		if(nvar == 0)
		{
			std::cerr << "Error: No information in reference and frequency files: \n" << \
					"reference file: " << par.ld_file[i] << "_ref.bim\n" << \
					"frequency file: " << par.ld_file[i] << "_frq.frq" << std::endl;
			exit(1);
		}
		for(int j = 1 ; j < nvar; ++j)
			if(curVars[j - 1] == curVars[j])
			{
				std::cerr << "Error: Duplicated information in reference and frequency files: \n" << \
						"reference file: " << par.ld_file[i] << "_ref.bim\n" << \
						"frequency file: " << par.ld_file[i] << "_frq.frq\n" << \
						"Line " << curVars[j - 1].idx + 1 << ":\n\t" << curVars[j - 1].coord.id \
						<< '\t' << curVars[j - 1].coord.pos << '\t' << \
						curVars[j - 1].coord.a1 << '\t' << curVars[j - 1].coord.a2 << '\n' << \
						"Line " << curVars[j].idx + 1 << ":\n" << \
						'\t' << curVars[j].coord.id << '\t' << curVars[j].coord.pos << '\t' << \
						curVars[j].coord.a1 << '\t' << curVars[j].coord.a2 << std::endl;
				exit(1);
			}

		std::vector<sumstats> &curSum(sum[i]);
		std::ifstream fpsum(par.sst_file[i].c_str());
		getline(fpsum, line);
		std::vector<std::string> buff;
		std::string info;
		int start(0);
		for(int j = 1 ; j < line.size(); ++j)
			if(line[j] == '\t')
			{
				buff.push_back(line.substr(start, j - start));
				start = j + 1;
			}
		buff.push_back(line.substr(start));
		int cchr(par.chr_col[i] - 1), csnp(par.snp_col[i] - 1), cpos(par.bp_col[i] - 1), \
				ca1(par.a1_col[i] - 1), ca2(par.a2_col[i] - 1), ceff(par.eff_col[i] - 1), \
				cse(par.se_col[i] - 1), cp(par.pval_col[i] - 1);
		int tot(buff.size());
		if(cchr >= tot || csnp >= tot || cpos >= tot || ca1 >= tot || ca2 >= tot || ceff >= tot \
				|| cse >= tot || cp >= tot)
		{
			std::cerr << "Error: No sufficient tokens in summary statistics file \n\t\"" << \
					par.sst_file[i] << '\"' << std::endl;
			std::cerr << tot << " tokens in header." << std::endl;
			std::cerr << line << std::endl;
			if(cchr >= tot)
				std::cerr << "chr: " << cchr + 1 << std::endl;
			if(csnp >= tot)
				std::cerr << "snp: " << csnp + 1 << std::endl;
			if(cpos >= tot)
				std::cerr << "bp: " << cpos + 1 << std::endl;
			if(ca1 >= tot)
				std::cerr << "a1: " << ca1 + 1 << std::endl;
			if(ca2 >= tot)
				std::cerr << "a2: " << ca2 + 1 << std::endl;
			if(ceff >= tot)
				std::cerr << "eff: " << ceff + 1 << std::endl;
			if(cse >= tot)
				std::cerr << "se: " << cse + 1 << std::endl;
			if(cp >= tot)
				std::cerr << "pval: " << cp + 1 << std::endl;
			exit(1);
		}
		bool isbeta(false);
		std::transform(buff[ceff].begin(), buff[ceff].end(), buff[ceff].begin(), ::toupper);
		if(buff[ceff].substr(0,4) == "BETA")
			isbeta = true;
		else if(buff[ceff].substr(0,2) == "OR")
			isbeta = false;
		else
		{
			std::cerr << "Error: The effect size should be started either with \"BETA\" or \"OR\".\n" \
					<< "\teff: " << ceff + 1 << " \"" << buff[ceff] << '\"' << std::endl;
			exit(1);
		}
		double nsqrt(sqrt(par.n_gwas[i]));
		int count(0);
		while(getline(fpsum, line))
		{
			start = 0;
			int n(0);
			for(int j = 1 ; j < line.size(); ++j)
				if(line[j] == '\t')
				{
					buff[n] = line.substr(start, j - start);
					start = j + 1;
					++n;
					if(n == tot)
					{
						std::cerr << "Error: More than " << tot << " token were found in Line "\
								<< count + 1<< " in summary statistics file \n\t\"" << \
								par.sst_file[i] << '\"' << std::endl;
						std::cerr << line << std::endl;
						exit(1);
					}
				}
			buff[n] = line.substr(start);
			if(n != tot - 1)
			{
				std::cerr << "Error: Less than " << tot << " token were found in Line "\
						<< count + 1 << " in summary statistics file \n\t\"" << \
						par.sst_file[i] << '\"' << std::endl;
				std::cerr << line << std::endl;
				exit(1);
			}
			sumstats cur;
			cur.coord.pos = atol(buff[cpos].c_str());
			if(buff[cchr] == par.chr && cur.coord.pos <= par.end && cur.coord.pos >= par.start)
			{
				cur.coord.id = buff[csnp];
				cur.coord.a1 = buff[ca1];
				cur.coord.a2 = buff[ca2];
				cur.pval = atof(buff[cp].c_str());
				cur.idx = count;
				double eff(atof(buff[ceff].c_str())), se(atof(buff[cse].c_str()));
				if(isbeta)
				{
					cur.beta = eff / se / nsqrt;
					cur.log10p = -log10(erfc(fabs(eff/se) * sqrt(0.5)));
				}
				else
				{
					cur.beta = log(eff) / se / nsqrt;
					cur.log10p = -log10(erfc(fabs(log(eff)/se) * sqrt(0.5)));
				}
				if(std::isnan(cur.beta))
				{
					std::cerr << "Error: Effect size of Line " << count + 1<< \
							" is not a number (NAN)" << std::endl;
					std::cerr << line << std::endl;
					std::cerr << "BETA/OR: " << buff[ceff] << std::endl;
					std::cerr << "SE: " << buff[cse] << std::endl;
					exit(1);
				}
				curSum.push_back(cur);
			}
			++count;
		}
	    std::cout << "... " << curSum.size() << " SNPs in the fine-mapping region read from " << \
	    		par.sst_file[i] << " ..." << std::endl;

		//Sort and check redundancy
		sort(curSum.begin(), curSum.end());
		int nsum(curSum.size());
		if(nsum == 0)
		{
			std::cerr << "Error: No information in summary statistics file: \n" << \
					par.sst_file[i] << std::endl;
			exit(1);
		}
		for(int j = 1 ; j < nsum; ++j)
		{
			if(curSum[j - 1] == curSum[j])
			{
				std::cerr << "Error: Duplicated information in summary statistics file: \n" << \
						par.sst_file[i] << '\n' << \
						"Line " << curSum[j - 1].idx << ":\n" << \
						'\t' << curSum[j - 1].coord.id << '\t' << curSum[j - 1].coord.pos << '\t' << \
						curSum[j - 1].coord.a1 << '\t' << curSum[j - 1].coord.a2 << '\n' << \
						"Line " << curSum[j].idx << ":\n" << \
						'\t' << curSum[j].coord.id << '\t' << curSum[j].coord.pos << '\t' << \
						curSum[j].coord.a1 << '\t' << curSum[j].coord.a2 << std::endl;
				exit(1);
			}
		}
	}


	std::vector<int> idx, idxSum, size, sizeSum;
	idx.resize(npop, 0);
	idxSum.resize(npop, 0);
	size.resize(npop, 0);
	sizeSum.resize(npop, 0);
	for(int i = 0 ; i < npop; ++i)
	{
		size[i] = vars[i].size();
		sizeSum[i] = sum[i].size();
	}

	std::vector<std::vector<int> > align;
	align.resize(npop);
	for(int i = 0 ; i < npop ; ++i)
		align[i].resize(vars[i].size(), -1);

	nsnp = 0;
	while(1)
	{
		int curRef(-1);
		for(int i = 0 ; i < npop; ++i)
			if(idx[i] < size[i])
			{
				if(curRef < 0)
					curRef = i;
				else if(vars[i][idx[i]] < vars[curRef][idx[curRef]])
					curRef = i;
			}
		//std::cout << idx[0] << '\t' << size[0] << std::endl;
		//std::cout << curRef << std::endl;
		if(curRef < 0)
			break;

		ldref & ref(vars[curRef][idx[curRef]]);
		//std::cout << "alignsize\t";
		//std::cout << align[curRef].size() << std::endl;
		//std::cout << idxMks << std::endl;
		align[curRef][vars[curRef][idx[curRef]].idx] = nsnp;
		//std::cout << "After align" << std::endl;
		++idx[curRef];

		mks.push_back(snp());
		snp & t(mks.back());
		//std::cout << "Before set" << std::endl;

		t.ref_set(curRef, ref);
		//std::cout << "After set" << std::endl;
		for(int i = 0 ; i < npop; ++i)
		{
			//std::cout << "idx\t" << idx[i] << std::endl;
			if(i != curRef && idx[i] < size[i] && var_eql_func(vars[i][idx[i]].coord, t.coord))
			{
				t.set(i, vars[i][idx[i]]);
				align[i][vars[i][idx[i]].idx] = nsnp;
				++idx[i];
			}
			//std::cout << "idx\t" << idx[i] << std::endl;
			//std::cout << "idxSum\t" << idxSum[i] << std::endl;
			while(idxSum[i] < sizeSum[i] && var_cmp_func(sum[i][idxSum[i]].coord, t.coord))
				++idxSum[i];
			//std::cout << "idxSum\t" << idxSum[i] << std::endl;
			if(idxSum[i] < sizeSum[i])
			{
				if(var_eql_func(sum[i][idxSum[i]].coord, t.coord))
					t.set_beta(i, sum[i][idxSum[i]]);
				else if(idxSum[i] > 0)
					--idxSum[i];
			}
			//std::cout << sum[i][idxSum[i]].coord.id << '\t';
			//std::cout << idxSum[i] << std::endl;
		}
		//std::cout << t.coord.id << std::endl;
		//std::cout << t.var() << std::endl;

		if(t.var())
			++nsnp;
		else
		{
			for(int i = 0 ; i < npop; ++i)
			{
				if(t.idxs[i * 3] & 1)
				{
					align[i][vars[i][idx[i] - 1].idx] = -1;
				}
			}
		}
	}
//	for(int i = 0 ; i < npop; ++i)
//		for(int j = 0 ; j < vars[i].size(); ++j)
//			if(align[i][j] == 133)
//				std::cout << i << '\t' << j << '\t' << align[i][j] << std::endl;
//	for(int i = 60 ; i < 70; ++i)
//		std::cout << i << '\t' << align[1][i] << std::endl;
	std::cout << "... " << nsnp << " common SNPs in the reference and sumstats ..." << std::endl;
	beta = new double[npop * nsnp];
	pval = new double[npop * nsnp];
	mkIdx = new int[nsnp];
	ind = new char[npop * nsnp];
	double *curBeta(beta), *curP(pval);
	char *curInd(ind);
	for(int i = 0 ; i < npop; ++i)
	{
		tau_sq[i] = 0;
		for(int j = 0 ; j < mks.size(); ++j)
		{
			if(mks[j].var())
			{
				if(mks[j].var(i))
				{
					*curBeta = mks[j].stats[i * 4 + 1];
					*curP = mks[j].stats[i * 4 + 2];
					*curInd = 1;
				}
				else
				{
					*curBeta = 0;
					*curP = 1;
					*curInd = 0;
				}
				double beta2((*curBeta) * (*curBeta));
				if(beta2 > tau_sq[i])
					tau_sq[i] = beta2;
				++curBeta;
				++curP;
				++curInd;
			}
		}
	}
	int *curMkIdx(mkIdx);
	for(int i = 0 ; i < mks.size(); ++i)
		if(mks[i].var())
		{
			*curMkIdx = i;
			++curMkIdx;
		}
	//std::cout << count << '\t' << mks.size() << '\t' << nsnp << std::endl;

//	for(int i = 0 ; i < mks.size(); ++i)
//	{
//		std::cout << mks[i].coord.id << '\t' << mks[i].idxs[0] << '\t' << mks[i].idxs[3] << '\t' << mks[i].var(0) << '\t' << mks[i].var(1) << '\t' << mks[i].ldflag(0) << '\t' << mks[i].ldflag(1) << std::endl;
//		std::cout << '\t' << mks[i].stats[1] << '\t' << mks[i].stats[5] << std::endl;
//	}

	const int MAX = 4096 * 256;
	const int LEN = MAX / sizeof(LDTYPE);
	char buff[MAX];
	//std::cout << "Done" << std::endl;
	//std::ofstream fpflag("Flag.txt");
	for(int i = 0 ; i < npop; ++i)
	{
		std::string path;
		path = par.ld_file[i] + ".ld.bin";
		std::cout << "... parse LD file: " << path << " ..." << std::endl;
		int fd;
		fd = open(path.c_str(), O_RDONLY);


		std::vector<int> & curAlign(align[i]);
		long long nvar(curAlign.size());
		long long len = FdGetFileSize(fd);
		if(nvar * nvar * 4 != len)
		{
			std::cerr << "Error: Wrong LD file size: \n" << \
					path << '\n' << \
					"Expect: " << nvar * nvar * 4 << "Byte\n" << \
					"Detect: " << len << "Byte" << std::endl;
			exit(1);
		}
		ld[i] = new LDTYPE*[nsnp];
		for(int j = 0 ; j < nsnp ; ++j)
		{
			ld[i][j] = new LDTYPE [nsnp];
			LDTYPE *curLD(ld[i][j]);
			memset(curLD, 0, sizeof(LDTYPE) * nsnp);
			curLD[j] = 1;
		}
		//std::cout << "nsnp " << nsnp << std::endl;
		//std::cout << "nvar " << nvar << std::endl;
		int p(0);
		read(fd, buff, MAX);
		LDTYPE * ldbuff;
		ldbuff = (LDTYPE*)buff;
		//std::vector<int> ldflag;
		//ldflag.resize(nsnp, 0);
		//for(int j = 0 ; j < nvar; ++j)
		//{
			//std::cout << j << std::endl;
			//std::cout << i << '\t' << j << '\t' << curAlign[j] << '\t' << mks[mkIdx[j]].ldflag(i) << std::endl;
			//std::cout << nvar << '\t' << curAlign.size() << '\t' << mkIdx[j] << '\t' << mks.size() << std::endl;
			//std::cout << nsnp << std::endl;
			//if(curAlign[j] >= 0)
			//	ldflag[curAlign[j]] = mks[mkIdx[curAlign[j]]].ldflag(i);
		//}
		//for(int j = 0 ; j < nsnp; ++j)
		//	fpflag << i << '\t' << j << '\t' << mks[mkIdx[j]].coord.id << '\t' << ldflag[j] << std::endl;
		//std::cout << ldflag[416] << std::endl;
		//std::cout << "done" << std::endl;
		for(int j = 0 ; j < nvar; ++j)
		{
			//std::cout << "Align:\t" <<  j << '\t' << curAlign[j] << std::endl;
			//if(curAlign[j] == 133)
			//	std::cout << "Found" << std::endl;
			if(curAlign[j] < 0)
			{
				p += nvar;
				//std::cout << p << std::endl;
				while(p >= LEN)
				{
					read(fd, buff, MAX);
					p -= LEN;
					//std::cout << p << std::endl;
				}
				//std::cout << curAlign[j] << '\t' << p << std::endl;
			}
			else
			{
				LDTYPE *curLD(ld[i][curAlign[j]]);
				for(int k = 0 ; k < nvar; ++k)
				{
					//std::cout << j << '\t' << curAlign[j] << '\t' << k << '\t' << curAlign[k] << std::endl;
					if(curAlign[k] >= 0)
					{
						//std::cout << p << std::endl;
						//std::cout << '\t' << dat[p] << std::endl;
						//std::cout << ld[0][0][0] << std::endl;
						//std::cout << curLD[0] << std::endl;
						//std::cout << curLD << std::endl;
						if(std::isnan(ldbuff[p]))
							curLD[curAlign[k]] = 0;
						else
							curLD[curAlign[k]] = ldbuff[p];
						//curLD[curAlign[k]] = ldflag[curAlign[j]] * ldflag[curAlign[k]] * ldbuff[p];
						//if(i == 1 && curAlign[j] == 133)
						//	std::cout << curAlign[k] << '\t' << ldbuff[p] << '\t' << p << std::endl;
					}
					//std::cout << p << std::endl;
					++p;
					if(p >= LEN)
					{
						read(fd, buff, MAX);
						p -= LEN;
					}
				}
				curLD[curAlign[j]] = 1;
			}
		}
		std::cout << "... (" << vars[i].size() << ',' << vars[i].size() << ") LD matrix read from " \
				<< path << " ..." << std::endl;
	}
//	std::cout << "LD load" << std::endl;

//	for(int i = 0 ; i < nsnp; ++i)
//		std::cout << mkIdx[i] << '\t' << beta[i] << '\t' << beta[i + nsnp] << std::endl;
/*
	std::ofstream fpold("LD.txt");
	std::ofstream fpob("Beta.txt");
	for(int i = 0 ; i < npop; ++i)
	{
		fpold << i << std::endl;
		for(int j = 0 ; j < nsnp; ++j)
		{
			for(int k = 0 ; k < nsnp; ++k)
			{
				fpold << ld[i][j][k] << '\t';
			}
			fpold << std::endl;
			fpob << mks[mkIdx[j]].coord.id << '\t' << beta[i * nsnp + j] << std::endl;
			//if(j == 416)
			{
			//	std::cout << mks[mkIdx[j]].ldflag(i) << std::endl;
			//	std::cout << mks[mkIdx[j]].coord.id << '\t' << mks[mkIdx[j]].var(i) << '\t' << mks[mkIdx[j]].idxs[i * 3] << std::endl;
			}
		}
	}
	*/
}


