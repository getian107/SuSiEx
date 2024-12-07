/*
 * model.cpp
 *
 *  Created on: Dec 27, 2021
 *      Author: kyuan
 */


# include "model.hpp"

# include <cstring>
# include <algorithm>
# include <fstream>
# include <map>
# include <stdint.h>

int csVar::key = 0;
int csVar::ncs = 5;

double csVar::pip_overall(const std::vector<cs>& csset) const
{
	double pip(1);
	for(int i = 0 ; i < ncs; ++i)
		if(!csset[i].fltOut)
			pip *= 1 - alpha[i];
	return 1 - pip;
}

int susiex::susie_sst_xethn()
{
	int flag;
	if(par.mult_step)
	{
		std::cout << "* Mult-step model fitting *" << std::endl;
		init(5);
		for(int i = 1 ; i <= par.max_iter; ++i)
		{
			std::cout << "--- iter-" << i << " ---" << std::endl;
			flag = iter();
			if(flag != 0)
				break;
		}
		if(flag == 1)
			ncs = cal_pip();

		if(flag == 1 && ncs == 5)
		{
			nsig = 11;
			flag = 0;
		}
		else if(flag != 1)
			nsig = 5;

		while(flag != 1)
		{
			--nsig;
			init(nsig);
			for(int i = 1 ; i <= par.max_iter; ++i)
			{
				std::cout << "--- iter-" << i << " ---" << std::endl;
				flag = iter();
				if(flag != 0)
					break;
			}
			if(flag == 1)
				ncs = cal_pip();
		}
	}
	else
	{
		init(par.n_sig);
		for(int i = 1 ; i <= par.max_iter; ++i)
		{
			std::cout << "--- iter-" << i << " ---" << std::endl;
			flag = iter();
			if(flag != 0)
				break;
		}
		if(flag == 1)
			ncs = cal_pip();
	}
	write_cs(flag);
	return 0;
}

const double PI = 3.141592653589793;

double dotSum(float* a, double* b, int len)
{
	double sum(0);
	for(int i = 0 ; i < len ; ++i)
	{
		sum += (*a) * (*b);
		++a;
		++b;
	}
	return sum;
}

double dotSum(double* a, double* b, int len)
{
	double sum(0);
	for(int i = 0 ; i < len ; ++i)
	{
		sum += (*a) * (*b);
		++a;
		++b;
	}
	return sum;
}

int susiex::init(int max_sig)
{
	std::cout << "... fit SuSiEx model with " << max_sig << " signals..." << std::endl;
	if(max_sig > _memsig)
	{
		clean();
		nsig = max_sig;
		alloc_mem();
	}

	nsig = max_sig;

	memset(alpha, 0, sizeof(double) * nsnp * nsig);
	memset(b, 0, sizeof(double) * nsnp * nsig * npop);
	memset(b_sq, 0, sizeof(double) * nsnp * nsig * npop);
	memset(loglik, 0, sizeof(double) * nsig);
	memset(pst_loglik, 0, sizeof(double) * nsig);
	memset(bhat, 0, sizeof(double) * nsnp * npop);
	memset(beta_ll, 0, sizeof(double) * nsnp * npop);

	for(int i = 0 ; i < npop; ++i)
		sigma_sq[i] = 1;

	elbo_old = -1e20;
	first = true;
	return 0;
}

//-1 err, 0 inter, 1 convergent
int susiex::iter()
{
	int snp_pop(npop * nsnp);
	for(int i = 0 ; i < nsig; ++i)
	{
#pragma omp parallel for
		for(int j = 0 ; j < snp_pop; ++j)
		{
			int snpidx(j % nsnp), popidx(j / nsnp);
			int bidx(i * snp_pop + popidx * nsnp);
			LDTYPE *curLD(dat.ld[popidx][snpidx]);
			bhat[j] -= dotSum(curLD, b + bidx, nsnp);
			beta_ll[j] = beta[j] - bhat[j];
			//std::cout << j / nsnp << '\t' << j % nsnp << '\t' << beta[j] << '\t' << bhat[j] << '\t' << beta_ll[j] << std::endl;
		}
		ser(i);
#pragma omp parallel for
		for(int j = 0 ; j < snp_pop; ++j)
		{
			int snpidx(j % nsnp), popidx(j / nsnp);
			int bidx(i * snp_pop + popidx * nsnp);
			LDTYPE *curLD(dat.ld[popidx][snpidx]);
			bhat[j] += dotSum(curLD, b + bidx, nsnp);
			//std::cout << j / nsnp << '\t' << j % nsnp << '\t' << bhat[j] << std::endl;
		}
	}
	memset(b_sq_sum, 0, sizeof(double) * npop);
	for(int i = 0 ; i < nsig; ++i)
	{
		for(int j = 0 ; j < npop; ++j)
		{
			double *curB_sq(b_sq + i * snp_pop + j * nsnp);
			for(int k = 0 ; k < nsnp; ++k)
			{
				b_sq_sum[j] += *curB_sq;
				++curB_sq;
			}
		}
	}
	for(int i = 0 ; i < npop; ++i)
	{
		double sum_b_beta(0), sum_bDb(0), sum_dig_bDb(0);
		LDTYPE **popLD(dat.ld[i]);
#pragma omp parallel for
		for(int j = 0 ; j < nsnp; ++j)
		{
			LDTYPE *curLD(popLD[j]);
			for(int k = 0 ; k < nsig; ++k)
				bDb[k * nsnp + j] = dotSum(curLD, b + k * snp_pop + i * nsnp, nsnp);
		}
		for(int j = 0; j < nsig; ++j)
		{
			double *curbDb(bDb + j * nsnp);
			for(int k = 0 ; k < nsig; ++k)
			{
				double prd(dotSum(curbDb, b + k * snp_pop + i * nsnp, nsnp));
				if(j == k)
					sum_dig_bDb += prd;
				sum_bDb += prd;
			}
		}
		double *curBeta(beta + i * nsnp);
		for(int j = 0 ; j < nsig; ++j)
			sum_b_beta += dotSum(b + j * snp_pop + i * nsnp, curBeta, nsnp);
		erss[i] = ngwas[i] - 2.0 * ngwas[i] * sum_b_beta + ngwas[i] * sum_bDb - ngwas[i] * sum_dig_bDb + ngwas[i] * b_sq_sum[i];
		//std::cout << ngwas[i] << std::endl;
		//std::cout << b_sq_sum[i] << std::endl;
		//std::cout << "Sum b beta: " << sum_b_beta << std::endl;
		//std::cout << "Sum bDb: " << sum_bDb << std::endl;
		//std::cout << "Sum dig bDb: " << sum_dig_bDb << std::endl;
		//std::cout << "ERSS: " << i << ' ' << erss[i] << std::endl;
		sigma_sq[i] = 1.0 / ngwas[i] * erss[i];
	}
	double elbo_new(0);
	for(int i = 0 ; i < nsig; ++i)
		elbo_new += loglik[i] - pst_loglik[i];
	for(int i = 0 ; i < npop; ++i)
		elbo_new = elbo_new - 0.5 * ngwas[i] * log(2 * PI * sigma_sq[i]) - 1.0/(2.0 * sigma_sq[i]) * erss[i];

	//std::cout << elbo_new << '\t' << elbo_old << '\t' << fabs(elbo_new - elbo_old) << std::endl;
	//for(int i = 0 ; i < nsig ; ++i)
	//	for(int j = 0 ; j < nsnp; ++j)
	//		std::cout << i << '\t' << j << '\t' << b[i * nsnp + j] << std::endl;
	if(first)
	{
		first = false;
		elbo_old = elbo_new;
		return 0;
	}

	if(std::isnan(elbo_new))
		return -1;
	else if(fabs(elbo_new - elbo_old) < par.tol)
		return 1;
	else
	{
		elbo_old = elbo_new;
		return 0;
	}
}

int susiex::ser(int cur_sig)
{
	int pop_snp(npop * nsnp);
	memset(v_sq, 0, sizeof(double) * pop_snp);
	memset(z, 0, sizeof(double) * pop_snp);
	memset(mu, 0, sizeof(double) * pop_snp);
	memset(mu_sq, 0, sizeof(double) * pop_snp);
	memset(phi_sq, 0, sizeof(double) * pop_snp);

	memset(maxPopBF, 0, sizeof(double) * npop);
	memset(w, 0, sizeof(double) * nsnp);

	double *curLogBF(logBF + cur_sig * pop_snp);
#pragma omp parallel for
	for(int i = 0 ; i < pop_snp; ++i)
	{
		int s(i / nsnp), p(i % nsnp);
		if(ind[i])
			v_sq[i] = sigma_sq[s] / ngwas[s];
		else
			v_sq[i] = 1e6;
		z[i] = beta_ll[i] / sqrt(v_sq[i]);
		curLogBF[i] = 0.5 * log(v_sq[i] / (v_sq[i] + tau_sq[s])) + 0.5 * z[i] * z[i] * tau_sq[s] / (v_sq[i] + tau_sq[s]);
		//std::cout << s << '\t' << p << '\t' << v_sq[i] << '\t' << tau_sq[s] << '\t' << z[i] << '\t' << beta_ll[i] << std::endl;

		phi_sq[i] = 1.0 / (1.0 / v_sq[i] + 1.0 / tau_sq[s]);
		mu[i] = phi_sq[i] / v_sq[i] * beta_ll[i];
		mu_sq[i] = phi_sq[i] + mu[i] * mu[i];
		//std::cout << i << '\t' << phi_sq[i] << '\t' << mu[i] << '\t' << mu_sq[i] << std::endl;
#pragma omp critical
		w[p] += curLogBF[i];
	}
	double sumMaxBF(0);
	for(int i = 0 ; i < npop; ++i)
	{
		maxPopBF[i] = curLogBF[i * nsnp];
		for(int j = i * nsnp; j < (i + 1) * nsnp; ++j)
			if(maxPopBF[i] < curLogBF[j])
				maxPopBF[i] = curLogBF[j];
		//std::cout << i << '\t' << maxPopBF[i] << std::endl;
		sumMaxBF += maxPopBF[i];
	}
	//std::cout << sumMaxBF << std::endl;
	double factor(1.0 / nsnp), sumW(0);
#pragma omp parallel for
	for(int i = 0 ; i < nsnp; ++i)
		w[i] = factor * exp(w[i] - sumMaxBF);

	for(int i = 0 ; i < nsnp; ++i)
		sumW += w[i];
		//std::cout << i << '\t' << w[i] << std::endl;

	double *curAlpha(alpha + cur_sig * nsnp);
	memset(curAlpha, 0, sizeof(double) * nsnp);
	for(int i = 0 ; i < nsnp; ++i)
		curAlpha[i] = w[i] / sumW;

	memset(b_sq_sum, 0, sizeof(double) * npop);
	double *curB(b + cur_sig * pop_snp), *curB_sq(b_sq + cur_sig * pop_snp);
#pragma omp parallel for
	for(int i = 0; i < pop_snp; ++i)
	{
		int p(i % nsnp);
		curB[i] = curAlpha[p] * mu[i];
		//std::cout << s << '\t' << p << '\t' << curAlpha[p] << '\t' << curB[i] << std::endl;
		curB_sq[i] = curAlpha[p] * mu_sq[i];
	}
	for(int i = 0; i < pop_snp; ++i)
		b_sq_sum[i / nsnp] += curB_sq[i];

	loglik[cur_sig] = sumMaxBF + log(sumW);
	//std::cout << "Loglik" << std::endl;
	//std::cout << loglik[cur_sig] << '\t' << sumMaxBF << '\t' << log(sumW) << std::endl;
	pst_loglik[cur_sig] = 0;
	for(int i = 0 ; i < npop; ++i)
	{
		loglik[cur_sig] = loglik[cur_sig] - 0.5 * ngwas[i] * log(2 * PI * sigma_sq[i]) - 1.0 / (2.0 * sigma_sq[i]) * ngwas[i];
		//std::cout << log(2 * PI * sigma_sq[i]) << '\t' << 1.0 / (2.0 * sigma_sq[i]) * ngwas[i] << std::endl;
		pst_loglik[cur_sig] = pst_loglik[cur_sig] - 0.5 * ngwas[i] * log(2 * PI * sigma_sq[i]) - 1.0 / (2.0 * sigma_sq[i]) * (ngwas[i] - 2.0 * ngwas[i] * dotSum(curB + i * nsnp, beta_ll + i * nsnp, nsnp) + ngwas[i] * b_sq_sum[i]);
		//std::cout << loglik[cur_sig] << '\t' << pst_loglik[cur_sig] << std::endl;
	}

	return 0;
}

int susiex::cal_pip()
{
	int snp_pop(nsnp * npop), snp_sig(nsnp * nsig);
	for(int i = 0 ; i < snp_sig; ++i) // Thanks Siru for catching this bug. for(int i = 0 ; i < snp_pop; ++i)
		if(std::isnan(alpha[i]))
		{
			ncs = 0;
			nsig = 0;
			return ncs;
		}

	csVar::ncs = nsig;
	csDetails.clear();
	csDetails.resize(nsnp);
	std::vector<double> sumPIP;
	sumPIP.resize(nsig, 0);
	for(int i = 0 ; i < nsig ; ++i)
	{
		for(int j = 0 ; j < nsnp; ++j)
		{
			csDetails[j].alpha[i] = alpha[i * nsnp + j];
			sumPIP[i] += csDetails[j].alpha[i];
		}
	}
	for(int i = 0 ; i < nsnp; ++i)
		csDetails[i].idx = i;
	csset.clear();
	csset.resize(nsig);
	//std::cout << par.level << std::endl;
	for(int i = 0 ; i < nsig; ++i)
	{
		csVar::key = i;
		sort(csDetails.begin(), csDetails.end(), cmp_csVar_alpha);
		double acc(0);
		for(int j = 0 ; j < nsnp; ++j)
			csDetails[j].alpha[i] /= sumPIP[i];
		for(int j = 0 ; j < nsnp; ++j)
		{
			double cur(csDetails[j].alpha[i]);
			acc += cur;
			csset[i].idx.push_back(csDetails[j].idx);
			csset[i].alpha.push_back(cur);
			if(acc > par.level)
				break;
		}
		std::vector<double> abs_corr;
		abs_corr.resize(npop, 1);
		int cssize(csset[i].idx.size());
		csset[i].minP = 1;
		csset[i].minPalPop.resize(npop, 1);
		for(int j = 0 ; j < npop; ++j)
		{
			for(int k = 0 ; k < cssize; ++k)
			{
				LDTYPE *curLD(dat.ld[j][csset[i].idx[k]]);
				for(int l = 0 ; l < cssize; ++l)
				{
					double d(fabs(curLD[csset[i].idx[l]]));
					if(abs_corr[j] > d)
						abs_corr[j] = d;
				}
				double curP(dat.pval[j * nsnp + csset[i].idx[k]]);
				if(curP < csset[i].minP)
					csset[i].minP = curP;
				if(curP < csset[i].minPalPop[j])
					csset[i].minPalPop[j] = curP;
			}
		}
		csset[i].purity = abs_corr[0];
		for(int j = 1 ; j < npop; ++j)
			if(csset[i].purity < abs_corr[j])
				csset[i].purity = abs_corr[j];
	}
	int ncs(0);
	for(int i = 0 ; i < nsig; ++i)
		if(csset[i].purity < par.min_purity || csset[i].minP > par.pth)
			csset[i].fltOut = true;

	unique(csset);

	for(int i = 0 ; i < nsig; ++i)
		if(!csset[i].fltOut)
			++ncs;

	memset(pip, 0, sizeof(double) * nsnp);
	if(ncs)
		for(int i = 0 ; i < nsnp; ++i)
			pip[csDetails[i].idx] = csDetails[i].pip_overall(csset);

	sort(csDetails.begin(), csDetails.end(), cmp_csVar_key);
	return ncs;
}

int susiex::write_cs(int convergent)
{
	std::vector<int> cs2id;
	if(convergent == 1)
	{
		int id(1);
		for(int i = 0 ; i < nsig ; ++i)
			if(!csset[i].fltOut)
			{
				cs2id.push_back(id);
				++id;
			}
			else
				cs2id.push_back(0);
	}
	else
	{
		nsig = 0;
		ncs = 0;
	}
	std::cout << "... write credible sets to file: " << par.out_dir << '/' << par.out_name << ".summary/.cs/.snp ..." << std::endl;
	std::string path;
	path = par.out_dir + "/" + par.out_name + ".snp";
	std::ofstream fpsnp(path.c_str());
	int pop_snp(npop * nsnp);
	if(par.key_by == 0)
	{
		fpsnp << "SNP";
        for(int i = 1 ; i <= ncs; ++i)
        {
        	fpsnp << "\tPIP(CS" << i << ")";
        	for(int j = 1; j <= npop; ++j)
        		fpsnp << "\tLogBF(CS" << i << ",Pop" << j << ")";
        }
		fpsnp << std::endl;
        for(int i = 0 ; i < nsnp; ++i)
        {
        	snp & cur(dat.mks[dat.mkIdx[i]]);
			if(cur.var())
			{
				fpsnp << cur.coord.id;
				for(int j = 0 ; j < nsig; ++j)
				{
					if(!csset[j].fltOut)
					{
						fpsnp << '\t' << csDetails[i].alpha[j];
						for(int k = 0 ; k < npop; ++k)
							fpsnp << '\t' << logBF[j * pop_snp + k * nsnp + i];
					}
				}
				fpsnp << std::endl;
			}
        }
	}
	else if(par.key_by == 1)
	{
		fpsnp << "BP\tA1\tA2";
        for(int i = 1 ; i <= ncs; ++i)
        {
        	fpsnp << "\tPIP(CS" << i << ")";
        	for(int j = 1; j <= npop; ++j)
        		fpsnp << "\tLogBF(CS" << i << ",Pop" << j << ")";
        }
		fpsnp << std::endl;
        for(int i = 0 ; i < nsnp; ++i)
        {
        	snp & cur(dat.mks[dat.mkIdx[i]]);
			if(cur.var())
			{
				fpsnp << cur.coord.pos << '\t' << cur.coord.a1 << '\t' << cur.coord.a2;
				for(int j = 0 ; j < nsig; ++j)
				{
					if(!csset[j].fltOut)
					{
						fpsnp << '\t' << csDetails[i].alpha[j];
						for(int k = 0 ; k < npop; ++k)
							fpsnp << '\t' << logBF[j * pop_snp + k * nsnp + i];
					}
				}
				fpsnp << std::endl;
			}
        }
	}
	else if(par.key_by == 2)
	{
        fpsnp << "BP\tSNP";
        for(int i = 1 ; i <= ncs; ++i)
        {
        	fpsnp << "\tPIP(CS" << i << ")";
        	for(int j = 1; j <= npop; ++j)
        		fpsnp << "\tLogBF(CS" << i << ",Pop" << j << ")";
        }
        fpsnp << std::endl;
        for(int i = 0 ; i < nsnp; ++i)
        {
        	snp & cur(dat.mks[dat.mkIdx[i]]);
        	if(cur.var())
        	{
        		fpsnp << cur.coord.pos << '\t' << cur.coord.id;
        		for(int j = 0 ; j < nsig; ++j)
        		{
					if(!csset[j].fltOut)
					{
						fpsnp << '\t' << csDetails[i].alpha[j];
						for(int k = 0 ; k < npop; ++k)
							fpsnp << '\t' << logBF[j * pop_snp + k * nsnp + i];
					}
        		}
        		fpsnp << std::endl;
        	}
        }
	}
	path = par.out_dir + "/" + par.out_name + ".cs";
	std::ofstream fpcs(path.c_str());
	path = par.out_dir + "/" + par.out_name + ".summary";
	std::ofstream fpsum(path.c_str());

	fpsum << "# chr" << par.chr << ':' << par.start << '-' << par.end << std::endl;

	snp::par = par;
	const int LIMT = 500;

	if(convergent == 1)
	{
		if(ncs == 0)
		{
			fpcs << "NULL" << std::endl;
			fpsum << "NULL" << std::endl;
		}
		else
		{
			uint64_t ncon(1);
			ncon = ncon << npop;
			fpcs << "CS_ID\tSNP\tBP\tREF_ALLELE\tALT_ALLELE\tREF_FRQ\tBETA\tSE\t-LOG10P\tCS_PIP\tOVRL_PIP" << std::endl;
			fpsum << "CS_ID\tCS_LENGTH\tCS_PURITY\tMAX_PIP_SNP\tBP\tREF_ALLELE\tALT_ALLELE\tREF_FRQ\tBETA\tSE\t-LOG10P\tMAX_PIP";
			if(npop > 1 && npop <= 64)
			{
				for(int i = 1 ; i <= npop; ++i)
					fpsum << "\tPOST-HOC_PROB_POP" << i;
				//for(int i = 1 ; i <= npop; ++i)
				//	fpsum << "\tBestPvalPop" << i;
				//for(uint64_t i = 0 ; i < ncon; ++i)
				//	fpsum << "\tConf" << i;
			}
			fpsum << std::endl;
			if(npop > 64)
				std::cout << "Warning: the post-hoc configuration analysis will not be performed for more than 64 populations." << std::endl;

			std::vector<uint64_t> popIdx;
			uint64_t t(1);
			for(int i = 0 ; i < npop; ++i)
			{
				popIdx.push_back(t);
				t = t << 1;
			}
			std::vector<std::vector<double> > tmpLogBF;
			tmpLogBF.resize(ncon);
			for(int i = 0 ; i < ncon; ++i)
				tmpLogBF[i].resize(nsnp, 0);
			std::vector<long double> prob;
			prob.resize(ncon, 0);

			for(int i = 0 ; i < nsig ; ++i)
			{
				csset[i].get_max();
				if(!csset[i].fltOut)
				{
					for(int j = 0 ; j < csset[i].idx.size(); ++j)
					{
						snp & cur(dat.mks[dat.mkIdx[csset[i].idx[j]]]);
						fpcs << cs2id[i] << '\t' << cur << '\t' << csset[i].alpha[j] << '\t' << pip[csset[i].idx[j]] << std::endl;
					}
					fpsum << cs2id[i] << '\t' << csset[i].idx.size() << '\t' << csset[i].purity << '\t' << dat.mks[dat.mkIdx[csset[i].idx[csset[i].max]]] << '\t' << csset[i].alpha[csset[i].max];
					if(npop > 1 && npop <= 64)
					{
						for(uint64_t curCon = 0; curCon < ncon; ++curCon)
						{
							std::vector<double>& curTmpLogBF(tmpLogBF[curCon]);
							std::fill(curTmpLogBF.begin(), curTmpLogBF.end(), 0);
							for(int j = 0 ; j < npop; ++j)
							{
								if(curCon & popIdx[j])
								{
									double *curLogBF(logBF + i * pop_snp + j * nsnp);
									for(int k = 0 ; k < nsnp; ++k)
										curTmpLogBF[k] += curLogBF[k];
								}
							}
							sort(curTmpLogBF.begin(), curTmpLogBF.end());
						}
						int max(tmpLogBF[0].back());
						for(int j = 1; j < ncon; ++j)
							if(max < tmpLogBF[j].back())
								max = tmpLogBF[j].back();
						double correct (max - LIMT);

						long double sumProb(0);
						for(uint64_t curCon = 0; curCon < ncon; ++curCon)
						{
							std::vector<double>& curTmpLogBF(tmpLogBF[curCon]);
							long double curProb(0);
							for(int k = 0 ; k < nsnp; ++k)
								curProb += exp(curTmpLogBF[k] - correct);
							prob[curCon] = curProb;
							sumProb += curProb;
						}
						std::vector<double> popProb;
						popProb.resize(npop, 0);
						for(uint64_t curCon = 0; curCon < ncon; ++curCon)
						{
							for(int j = 0 ; j < npop; ++j)
							{
								if(curCon & popIdx[j])
									popProb[j] += prob[curCon] / sumProb;
							}
						}
						for(int k = 0 ; k < npop; ++k)
							fpsum << '\t' << popProb[k];
						//for(int k = 0 ; k < npop; ++k)
						//	if(csset[i].minPalPop[i] == 1)
						//		fpsum << "\tNA";
						//	else
						//		fpsum << '\t' << -log10(csset[i].minPalPop[k]);
						//for(uint64_t k = 0 ; k < ncon; ++k)
						//	fpsum << '\t' << prob[k] / sumProb;
					}
					fpsum << std::endl;
				}
			}
		}
	}
	else
	{
		fpcs << "FAIL" << std::endl;
		fpsum << "FAIL" << std::endl;
	}
	return 0;
}

void unique(std::vector<cs>& csset)
{
	int ncs(csset.size());
	for(int i = 0 ; i < ncs; ++i)
	{
		if(!csset[i].fltOut)
			for(int j = i + 1 ; j < ncs; ++j)
				if(!csset[j].fltOut && csset[i] == csset[j])
					csset[j].fltOut = true;
	}
}
