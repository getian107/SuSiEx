/*
 * model.hpp
 *
 *  Created on: Dec 27, 2021
 *      Author: kyuan
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

# include "data.hpp"

# include <limits>
# include <cstring>
# include <vector>
# include <algorithm>
# include "omp.h"

class cs;

class csVar
{
public:
	csVar() : idx(0)
	{
		alpha.resize(ncs, 0);
	}

	double pip_overall(const std::vector<cs>& csset) const;

	int idx;
	std::vector<double> alpha;
	static int ncs, key;
};

inline bool cmp_csVar_alpha(const csVar& a, const csVar& b)
{
	return a.alpha[csVar::key] > b.alpha[csVar::key];
}

inline bool cmp_csVar_key(const csVar& a, const csVar&b)
{
	return a.idx < b.idx;
}

class cs
{
public:
	cs() : purity(1), minP(1), max(0), fltOut(false) {}
	std::vector<int> idx;
	std::vector<double> alpha, minPalPop;
	double purity, minP;
	int max;
	bool fltOut;

	void get_max()
	{
		if(idx.size() == 0)
		{
			fltOut = true;
			return;
		}
		max = 0;
		double maxA(alpha[0]);
		for(int i = 1 ; i < alpha.size(); ++i)
			if(alpha[i] > maxA)
			{
				maxA = alpha[i];
				max = i;
			}
	}
};

inline bool operator==(const cs& a, const cs& b)
{
	if(a.idx.size() != b.idx.size())
		return false;
	int n(a.idx.size());
	std::vector<int> ia(a.idx), ib(b.idx);
	sort(ia.begin(), ia.end());
	sort(ib.begin(), ib.end());
	for(int i = 0 ; i < n ; ++i)
		if(ia[i] != ib[i])
			return false;
	return true;
}

class susiex
{
public:
	susiex(int _npop, int _nsnp, dataset& _dat, const softpar& _par) :\
			dat(_dat), par(_par), npop(_npop), nsnp(_nsnp), \
			elbo_old(-std::numeric_limits<double>::infinity()), ncs(0), first(false)
	{
		if(par.mult_step)
			nsig = 5;
		else
			nsig = par.n_sig;

		alloc_mem();

		ngwas = new int[npop];
		for(int i = 0 ; i < npop; ++i)
			ngwas[i] = par.n_gwas[i];

		ind = dat.ind;
		beta = dat.beta;
		tau_sq = dat.tau_sq;

		sigma_sq = new double[npop];
		bhat = new double[npop * nsnp];
		beta_ll = new double[npop * nsnp];
		erss = new double[npop];

		v_sq = new double[npop * nsnp];
		z = new double[npop * nsnp];
		//logBF = new double[npop * nsnp];
		mu = new double[npop * nsnp];
		mu_sq = new double[npop * nsnp];
		phi_sq = new double[npop * nsnp];
		w = new double[nsnp];
		maxPopBF = new double[npop];
		b_sq_sum = new double[npop];
		pip = new double[nsnp];
	}

	~susiex()
	{
		clean();

		delete []ngwas;

		delete []sigma_sq;
		delete []bhat;
		delete []beta_ll;

		delete []v_sq;
		delete []z;
		//delete []logBF;
		delete []mu;
		delete []mu_sq;
		delete []phi_sq;
		delete []w;
		delete []maxPopBF;
		delete []b_sq_sum;
		delete []pip;
	}

	dataset& dat;
	const softpar& par;

	void alloc_mem()
	{
		alpha = new double[nsnp * nsig];
		b = new double[npop * nsig * nsnp];
		b_sq = new double[npop * nsig * nsnp];
		loglik = new double[nsig];
		pst_loglik = new double[nsig];
		bDb = new double[nsig * nsnp];

		logBF = new double[nsig * npop * nsnp];

		cs_purity = new double[nsig];

		_memsig = nsig;
	}

	void clean()
	{
		delete []alpha;
		delete []b;
		delete []b_sq;
		delete []loglik;
		delete []pst_loglik;
		delete []bDb;

		delete []logBF;

		delete []cs_purity;
	}

	int susie_sst_xethn();

	int npop, nsnp, nsig, _memsig;

	double elbo_old;

	//par
	int *ngwas;

	//dataset
	char *ind;
	double *beta, *tau_sq;

	//SuSiEx
	double *alpha, *b, *b_sq, *loglik, *pst_loglik, *sigma_sq, *bhat, *beta_ll, *erss, *bDb;

	//SER
	double *v_sq, *z, *logBF, *mu, *mu_sq, *phi_sq, *w, *maxPopBF, *b_sq_sum;

	int ncs;
	//CS
	double *cs_purity, *pip;
	char *cs_bin;

	bool first;

	std::vector<cs> csset;

	std::vector<csVar> csDetails;

private:
	int init(int max_sig);
	int iter();
	int ser(int cur_sig);
	int cal_pip();
	int write_cs(int convergent);
};

void unique(std::vector<cs>& csset);

#endif /* MODEL_HPP_ */
