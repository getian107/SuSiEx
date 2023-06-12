/*
 * data.hpp
 *
 *  Created on: Dec 27, 2021
 *      Author: kyuan
 */

#ifndef DATA_HPP_
#define DATA_HPP_

# include <iostream>
# include <vector>
# include <cmath>
# include <cstdlib>
# include <cstring>

# define LDTYPE float

extern int npop;

class softpar
{
public:
	softpar() : out_dir(""), out_name(""), chr(""), start(0), end(0), plink(""), \
	keep_ambig(false), mult_step(false), precmp(false), maf(0.005), level(0.95), \
	min_purity(0.5), pth(1e-5), tol(1e-4), n_sig(5), max_iter(100), nthreads(1), \
	key_by(2) {}

	void split_str(const std::string& info, const std::string& label)
	{
		int st(0);
		std::vector<std::string> *str;
		if(label == "sst")
			str = &sst_file;
		else if(label == "ref")
			str = &ref_file;
		else if(label == "ld")
			str = &ld_file;
		else
		{
			std::cerr << "Cannot identify " << label << std::endl;
			exit(1);
		}
		for(int i = 0 ; i < info.size(); ++i)
			if(info[i] == ',')
			{
				str->push_back(info.substr(st, i - st));
				st = i + 1;
			}
		str->push_back(info.substr(st));
	}

	void split_int(const std::string& info, const std::string& label)
	{
		int st(0);
		std::vector<int> *num;
		if(label == "n")
			num = &n_gwas;
		else if(label == "chr")
			num = &chr_col;
		else if(label == "snp")
			num = &snp_col;
		else if(label == "bp")
			num = &bp_col;
		else if(label == "a1")
			num = &a1_col;
		else if(label == "a2")
			num = &a2_col;
		else if(label == "eff")
			num = &eff_col;
		else if(label == "se")
			num = &se_col;
		else if(label == "pval")
			num = &pval_col;
		else if(label == "stend")
		{
			for(int i = 0 ; i < info.size(); ++i)
				if(info[i] == ',')
				{
					start = atol(info.substr(0, i).c_str());
					end = atol(info.substr(i + 1).c_str());
				}
			return;
		}
		else
		{
			std::cerr << "Cannot identify " << label << std::endl;
			exit(1);
		}
		for(int i = 0 ; i < info.size(); ++i)
			if(info[i] == ',')
			{
				num->push_back(atol(info.substr(st, i - st).c_str()));
				st = i + 1;
			}
		num->push_back(atol(info.substr(st).c_str()));
	}

	void show() const
	{
		std::cout << "Software parameters:" << std::endl;
		//std::cout << "Number of populations: " << npop << std::endl;
		std::cout << "--sst_file = " << sst_file[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << sst_file[i];
		std::cout << std::endl;
		if(!precmp)
		{
			std::cout << "--ref_file = " << ref_file[0];
			for(int i = 1 ; i < npop; ++i)
				std::cout << ',' << ref_file[i];
			std::cout << std::endl;
		}
		std::cout << "--ld_file = " << ld_file[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << ld_file[i];
		std::cout << std::endl;
		std::cout << "--n_gwas = " << n_gwas[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << n_gwas[i];
		std::cout << std::endl;
		std::cout << "--out_dir = " << out_dir << std::endl;
		std::cout << "--out_name = " << out_name << std::endl;
		std::cout << "--chr = " << chr << std::endl;
		std::cout << "--bp = " << start << ',' << end << std::endl;
		std::cout << "--chr_col = " << chr_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << chr_col[i];
		std::cout << std::endl;
		std::cout << "--snp_col = " << snp_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << snp_col[i];
		std::cout << std::endl;
		std::cout << "--bp_col = " << bp_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << bp_col[i];
		std::cout << std::endl;
		std::cout << "--a1_col = " << a1_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << a1_col[i];
		std::cout << std::endl;
		std::cout << "--a2_col = " << a2_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << a2_col[i];
		std::cout << std::endl;
		std::cout << "--eff_col = " << eff_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << eff_col[i];
		std::cout << std::endl;
		std::cout << "--se_col = " << se_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << se_col[i];
		std::cout << std::endl;
		std::cout << "--pval_col = " << pval_col[0];
		for(int i = 1 ; i < npop; ++i)
			std::cout << ',' << pval_col[i];
		std::cout << std::endl;
		std::cout << "--plink = " << plink << std::endl;
		std::cout << "--keep-ambig = " << std::boolalpha << keep_ambig << std::noboolalpha << std::endl;
		std::cout << "--mult-step = " << std::boolalpha << mult_step << std::noboolalpha << std::endl;
		std::cout << "--precmp = " << std::boolalpha << precmp << std::noboolalpha << std::endl;
		std::cout << "--maf = " << maf << std::endl;
		std::cout << "--level = " << level << std::endl;
		std::cout << "--min_purity = " << min_purity << std::endl;
		std::cout << "--pval_thresh = " << pth << std::endl;
		std::cout << "--tol = " << tol << std::endl;
		std::cout << "--n_sig = " << n_sig << std::endl;
		std::cout << "--max_iter = " << max_iter << std::endl;
		std::cout << "--threads = " << nthreads << std::endl;

	}

	std::vector<std::string> sst_file, ref_file, ld_file;
	std::vector<int> n_gwas;
	std::string out_dir, out_name;
	std::string chr;
	long start, end;
	std::vector<int> chr_col, snp_col, bp_col, a1_col, a2_col, eff_col, se_col, pval_col;
	std::string plink;
	bool keep_ambig, mult_step, precmp;
	double maf, level, min_purity, pth, tol;
	int n_sig, max_iter, nthreads;

	/*
	 * 0 rsID, 1 position ref / alt, 2 position + rsID
	 * 2 as default
	 */
	int key_by;
};

class varCoord
{
public:
	std::string id, a1, a2;
	long pos;
};

inline bool cmp_rs(const varCoord& a, const varCoord& b)
{
	return a.id < b.id;
}

inline bool eql_rs(const varCoord& a, const varCoord& b)
{
	return a.id == b.id;
}

inline bool cmp_pos_ref(const varCoord& a, const varCoord& b)
{
	if(a.pos != b.pos)
		return a.pos < b.pos;
	if(a.a1 != b.a1)
		return a.a1 < b.a1;
	return a.a2 < b.a2;
}

inline bool eql_pos_ref(const varCoord& a, const varCoord& b)
{
	return a.pos == b.pos && a.a1 == b.a1 && a.a2 == b.a2;
}

inline bool cmp_pos_rs(const varCoord& a, const varCoord& b)
{
	if(a.pos != b.pos)
		return a.pos < b.pos;
	return a.id < b.id;
}

inline bool eql_pos_rs(const varCoord& a, const varCoord& b)
{
	return a.pos == b.pos && a.id == b.id;
}

static bool (*var_cmp_func)(const varCoord& a, const varCoord& b) = cmp_pos_rs;

static bool (*var_eql_func)(const varCoord& a, const varCoord& b) = eql_pos_rs;

class sumstats
{
public:
	varCoord coord;
	int idx;
	double beta;
	long double pval;
};

inline bool operator<(const sumstats& a, const sumstats& b)
{
	return var_cmp_func(a.coord, b.coord);
}

inline bool operator==(const sumstats& a, const sumstats& b)
{
	return var_eql_func(a.coord, b.coord);
}

class ldref
{
public:
	varCoord coord;
	int idx;
	double frq;
};

inline bool operator<(const ldref& a, const ldref& b)
{
	return var_cmp_func(a.coord, b.coord);
}

inline bool operator==(const ldref& a, const ldref& b)
{
	return var_eql_func(a.coord, b.coord);
}

inline int stat_amb(const varCoord& ref, const varCoord& q)
{
	if(ref.pos != q.pos)
		return 0;
	if(ref.a1 == q.a1 && ref.a2 == q.a2)
		return 3;
	else if(ref.a1 == q.a2 && ref.a2 == q.a1)
		return 1;
	return 0;
}

inline int stat_con(const varCoord& ref, const varCoord& q)
{
	if(ref.pos != q.pos)
		return 0;
	if(ref.a1 == q.a1 && ref.a2 == q.a2)
		return 3;
	return 0;
}

static int (*stat_func)(const varCoord& ref, const varCoord& q) = stat_amb;

class snp
{
public:
	snp()
	{
		idxs = new int[npop * 3];
		memset(idxs, 0, sizeof(int) * 3 * npop);
		stats = new double[npop * 4];
		memset(stats, 0, sizeof(double) * 4 * npop);
	}
	~snp()
	{
		delete []idxs;
		delete []stats;
	}

	snp(const snp& dat)
	{
		if(this != & dat)
		{
			coord = dat.coord;
			idxs = new int[npop * 3];
			memcpy(idxs, dat.idxs, sizeof(int) * npop * 3);
			stats = new double[npop * 4];
			memcpy(stats, dat.stats, sizeof(double) * npop * 4);
		}
	}

	snp & operator=(const snp& dat)
	{
		if(this != & dat)
		{
			coord = dat.coord;
			memcpy(idxs, dat.idxs, sizeof(int) * npop * 3);
			memcpy(stats, dat.stats, sizeof(double) * npop * 4);
		}
		return *this;
	}

	void ref_set(int popIdx, const ldref& s)
	{
		coord = s.coord;
		idxs[popIdx * 3] = 3;
		idxs[popIdx * 3 + 1] = s.idx;
		stats[popIdx * 4] = s.frq;
	}
	void set(int popIdx, const ldref& s)
	{
		idxs[popIdx * 3] = stat_func(coord, s.coord);
		idxs[popIdx * 3 + 1] = s.idx;
		stats[popIdx * 4] = s.frq;
	}
	void set_beta(int popIdx, const sumstats& s)
	{
		int st = stat_func(coord, s.coord);

		if(st)
		{
			if(st == idxs[popIdx * 3])
				stats[popIdx * 4 + 1] = s.beta;
			else
				stats[popIdx * 4 + 1] = - s.beta;
			idxs[popIdx * 3] += 4;
			idxs[popIdx * 3 + 2] = s.idx;
			stats[popIdx * 4 + 2] = s.pval;
			stats[popIdx * 4 + 3] = - log10(s.pval);
		}
	}
	bool var() const
	{
		for(int i = 0 ; i < npop; ++i)
			if((idxs[i * 3] & 5) == 5)
				return true;
		return false;
	}
	bool var(int idxPop) const
	{
		return ((idxs[idxPop * 3] & 5) == 5);
	}
	/*
	int ldflag(int idxPop) const
	{
		return 1;
		if((idxs[idxPop * 3] & 3) == 3)
			return 1;
		else if((idxs[idxPop * 3] & 3) == 1)
			return -1;
		else
			return 0;
	}
	*/

	varCoord coord;

	//[[1 in LD ref file, 2 non-rev in LD ref file, 4 in gwas file]
	//[index of LD ref file] [index of gwas file]]
	int *idxs;

	//[[frq] [beta] [p] [logp]]
	double *stats;

	static softpar par;
};

inline std::ostream & operator << (std::ostream & os, const snp& var)
{
	os << var.coord.id << '\t' << var.coord.pos;
	if((var.idxs[0] & 7) == 7)
		os << '\t' << var.coord.a1;
	else if((var.idxs[0] & 7) == 5)
		os << '\t' << var.coord.a2;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 7) == 7)
			os << ',' << var.coord.a1;
		else if((var.idxs[i * 3] & 7) == 5)
			os << ',' << var.coord.a2;
		else
			os << ",NA";

	if((var.idxs[0] & 7) == 7)
		os << '\t' << var.coord.a2;
	else if((var.idxs[0] & 7) == 5)
		os << '\t' << var.coord.a1;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 7) == 7)
			os << ',' << var.coord.a2;
		else if((var.idxs[i * 3] & 7) == 5)
			os << ',' << var.coord.a1;
		else
			os << ",NA";

	if((var.idxs[0] & 5) == 5)
		os << '\t' << var.stats[0];
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << var.stats[i * 4];
		else
			os << ",NA";

	double sqrtfrq(sqrt(2 * var.stats[0] * (1 - var.stats[0])));
	if((var.idxs[0] & 7) == 7)
		os << '\t' << var.stats[1] / sqrtfrq;
	else if((var.idxs[0] & 7) == 5)
		os << '\t' << - var.stats[1] / sqrtfrq;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
	{
		sqrtfrq = sqrt(2 * var.stats[i * 4] * (1 - var.stats[i * 4]));
		if((var.idxs[i * 3] & 7) == 7)
			os << ',' << var.stats[i * 4 + 1] / sqrtfrq;
		else if((var.idxs[i * 3] & 7) == 5)
			os << ',' << - var.stats[i * 4 + 1] / sqrtfrq;
		else
			os << ",NA";
	}

	sqrtfrq = sqrt(2 * var.stats[0] * (1 - var.stats[0]));
	if((var.idxs[0] & 5) == 5)
		os << '\t' << 1.0 / sqrt(snp::par.n_gwas[0]) / sqrtfrq;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
	{
		sqrtfrq = sqrt(2 * var.stats[i * 4] * (1 - var.stats[i * 4]));
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << 1.0 / sqrt(snp::par.n_gwas[i]) / sqrtfrq;
		else
			os << ",NA";
	}

	if((var.idxs[0] & 5) == 5)
		os << '\t' << var.stats[3];
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << var.stats[i * 4 + 3];
		else
			os << ",NA";

    return os;
}

class dataset
{
public:
	dataset() : beta(NULL), pval(NULL), mkIdx(NULL), ind(NULL), nsnp(0)
	{
		tau_sq = new double[npop];
		ld = new LDTYPE**[npop];
	}

	~dataset()
	{
		if(npop)
		{
			if(nsnp)
				for(int i = 0 ; i < npop; ++i)
					for(int j = 0 ; j < nsnp; ++j)
						delete []ld[i][j];
			delete []tau_sq;
			delete []beta;
			delete []pval;
			delete []mkIdx;
			delete []ind;
			delete []ld;
		}
	}

	void load(const softpar& par);

	//pop
	double *tau_sq;

	//pop X snp
	double *beta, *pval;
	//snp
	int *mkIdx;
	//pop X snp
	char *ind;

	//pop snp snp
	LDTYPE ***ld;
	std::vector<snp> mks;

	int nsnp;
};


#endif /* DATA_HPP_ */
