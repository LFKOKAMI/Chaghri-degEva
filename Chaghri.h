#ifndef _CHAGHRI_H_
#define _CHAGHRI_H_
#include "gurobi_c++.h"
#include <vector>
using namespace std;

typedef unsigned long long u64;
const u64 MAX = 0x8000000000000000;

class Chaghri {
public:
	void addition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& d,int fs);
	void modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int fs);
	void compare(GRBModel& model, vector<GRBVar>& a, u64 m, int fs);
	void computeCoefficientFrequency(int r, vector<vector<u64> >& n);

	int multivariatePoly(GRBEnv& env, int r, bool isRefined, bool isReduceUsed);
	int univariatePoly(GRBEnv& env, int r,bool isRefined, bool isReduceUsed);
	void initializeVec(vector<vector<GRBVar> >& vec,int size);
	void initializeVar(GRBModel& model, vector<vector<GRBVar> >& vec);

	void reduce(vector<u64>& vec, int n, int t);
	void reduce(u64 vec[], int n, int t);
	void computeCoefficientFrequencyWithReduction(int r, vector<vector<u64> >& n, int t, int k0, int k1, int size);
	int constructMiMCModel(int r, GRBEnv& env,bool isRefined);

};

#endif
