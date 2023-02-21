#include "Chaghri.h"
#include "gurobi_c++.h"
#include <vector>
using namespace std;

//find equivalent sets (reduction algorithm)
void Chaghri::reduce(vector<u64>& vec, int n, int t) {
	int times = 2;
	while (times > 0) {
		for (int i = 0; i < n; i++) {
			if (vec[i]>=t && (vec[i] - t) % 2 == 0) {
				vec[(i + 1) % n] = vec[(i + 1) % n] + (vec[i] - t) / 2;
				vec[i] = t;
			}
			if (vec[i]>=t && (vec[i] - t) % 2 == 1) {
				vec[(i + 1) % n] = vec[(i + 1) % n] + (vec[i] - t - 1) / 2;
				vec[i] = t + 1;
			}
		}
		times--;
	}
}


void Chaghri::addition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& d,int fs) {
	for (int i = 0; i < fs; i++) {
		model.addConstr(a[(i - sa + fs) % fs] + b[(i - sb + fs) % fs] + c[i] == d[i] + 2 * c[i + 1]);
	}
}

void Chaghri::modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int fs) {
	model.addConstr(c[0] == 0);
	for (int i = 0; i < fs; i++) {
		model.addConstr(a[(i - sa + fs) % fs] + b[(i - sb + fs) % fs] + c[i] == q[i] + 2 * c[i + 1]);
	}
	model.addConstr(c1[0] == c[fs]);
	for (int i = 0; i < fs; i++) {
		model.addConstr(2 * c1[i + 1] + d[i] == q[i] + c1[i]);
	}
}

void Chaghri::compare(GRBModel& model, vector<GRBVar>& a, u64 m, int fs) {
	//get the binary representation of m
	vector<bool> mvec(fs);
	for (int i = 0; i < fs; i++) {
		mvec[i] = (m >> i) & 0x1;
	}
	//find the set i_1, i_2, ..., i_h
	vector<int> p;
	p.clear();
	for (int i = 62; i >= 0; i--) {
		if (mvec[i]) {
			p.push_back(i);
		}
	}
	//deduce the 63-p.size() linear inequalities
	//deduce the equalities
	for (int i = p[0]+1; i <=62; i++) {
		model.addConstr(a[i] == 0);
	}
	//deduce the inequalities
	for (int len = 1; len <= p.size(); len++) {
		//until p[len-1]
		GRBLinExpr sum = 0;
		for (int t = 0; t < len; t++) {
			sum += (1 - a[p[t]]);
		}
		if (len < p.size()) {
			for (int j = p[len] + 1; j <= p[len - 1]-1; j++) {
				model.addConstr(sum - a[j] >= 0);
			}
		}
		else {
			for (int j = 0; j <= p[len - 1] - 1; j++) {
				model.addConstr(sum - a[j] >= 0);
			}
		}
	}
}

void Chaghri::computeCoefficientFrequency(int r, vector<vector<u64> >& n) {
	int k = 32;
	bool flag = false;
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < 63; j++) {
			if (n[i][(j - 3 + 63) % 63] >= MAX && n[i][(j - 3 - k + 63) % 63] >= MAX) {
				cout << "maximal degree is reached" << endl;
				flag = true;
			}
			else {
				n[i + 1][j] = n[i][(j - 3 + 63) % 63] + n[i][(j - 3 - k + 63) % 63];
			}
		}
		if (flag) {
			break;
		}
	}
}

//using the equivalence theorem, we can use the reduction after each step update for the vector
void Chaghri::computeCoefficientFrequencyWithReduction(int r, vector<vector<u64> >& n,int t, int k0, int k1,int size) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < size; j++) {
			n[i + 1][j] = n[i][(j - k1 + size) % size] + n[i][(j - k1 - k0 + 2*size) % size];	
		}
		/*for (int j = 0; j < size; j++) {
			cout << n[i + 1][j] << " ";
		}
		cout << endl;*/
		reduce(n[i + 1], size, t);
		/*for (int j = 0; j < size; j++) {
			cout << n[i + 1][j] << " ";
		}
		cout << endl;
		system("pause");*/

	}
}


void Chaghri::initializeVec(vector<vector<GRBVar> >& vec, int size) {
	int len = vec.size();
	for (int i = 0; i < len; i++)
		vec[i].resize(size);
}

void Chaghri::initializeVar(GRBModel& model, vector<vector<GRBVar> >& vec) {
	for (int i = 0; i < vec.size(); i++) {
		for (int j = 0; j < vec[i].size(); j++) {
			vec[i][j]=model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}

int Chaghri::univariatePoly(GRBEnv& env, int r,bool isRefined,bool isReduceUsed) {
	vector<vector<u64> > n;
	n.resize(r + 1);
	for (int i = 0; i < r + 1; i++) {
		n[i].resize(63);
		for (int j = 0; j < 63; j++) {
			n[i][j] = 0;
		}
	}
	n[0][0] = 1;


	if (isReduceUsed) {//can improve the performance
		computeCoefficientFrequencyWithReduction(r, n, 1, 32, 3, 63);//not big integer occurs, 1 variable
	}
	else {
		computeCoefficientFrequency(r, n);//be careful of the overflow of big integer addition
	}

	//get the nonempty group
	vector<int> pos;
	pos.clear();
	for (int i = 0; i < 63; i++) {
		if (n[r][i] != 0) {
			pos.push_back(i);
		}
	}
	int size = pos.size();

	GRBModel model = GRBModel(env);

	int fs = 63;//field size
	vector<GRBVar> zero(fs);
	for (int i = 0; i < fs; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
	}

	vector<vector<GRBVar> > alpha;
	alpha.resize(size);
	initializeVec(alpha, fs);
	initializeVar(model, alpha);
	for (int i = 0; i < size; i++) {//alpha_i <= n_i
		compare(model, alpha[i], n[r][pos[i]], fs);
	}

	vector<vector<GRBVar> > c, d;
	c.resize(size * 4);
	d.resize(size * 4 + 1);
	initializeVec(c, fs + 1);
	initializeVar(model, c);
	initializeVec(d, fs);
	initializeVar(model, d);


	//sum(2^pos[i]*alpha_i)
	for (int i = 0; i < 63; i++) {
		model.addConstr(d[0][i] == 0);
	}

	for (int i = 0; i < size; i++) {
		if (isRefined == false) {
			modAddition(model, d[i * 2], 0, alpha[i], pos[i], c[i * 2], c[i * 2 + 1], d[i * 2 + 1], d[i * 2 + 2], fs);
		}
		else {
			modAddition(model, d[i * 4], 0, alpha[i], pos[i], c[i * 4], c[i * 4 + 1], d[i * 4 + 1], d[i * 4 + 2], fs);
			modAddition(model, d[i * 4+2], 0, alpha[i], (pos[i] + 32) % 63, c[i * 4+2], c[i * 4 + 3], d[i * 4 + 3], d[i * 4 + 4], fs);
		}
	}

	GRBLinExpr sum = 0;
	for (int i = 0; i < fs; i++) {
		if (isRefined == false) {
			sum += d[2 * size][i];
		}
		else {
			sum += d[4 * size][i];
		}
	}

	model.setObjective(sum, GRB_MAXIMIZE);
	model.optimize();

	int res = model.getObjective().getValue();

	/*u64 one = 1;
	for (int i = 0; i < size; i++) {
		u64 val = 0;
		for (int j = 0; j < 63; j++) {
			if (alpha[i][j].get(GRB_DoubleAttr_X) == 1) {
				val = val | (one << j);
			}
		}
		cout << hex << val << " <= " << n[r][pos[i]] << endl;
	}

	cout << "hamming weight:" << dec<<res << endl;*/
	return res;
}


int Chaghri::multivariatePoly(GRBEnv& env, int r, bool isRefined,bool isReduceUsed) {
	vector<vector<u64> > n;
	n.resize(r + 1);
	for (int i = 0; i < r + 1; i++) {
		n[i].resize(63);
		for (int j = 0; j < 63; j++) {
			n[i][j] = 0;
		}
	}
	n[0][0] = 1;
	
	if (isReduceUsed) {//can improve the performance
		computeCoefficientFrequencyWithReduction(r, n, 2, 32, 3, 63);//not big integer occurs, 2 variable
	}
	else {
		computeCoefficientFrequency(r, n);//be careful of the overflow of big integer addition
	}



	//get the nonempty group
	vector<int> pos;
	pos.clear();
	for (int i = 0; i < 63; i++) {
		if (n[r][i] != 0) {
			pos.push_back(i);
		}
	}
	int size = pos.size();

	GRBModel model = GRBModel(env);

	int fs = 63;//field size
	vector<GRBVar> zero(fs), zero1(fs);
	for (int i = 0; i < fs; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		zero1[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
		model.addConstr(zero1[i] == 0);
	}

	vector<vector<GRBVar> > alpha, beta, zeta, carry;
	alpha.resize(size);
	beta.resize(size);
	zeta.resize(size);
	carry.resize(size);

	initializeVec(alpha, fs);
	initializeVar(model, alpha);
	initializeVec(beta, fs);
	initializeVar(model, beta);
	initializeVec(zeta, fs);
	initializeVar(model, zeta);
	initializeVec(carry, fs + 1);
	initializeVar(model, carry);
	for (int i = 0; i < size; i++) {
		model.addConstr(carry[i][0] == 0);
		addition(model, alpha[i], 0, beta[i], 0, carry[i], zeta[i], fs);
		model.addConstr(carry[i][fs] == 0);//directly exit because the maximal degree is reached
		compare(model, zeta[i], n[r][pos[i]], fs);
	}

	vector<vector<GRBVar> > c, d;
	c.resize(size * 4);
	d.resize(size * 4 + 1);
	initializeVec(c, fs + 1);
	initializeVar(model, c);
	initializeVec(d, fs);
	initializeVar(model, d);

	vector<vector<GRBVar> > cbeta, dbeta;
	cbeta.resize(size * 4);
	dbeta.resize(size * 4 + 1);
	initializeVec(cbeta, fs + 1);
	initializeVar(model, cbeta);
	initializeVec(dbeta, fs);
	initializeVar(model, dbeta);


	//sum(2^pos[i]*alpha_i)
	for (int i = 0; i < 63; i++) {
		model.addConstr(d[0][i] == 0);
		model.addConstr(dbeta[0][i] == 0);
	}

	for (int i = 0; i < size; i++) {
		if (isRefined == false) {
			//alpha
			modAddition(model, d[i * 2], 0, alpha[i], pos[i], c[i * 2], c[i * 2 + 1], d[i * 2 + 1], d[i * 2 + 2], fs);
			//beta
			modAddition(model, dbeta[i * 2], 0, beta[i], pos[i], cbeta[i * 2], cbeta[i * 2 + 1], dbeta[i * 2 + 1], dbeta[i * 2 + 2], fs);
		}
		else {
			//alpha
			modAddition(model, d[i * 2], 0, alpha[i], pos[i], c[i * 2], c[i * 2 + 1], d[i * 2 + 1], d[i * 2 + 2], fs);
			//beta
			modAddition(model, dbeta[i * 4], 0, beta[i], pos[i], cbeta[i * 4], cbeta[i * 4 + 1], dbeta[i * 4 + 1], dbeta[i * 4 + 2], fs);
			modAddition(model, dbeta[i * 4 + 2], 0, beta[i], (pos[i] + 32) % 63, cbeta[i * 4 + 2], cbeta[i * 4 + 3], dbeta[i * 4 + 3], dbeta[i * 4 + 4], fs);
		}
	}

	GRBLinExpr sum0 = 0, sum1 = 0;
	for (int i = 0; i < fs; i++) {
		if (isRefined == false) {
			sum0 = sum0 + d[2 * size][i];
			sum1 = sum1 + dbeta[2 * size][i];
		}

		else {
			sum0 = sum0 + d[2 * size][i];
			sum1 = sum1 + dbeta[4 * size][i];
		}
	}

	if(isRefined)
		model.addConstr(sum0 == 63);
	model.setObjective(sum0+sum1, GRB_MAXIMIZE);
	model.optimize();

	int res = model.getObjective().getValue();

	return res;
}


int Chaghri::constructMiMCModel(int r, GRBEnv& env,bool isRefined) {
	int n = 129;
	int t = 1;
	int end = 2;
	vector<vector<u64> > vec(r + 1);
	for (int i = 0; i < r + 1; i++) {
		vec[i].resize(n);
		for (int j = 0; j < n; j++) {
			vec[i][j] = 0;
		}
	}
	vec[0][0] = 1;

	for (int i = 1; i <= r; i++) {
		for (int j = 0; j < n; j++) {
			vec[i][j] = vec[i - 1][j] + vec[i - 1][(j - 1 + n) % n];
		}
		reduce(vec[i], n, t);
	}

	int size = 0;
	vector<int> pos;
	pos.clear();
	for (int i = 0; i < n; i++) {
		if (vec[r][i] > 0) {
			size++;
			pos.push_back(i);
		}
	}

	GRBModel model = GRBModel(env);

	vector<GRBVar> zero(n);
	for (int i = 0; i < n; i++) {
		zero[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		model.addConstr(zero[i] == 0);
	}

	vector<vector<GRBVar> > a;
	a.resize(size);
	initializeVec(a, n);
	initializeVar(model, a);

	for (int i = 0; i < size; i++) {//a_i <= vec[i]
		for (int j = n - 1; j >= end; j--) {
			model.addConstr(a[i][j] == 0);
		}
		if (vec[r][pos[i]] == 2) {
			model.addConstr(a[i][1] + a[i][0] <= 1);
		}
		else if (vec[r][pos[i]] == 1) {
			model.addConstr(a[i][1] == 0);
		}
	}

	vector<vector<GRBVar> > c, d;
	c.resize(size * 4);
	d.resize(size * 4 + 1);
	initializeVec(c, n + 1);
	initializeVar(model, c);
	initializeVec(d, n);
	initializeVar(model, d);

	//sum(2^pos[i]*alpha_i)
	for (int i = 0; i < n; i++) {
		model.addConstr(d[0][i] == 0);
	}

	for (int i = 0; i < size; i++) {
		//model.addConstr(c[i * 4][0] == 0);
		//addition(model, d[i * 4], 0, a[i], pos[i], c[i * 4], d[i * 4 + 1], n);
		//model.addConstr(c[i * 4][n] == c[i * 4 + 1][0]);
		//addition(model, d[i * 4 + 1], 0, zero, 0, c[i * 4 + 1], d[i * 4 + 2], n);
		if (isRefined) {
			modAddition(model, d[i * 4], 0, a[i], pos[i], c[i * 4], c[i * 4 + 1], d[i * 4 + 1], d[i * 4 + 2], n);
			modAddition(model, d[i * 4 + 2], 0, a[i], (pos[i] + 1) % n, c[i * 4 + 2], c[i * 4 + 3], d[i * 4 + 3], d[i * 4 + 4], n);
		}
		else {
			modAddition(model, d[i * 2], 0, a[i], pos[i], c[i * 2], c[i * 2 + 1], d[i * 2 + 1], d[i * 2 + 2], n);
		}
		//model.addConstr(c[i * 4 + 2][0] == 0);
		//addition(model, d[i * 4 + 2], 0, a[i], (pos[i] + 1) % n, c[i * 4 + 2], d[i * 4 + 3], n);
		//model.addConstr(c[i * 4 + 2][n] == c[i * 4 + 3][0]);
		//addition(model, d[i * 4 + 3], 0, zero, 0, c[i * 4 + 3], d[i * 4 + 4], n);

	}

	GRBLinExpr sum = 0;
	for (int i = 0; i < n; i++) {
		if(isRefined)
			sum += d[4 * size][i];
		else {
			sum += d[2 * size][i];
		}
	}

	model.setObjective(sum, GRB_MAXIMIZE);
	model.optimize();

	int res = model.getObjective().getValue();
	//cout << r << ": " << res << endl;

	for (int i = 0; i < r + 1; i++)
		vec[i].clear();
	vec.clear();

	return res;
}