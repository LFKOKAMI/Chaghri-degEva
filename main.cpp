#include "Chaghri.h"
#include <iostream>
#include <ctime>
using namespace std;

int main() {
	Chaghri ch;
	GRBEnv env;
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 4);

	cout << "=======================================" << endl;
	cout << "unrefined model (univariate, using reduction algorithm):" << endl;
	for (int i = 1; i < 28; i++) {
		cout << "degree at step " << i << ":" << ch.univariatePoly(env, i, 0, 0) << endl;
	}

	cout << "unrefined model (univariate, not using reduction algorithm):" << endl;
	for (int i = 1; i < 28; i++) {
		cout << "degree at step " << i << ":" << ch.univariatePoly(env, i, 0, 0) << endl;
	}

	cout << "=======================================" << endl;
	cout << "refined model (univariate, using reduction algorithm):" << endl;
	for (int i = 2; i < 28; i++) {
		cout << "degree at step " << i << ":" << ch.univariatePoly(env, i - 1, 1, 1) << endl;
	}

	cout << "=======================================" << endl;
	cout << "refined model (univariate, not using reduction algorithm, slightly slow):" << endl;
	for (int i = 2; i < 28; i++) {
		cout << "degree at step " << i << ":" << ch.univariatePoly(env, i - 1, 1, 0) << endl;
	}

	cout << "unrefined model (multivariate,using the reduction algorithm):" << endl;
	for (int i = 1; i < 28; i++) {
		cout << "degree at step " << i << ":" << ch.multivariatePoly(env, i, 0,1) << endl;
	}

	cout << "=======================================" << endl;
	cout << "refined model (multivariate, using reduction algorithm):" << endl;
	for (int i = 27; i < 29; i++) {
		cout << "degree at step " << i << ":" << ch.multivariatePoly(env, i - 1, 1,1) << endl;
	}

	cout << "=======================================" << endl;
	cout << "unrefined model for MiMC (using reduction algorithm):" << endl;
	for (int i = 1; i <= 80; i++) {
		cout << "degree at step " << i << ":" << ch.constructMiMCModel(i, env, 0) << endl;
	}

	cout << "=======================================" << endl;
	cout << "refined model for MiMC (using reduction algorithm):" << endl;
	for (int i = 2; i <= 80; i++) {
		cout << "degree at step " << i << ":" << ch.constructMiMCModel(i - 1, env, 1) << endl;
	}

	return 1;
}