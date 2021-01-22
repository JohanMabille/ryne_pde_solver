#include "utils.h"
#include <iostream>

using namespace std;

vector<double> apply1D(const function<double(double)>& f, const vector<double>& x) {
	vector<double> res(x.size());
	for (int i = 0; i < x.size(); ++i) {
		res[i] = f(x[i]);
	}
	return res;
}

vector<double> apply2D(const function<double(double, double)>& f, const vector<double>& x, double t) {
	vector<double> res(x.size());
	for (int i = 0; i < x.size(); ++i) {
		res[i] = f(x[i], t);
	}
	return res;
}

vector<vector<double>> trigMat(const vector<vector<double>>& A,const vector<double>& b) {
	int n = A[1].size();

	vector<double> newA(n - 1);
	newA[0] = A[0][0] / A[1][0];
	for (int i = 1; i < n - 1; ++i) {
		double val = A[0][i] / (A[1][i] - A[0][i] * A[2][i]);
		newA[i] = val;
	}

	vector<double> newB(n);
	newB[0] = b[0] / A[1][0];

	for (int i = 1; i < n; ++i) {
		double val = b[i] - newB[i - 1] * A[2][i - 1];
		val /= A[1][i] - newA[i - 1] * A[2][i - 1];
		newB[i] = val;
	}
	return { newA,newB };
}


