#include "utils.h"
#include <cassert>
#include <fstream>
#include <string>
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
        double val = A[0][i] / (A[1][i] - A[0][i - 1] * A[2][i - 1]);
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

vector<double> trigMult(const vector<vector<double>>& tridiag, const vector<double>& f) {
    int n = f.size();
    vector<double> res(n);
    res[0] = f[0] * tridiag[1][0] + f[1] * tridiag[0][0];
    for (int i = 1; i < n - 1; ++i) {
        res[i] = f[i - 1] * tridiag[2][i-1] + f[i] * tridiag[1][i] + f[i+1] * tridiag[0][i];
    }
    res[n - 1] = f[n - 2] * tridiag[2][n - 2] + f[n - 1] * tridiag[1][n - 1];
    return res;
}

vector<double> addVec(const vector<double>& a, const vector<double>& b) {
    assert(a.size() == b.size());
    vector<double> res(a.size());
    for (int i = 0; i < a.size(); ++i) {
        res[i] = a[i] + b[i];
    }
    return res;
}

vector<double> subVec(const vector<double>& a, const vector<double>& b) {
    assert(a.size() == b.size());
    vector<double> res(a.size());
    for (int i = 0; i < a.size(); ++i) {
        res[i] = a[i] - b[i];
    }
    return res;
}

void write_csv(string fname,vector<double> meshX,vector<double> meshT, function<vector<double>(int)> values){
	ofstream file;
	file.open(fname);
	// add header
	string header = "T,";
	for (double i : meshX) {
		header += to_string(i) + ",";
	}
	file << header.substr(0, header.size() - 1) << endl;
	for (int pt = 0; pt < meshT.size(); ++pt) {
		vector<double> vect = values(pt);
		double T = meshT[pt];
		string line = to_string(T) + ",";
		for (double val: vect) {
			line += to_string(val) + ",";
		}
		file << line.substr(0, line.size() - 1) << endl;
	}

	file.close();
}