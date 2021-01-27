#include "greeks.h"
#include <iostream>

using namespace std;

vector<double> delta(const PDE& pde, int t_idx, bool convert) {
	int n = pde.values[0].size();
	const vector<double>& price = pde.values[t_idx];
	double t = pde.meshT[0];
	const vector<double>& meshX = pde.meshX;
	const double dx = pde.dx;
	const PDEBounds& bound = pde.bound;

	vector<double> res(n);

	res[0] = (price[1] - pde.bound.xmin(meshX[0], t)) / (2 * dx);
	if (convert) {
		res[0] /= exp(meshX[1]);
	}
	for (int i = 1; i < n - 1; ++i) {
		res[i] = (price[i + 1] - price[i - 1]) / (2 * dx);
		if (convert) {
			res[i] /= exp(meshX[i + 1]);
		}
	}
	res[n - 1] = (bound.xmax(meshX[n + 1], t) - price[n - 2]) / (2 * dx);
	if (convert) {
		res[n - 1] /= exp(meshX[n]);
	}
	return res;
}

void delta_csv(const PDE& pde, string fname, bool convert) {

	const vector<double>& meshX = pde.meshX;
	const vector<double>& meshT = pde.meshT;

	cout << "delta: ";
	vector<double> delta_0 = delta(pde, 0, convert);
	int nX = (meshX.size() - 1) / 2;
	cout << delta_0[nX-1] << endl;

	vector<double> mesh_X(meshX.begin() + 1, meshX.end() - 1);
	mesh_X = (convert) ? apply1D([](double x) {return exp(x); }, mesh_X) : mesh_X;
	write_csv(fname, mesh_X, meshT, [pde, convert](int idx) {
		return delta(pde, idx, convert);
		});
}