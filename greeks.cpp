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

vector<double> gamma(const PDE& pde,int t_idx, bool convert){
    int n = pde.values[0].size();
    const vector<double>& price = pde.values[t_idx];
    double t = pde.meshT[t_idx];
    vector<double> firstOrder;
    if (convert) {
        firstOrder = delta(pde,t_idx, false);
    }
    const vector<double>& meshX = pde.meshX;
    const double dx = pde.dx;
    const PDEBounds& bound = pde.bound;
    vector<double> res(n);
    res[0] = (price[1] - 2 * price[0] + bound.xmin(meshX[0], t)) / (dx * dx);
    if (convert) {
        double s = exp(meshX[1]);
        res[0] = (res[0] - firstOrder[0]) / (s * s);
    }
    for (int i = 1; i < n - 1; ++i) {
        res[i] = (price[i + 1] - 2 * price[i] + price[i - 1]) / (dx * dx);
        if (convert) {
            double s = exp(meshX[i + 1]);
            res[i] = (res[i] - firstOrder[i]) / (s * s);
        }
    }
    res[n - 1] = (bound.xmax(meshX[n + 1], t) - 2 * price[n - 1] + price[n - 2]) / (dx * dx);
    if (convert) {
        double s = exp(meshX[n]);
        res[n - 1] = (res[n - 1] - firstOrder[n - 1]) / (s * s);
    }
    return res;
}

void gamma_csv(const PDE& pde,string fname, bool convert) {
    const vector<double>& meshX = pde.meshX;
    const vector<double>& meshT = pde.meshT;
    
    cout << "gamma: ";
    vector<double> gamma_0 = gamma(pde, 0, convert);
    int nX = (meshX.size() - 1) / 2;
    cout << gamma_0[nX] << endl;

    vector<double> mesh_X(meshX.begin() + 1, meshX.end() - 1);
    mesh_X = (convert) ? apply1D([](double x) {return exp(x); }, mesh_X) : mesh_X;
    write_csv(fname, mesh_X, meshT, [pde, convert](int idx) {
        return gamma(pde,idx, convert);
    });
}
