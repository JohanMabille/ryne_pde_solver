#include "blackscholes.h"


pdefunc generateCallPayoff(double K) {
    return [K](double x, double t)->double {
        double res = exp(x) - K;
        return (res > 0) ? res : 0;
    };
}

pdefunc generatePutPayoff(double K) {
	return [K](double x, double t)->double {
		double res = K - exp(x);
		return (res > 0) ? res : 0;
	};
}

pdefunc generateCallBound(double r, double K, double T) {
    return [r, K, T](double x, double t)->double {
        double res = exp(x) - K * exp(-r * (T - t));
        return (res >= 0) ? res : 0;
    };
}

pdefunc generatePutBound(double r, double K, double T) {
	return [r, K, T](double x, double t)->double {
		double res = K * exp(-r * (T - t)) - exp(x);
		return (res >= 0) ? res : 0;
	};
}

pdefunc generateConstant(double c) {
	return [c](double x, double t)->double {
		return c;
	};
}

vector<double> generateMesh(double start, double end, int size) {
    double step = (end - start) / (size - 1);
    vector<double> res(size);
    for (int i = 0; i < size; ++i) {
        res[i] = start + i * step;
    }
    return res;
}

vector<double> generateCenteredMesh(double center, double halfLen, double halfN) {
	vector<double> res(2 * halfN + 1);
	double step = halfLen / halfN;
	for (int i = 0; i < 2 * halfN + 1; ++i) {
		res[i] = center - halfLen + i * step;
	}
	return res;
}

vector<double> generateSpotMesh(double spot, double std, int n) {
	return generateCenteredMesh(log(spot), 5 * std, n);
}

PDE generate_BS_PDE(pdefunc r, pdefunc sigma, vector<double> meshX, vector<double> meshT, PDEBounds bounds, double theta, bool constCoef, bool constBound) {
	PDECoefs coefs = {
		[sigma](double x,double t) {return -.5 * sigma(x,t) * sigma(x,t); },
		[sigma,r](double x,double t) {return .5 * sigma(x,t) * sigma(x,t) - r(x,t); },
		r,
		generateConstant(0),
	};
	return PDE(coefs, meshX, meshT, bounds, theta, constCoef, constBound);
}
