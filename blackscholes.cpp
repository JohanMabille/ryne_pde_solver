#include "blackscholes.h"


pdefunc generateConstant(double c) {
	return [c](double x, double t)->double {
		return c;
	};
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