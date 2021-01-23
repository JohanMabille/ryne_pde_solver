#include "PDE.h"

using namespace std;

double PDE::alpha(double x, double t, double theta) const {
	return (dt * theta) * (coef.a(x, t) / dx - coef.b(x, t) / 2) / dx;
}

double PDE::beta(double x, double t, double theta) const {
	return (dt * theta) * (coef.c(x, t) - 2 * coef.a(x, t) / (dx * dx)) + 1;
}

double PDE::gamma(double x, double t, double theta) const {
	return (dt * theta) * (coef.a(x, t) / dx + coef.b(x, t) / 2) / dx;
}

PDE::PDE(PDECoefs coef, vector<double> meshX, vector<double> meshT, PDEBounds bound, double theta) :
	coef(coef), meshX(meshX), meshT(meshT), bound(bound), theta(theta) {
	this->dx = meshX[1] - meshX[0];
	this->dt = meshT[1] - meshT[0];
	this->meshX_center = vector<double>(meshX.begin() + 1, meshX.end() - 1);
	this->meshX_alpha = vector<double>(meshX.begin() + 1, meshX.end() - 2);
	this->meshX_gamma = vector<double>(meshX.begin() + 2, meshX.end() - 1);
	this->values = vector<vector<double>>(meshT.size());
	this->values[meshT.size() - 1] = apply2D(bound.tbound, this->meshX_center, *(meshT.end() - 1));
	current = meshT.size() - 2;
}

vector<vector<double>> PDE::getWeight(int n, double theta) const {
	// Af(n+1)+B = Cf(n)+D
	//instead of using sparse weight matrixes (3 diagonals only), we represent a tridiagonal matrix with three diagonal vectors
	//this method computes either A or C depending on the n and theta inputs
	double t = this->meshT[n];

	vector<double> betaVec = apply2D([this, theta](double x, double t) {return this->beta(x, t, theta); }, this->meshX_center, t);
	vector<double> alphaVec = apply2D([this, theta](double x, double t) {return this->alpha(x, t, theta); }, this->meshX_alpha, t);
	vector<double> gammaVec = apply2D([this, theta](double x, double t) {return this->gamma(x, t, theta); }, this->meshX_gamma, t);

	return { gammaVec,betaVec, alphaVec };
}

vector<double> PDE::getBias(int n, double theta) const {
	// Af(n+1)+B = Cf(n)+D
	//this function computes either B or D depending on the n and theta inputs
	double t = this->meshT[n];
	int size = this->meshX.size();
	vector<double> dVec = apply2D([theta, this](double x, double t) {return (dt * theta) * this->coef.d(x, t); }, this->meshX_center, t);
	dVec[0] += this->alpha(this->meshX[0], t, theta) * this->bound.xmin(this->meshX[0], t);
	dVec[size - 3] += this->gamma(this->meshX[size - 1], t, theta) * this->bound.xmax(this->meshX[size - 1], t);
	return dVec;
}