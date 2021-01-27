#pragma once
#include <iostream>
#include <algorithm>
#include "utils.h"

using std::string;

typedef function<double(double, double)> pdefunc;

/*
*	2nd degree PDE df/dt = af''+bf'+cf+d
*	PDE coefs a,b,c,d are functions of  (x,t)
*	In a matrix form, the linear combination 
	of the implicit and explict schemes can be written
*	as follows:
*   Af(n+1)+B = Cf(n)+D
*	with A and C called weight matrixes and
*	B and D vectors called bias matrixes
*	Moreover, since A and C are sparse tridiagonal matrixes,
*	they will be represented as three vectors called 
*	gammaVec,betaVec, alphaVec in the getWeights() method
*/

typedef struct PDECoefs {
	pdefunc a;
	pdefunc b;
	pdefunc c;
	pdefunc d;
} PDECoefs;

typedef struct PDEBounds {
	pdefunc xmin;
	pdefunc xmax;
	pdefunc tbound;
} PDEBounds;

typedef struct stepMat {
	vector<vector<double>> A;
	vector<double> B;
	vector<vector<double>> C;
	vector<double> D;
}stepMat;



class PDE
{
private:
	PDECoefs coef;
	vector<double> meshX;
	vector<double> meshT;
	vector<double> meshX_center;
	vector<double> meshX_alpha;
	vector<double> meshX_gamma;
	PDEBounds bound;
	double dx;
	double dt;
	double theta;
	int current;
	bool isConst;
	bool constBound;
	vector<vector<double>> cacheA;
	vector<double> cacheB;
	vector<vector<double>> cacheC;
	vector<double> cacheD;
	void step();


public:
	vector<vector<double>> values;
	/**
	 * Generic PDE solver for 2nd order PDEs
	 *
	 * @param coef the pde coefficients ( as function of (x,t) ) ( df/dt = af''+bf'+cf+d)
	 * @param meshX the X mesh ( assumes the mesh has a constant dx step )
	 * @param meshT the T mesh ( assumes the mesh has a constant dt step )
	 * @param bound the boundary conditions and payoff
	 * @param theta theta scheme coefficient
	 * @param isConst whether pde coefficients are constant. if true caches the weight matrix
	 * @param constBound whether boundary conditions are constant. for further caching
	 */

	PDE(PDECoefs coef, vector<double> meshX, vector<double> meshT, PDEBounds bound, double theta, 
	bool isConst = false,bool constBound = false));
	double alpha(double x, double t, double theta) const;
	double beta(double x, double t, double theta) const;
	double gamma(double x, double t, double theta) const;
	vector<vector<double>> getWeight(int n, double theta) const;
	vector<double> getBias(int n, double theta) const;
	vector<double> getBiasfromCache(int n, bool current) const;
	stepMat getStepMatrices(int n) const;
	void solve();
	void to_csv(string fname) const;
};

