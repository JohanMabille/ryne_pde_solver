#pragma once
#include <iostream>
#include <algorithm>
#include "utils.h"

using std::string;

// More C++14: using pdefunc = std::function<double(double, double)>;
typedef function<double(double, double)> pdefunc; // pdefunc are functions of x and t 

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

// This is C style coding
// a simple struc PDECoefs {}; is enough in C++
// A hierachy of classes with virtual methods for computing
// the coefficients would be more appropriate: it allows you
// some optimizations such as precomputing some coefficients,
typedef struct PDECoefs { // PDECoefs is a strcuture of 4 functions corresponding to coefficients a, b, c and d
	pdefunc a;
	pdefunc b;
	pdefunc c;
	pdefunc d;
} PDECoefs;

// A hierarchy of classes would have been more expressive
// (you could have classes for DIrichlet and NeuMann inheriting
// from an abstract boundary condition class)
typedef struct PDEBounds { // PDEBounds is a structure of 3 functions of x and t 
	pdefunc xmin;  // xmin is a function of t and x which applies to the lower bound
	pdefunc xmax; // xmax is a function of t and x which applies to the upper bound
	pdefunc tbound; // tbound is a function of x which applies at time T (payoff of the option)
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
    // An improvement would be to split this class
    // into a pure solver class and a PDE Grid class
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
        // This should be a private data member. A public const method should
        // be provided here for users that need ot read the values
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
	bool isConst = false,bool constBound = false);
	double alpha(double x, double t, double theta) const;
	double beta(double x, double t, double theta) const;
	double gamma(double x, double t, double theta) const;
	vector<vector<double>> getWeight(int n, double theta) const;
	vector<double> getBias(int n, double theta) const;
	vector<double> getBiasfromCache(int n, bool current) const;
	stepMat getStepMatrices(int n) const;
	void solve();
	void to_csv(string fname) const;
        // PDE should provide a public API that avoids all these friend functions
	friend void delta_csv(const PDE& pde, string fname, bool convert);
	friend vector<double> delta(const PDE& pde, int t_idx, bool convert);
        friend void gamma_csv(const PDE& pde, string fname, bool convert);
        friend vector<double> gamma(const PDE& pde, int t_idx, bool convert);
	friend void theta_csv(const PDE& pde, string fname, bool convert);
	friend vector<double> theta_greek(const PDE& pde, int x_idx);
	friend vector<vector<vector<double>>> sigmaIter(const PDE& pde, vector<double> meshSigma);
        friend void vega_t_csv(const PDE& pde, string fname, vector<double> meshSigma, const vector<vector<vector<double>>>& sigmaIt, int t_idx);
        friend void vega_x_csv(const PDE& pde, string fname, vector<double> meshSigma, const vector<vector<vector<double>>>& sigmaIt, int x_idx);
};

