#include "closed_form.hpp"
#include <iostream>
#include "PDE.h"
#include <functional>
#include <cassert>
#include "blackscholes.h"
#include "greeks.h"

// Guidelines:
//
// 1] Equation
// Even if the volatility and the rate are constant in the BS model,
// we can prove that, under certain assumptions, we can use the same
// PDE equation with volatility surface and rate curves. We could also
// want to take into account the repo (even if it could theoretically
// be part of the r factor). Therefore, it is more generic to solve
// the following generic equation;
// df/dt = a(x,t)d2f/dx2 + b(x,t)df/dx + c(x, t)f + d(x, t).
// The pricer should ask the coefficients a, b, c and d to an
// abstract class that will be inherited by classes implementing
// different models.
// 
// 2] Payoff
// The pricer should be able to price exotic options, such as
// barriers or asian options. To do so, the pricer must be able
// to call an external function between each step. Define an API
// that allows to register a function or an abstract class modeling
// a payoff.

using namespace std;

int main(int argc, const char * argv[])
{
	double spot = 100;

	double sigma = .2;
	double r = .0;

	double T = 1.;
	double nT = 10000;
	double K = 100;
	double nX = 50;

	//there is a condition on nX^2 and nT for stability purposes

	double theta = .5;
	bool isCall = true;

	vector<double> meshX = generateSpotMesh(spot, sigma * sqrt(T), nX);
	vector<double> meshT = generateMesh(0, T, nT * T);

	PDEBounds CallBounds = {
		generateCallBound(r,K,T),
		generateCallBound(r,K,T),
		generateCallPayoff(K)
	};

	PDEBounds PutBounds = {
		generatePutBound(r,K,T),
		generatePutBound(r,K,T),
		generatePutPayoff(K)
	};

	PDEBounds bounds = (isCall) ? CallBounds : PutBounds;

	PDE solution = generate_BS_PDE(
		generateConstant(r),
		generateConstant(sigma),
		meshX,
		meshT,
		bounds,
		theta,
		true,
		false
	);

	solution.solve();
	cout << "price: " << solution.values[0][nX-1] << endl;
	cout << "price: " << dauphine::bs_price(spot, K, sigma, T, isCall) << endl;
	
	
	solution.to_csv("test.csv");
	delta_csv(solution, "delta.csv");
    gamma_csv(solution,"gamma.csv");
	theta_csv(solution, "theta.csv");
	vector<double> meshSigma = generateMesh(0.1, 0.5, 20);
	vector<vector<vector<double>>> sigmaRes = sigmaIter(solution, meshSigma);
    vega_t_csv(solution,"vega.csv", meshSigma, sigmaRes, 0);
	cin.get();
    return 0;
}
