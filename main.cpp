#include "closed_form.hpp"
#include <iostream>
#include "PDE.h"
#include <functional>
#include <cassert>
#include "blackscholes.h"
#include "greeks.h"
#include <cmath>

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


int main(int argc, const char * argv[])
{
	double spot = 100;	    // Spot price 
	double sigma = .2;	   // Volatility 
	double r = .0;		  // Interest rate
	double T = 1.;		 // Maturity in years
	double nT = 10000;  // Size of the time mesh
	double K = 100;    // Strike Price 
	double nX = 50;   // will result in 2*Nx +1 values in the space mesh centered in log(spot)

	// For the method to converge, we should have nT that is greater than nX^2 

	double theta = .5;		// theta of the scheme
	bool isCall = true;		

	vector<double> meshX = generateSpotMesh(spot, sigma * std::sqrt(T), nX);
	vector<double> meshT = generateMesh(0, T, nT * T);

	PDEBounds CallBounds = {	//structure explained in PDE.h
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
	std::cout << "Price: " << solution.values[0][nX-1] << std::endl;
	std::cout << "Price with closed form: " << dauphine::bs_price(spot, K, sigma, T, isCall) << std::endl;
	
	
	solution.to_csv("test.csv");
	delta_csv(solution, "delta.csv");
    gamma_csv(solution,"gamma.csv");
	theta_csv(solution, "theta.csv");

	int nS = 3;		// like for nX, will result in a volatility mesh with 2*NS +1 values, centered in sigma
	double delta_vol = 0.1;   // will result in a volatility mesh with bounds = sigma +/- delta_vol
	vector<double> meshSigma = generateCenteredMesh(sigma, delta_vol, nS);
	vector<vector<vector<double>>> sigmaRes = sigmaIter(solution, meshSigma);

	vega_t_csv(solution, "vega.csv", meshSigma, sigmaRes, 0);
	vector<double> vega0 = vega_xt(sigmaRes, meshSigma, nX - 1, 0);
	std::cout << "Vega (for 1% vol move):" << vega0[nS] / 100 << std::endl;
	
	std::cin.get();
    return 0;
}
