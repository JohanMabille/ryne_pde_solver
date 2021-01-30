#pragma once
#include "PDE.h"

void delta_csv(const PDE& pde, string fname, bool convert);
vector<double> delta(const PDE& pde, int t_idx, bool convert);
void gamma_csv(const PDE& pde, string fname, bool convert);
vector<double> gamma(const PDE& pde, int t_idx, bool convert);
void theta_csv(const PDE& pde, string fname, bool convert);
vector<double> theta_greek(const PDE& pde, int x_idx);
vector<vector<vector<double>>> sigmaIter(const PDE& pde, const PDE& solution, vector<double> meshSigma);
vector<double> vega_xt(const PDE& pde, const vector<vector<vector<double>>>& sigmaIt, vector<double> meshSigma, int x_idx, int t_idx) ;
void vega_t_csv(const PDE& pde, string fname, vector<double> meshSigma, const vector<vector<vector<double>>>& sigmaIt, int t_idx) ;
void vega_x_csv(const PDE& pde, string fname, vector<double> meshSigma, const vector<vector<vector<double>>>& sigmaIt, int x_idx) ;
