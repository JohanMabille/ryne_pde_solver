#pragma once
#include "PDE.h"

pdefunc generateCallPayoff(double K);
pdefunc generatePutPayoff(double K);
pdefunc generateCallBound(double r, double K, double T);
pdefunc generatePutBound(double r, double K, double T);
pdefunc generateConstant(double c);
vector<double> generateMesh(double start, double end, int size);
vector<double> generateSpotMesh(double spot, double std, int n);

PDE generate_BS_PDE(pdefunc r, pdefunc sigma, vector<double> meshX, vector<double> meshT, PDEBounds bounds, double theta, bool constCoef = false, bool constBound = false);
