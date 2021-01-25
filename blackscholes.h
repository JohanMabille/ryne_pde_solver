#pragma once
#include "PDE.h"

pdefunc generateConstant(double c);

PDE generate_BS_PDE(pdefunc r, pdefunc sigma, vector<double> meshX, vector<double> meshT, PDEBounds bounds, double theta, bool constCoef = false, bool constBound = false);