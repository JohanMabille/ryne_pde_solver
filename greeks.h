#pragma once
#include "PDE.h"

void delta_csv(const PDE& pde, string fname, bool convert);
vector<double> delta(const PDE& pde, int t_idx, bool convert);
