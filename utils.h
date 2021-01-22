#pragma once
#include <functional>
#include <vector>

using std::function;
using std::vector;

vector<double> apply1D(const function<double(double)>& f, const vector<double>& x);
vector<double> apply2D(const function<double(double, double)>& f, const vector<double>& x, double t);
vector<vector<double>> trigMat(const vector<vector<double>>& A, const vector<double>& b);
