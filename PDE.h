#pragma once
#include <iostream>
#include "utils.h"

using namespace std;

typedef function<double(double, double)> pdefunc;

/*
*	2nd degree PDE df/dt = af''+bf'+cf+d
*	PDE coefs a,b,c,d are functions of  (x,t)
*/

typedef struct PDECoefs {
	pdefunc a;
	pdefunc b;
	pdefunc c;
	pdefunc d;
} PDECoefs;


class PDE
{
private:
	PDECoefs coef;
	vector<double> meshX;
	vector<double> meshT;
	double dx;
	double dt;
	double theta;
public:
	double alpha(double x, double t, double theta) const;
	double beta(double x, double t, double theta) const;
	double gamma(double x, double t, double theta) const;
};

