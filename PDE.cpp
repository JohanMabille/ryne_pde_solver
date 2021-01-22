#include "PDE.h"

double PDE::alpha(double x, double t, double theta) const {
	return (dt * theta) * (coef.a(x, t) / dx - coef.b(x, t) / 2) / dx;
}

double PDE::beta(double x, double t, double theta) const {
	return (dt * theta) * (coef.c(x, t) - 2 * coef.a(x, t) / (dx * dx)) + 1;
}

double PDE::gamma(double x, double t, double theta) const {
	return (dt * theta) * (coef.a(x, t) / dx + coef.b(x, t) / 2) / dx;
}