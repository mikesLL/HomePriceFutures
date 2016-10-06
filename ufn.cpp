/*
ufc: felicitous utility function 
accepts double c_in (non-durable consumption), double hu_in (housing consumption)
returns double cobb-douglas evaluated at those values

_cd, calpha_cd, rho, are parameters initialized in the calibration.h

Copyright A. Michael Sharifi, 2016
*/

#include "headers.h"

double ufn( double c_in, double hu_in, int pref_in){
	int pref = pref_in;
	double c_comp;
	double uc;

	if (c_in > 0.0) {                                                 // case: almost negative consumption

		if (pref == 0) {
			c_comp = pow(c_in, alpha_cd) * pow(hu_in, calpha_cd);     // composite utility in the current period is Cobb-Douglas
			uc = 1.0 / (1.0 - rho)* pow(c_comp, 1.0 - rho);            // power utility in composite
		}
		else {
			c_comp = pow(c_in, alpha_ces) + gammad*pow(hu_in*1.0/100.0, alpha_ces);
			uc = 1.0 / (1.0 - rho) * pow(c_comp, (1.0 - rho) / alpha_ces);
		}
	}
	else {
		uc = -1.0e8*(1.0 + pow(c_in, 2.0));
	}

	return uc;
}
