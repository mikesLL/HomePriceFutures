// 
// Copyright A. Michael Sharifi
#include "headers.h"

snodes::snodes() {

	// initialize gammat state-state transition matrix
	vector<vector<vector<double>>> zeros_T_NS_NS(T_max + 1, vector<vector<double>>(n_s, vector<double>(n_s, 0.0)));
	gammat = zeros_T_NS_NS;

	// initialize i2s and s2i mappings
	vector<vector<vector<int>>> zeros_NPH_NRENT_NYI(n_ph, vector<vector<int>>(n_rent, vector<int>(n_yi, 0)));
	vector<int> zeros_NS(n_s, 0);

	i2s_map = zeros_NPH_NRENT_NYI;
	s2i_ph = zeros_NS;
	s2i_rent = zeros_NS;
	s2i_yi = zeros_NS;

	i_s = 0;

	for (i_ph = 0; i_ph < n_ph; i_ph++) {
		for (i_rent = 0; i_rent < n_rent; i_rent++) {
			for (i_yi = 0; i_yi < n_yi; i_yi++) {
				i2s_map[i_ph][i_rent][i_yi] = i_s;   // initialize i2s map
				s2i_ph[i_s] = i_ph;                 // initialize s2_i maps
				s2i_rent[i_s] = i_rent;
				s2i_yi[i_s] = i_yi;
				i_s++;
			}
		}
	}

	// here: compute i_s_mid;
	int i_ph_mid = (int) floor( .5*n_ph );
	int i_rent_mid = (int) floor( .5*n_rent );
	int i_yi_mid = (int) floor( .5*n_yi );

	i_s_mid = i2s_map[i_ph_mid][i_rent_mid][i_yi_mid];

    for( i_ph = 0; i_ph < n_ph; i_ph++ ){
		s_ph_midry[i_ph] = i2s_map[i_ph][i_rent_mid][i_yi_mid];
    }

}


void snodes::adj_tax() {

	// variable of interest: 
	// double yi_gridt[T_max + 1][n_yi]; 

	int i_t2 = 0;
	int i_y2 = 0;

	double tax_brack[] = {0.0, 0.4385, 1.0595, 1.61, 2.88 };
	double tax_rate[] = { 0.15, 0.28, 0.31, 0.36, 0.396 };

	double y_btax = 0.0;
	double y_atax = 0.0;
	double y_diff = 0.0;

	for (i_t2 = 0; i_t2 < (T_max + 1); i_t2++) {
		for (i_y2 = 0; i_y2 < n_yi; i_y2++){ 
			
			y_btax = yi_gridt[i_t2][i_y2];           // load pre-tax income
			y_atax = 0.0;

			if (y_btax >= tax_brack[0]) {
				y_diff = min(y_btax - tax_brack[0], tax_brack[1] - tax_brack[0] );
				y_atax = y_atax + tax_rate[0] * y_diff;				
			}

			if (y_btax >= tax_brack[1]) {
				y_diff = min( y_btax - tax_brack[1], tax_brack[2] - tax_brack[1] );
				y_atax = y_atax + tax_rate[1] * y_diff;
			}

			if (y_btax >= tax_brack[2]) {
				y_diff = min(y_btax - tax_brack[2], tax_brack[3] - tax_brack[2]);
				y_atax = y_atax + tax_rate[2] * y_diff;
			}

			if (y_btax >= tax_brack[3]) {
				y_diff = min(y_btax - tax_brack[3], tax_brack[4] - tax_brack[3]);
				y_atax = y_atax + tax_rate[3] * y_diff;
			}

			if (y_btax >= tax_brack[4]) {
				y_diff = y_btax - tax_brack[4];
				y_atax = y_atax + tax_rate[4] * y_diff;
			}

			yi_gridt[i_t2][i_y2] = y_tax;

		}
	}

}
