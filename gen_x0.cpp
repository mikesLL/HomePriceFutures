/*
gen_x0.cpp

Given coh (cash on hand), b_min (lower bound for bonds / mortgage borrowing), ph_i1 (current price), t_i2 (tenure in next period),
get_x0_loop loops through choice space in C,B,X,CSFp,CSFn to find an initial guess for the policy that maximizes the value fn
then get_x0_loop runs a hill-climbing algorithm to find the optimal choice values and returns these as a vector

Copyright A. Michael Sharifi
*/

#include "headers.h"

vector<double> gen_x0(double coh_in, double b_min, void *vf1_in, void *vf2_in, void *ufnEV21_in, vector<double> x0_in) {
	
	vfn *vf1 = (vfn *)vf1_in;             // pointer to vf1 or V1 value fn
	vfn *vf2 = (vfn *)vf2_in;             // pointer to vf2 or V2 value fn

	ufnEV2 * ufnEV21 = (ufnEV2*)ufnEV21_in;

	double coh = coh_in;
	double csf_min = 0.0;

	double h_tol = 0.001; //.0005; //  50 dollar accuracy

	int opt_flag = 1, nds = 4, it_max = 100000;
	int N_control2 = N_control - 3;  // For part without csf
	int N_control3 = N_control - 1;  // for part with csf
	int N_controlh;
	double h_step0 = 0.2;
	double h_step_mult = 0.25;
	double h_step = h_step0;

	vector<double> x0 = x0_in;    
	vector<double> x0_nocsf = x0_in;
	vector<double> x0_h(5, 0.0);
	vector<double> x1(5, 0.0);
	
	double lb[] = { 0.01, b_min, 0.0, 0.0, 0.0 };
	
	int it = 0;
	double v0, v0_h, v1, vi_max, vi_min;
	double v0_nocsf;

	int i_max, i_min, i_min_flag, i, j, i_inner_loop;
	//int i_csf;
	//double vi_csf;

	// starting guess for v0, adjusted to current coh, b_min
	csf_min = min(x0[3], x0[4]);
	x0[3] = x0[3] - csf_min;
	x0[4] = x0[4] - csf_min;
	x0[1] = max(x0[1], b_min);                       
	x0[0] = coh - x0[1] - x0[2] - x0[3] - x0[4];
	v0 = (*ufnEV21).eval(x0);       

	x0_nocsf[3] = 0.0;
	x0_nocsf[4] = 0.0;
	x0_nocsf[1] = max(x0_nocsf[1], b_min);
	x0_nocsf[0] = coh - x0_nocsf[1] - x0_nocsf[2];
	v0_nocsf = (*ufnEV21).eval(x0_nocsf);
	
	double h_diff = 0.001;
	//
	// First Loop
	//
	opt_flag = 0;
	while (opt_flag) {

		i_max = 0;
		i_min = 0;
		i_min_flag = 0;
		
		vi_max = -1.0e6;
		vi_min = 1.0e6;

		x1 = x0_nocsf;
		csf_min = min(x1[3], x1[4]);
		x1[3] = x1[3] - csf_min;
		x1[4] = x1[4] - csf_min;
		x1[1] = max(x1[1], b_min);
		x1[0] = coh - x1[1] - x1[2] - x1[3] - x1[4];
		v1 = (*ufnEV21).eval(x1);

		if (v1 > v0_nocsf) {
			x0_nocsf = x1;
			v0_nocsf = v1;
		}
		
		// find min, max partials; restrict to x: 0,1,2
		for (i = 0; i < N_control2; i++) {
			x0_h = x0_nocsf;
			x0_h[i] = x0_nocsf[i] + h_step;

			v0_h = (*ufnEV21).eval(x0_h);

			if (v0_h > vi_max) {
				i_max = i;
				vi_max = v0_h;
			}

			if ( (v0_h < vi_min) && ((x0_nocsf[i] - h_step) >= lb[i])) {
				i_min = i;
				i_min_flag = 1;
				vi_min = v0_h;
			}
		}
		
		if (i_min_flag >= 1) {
			x1 = x0_nocsf;
			x1[i_max] = x0_nocsf[i_max] + h_step;
			x1[i_min] = x0_nocsf[i_min] - h_step;
			v1 = (*ufnEV21).eval(x1);
		} 
		
		if (v1 > v0_nocsf) {
			v0_nocsf = v1;
			x0_nocsf = x1;
		} else {
			h_step = h_step * h_step_mult;
		}

		it++;

		if ( (it > it_max) || (h_step < h_tol) ) {
			opt_flag = 0;
		}
	}


	if (csfLev > 0.0) {
		h_step = h_step0;
		it = 0;
		opt_flag = 1;
	}

	//
	// Second Loop
	//
	opt_flag = 1;
	N_controlh = N_control2; 
	while (opt_flag) {

		i_max = 0;
		i_min = 0;
		i_min_flag = 0;

		vi_max = -1.0e6;
		vi_min = 1.0e6;

		//i_csf = 3;
		//vi_csf = -1.0e6;

		x1 = x0;
		csf_min = min(x1[3], x1[4]);
		x1[3] = x1[3] - csf_min;
		x1[4] = x1[4] - csf_min;
		x1[1] = max(x1[1], b_min);
		x1[0] = coh - x1[1] - x1[2] - x1[3] - x1[4];
		v1 = (*ufnEV21).eval(x1);

		if (v1 > v0) {
			x0 = x1;
			v0 = v1;
		}

		x1 = x0;
		x1[3] = 0.0;
		x1[4] = 0.0;
		x1[1] = max(x1[1], b_min);
		x1[0] = coh - x1[1] - x1[2] - x1[3] - x1[4];
		v1 = (*ufnEV21).eval(x1);

		if (v1 > v0) {
			x0 = x1;
			v0 = v1;
		}

		// find min, max partials; restrict to x: 0,1,2
		if (h_step <= 0.01) {
			N_controlh = N_control3; 
		}

		for (i = 0; i < N_controlh; i++) {
			x0_h = x0;
			x0_h[i] = x0[i] + h_step;

			v0_h = (*ufnEV21).eval(x0_h);

			if (v0_h > vi_max) {
				i_max = i;
				vi_max = v0_h;
			}

			if ((v0_h < vi_min) && ((x0[i] - h_step) >= lb[i])) {
				i_min = i;
				i_min_flag = 1;
				vi_min = v0_h;
			}
		}

		if (i_min_flag >= 1) {
			x1 = x0;
			x1[i_max] = x0[i_max] + h_step;
			x1[i_min] = x0[i_min] - h_step;
			v1 = (*ufnEV21).eval(x1);
		}

		if (v1 > v0) {
			v0 = v1;
			x0 = x1;
		}
		else {
			h_step = h_step * h_step_mult;
		}

		it++;

		if ((it > it_max) || (h_step < h_tol)) {
			opt_flag = 0;
		}
	}

	if (v0_nocsf > v0) {
		v0 = v0_nocsf;
		x0 = x0_nocsf;
	}

	return x0;
}




/*
if (i_inner_loop) {
v1 = v0_nocsf; // reset guess for x1
x1 = x0_nocsf;

//for (i = 0; i < N_control3; i++) {
for (i = 0; i < N_control2; i++) {
for (j = 0; j < N_control2; j++) {

if (i != j) {
// confirm guess is valid with lower bounds
if ( (x0_nocsf[j] - h_step) >= lb[j]) {
x0_h = x0_nocsf;
x0_h[i] = x0_nocsf[i] + h_step;
x0_h[j] = x0_nocsf[j] - h_step;

v0_h = (*ufnEV21).eval(x0_h);

if (v0_h > v1) {
v1 = v0_h;
x1 = x0_h;
}
}
}
}
}

if (v1 > v0_nocsf) {
v0_nocsf = v1;
x0_nocsf = x1;
} else {
h_step = h_step * h_step_mult;
}
}
*/