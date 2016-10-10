// vfn.cpp
// Copyright A. Michael Sharifi, 2016

#include "headers.h"

void vfn::enter_data(void *snodes_in, double phr_in, int t_id_in, int t_num_in, double csf_1yr_in, int pref_in, int T_max_in) {

	T_max = T_max_in; 
	snodes1 = (snodes *)snodes_in;

	csf_1yr = csf_1yr_in;
	t_id = t_id_in;                 // current year 
	phr = phr_in;                   // rent 
	t_num = t_num_in;               // planning horizon
	pref = pref_in;
	lcount = 0;
	def_flag = 0;

	cout << "t_num = " << t_num << endl;
	cout << "t_id = " << t_id << endl;

	vector<vector<vector<int>>> zeros_int_TN_NS_WN(t_n, vector<vector<int>>(n_s, vector<int>(w_n, 0)));
	vector<vector<vector<double>>> zeros_TN_NS_WN(t_n, vector<vector<double>>(n_s, vector<double>(w_n, 0.0)));
	vector<vector<vector<double>>> neg_TN_NS_WN(t_n, vector<vector<double>>(n_s, vector<double>(w_n, num_small)));

	vector<vector<double>> neg_NS_WN(n_s, vector<double>(w_n, -1.0e6));

	vector<double> zeros_WN(w_n, 0.0);

	xt_grid = zeros_int_TN_NS_WN;         // tenure control grid; ordered: tenure, state, wealth
	x1_grid = zeros_TN_NS_WN;             // consumption control grid; ordered: tenure, state, wealth
	x2_grid = zeros_TN_NS_WN;             // bond control grid; ordered: tenure, state, wealth
	x3_grid = zeros_TN_NS_WN;             // stock control grid; ordered: tenure, state, wealth
	x4_grid = zeros_TN_NS_WN;             // futures:+ control grid; ordered: tenure, state, wealth
	x5_grid = zeros_TN_NS_WN;             // futures:- control grid; ordered: tenure, state, wealth
	
	vw3_grid = neg_TN_NS_WN;            // value function; ordered: tenure, state, wealth
	vw3_def_grid = neg_NS_WN;            // value function; ordered: tenure, state, wealth

	vw3_d_grid = zeros_TN_NS_WN;            // value function: 1st deriv
	vw3_dd_grid = zeros_TN_NS_WN;            // value function: 2nd deriv

	lambda_grid = zeros_TN_NS_WN;         // value function interpolation parameters; ordered: tenure, state, wealth

	v_move = zeros_WN;
	vw3_grid_norm = zeros_WN;

	// wealth grid
	for (int w_i = 0; w_i < w_n; w_i++) {
		w_grid[w_i] = w_min + (w_max - w_min)* double(w_i) / (double(w_n) - 1.0);
	}

}


void vfn::get_pol(int i_t_in, int i_s_in, int i_w_in, vector<double> &x_pol) {

	if (i_w_in >= 0) {
		x_pol[0] = x1_grid[i_t_in][i_s_in][i_w_in];
		x_pol[1] = x2_grid[i_t_in][i_s_in][i_w_in];
		x_pol[2] = x3_grid[i_t_in][i_s_in][i_w_in];
		x_pol[3] = x4_grid[i_t_in][i_s_in][i_w_in];
		x_pol[4] = x5_grid[i_t_in][i_s_in][i_w_in];
	}
	else {
		x_pol[0] = c_fs; x_pol[1] = 0.0; x_pol[2] = 0.0; x_pol[3] = 0.0; x_pol[4] = 0.0;
	}
}


void vfn::set_pol_ten_v(int i_t_in, int i_s_in, int i_w_in, vector<double> &x_pol, int t_i2_in, double v0_in) {
	x1_grid[i_t_in][i_s_in][i_w_in] = x_pol[0];
	x2_grid[i_t_in][i_s_in][i_w_in] = x_pol[1];
	x3_grid[i_t_in][i_s_in][i_w_in] = x_pol[2];
	x4_grid[i_t_in][i_s_in][i_w_in] = x_pol[3];
	x5_grid[i_t_in][i_s_in][i_w_in] = x_pol[4];
	xt_grid[i_t_in][i_s_in][i_w_in] = t_i2_in;       // optimal tenure choice
	vw3_grid[i_t_in][i_s_in][i_w_in] = v0_in;
}

// input: tenure state t_i, employment state y_i, home price state ph_i, wealth w_in
// output: value fn evaluated at the above state
eval_res vfn::eval_v(int i_t_in, int i_s_in, double w_in) {
	eval_res res1;                                       // v_out structure
	int w_i_low;
	double w_i_d;
	double alpha1;
	double lambda1;
	double v_tlower;
	double v_tupper;
	double num1, den1, w_diff1, w_diff2;

	if ( (w_in >= w_min ) && (w_in < w_max) ) {
		w_i_d = (double) (w_n - 1.0)*(w_in - w_min) / (w_max - w_min);    // map w_in to w_i
		w_i_low = (int) floor(w_i_d);
		
		w_diff2 = w_grid[1] - w_grid[0];

		alpha1 = ( w_in - w_grid[w_i_low] ) / w_diff2;

		v_tlower = vw3_grid[i_t_in][i_s_in][w_i_low];
		

		v_tupper = vw3_grid[i_t_in][i_s_in][w_i_low + 1];
		
		res1.v_out = v_tlower + alpha1 * (v_tupper - v_tlower);
		

		if (res1.v_out != res1.v_out) {
			res1.v_out = v_tlower;
		}
		res1.w_i_floor = w_i_low;
		res1.v_i_floor = vw3_grid[i_t_in][i_s_in][w_i_low];  // TODO: double check; this may be an issue when considering tenure changes
	}
	else {
		if (w_in >= w_max) {
			w_i_low = w_n - 1;
			v_tlower = vw3_grid[i_t_in][i_s_in][w_i_low];
			

			res1.v_out = max(vw3_grid[i_t_in][i_s_in][w_i_low], v_tlower);

		} else {
			res1.v_out = vw3_grid[i_t_in][i_s_in][0] - 1.0e6*pow(w_in - w_min, 2);
			res1.w_i_floor = 0;
			res1.v_i_floor = res1.v_out;
		}
	}


	return res1;
}


eval_res vfn::eval_v_def(int i_s_in, double w_in) {
	eval_res res1;                                       // v_out structure
	int w_i_low;
	double w_i_d;
	double alpha1;
	double lambda1;
	double v_tlower;
	double v_tupper;
	double num1, den1, w_diff1, w_diff2;

	if ((w_in >= w_min) && (w_in < w_max)) {
		w_i_d = (double)(w_n - 1.0)*(w_in - w_min) / (w_max - w_min);    // map w_in to w_i
		w_i_low = (int)floor(w_i_d);

		w_diff2 = w_grid[1] - w_grid[0];

		alpha1 = (w_in - w_grid[w_i_low]) / w_diff2;

		v_tlower = vw3_def_grid[i_s_in][w_i_low];
		v_tupper = vw3_def_grid[i_s_in][w_i_low + 1];

		res1.v_out = v_tlower + alpha1 * (v_tupper - v_tlower);

		if (res1.v_out != res1.v_out) {
			res1.v_out = v_tlower;
		}
		res1.w_i_floor = w_i_low;
		res1.v_i_floor = vw3_def_grid[i_s_in][w_i_low];  // TODO: double check; this may be an issue when considering tenure changes
	}
	else {
		if (w_in >= w_max) {
			w_i_low = w_n - 1;
			v_tlower = vw3_def_grid[i_s_in][w_i_low];
			res1.v_out = max(vw3_def_grid[i_s_in][w_i_low], v_tlower);

		}
		else {
			res1.v_out = vw3_def_grid[i_s_in][0] - 1.0e6*pow(w_in - w_min, 2);
			res1.w_i_floor = 0;
			res1.v_i_floor = res1.v_out;
		}
	}
	return res1;
}

eval_res vfn::eval_v_norm( double w_in ) {
	eval_res res1;                                       // v_out structure
	int w_i_low;
	double w_i_d;
	double alpha1;
	double lambda1;
	double v_tlower;
	double v_tupper;
	double num1, den1, w_diff1, w_diff2;

	if ((w_in >= w_min) && (w_in < w_max)) {
		w_i_d = (double)(w_n - 1.0)*(w_in - w_min) / (w_max - w_min);    // map w_in to w_i
		w_i_low = (int)floor(w_i_d);

		w_diff2 = w_grid[1] - w_grid[0];

		alpha1 = (w_in - w_grid[w_i_low]) / w_diff2;

		v_tlower = vw3_grid_norm[w_i_low];


		v_tupper = vw3_grid_norm[w_i_low + 1];

		res1.v_out = v_tlower + alpha1 * (v_tupper - v_tlower);


		if (res1.v_out != res1.v_out) {
			res1.v_out = v_tlower;
		}
		res1.w_i_floor = w_i_low;
		res1.v_i_floor = vw3_grid_norm[w_i_low];  // TODO: double check; this may be an issue when considering tenure changes
	}
	else {
		if (w_in >= w_max) {
			w_i_low = w_n - 1;
			v_tlower = vw3_grid_norm[w_i_low];


			res1.v_out = max(vw3_grid_norm[w_i_low], v_tlower);

		}
		else {
			res1.v_out = vw3_grid_norm[0] - 1.0e6*pow(w_in - w_min, 2);
			res1.w_i_floor = 0;
			res1.v_i_floor = res1.v_out;
		}
	}


	return res1;
}


// value function in the event agent moves to another area
eval_res vfn::eval_v_move(int i_t_in, int i_s_in, double w_in, int t_left) {
	eval_res res1;
	res1.v_out = -1e10;
	res1.v_i_floor = -1e10;
	res1.w_i_floor = 0;

	double v_out = -1e10;
	w_in = w_in + 1.0;

	double v_max_slope = 1.0 / (w_grid[w_n - 1] - w_grid[w_n - 2])
		*(v_move[w_n - 1] - v_move[w_n - 2]);

	double v_max_slope2 = 1.0 / (w_grid[w_n - 2] - w_grid[w_n - 3])
		*(v_move[w_n - 2] - v_move[w_n - 3]);

	double v_max_d2 = (v_max_slope - v_max_slope2) / (w_grid[w_n - 1] - w_grid[w_n - 2]);

	if ((w_in >= w_min) && (w_in < w_max)) {
		double w_i_dirty = (double)(w_n - 1)*(w_in - w_min) / (w_max - w_min);    // work here ;
		int w_i_low = int(floor(w_i_dirty));
		double a_h = w_i_dirty - double(w_i_low);

		v_out = v_move[w_i_low] +
			a_h*(v_move[min(w_i_low + 1, w_n - 1)] - v_move[w_i_low]);

		res1.w_i_floor = w_i_low;

		res1.v_i_floor = v_move[w_i_low];
	}
	else if (w_in < w_min) {
		v_out = v_move[0] - abs(w_in - w_min)*1.0e6;  // trying this

		res1.w_i_floor = 0;
		res1.v_i_floor = v_move[0] - abs(w_in - w_min) * 1.0e6;
	}
	else {
		v_out = v_move[w_n - 1];
		v_out = min(v_out, 0.0);
		res1.w_i_floor = (w_n - 1);
		res1.v_i_floor = v_move[w_n - 1];
	}

	res1.v_out = v_out;
	return res1;
}

// set value function in terminal case
void vfn::set_terminal(double phr_in) {
	double coh_perm, coh_perm2;
	double c_perm;

	double w_adj = 0.0;
	double c_comp_perm;
	double uc_perm;
	double V_perm;
	double V_fs = (1.0 - pow(beta, 20.0)) / (1.0 - beta) * ufn(c_fs, hu_ten_def, pref);

	int i_ph3;
	int i_t, i_s, i_w, i_yi, i_rent;

	for (i_t = 0; i_t < t_n; i_t++) {
		for (i_s = 0; i_s < n_s; i_s++) {
			for (i_w = 0; i_w < w_n; i_w++) {
				
				w_adj = c_fs + max(w_grid[i_w] , 0.0);   // terminal wealth

				V_perm = 1.0 / (1.0 - rho)*pow(w_adj, 1.0 - rho);
				
				vw3_grid[i_t][i_s][i_w] = b_motive*V_perm;        // bequest value

				if (i_t == 0) {
					vw3_def_grid[i_s][i_w] = V_perm;
				}

				v_move[i_w] = b_motive*V_perm;
			}

			interp_vw3(i_t, i_s);
		}
	}
}


// method checks and corrects for monotonicity and then estimates
// lambda interpolation parameters
void vfn::interp_vw3(int i_t_in, int i_s_in) {

	//vector<double> x0_default = { c_fs, 0.0, 0.0, 0.0, 0.0 };
	int w_i4;
	//double v0_default = 0.0;

	for (w_i4 = (w_n - 2); w_i4 >= 0; w_i4--) {
		vw3_grid[i_t_in][i_s_in][w_i4] = min(vw3_grid[i_t_in][i_s_in][w_i4], vw3_grid[i_t_in][i_s_in][w_i4 + 1]);
	}

}


// check for monotonicity
void vfn::clean_vw3_grid(int i_t_in, int i_s_in) {
	int w_i0;
	double vd1_min_tmp = 1e10, vd2_max_tmp = -1e10;
	double h_tmp = (w_grid[1] - w_grid[0]);

	double vd1_1;
	double vd1_0;
	double vd2;
	double vd2_0;

	for (w_i0 = (w_n - 2); w_i0 >= 0; w_i0--) {
		vw3_grid[i_t_in][i_s_in][w_i0] = min(vw3_grid[i_t_in][i_s_in][w_i0], vw3_grid[i_t_in][i_s_in][w_i0 + 1]);
	}

	vd1_1 = (vw3_grid[i_t_in][i_s_in][w_n-1] - vw3_grid[i_t_in][i_s_in][w_n-2]) / h_tmp;
	vd1_0 = (vw3_grid[i_t_in][i_s_in][w_n - 2] - vw3_grid[i_t_in][i_s_in][w_n - 3]) / h_tmp;
	
	vd2 = (vd1_1 - vd1_0) / h_tmp;

}


double vfn::get_h_step(int i_t_in, int i_s_in, int i_w_in) {
	int i_diff, i_max;
	double h_step = 0.2, d_max = 0.0;
	vector<double> diff(5, 0.0);

	if (i_w_in >= 2) {
		diff[0] = abs(x1_grid[i_t_in][i_s_in][i_w_in - 1] - x1_grid[i_t_in][i_s_in][i_w_in - 2]);
		diff[1] = abs(x2_grid[i_t_in][i_s_in][i_w_in - 1] - x2_grid[i_t_in][i_s_in][i_w_in - 2]);
		diff[2] = abs(x3_grid[i_t_in][i_s_in][i_w_in - 1] - x3_grid[i_t_in][i_s_in][i_w_in - 2]);
		diff[3] = abs(x4_grid[i_t_in][i_s_in][i_w_in - 1] - x4_grid[i_t_in][i_s_in][i_w_in - 2]);
		diff[4] = abs(x5_grid[i_t_in][i_s_in][i_w_in - 1] - x5_grid[i_t_in][i_s_in][i_w_in - 2]);

		for (i_diff = 0; i_diff < 5; i_diff++) {
			if (diff[i_diff] > d_max) {
				i_max = i_diff;
				d_max = diff[i_diff];
			}			
		}
		
		h_step = 10.0*d_max;
		//cout << "vfn.cpp: get_h_step: h_step = " << h_step << endl;
	}
	return h_step;
}
