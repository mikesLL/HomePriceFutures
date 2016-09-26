// vfn.h
// vfn is main data structure in project
// Copyright A. Michael Sharifi, 2016

using namespace std;

// eval_res: value associated with fn at that point
typedef struct{
	double v_out;              // value function
	int w_i_floor;       
	double v_i_floor;
} eval_res;


#ifndef VFN_H
#define VFN_H

class snodes;

class vfn {

public:

	int t_i2, t_i1, w_i1, i_s1;
	int def_flag; 
	int lcount;
	double ph1, phr, w_comp;
	double w_grid[w_n];
	double v_opt_loop;
	int t_num;
	int t_id;
	int pref;
	double csf_1yr;
	int T_max;

	snodes *snodes1;    // pointer to 

	vector<double> v_move;

	vector<vector<vector<int>>> xt_grid;
	vector<vector<vector<double>>> x1_grid, x2_grid, x3_grid, x4_grid, x5_grid; // control variables
	vector<vector<vector<double>>> vw3_grid;                                    // value function
	vector<vector<double>> vw3_def_grid;                                // value function (default)

	vector<double> vw3_grid_norm;                               // value function (normalized)
	vector<vector<vector<double>>> lambda_grid;                                 // lambda interpolation parameters

	vector<vector<vector<double>>> vw3_d_grid;                                 // first derivative
	vector<vector<vector<double>>> vw3_dd_grid;                                 // second derivative

	void enter_data( void *snodes_in, double phr_in,int t_id, int t_num_in, double csf_1yr_in, int pref_in, int T_max_in );	
	
	void get_pol(int i_t_in, int i_s_in, int i_w_in, vector<double> &x_pol);

	void set_pol_ten_v(int i_t_in, int i_s_in, int i_w_in, vector<double> &x_pol, int t_i2_in, double v0_in);
	

	eval_res eval_v(int i_t_in, int i_s_in, double w_in);
	eval_res eval_v_def( int i_s_in, double w_in);


	eval_res eval_v_norm(double w_in);

	eval_res eval_v_move(int i_t_in, int i_s_in, double w_in, int t_left);


	void set_terminal( double phr_in );

	void interp_vw3(int i_t_in, int i_s_in);

	void clean_vw3_grid(int i_t_in, int i_s_in);

	double get_h_step(int i_t_in, int i_s_in, int w_i_in );
};

#endif
