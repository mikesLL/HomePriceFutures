/*
Given the current value fn VFN_3d_2, this program solve for the preceeding value fn VFN_3d_1

Start with the renter problem: t_i= 0 and loop through home price states ph_i and wealth states w_i
Given each state ph_i and w_i, cycle through next period tenure t_i2 = 0,1,2 
Solve optimization problem for that inidividual state
A solution includes policies x[5] = {C, B, X, CSFp, CSFn}, t_i2 in {0,1,2}
and a value of the current state V1()

C: Consumption 
B: Bond (Mortgage) holdings
X: Equities
CSFp: Case-Shiller Index Futures, Long Position, Margin
CSFn: Case-Shiller Index Futures, Short Position, Margin
t_i, t_i2 in {0: rent, 1: starter home, 2: full-size home}

Copyright A. Michael Sharifi, 2016

TODO: add an algorithm to use the concavity of the
value function to generate quick guesses at the begining of the state
also exploits fact that policies / strategies should be smooth
Would be kind of like the MCMC analog to VFI
*/

#include "headers.h"

void gen_VP(void *snodes_in, void *VFN_3d_1, void *VFN_3d_2 ){
	
    double duration;
	clock_t start;

	snodes *snodes1 = (snodes *)snodes_in;
	vfn *rr2 = (vfn *) VFN_3d_2;          // address to initialized V2
	vfn *rr1 = (vfn *) VFN_3d_1;          // address to initialized V1

    gen_res res1;
	eval_res res0;
	eval_res res_t_lag;
	eval_res res_t_0;
	
	double b_min;
	
	vector<double> theta1;
	
	int y_i; 
	int t_i, t_i2 , ph_i, w_i, t_adj;
	int i_rent, i_ph, i_yi;

	int i_s;
	int t_hor = (*snodes1).t_hor;

	if (t_hor == 0) {
		cout << "first horizon" << endl; 
	}

	double coh, w_adj, v_adj, v_i_floor, v0_opt, v0, beg_equity, mpmt, b_min2;
	double v1;
	double v_lag_w; // value function guess from previous wealth level
	double v_lag_t; // value function guess from previous tenure computation

	vector<double> x(5,0.0);                     // policy / strategy
	vector<double> x_lag_w(5, 0.0);
	vector<double> x0(5, 0.0);
	vector<double> x1(5, 0.0);

	vector<double> x_guess(5, 0.0);

	vector<double> x_ti1(5, 0.0);       // policy: owns smallest home
	double v_ti1 = -1.0e20;              // value: owns smallest home
	//x1 = x_ti1 * 1.2;
	//v1 = v_ti1 * (1.2);

	double h_mult = 0.0;
	double bond_est = 0.0;
	
	vector<vector<double>> x_lag_wt(t_n, vector<double> (5, 0.0) );

	for (i_ph = 0; i_ph < n_ph; i_ph++) {
		cout << (*snodes1).p_gridt[t_hor][i_ph] << "..." << endl;
	}

	t_i = 0; // consider t_i = 0 first; case: begin with renter 
	int res1_flag = 1;
	int t_i2_lag_w = 0;

	for (i_s = 0; i_s < n_s; i_s++) {
		i_yi = (*snodes1).s2i_yi[i_s];
		i_rent = (*snodes1).s2i_rent[i_s];
		i_ph = (*snodes1).s2i_ph[i_s];

		start = clock();
		cout << "i_s = " << i_s << "  t_i = " << t_i << "  i_yi = " << i_yi 
			 << "  i_ph = " << i_ph << "  begin renter problem" << endl;

		for (w_i = 0; w_i < w_n; w_i++) {
			for (t_i2 = 0; t_i2 < t_n; t_i2 ++){
				
				if (w_i % 100 == 0) {
					cout << "w_i = " << w_i << endl;
				}
	
				if (t_i2 == 0) {
					b_min = b_min_unsec;
					beg_equity = -1.0e6;
					mpmt = -1.0e6;
				}
				else {
					b_min =  - max_ltv*(*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph] + b_min_unsec;          // lower bound on bond / mortgage
					beg_equity = min_dpmt * (*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph];
				    //mpmt = ( rb + mort_spread - 1.0 ) * (*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph];
					mpmt = -1.0e6;
					b_min2 = -max_lti * (*snodes1).yi_gridt_btax[t_hor][i_yi] / (rb + mort_spread - 1.0);

					b_min = max(b_min, b_min2);
				}

				// load previous w_i policy as a benchmark
				(*rr1).get_pol(t_i, i_s, w_i - 1, x_lag_w);   					 
				t_i2_lag_w = (*rr1).xt_grid[t_i][i_s][max(w_i - 1, 0)];
				v_lag_w = (*rr1).vw3_grid[t_i][i_s][max(w_i - 1, 0)];

				// so now have:
				// x_lag_w, t_i2_lag_w, v_lag_w;
				
				// load value given previous t_i as a benchmark
				v_lag_t = (*rr1).vw3_grid[t_i][i_s][w_i];

				// load t_i2-restricted guess as an initial starting point
				x_guess = x_lag_wt[t_i2];
						
				coh = (*rr1).w_grid[w_i] + y_atax*(*snodes1).yi_gridt[t_hor][i_yi] -
					(1.0 + phi_buy) * (*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph];

				if (t_i2 == 0) {
					coh = coh - min( (*snodes1).rent_gridt[t_hor][i_rent] * (*snodes1).rent_adj,
						1.0 / 3.0*(*snodes1).yi_gridt[t_hor][i_yi] );
				}
				
				(*rr2).i_s1 = i_s;  // pass in current state to next-period value function
				(*rr2).t_i1 = t_i;
				(*rr2).w_i1 = w_i;                      
				(*rr2).t_i2 = t_i2; //  t_i2 is a choice variable, so adding it to fn pointer 

				if ( t_i2 == 0 ) {
					(*rr2).def_flag = 1;
					res1 = gen_VPw(snodes1, rr1, rr2, coh, x_guess, b_min, beg_equity, mpmt);
					(*rr1).vw3_def_grid[i_s][w_i] = res1.v_opt; 
					//cout << res1.v_opt << endl; 
					(*rr2).def_flag = 0;
				}
				
				res1 = gen_VPw(snodes1, rr1, rr2, coh, x_guess, b_min, beg_equity, mpmt);                               // pass problem into gen_V1_w
				
				if (t_i2 == 1) {
					v_ti1 = res1.v_opt;
					x_ti1 = res1.x_opt;
				}

				v1 = res1.v_opt;   // guess for the current vfn given current t_i2
				x1 = res1.x_opt;    // guess for current policy given current t_i2

				if (t_i2 == 0) {
					(*rr1).v_move[w_i] = v1;
				}
				
				// store the just retrieved x_opt as a guess for the next period under same tenure
				x_lag_wt[t_i2] = x1;

				if ((v1 > v_lag_t) || (t_i2 == 0)) {                      // current result is better than value given previous tenure
					(*rr1).set_pol_ten_v(t_i, i_s, w_i, x1, t_i2, v1);  // set x, t_i2, v0 in
					
					if (res1.valid_flag == 0) {
						(*rr1).set_pol_ten_v(t_i, i_s, w_i, x1, 0, v1);
					}
				}
			}
		}
		(*rr1).interp_vw3(t_i, i_s);  // clean and interpolate grid			

		duration = (clock() - start) / (double)CLOCKS_PER_SEC;
		cout << "time elapsed: " << duration << '\n';
	}
	

	// now that the t_i = 0 case has been solved, t_i = 1,2 cases when sale or trade-up are equivalent
	//to t_i = 0 with the correct downward wealth adjustment; 
	//do not need to cycle through t_i2 = {0, 1,2}, only t_i2 = t_i
	
	cout << "gen_VP.cpp: begin homeowner problem" << endl; 
	for (t_i = 1; t_i < t_n; t_i++) {                        // consider t_i = 0 first; case: begin with renter 
		for (i_s = 0; i_s < n_s; i_s++) {
			i_yi = (*snodes1).s2i_yi[i_s];
			i_rent = (*snodes1).s2i_rent[i_s];
			i_ph = (*snodes1).s2i_ph[i_s];

			start = clock();
			cout << "i_s = " << i_s << "   t_i = " << t_i << "i_yi = " << i_yi 
				 << " i_ph = " << i_ph << " begin homeowner problem " << endl;

			for (w_i = 0; w_i < w_n; w_i++) {
				t_i2 = t_i;
				if (t_i2 == 0) {
					b_min = b_min_unsec;
				}
				else {
					b_min = -max_ltv*(*snodes1).ten_w[t_i2] * (*snodes1).p_gridt[t_hor][i_ph] + b_min_unsec;
					b_min2 = - max_lti * (*snodes1).yi_gridt_btax[t_hor][i_yi] / (rb + mort_spread - 1.0);
					//b_min = max(b_min, b_min2);
				}

				// load previous w_i policy as an initial guess
				(*rr1).get_pol(t_i, i_s, w_i - 1, x_lag_w);                       // get x pol sol from previous w_i and assign to x
				t_i2_lag_w = (*rr1).xt_grid[t_i][i_s][max(w_i - 1, 0)];
				v_lag_w = (*rr1).vw3_grid[t_i][i_s][max(w_i - 1, 0)];
      
				coh = (*rr1).w_grid[w_i] + y_atax*(*snodes1).yi_gridt[t_hor][i_yi] - (*snodes1).ten_w[t_i] * maint_mult * (*snodes1).p_gridt[t_hor][i_ph];       // subtract housing wealth from cash on hand  			
				
				mpmt = -1.0e6; 
				beg_equity = -1.0e6;

				(*rr2).i_s1 = i_s;
				(*rr2).t_i1 = t_i;
				(*rr2).w_i1 = w_i;
				(*rr2).t_i2 = t_i2;                                                  // impose t_i2 = t_i    	      

				res1 = gen_VPw(snodes1, rr1, rr2, coh, x_lag_w, b_min, beg_equity, mpmt);                  // solve optimization problem
				v1 = res1.v_opt;  // guess for current solution: vfn
				x1 = res1.x_opt;  // guess for current policy

				(*rr1).set_pol_ten_v(t_i, i_s, w_i, x1, t_i2, v1);

				w_adj = (*rr1).w_grid[w_i] - phi_sell*(*snodes1).ten_w[t_i] * (*snodes1).p_gridt[t_hor][i_ph];    //calc costs if agent sells and converts to renter
				res_t_0 = (*rr1).eval_v(0, i_s, w_adj); // evaluate value fn if agent sells and converts to renter

				if ( (w_adj >= 0.0 ) && (res_t_0.v_i_floor > v1) && (res_t_0.w_i_floor >= 0)) {

					(*rr1).get_pol(0, i_s, res_t_0.w_i_floor, x);                         // submit x as reference and load in x pol from t1 = 0
					t_adj = (*rr1).xt_grid[0][i_s][res_t_0.w_i_floor];                    // get t2 pol from t1 = 0; simulated sale
					(*rr1).set_pol_ten_v(t_i, i_s, w_i, x, t_adj, res_t_0.v_i_floor);     // first arguments are current state variables, x containts updated policy
				}

				//cout << "gen_Vp: (*rr1).vw3_def_grid[i_s][w_i_zero] " << (*rr1).vw3_def_grid[i_s][w_i_zero] << endl;
				//double vw3_def_test = (*rr2).vw3_def_grid[i_s][w_i_zero]; 

				/*
				if (  (*rr1).vw3_def_grid[i_s][w_i_zero] > max(v1, res_t_0.v_i_floor ) ) {

					(*rr2).def_flag = 1;
					coh = 0.0 + y_atax*(*snodes1).yi_gridt[t_hor][i_yi] -
						min((*snodes1).rent_gridt[t_hor][i_rent] * (*snodes1).rent_adj,
							1.0 / 3.0 * (*snodes1).yi_gridt[t_hor][i_yi]);

					b_min = b_min_unsec; // TODO: verify mod here
					beg_equity = -1.0e6;

					res1 = gen_VPw(snodes1, rr1, rr2, coh, x_guess, b_min, beg_equity, mpmt);

					(*rr1).set_pol_ten_v(0, i_s, w_i, res1.x_opt, 0, (*rr1).vw3_def_grid[i_s][w_i_zero] );
					(*rr2).def_flag = 0;
				}
				*/
				
				
			}
			(*rr1).interp_vw3(t_i, i_s);  // clean grid
			duration = (clock() - start) / (double)CLOCKS_PER_SEC;
			cout << "time elapsed: " << duration << "  lcount = " << (*rr2).lcount << endl;
		}
	}

	// clean grid again
	for (t_i = 0; t_i < t_n; t_i++) {                       
		for (i_s = 0; i_s < n_s; i_s++) {
			(*rr1).interp_vw3(t_i, i_s);  // clean grid
		}
	}

	cout << " gen" << to_string((*rr1).t_num) << "completed" << endl;
}
