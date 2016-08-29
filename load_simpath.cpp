/*
load_simpath.cpp
given intial parameters, this program simulates a home price path
and returns interpolating prices as well as a transition matrix

Copyright A. Michael Sharifi, 2016

*ps0: pointer to price structure
*gs0: pointer to transition matrix stucture
rent_in: initial rent
ph0_in: initial home price
ret0_in: most recently observed hpa
csf_1_yr_in: currently observed futures price
t_id: initial year
city_init_in: initial city

impose forecasted 1st year return equals that implied by csf
ph0 : current price
csf1_yr: forecasted 1yr price
ret_tn = (csf_1yr - ph0) / ph0 + sigma_ret*randn_draw;

In actual regression:
1. Estimate cointegrating relation ship with rent to price ratio
2. regress nominal returns on lagged nominal returns
3. Done; parameters can be entered here

Numeraire: non-durable consumption
Rent grows in real terms
Income grows in real terms as the agent ages and city-wide income increases
*/

#include "headers.h"

void load_simpath(void *snodes_in, double rent_in, double ph0_in, double ret0_in, double csf_1yr_in, int t_id, string city_init_in, int city_id, int age_begin_in) {

	int s1, s2, s_test;                             // state in current period, state in next period
	int N_print =  10000;                          // number of observations to print to file
	int N_sim = 2000000;                                            // number of simulations
	int i_ph, i_rent, i_yi, i_s;                                 // state and individual dimension indices

	cout << "load ppath: ph0_in:  " << ph0_in << endl;

	// Parameters for each city
	// year: 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014
	double inf_mult[] = {104.5, 106.9, 109.7, 113.4, 117.1, 120.4, 125.0, 124.6, 126.6, 130.6, 133.3, 135.3, 137.6}; 

	//determnistic (real) rent growth
	double g_rent_store[] = {
		0.00279512280912727,
		0.00372666252799066,
		0.00374066193750439,
		0.00158055419701987,
		0.00503271347985580,
		- 0.000102926026248194,
		0.00425295851646686,
		0.00612790147006233
	};

	double g_rent = g_rent_store[city_id];

	ph0_in = ph0_in / inf_mult[t_id] * 100.0;
	rent_in = rent_in / inf_mult[t_id] * 100.0;
	csf_1yr_in = csf_1yr_in / inf_mult[t_id] * 100.0; 

	double alpha_store[] = {
	0.0104607315801136,
	0.0191074838324645,
	0.00704236034135076,
	0.00381705100823080,
	-0.000844258143867715,
	0.0116455917169509,
	0.00635205965167907,
	0.000754767829009388
	};

double rhof_store[] = {
	0.498136281955623,
	0.214709512333240,
	0.586569953107669,
	0.617686906437158,
	0.525396917591809,
	0.367979365128841,
	0.660449824037308,
	0.685848074672189
};

double theta_store[] = {
		0.320604675683815,
		0.517775013381016,
		0.367766407618239,
		0.298417763209689,
		0.243452622039727,
		0.518388584994886,
		0.275703121600550,
		0.217163204014952
};

	
double gamma0_store[] = {
	-2.24523949617400,
	-2.06691403017790,
	-2.17662237312747,
	-2.33405568060629,
	-2.58660425485324,
	-2.55177401345514,
	-2.39513031942948,
	-2.36991753806891
};

	double gamma1_store[] = {
		0.568455966086471,
		0.544596925825744,
		0.496007575833706,
		0.586733373192217,
		0.731138993521001,
		0.611450534934516,
		0.548719284867125,
		0.600959631692874
	};

	double sigma_ret_store[] = {
		0.0842179541145515,
		0.0921669031549657,
		0.0850769610135847,
		0.0374637817567246,
		0.0392694194702332,
		0.0302220840642717,
		0.0717777366157635,
		0.0340633256106433
	};

	double y_city_mult[] = { 1.23, 1.54, 1.02, 1.34, 1.02, 1.00,  0.95, 1.39 };

	double alpha_hat = alpha_store[city_id];
	double rhof_hat = rhof_store[city_id];
	double theta_hat = theta_store[city_id];
	double sigma_ret = sigma_ret_store[city_id];
	double gamma0_hat = gamma0_store[city_id];
	double gamma1_hat = gamma1_store[city_id];
		 
	double p_min = 0.01;                                                         // lower bound on home prices										   
	
	// stdev centers (tauchen-style discretization) 
	double ph_nd_std[] = { -2.0, -1.5, -1.0, -0.5, 0.0, .5, 1.0, 1.5, 2.0 };     
	double rent_nd_std[] = { 0.0 };                                           // can also set { -1.0, 0.0, 1.0 };
	double yi_nd_std[] = {  -1.5, 0.0, 1.5 }; 
	
	string city_init = city_init_in;                                          // initial conditions for simulation
	    
	double rent0 = rent_in;

	cout << "load ppath; rent = " << rent0 << endl;
	
	double ph0 = ph0_in;
	double ret0 = ret0_in;
	double csf_1yr = csf_1yr_in;
	
	snodes *snodes1 = (snodes*)snodes_in;
	int T_max = (*snodes1).T_max;

	cout << "begin load_simpath" << endl ;
    double duration;
	clock_t start;                          // instantiate timer
	start = clock();
	

	random_device rd;                       // instantiate random generator
	mt19937 gen(rd());
	normal_distribution<> dist(0.0, 1.0);   // standard normal distribution

	int i = 0, j = 0;
	cout << " N_sim = " << N_sim << endl;
	
	int T_sim = T_max + 1;                                                              // T_max is a measure of the investor's horizon
	
	cout << "T_sim = " << T_sim << endl;
	
	double ret_tn, ret_lag, ecm;

	vector<vector<double>> ph_str(T_sim + 1, vector<double>(N_sim, 0.0) );              // store city-wide home prices
	vector<vector<double>> rent_str(T_sim + 1, vector<double>(N_sim, 0.0));             // store city-wide rent
	vector<vector<double>> yi_str(T_sim + 1, vector<double>(N_sim, 0.0));               // individual income

	vector<vector<double>> ph_str_nds(T_sim + 1, vector<double> (n_ph, 0.0) );          // home price nds
	vector<vector<double>> rent_str_nds(T_sim + 1, vector<double> (n_rent, 0.0) );      // rent nds
	vector<vector<double>> yi_str_nds(T_sim + 1, vector<double>(n_yi, 0.0) );         // individual-income nds

	double v_t0 = 0.0;                                   // initial permanent income shock
	double v_t = 0.0;                                    // permanent income state
	
	
	double log_fe = 2.3831 + 0.4831 - 0.0228*2.0;        //2.18 // fixed effects
	double y_mu_city = 0.005;                            // city-wide real income/rent growth
	double y_sig_city = 0.005;

	double sigma_u = pow(0.0106, 0.5);  //0.005;                              // permanent income shock stdev
	double sigma_e = pow(0.0738, 0.5); // 0.005;                              // transitory income shock stdev

	double mu_yc =  0.000;    // 0.0                             // avg real income / rent growth
	double sigma_yc = 0.0001; //0.01;   // 0.0                  // avg real income / rent growth stdev

	double mult_2005 =   .56 / .51; // 1.5;                              // multiplier: convert $1992 to $2005
	double mult_unit = 0.01;                             // multiplier: convert $1000s from parameters to $100000s in paper

	double log_y_rand = 0.0;                             // random componenet of income (log)
	double log_y_age;                                    // age-based (deterministic) component of income
    double g_yc_t = 0.0;                                 // city-wide income growth
    double u_t = 0.0;                                    // individual permanent income shock
    double e_t = 0.0;                                    // idiosyncratic income shock

	int t = 0, n = 0, n1 = 0, n2 = 0;                                          // time and simulation indices
	double res_sum, res_sq_sum, res_mean, res_std, row_sum;                    // sums for computing means and standard deviations

    //int i_age; // Track agent's age
	//int age_t = age_begin_in;  // track household age
	double age_td; // track household age as double

	// loop through simulations
	cout << "load_simpath.cpp: Begin Simulations" << endl;
	for (n = 0; n < N_sim; n++){                  
		if ( n % 10000 == 0) {
			cout << "n = " << n << endl;
		}
		
		// t = 0: initial year
		t = 0;                                 
		age_td = double(age_begin_in) + double(t); // initial age
		
		ph_str[t][n] = log(ph0);
		rent_str[t][n] = rent0;
	
		// compute age-related component of log-income (deterministic)
		
		log_y_age = -4.3148 + 0.3194*age_td - 0.0577*pow(age_td, 2.0) / 10.0 
			+ 0.0033*pow(age_td, 3.0) / 100.0;
		
		yi_str[t][n] = y_city_mult[city_id]*mult_unit*mult_2005*exp(log_fe + log_y_age);
		yi_str[t][n] = log(yi_str[t][n]);

		v_t0 = 0.0;                             // set permanent income shock

		// t = 1: first year
		t = 1;                                  // first year: draw shocks
		age_td = double(age_begin_in) + double(t);   // current age

		g_yc_t = mu_yc + sigma_yc*dist(gen);    // city-wide income / rent
		u_t = sigma_u*dist(gen);                // individual: permanent shock
		e_t = sigma_e*dist(gen);                // individual: transient shock

		v_t = v_t0 + u_t;  // upadate income permanent component
		v_t0 = v_t;                             // update transient component
		
		// compute age-related component of log-income (deterministic)
		log_y_age = -4.3148 + 0.3194*age_td - 0.0577*pow(age_td, 2.0) / 10.0 + 0.0033*pow(age_td, 3.0) / 100.0;

		// compute home price appreciation
		ret_tn = (csf_1yr - ph0) / ph0 + sigma_ret*dist(gen);  // impose first year expected return equals the futures-based forecast

		rent_str[t][n] = exp(g_rent)*rent_str[t-1][n]; // update rent path
		
		yi_str[t][n] = y_city_mult[city_id]*mult_unit*mult_2005*exp( log_fe + log_y_age + v_t + e_t );  // Cocco, Maenhoust, Gomes parameters are in 1000's; convert to 100k
		yi_str[t][n] = log(yi_str[t][n]);

		ph_str[t][n] = (ret_tn) + ph_str[t - 1][n];  // home prices in sim are in logs

		// simulate later time periods
		for (t = 2; t < (T_sim + 1); t++){
			age_td = double(age_begin_in) + double(t);  // compute age

			// compute age-related component of log-income (deterministic)
			log_y_age = -4.3148 + 0.3194*age_td - 0.0577*pow(age_td, 2.0) / 10.0 + 0.0033*pow(age_td, 3.0) / 100.0;

			// draw shocks
			g_yc_t = mu_yc + sigma_yc*dist(gen);  // city-wide income / rent
			u_t = sigma_e*dist(gen);              // individual: permanent shock
			e_t = sigma_e*dist(gen);              // individual: transient shock

			// update income path 
			v_t = v_t0 + u_t;
			v_t0 = v_t;                             // update transient component
			
			// compute home price
			ret_lag = ph_str[t - 1][n] - ph_str[t - 2][n];                                       // ph_str is in logs 
			ecm = log(rent_str[t - 1][n]) - gamma0_hat - gamma1_hat*(ph_str[t - 1][n]);          // cointegrate rents, prices
			
			ret_tn = alpha_hat + rhof_hat*ret_lag + theta_hat*ecm + sigma_ret*dist(gen);         // return series
		
			rent_str[t][n] = exp(g_rent)*rent_str[t-1][n];

			yi_str[t][n] = mult_unit*mult_2005*exp(log_fe + log_y_age + v_t + e_t);
			yi_str[t][n] = log(yi_str[t][n]);

			ph_str[t][n] = ret_tn + ph_str[t - 1][n];
		}
	}

	cout << "load_simpath.cpp: Finished Running Simulations" << endl;

	cout << "load_simpath.cpp: Compute Home Price Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(ph_str[t].begin(), ph_str[t].end(), 0.0);
		res_mean = res_sum / (double)ph_str[t].size();
		res_sq_sum = inner_product(ph_str[t].begin(), ph_str[t].end(), ph_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double) ph_str[t].size() - res_mean * res_mean);

		if (t == 0) {
			res_std = 0.01;  // initial period: manually set stdev 
		}

		for (n = 0; n < n_ph; n++) {
			ph_str_nds[t][n] = res_mean + ph_nd_std[n] * res_std;  // compute home price nodes (log scale)
			(*snodes1).p_gridt[t][n] = exp(ph_str_nds[t][n]);    // load home price nodes into p_gridt
			cout << (*snodes1).p_gridt[t][n] << "...";      // print out price grid entry
		}
		cout << endl;
	}

	cout << "load_simpath.cpp: Compute Rent Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(rent_str[t].begin(), rent_str[t].end(), 0.0);
		res_mean = res_sum / (double) rent_str[t].size();

		res_sq_sum = inner_product(rent_str[t].begin(), rent_str[t].end(), rent_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double) rent_str[t].size() - res_mean * res_mean);

		// check null stdev case
		if ( (res_std != res_std ) || (res_std <= 0.01) ) {
			res_std = 0.01;  // initial period: manually set stdev
		}

		for (n = 0; n < n_rent; n++) {
			rent_str_nds[t][n] = res_mean + rent_nd_std[n] * res_std;
			(*snodes1).rent_gridt[t][n] = rent_str_nds[t][n];    // load rent nodes into rent_gridt
		}
	}

	cout << "load_simpath.cpp: Compute Income Nodes " << endl;
	for (t = 0; t < T_sim; t++) {
		res_sum = accumulate(yi_str[t].begin(), yi_str[t].end(), 0.0);
		res_mean = res_sum / (double)yi_str[t].size();

		res_sq_sum = inner_product(yi_str[t].begin(), yi_str[t].end(), yi_str[t].begin(), 0.0);
		res_std = sqrt(res_sq_sum / (double)yi_str[t].size() - res_mean * res_mean);

		if (t == 0) {
			res_std = 0.01;
		}
		for (n = 0; n < n_yi; n++) {
			yi_str_nds[t][n] = res_mean + yi_nd_std[n] * res_std; // compute income nodes (log scale)
			(*snodes1).yi_gridt[t][n] = exp(yi_str_nds[t][n]);    // load individual-income nodes into income_gridt
		}
	}

	cout << "load_simpath.cpp: Begin Computing State-Space Transition Matrix" << endl;

	vector<vector<int>> zeros_NS_NS(n_s, vector<int>(n_s, 0)); // default matrix for state to state transition
	vector<vector<int>> gamma_store = zeros_NS_NS;             // gamma matrix tracks state to state transition

	t = 0;
	double ph_step;
	double rent_step;
	double yi_step;

	double ph_step1;
	double rent_step1;
	double yi_step1;

	cout << "load_simpath.cpp: Cycle Through Observations" << endl;
	for (t = 0; t < T_sim; t++) {
		gamma_store = zeros_NS_NS;                                  // initialize current-period transition matrix

		ph_step = ph_str_nds[t][1] - ph_str_nds[t][0];          // step-sizes
		rent_step = 0.001; //rent_str_nds[t][1] - rent_str_nds[t][0];
		yi_step = yi_str_nds[t][1] - yi_str_nds[t][0];

		ph_step1 = ph_str_nds[t + 1][1] - ph_str_nds[t + 1][0];
		rent_step1 = 0.001; //rent_str_nds[t+1][1] - rent_str_nds[t+1][0];
		yi_step1 = yi_str_nds[t + 1][1] - yi_str_nds[t + 1][0];


		for (n = 0; n < N_sim; n++) {

			// current period: map observation to closest node
			i_ph = (int) round( (ph_str[t][n] - ph_str_nds[t][0]) / ph_step);
			i_rent = (int) round( (rent_str[t][n] - rent_str_nds[t][0]) / rent_step);
			i_yi = (int) round( (yi_str[t][n] - yi_str_nds[t][0]) / yi_step);
			
			i_ph = min(max(i_ph, 0), n_ph - 1);              // bound extremes
			i_rent = min(max(i_rent, 0), n_rent - 1);
			i_yi = min(max(i_yi, 0), n_yi - 1);

			s1 = (*snodes1).i2s_map[i_ph][i_rent][i_yi];                // map node to state

			// next period: map observation to closest node
			i_ph = (int)round((ph_str[t + 1][n] - ph_str_nds[t + 1][0]) / ph_step1);
			i_rent = (int)round((rent_str[t + 1][n] - rent_str_nds[t + 1][0]) / rent_step1);
			i_yi = (int)round((yi_str[t + 1][n] - yi_str_nds[t + 1][0]) / yi_step1);

			i_ph = min(max(i_ph, 0), n_ph - 1);              // bound extremes
			i_rent = min(max(i_rent, 0), n_rent - 1);
			i_yi = min(max(i_yi, 0), n_yi - 1);
			
			s2 = (*snodes1).i2s_map[i_ph][i_rent][i_yi];                 // map node to state
			
			gamma_store[s1][s2] += 1;                                // store  in transition matrix
		}

		cout << "load_simpath.cpp: t = " << t  << "  load transition matrix into data structure "<< endl;

		// load transition matrix into data structure
		for (i = 0; i < n_s; i++) {
			row_sum = (double) accumulate(gamma_store[i].begin(), gamma_store[i].end(), 0);

			for (j = 0; j < n_s; j++) {
				if (row_sum >= 1.0) {
					(*snodes1).gammat[t][i][j] = (double)gamma_store[i][j] / row_sum;   // load transition matrix into structure
				} else {
					(*snodes1).gammat[t][i][j] = 0.0;
				}
			}
		}
	}

	cout << "load_simpath.cpp: Finished Computing Transition Matrix" << endl;

	cout << "load_simpath.cpp: gamma T max:" << endl;
	
	string fn_beg = "price_results/" + city_init + "_age" + to_string((*snodes1).age0) + "_yr" + to_string(t_id);

	// print sample price paths
	ofstream price_obs_file;                                        
	//price_obs_file.open("price_results/" + string("age") + to_string( (*snodes1).age0 ) + "/"  + city_init + "yr" + to_string(t_id) + "price_obs_file.csv", ios::out | ios::trunc);
	price_obs_file.open(fn_beg + "_price_obs_file.csv", ios::out | ios::trunc);

	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			price_obs_file << exp(ph_str[t][n]) << ",";
		}
		price_obs_file << endl;
	}
	price_obs_file.close();

	// print sample rent paths
	ofstream rent_obs_file;                                       
	//rent_obs_file.open("price_results/" + city_init + "age" + to_string((*snodes1).age0)  + "yr" + to_string(t_id) + "rent_obs_file.csv", ios::out | ios::trunc);
	rent_obs_file.open( fn_beg + "_rent_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			rent_obs_file << rent_str[t][n] << ",";
		}
		rent_obs_file << endl;
	}
	rent_obs_file.close();

	// print sample income paths
	ofstream yi_obs_file;                                     
	//yi_obs_file.open("price_results/" +  city_init + "_age" + to_string((*snodes1).age0) + "yr" + to_string(t_id) + "yi_obs_file.csv", ios::out | ios::trunc);
	yi_obs_file.open( fn_beg +  "_yi_obs_file.csv", ios::out | ios::trunc);
	for (n = 0; n < N_print; n++) {
		for (t = 0; t < T_sim; t++) {
			yi_obs_file << exp(yi_str[t][n]) << ",";
		}
		yi_obs_file << endl;
	}
	yi_obs_file.close();


	// print out price grid
	ofstream pstruct_file;                                                  
	//pstruct_file.open("price_results/" + city_init + "_age" + to_string((*snodes1).age0)  + "_yr" + to_string(t_id) + "pstruct_file.csv", ios::out | ios::trunc);
	pstruct_file.open( fn_beg + "_pstruct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_ph; n++) {
			pstruct_file << (*snodes1).p_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		pstruct_file << endl;
	}

	// print out income grid
	ofstream ystruct_file;                                                              
	//ystruct_file.open("price_results/" + city_init + "_age" + to_string((*snodes1).age0) + "_yr" + to_string(t_id) + "ystruct_file.csv", ios::out | ios::trunc);
	ystruct_file.open( fn_beg + "_ystruct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_yi; n++) {
			ystruct_file << (*snodes1).yi_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		ystruct_file << endl;
	}

	// print out rent grid
	ofstream rstruct_file;
	//rstruct_file.open("price_results/" +  city_init + "_age" + to_string((*snodes1).age0) +  "_yr" + to_string(t_id) + "rstruct_file.csv", ios::out | ios::trunc);
	rstruct_file.open(fn_beg + "_rstruct_file.csv", ios::out | ios::trunc);
	for (t = 0; t < T_sim; t++) {
		for (n = 0; n < n_rent; n++) {
			rstruct_file << (*snodes1).rent_gridt[t][n] << ",";         // alt: ph_str_nds[t][n];
		}
		rstruct_file << endl;
	}

	// print out transition matrix
	ofstream gstruct_file;                                                           
	//gstruct_file.open("price_results/" + city_init + "_age" + to_string((*snodes1).age0) + "_yr" + to_string(t_id) + "gstruct_file.csv", ios::out | ios::trunc);
	gstruct_file.open(fn_beg + "_gstruct_file.csv", ios::out | ios::trunc);

	for (t = 0; t < T_sim; t++) {
		for (i = 0; i < n_s; i++) {
			for (j = 0; j < n_s; j++) {
				gstruct_file << (*snodes1).gammat[t][i][j] << ",";
			}
			gstruct_file << endl;
		}
		gstruct_file << endl;
	}

	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout<<"Finished load_rand. Time elapsed: "<< duration <<'\n';
}

