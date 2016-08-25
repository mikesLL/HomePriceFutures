// store_data.cpp
// Writes value function structure to .csv
// Inputs: value function structure (*vfn_in), file name (file_name_in),
// starting year (year1_id_in), and a year horizon (year1t_id_in)
// which contains all results
// Copyright A. Michael Sharifi

#include "headers.h"

void store_data(void *snodes_in, void *vfn_in, string file_name_in, int year1_id_in, int year1t_id_in) {

	int i_s, i_ph, i_w, t_i, w_i, i_yi, i_rent;
	int y_i_def = 0;
	int year1_id = year1_id_in;                                   // read-in / curent year
	int year1t_id = year1t_id_in;                                 // year horizon

	vfn *vfnt = (vfn *)vfn_in;                              // address to initialized V2

	snodes *snodes1 = (snodes *)snodes_in;

	// this is the state of interest; it represents median income, median home price, median rent out of all states
	int i_s_mid = (*snodes1).i_s_mid;

	ofstream v1_file;                                             // open output file stream    

	string file_name_def =
		file_name_in + to_string(year1_id) + "yr1t" + to_string(year1t_id) + "pref" + to_string(pref) + "rho" + to_string(rhoi) + "gamma" + to_string(gammai) + "csfLev" + to_string(csfLevi) + "w_n" + to_string(w_n);

	string file_name_def2 =
		file_name_in + to_string(year1_id) + "pref" + to_string(pref) + "rho" + to_string(rhoi) + "gamma" + to_string(gammai) + "csfLev" + to_string(csfLevi) + "w_n" + to_string(w_n);

	//string file_name = file_name_def + ".csv";
	string file_name = "vfn_results/" + file_name_def + ".csv";


	cout << "store_data.cpp: list prices" << endl;
	for (i_ph = 0; i_ph < n_ph; i_ph++) {
		cout << "ph = " << (*snodes1).p_gridt[year1t_id][i_ph] << "..." << endl;
	}
	
	cout << "store_data.cpp: write v1_file" << endl;
	v1_file.open(file_name, ios::out | ios::trunc);                                 // outstream, truncate

	// print t_i, ph_i, ph and state headers
	// commas represent spaces for control variables
	for ( t_i = 0; t_i < t_n; t_i++) {                                               
		for (i_ph = 0; i_ph < n_ph; i_ph++) {
			for (i_rent = 0; i_rent < n_rent; i_rent++) {
				for (i_yi = 0; i_yi < n_yi; i_yi++) {
					// for i_yi = 0; i_yi < n_y ; cycle across labor income states
					//v1_file << "t_i = " << t_i << "," << "i_ph = " << i_ph << "," << "ph = " << (*snodes1).p_gridt[year1t_id][i_ph] << " , , , , , , ";
					v1_file << "t_i = " << t_i << "," << "i_ph = " << i_ph << "," << "ph = " << (*snodes1).p_gridt[year1t_id][i_ph]
						<< "," << "i_rent = " << i_rent << "," << "i_yi = " << i_yi << ","
						<< "y = " << (*snodes1).yi_gridt[year1t_id][i_yi] << ",";

					i_s = (*snodes1).i2s_map[i_ph][i_rent][i_yi]; 
					if (i_s == i_s_mid) {
						v1_file << "  i_s_mid" ;
					}
					v1_file << " , , ";

				}
			}
		}
	}
	v1_file << endl;

	// print policy / control variable headers
	for ( t_i = 0; t_i < t_n; t_i++) {                                               
		for ( i_ph = 0; i_ph < n_ph; i_ph++) {
			for (i_rent = 0; i_rent < n_rent; i_rent++) {
				for (i_yi = 0; i_yi < n_yi; i_yi++) {
					// for i_yi = 0, ..., 
					v1_file << "W,C,B,X,CSFp,CSFn,t_i2,V,";
				}
			}
		}
	}
	v1_file << endl;

	/* This part of file cycles through wealth and states and prints policy and value at that state
	   Main variable of interest is how policies change with wealth
	   If there are many states, it is not clear which one to prioritize
	   Suggestion: compute i_s_mid which represents median home price, rent, and income
	   then middle state is: i_s_mid
	*/

	//(*snodes1).s_ph_midry

	for ( w_i = 0; w_i < w_n; w_i++) {
		for ( t_i = 0; t_i < t_n; t_i++) {
			for (i_ph = 0; i_ph < n_ph; i_ph++) {
				for (i_rent = 0; i_rent < n_rent; i_rent++) {
					for (i_yi = 0; i_yi < n_yi; i_yi++) {
						// given i_ph, i_rent, i_yi, get state
						i_s = (*snodes1).i2s_map[i_ph][i_rent][i_yi];
						// for y_i = 0,..,n_y
						// i_s = (*snodes1).i2s[i_ph][y_i][rent_i]; check ordering
						//i_s = (*snodes1).s_ph_midry[i_ph];  // pass in price dimension of interest, receive state given median rent, income
						v1_file << (*vfnt).w_grid[w_i] << ","
							<< (*vfnt).x1_grid[t_i][i_s][w_i] << "," << (*vfnt).x2_grid[t_i][i_s][w_i] << "," << (*vfnt).x3_grid[t_i][i_s][w_i] << ","
							<< (*vfnt).x4_grid[t_i][i_s][w_i] << "," << (*vfnt).x5_grid[t_i][i_s][w_i] << "," << (*vfnt).xt_grid[t_i][i_s][w_i] << ","
							<< (*vfnt).vw3_grid[t_i][i_s][w_i] << ",";
					}
				}
			}
		}
		v1_file << endl;
	}

	// TODO: double check file works; only have 1 rent state now so added complexity shouldnt be too much
	v1_file.close();                                                                 // finished writing file


	cout << "store_data.cpp: begin writing flat file" << endl;

	string file_name_flat = "vfn_results/" + file_name_def2 + "_flat" + ".csv";               // flat file written and appended throughout project

	if (year1t_id == T_max) {
		ofstream v1_file_flat;
		v1_file_flat.open(file_name_flat, ios::out | ios::trunc);
		//v1_file_flat << "year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, w_i, W, C, B, X, CSFp, CSFn, t_i2, V";
		v1_file_flat << "year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, i_rent, i_yi, w_i, W, C, B, X, CSFp, CSFn, t_i2, V";
		v1_file_flat << endl;
		v1_file_flat.close();
	}

	ofstream v1_file_flat;
	v1_file_flat.open(file_name_flat, ios::out | ios::app);               // outstream, truncate

	for (w_i = 0; w_i < w_n; w_i++) {
		for (t_i = 0; t_i < t_n; t_i++) {
			for (i_ph = 0; i_ph < n_ph; i_ph++) {
				for (i_rent = 0; i_rent < n_rent; i_rent++) {
					for (i_yi = 0; i_yi < n_yi; i_yi++) {
						i_s = (*snodes1).s_ph_midry[i_ph];  // pass in price dimension of interest, receive state given median rent, income
						v1_file_flat << year1_id << "," << year1t_id << "," << rhoi << "," << gammai << "," << csfLev << "," << w_n << ","; // start working here; looking good so far
						//v1_file_flat << t_i << "," << ph_i << "," << w_i << ",";
						v1_file_flat << t_i << "," << i_ph << "," << i_rent << "," << i_yi << "," << w_i << ",";
						v1_file_flat << (*vfnt).w_grid[w_i] << ","
							<< (*vfnt).x1_grid[t_i][i_s][w_i] << "," << (*vfnt).x2_grid[t_i][i_s][w_i] << "," << (*vfnt).x3_grid[t_i][i_s][w_i] << ","
							<< (*vfnt).x4_grid[t_i][i_s][w_i] << "," << (*vfnt).x5_grid[t_i][i_s][w_i] << "," << (*vfnt).xt_grid[t_i][i_s][w_i] << ","
							<< (*vfnt).vw3_grid[t_i][i_s][w_i] << ",";

						v1_file_flat << endl;
					}
				}
			}
		}
	}
	v1_file_flat.close();
}



/*
for(int w_i = 0; w_i < w_n; w_i ++ ){
for(int t_i = 0; t_i < t_n; t_i ++){
for(int ph_i = 0; ph_i < ph_n; ph_i ++){
v1_file << (*vfnt).w_grid[w_i] << ","
<< (*vfnt).x1_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x2_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x3_grid[t_i][y_i_def][ph_i][w_i] << ","
<< (*vfnt).x4_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x5_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).xt_grid[t_i][y_i_def][ph_i][w_i] << ","
<< (*vfnt).vw3_grid[t_i][y_i_def][ph_i][w_i] << ",";
}
}
v1_file<< endl;
}
*/


/*

for (int w_i = 0; w_i < w_n; w_i++){
for (int t_i = 0; t_i < t_n; t_i++){
for (int ph_i = 0; ph_i < ph_n; ph_i++){
v1_file_flat << year1_id << "," << year1t_id << "," << rhoi << "," << gammai << "," << csfLev << "," << w_n << ","; // start working here; looking good so far
v1_file_flat << t_i << "," << ph_i << "," << w_i << "," ;
v1_file_flat << (*vfnt).w_grid[w_i] << ","
<< (*vfnt).x1_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x2_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x3_grid[t_i][y_i_def][ph_i][w_i] << ","
<< (*vfnt).x4_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).x5_grid[t_i][y_i_def][ph_i][w_i] << "," << (*vfnt).xt_grid[t_i][y_i_def][ph_i][w_i] << ","
<< (*vfnt).vw3_grid[t_i][y_i_def][ph_i][w_i] << ",";

v1_file_flat << endl;
}
}

}

*/

/*
want:
ID	year0	rho	gamma	csfLev	wn	vt	t_i	ph_i	ph	W	C	B	X	CSFp	CSFn	t_i2	V
ID: autonumber;
rho: const
gamma: const
csfLev: const
wn: const
year0_id:  go to main and look at t; ex: t = 0,1,2,...,
; stays same for most of program
vt_id:

*/