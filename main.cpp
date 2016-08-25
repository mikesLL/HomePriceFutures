/*
This program numerically solves a household's (HH's) lifecycle portfolio optimization 
problem when the HH is able to take S&P Case-Shiller Index Futures positions. The HH's value fn 
and status in t = 11 (11 year horizon) is modeled as vf_11. Given vf_11, program solves for 
vf_10,..., vf_1 by value function iteration using gen_v1.

An integral part of the portfolio problem is estimating the future 
of home prices. I simulate the home price path in load_rand.

pstruct contains possible home prices in given number of time periods
gstruct contains a transition matrix for each time period.
sd_read_in.csv, sf_read_in.csv, and lax_read_in.csv include housing related data for 
that particular year (median price, median rent, lagged home price appreciation, 
Case-Shiller Inex Futures Price, etc...)
(not included on GitHub due to proprietary data)
Currently works for on data from 2007-2014
Cities: San Diego, San Francisco, Los Angeles, Boston, Chicago, 
Denver, Miami, New York 

Copyright A. Michael Sharifi, 2016
*/

#include "headers.h"


int main(){
	string city_init, city_filename;                                   // city name
	hdata city_data;                                     // housing data structure stores previous rents and lagged returns
	
	const string city_init_vec[] = {"sd", "sf", "lax", "bos", "chi", "den", "mia", "nym"  };
	const string city_filename_vec[] = { "sd_read_in.csv", "sf_read_in.csv", "lax_read_in.csv",
		                                 "bos_read_in.csv", "chi_read_in.csv", "den_read_in.csv",
		                                 "mia_read_in.csv", "nym_read_in.csv" };

	int city_id;
	for (city_id = city_begin; city_id <= city_end; city_id++) {

		city_filename = city_filename_vec[city_id];
		city_init = city_init_vec[city_id];
		cout << "city_init = " << city_init << endl;
		load_csv(&city_data, city_filename);

		int t;

#pragma omp parallel for
		for (t = t_begin; t <= t_end; t++) {

			int t_hor = T_max; // optimization problem horizon;
			double duration;
			double phr_in = city_data.rent[t];     // load in current median rent

			clock_t start = clock();
			snodes snodes1;                        // discretized states including home prices, rents, incomes

			snodes1.city_id = city_id; 
			snodes1.hu_ten[0] = hu_ten_store[city_id][0];
			snodes1.hu_ten[1] = hu_ten_store[city_id][1];
			snodes1.hu_ten[2] = hu_ten_store[city_id][2];

			snodes1.ten_w[0] = 0.0;
			snodes1.ten_w[1] = snodes1.hu_ten[1] / hu_med[city_id]; 
			snodes1.ten_w[2] = snodes1.hu_ten[2] / hu_med[city_id];

			snodes1.rent_adj = snodes1.hu_ten[0] / hu_med[city_id];

			snodes1.adj_tax();

			cout << "housing tenure: " << endl;
			cout << snodes1.hu_ten[0] << "..." << snodes1.hu_ten[1] << "..." << snodes1.hu_ten[2] << "..." << endl;

			cout << "housing wealth weight: " << endl;
			cout << snodes1.ten_w[0] << "..." << snodes1.ten_w[1] << "..." << snodes1.ten_w[2] << "..." << endl;

			cout << "rent adj: " << endl;
			cout << snodes1.rent_adj << endl;

			vfn vf_F;
			vfn vf_P;

			// load city_data and into ps1 and gs1; include current rent, current home price,
			// lagged home price appreciation, Case-Shiller Futures Price, current time
			// load_simpath store discretized approximation in snodes1 structure 
			load_simpath(&snodes1, city_data.rent[t], city_data.price[t], city_data.ret_lag[t], city_data.csf_1yr[t], t, city_init, city_id);

			cout << "main.cpp: begin data" << endl;
			vf_F.enter_data(&snodes1, phr_in, t, t_hor, city_data.csf_1yr[t], pref);

			cout << "main.cpp: set terminal" << endl;
			vf_F.set_terminal(phr_in);

			cout << "main.cpp: store_data" << endl;
			store_data(&snodes1, &vf_F, city_init + "yr", t, t_hor);

			for (t_hor = (T_max - 1); t_hor >= 0; t_hor--) {
				snodes1.t_hor = t_hor;
				cout << "main.cpp: vf_P: enter data" << endl;
				vf_P.enter_data(&snodes1, phr_in, t, t_hor, city_data.csf_1yr[t], pref);

				cout << " main.cpp: gen_VP" << endl;
				gen_VP(&snodes1, &vf_P, &vf_F);

				cout << "main.cpp: begin store_data" << endl;
				store_data(&snodes1, &vf_P, city_init + "yr", t, t_hor);

				vf_F = vf_P;
			}

			t_hor = 0;

			cout << "\n Value Function Iteration Completed \n";
			duration = (clock() - start) / (double)CLOCKS_PER_SEC;
			cout << "total run time: " << duration << '\n';
		}
	}

	cin.get();
	return 0;
}


