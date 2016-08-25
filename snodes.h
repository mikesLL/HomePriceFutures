// class snodes.h
// contains nodes used to discretize home prices, rents, and income 
// as well as transition matrix                                                
// Copyright A. Michael Sharifi, 2016

#ifndef SNODES_H
#define SNODES_H

class snodes {
	int i_s, i_rent, i_ph, i_yi;
	
public:
	// adding these here
	int city_id;
	double hu_ten[t_n]; //= hu_ten_store[N_cities][t_n];
	double ten_w[t_n]; 
	double rent_adj;

	int i_s_mid;
	int t_hor;
	int s_ph_midry[n_ph];                               // states where home prices index from low to high, but rent and yi are always median
	double p_gridt[T_max + 1][n_ph];                    // p_grid arguments: time horizon, ph_n home price grid points
	double rent_gridt[T_max + 1][n_rent];               // rent_grid arguments: time horizon, ph_n home price grid points
	double yi_gridt[T_max + 1][n_yi];                 // yi_grid arguments: time horizon, ph_n home price grid points

	vector<vector<vector<double>>> gammat;              // state transition matrix for each time period

	vector<vector<vector<int>>> i2s_map;                // maps individual dimensions to state
	vector<int> s2i_ph;                                 // maps state to current price
	vector<int> s2i_rent;                               // maps state to current rent
	vector<int> s2i_yi;                                 // maps state to current income

	snodes();

	void adj_tax();
};

#endif
