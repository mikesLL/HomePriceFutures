// class snodes.h
// contains nodes used to discretize home prices, rents, and income 
// as well as transition matrix                                                
// Copyright A. Michael Sharifi, 2016

#ifndef SNODES_H
#define SNODES_H

class snodes {
	int i_s, i_rent, i_ph, i_yi;
	
public:

	int city_id;
	double hu_ten[t_n]; //= hu_ten_store[N_cities][t_n];
	double ten_w[t_n]; 
	double rent_adj;
	double csfLevSn;

	int age0;
	int T_max;
	int i_s_mid;
	int t_hor;
	int s_ph_midry[n_ph];                               // states where home prices index from low to high, but rent and yi are always median

	vector<vector<double>> p_gridt, rent_gridt, yi_gridt, yi_gridt_btax;

	vector<vector<vector<double>>> gammat;              // state transition matrix for each time period

	vector<vector<vector<int>>> i2s_map;                // maps individual dimensions to state
	vector<int> s2i_ph;                                 // maps state to current price
	vector<int> s2i_rent;                               // maps state to current rent
	vector<int> s2i_yi;                                 // maps state to current income

	snodes(int age0_in, int T_max_in, int city_id_in );
	
	void adj_tax();
};

#endif
