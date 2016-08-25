// hdata.h 
// Class stores information relevant to a particular metro housing market by year 
// including: year, rent, price, ret_lag, and csf_1yr 
// Is analagous to data in __read_in.csv 
// Copyright A. Michael Sharifi, 2016


#ifndef HDATA_H
#define HDATA_H

class hdata {

public:
	vector<double> year;           
	vector<double> rent; 
	vector<double> price; 
	vector<double> ret_lag; 
	vector<double> csf_1yr;
	int max_len = 20;
	
	hdata ( ){
		year.resize(max_len);
		rent.resize(max_len);
		price.resize(max_len);
		ret_lag.resize(max_len);
		csf_1yr.resize(max_len);
	}
};

#endif