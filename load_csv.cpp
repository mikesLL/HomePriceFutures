/*
load_csv.cpp
This program simulates a home price path and returns interpolating prices as well as a transition matrix
Inputs: *city_in refers to a housing data structure and string file_name refers to a .csv file which contains
metro-level housing data like median price, median rent, recently observed home price appreciation and others

reads from ../read_in

Copyright A. Michael Sharifi, 2016
*/

#include "headers.h"

void load_csv( void *city_in, string file_name ){

	hdata *city = (hdata *) city_in;
	
	ifstream read_in("read_in/"+file_name);            // open input file stream to .csv file
	string line;
	string line2;
	string parsed;
	int id, j, k, max_id = 0;
	int line_length, c1, c2, c3, c4, c5;
	const int max_len = 40, max_col = 5;

	double data_array[max_len][max_col]; // col of double data array stores info on:  year, rent, price, lag_ret, csf_1yr
   
	while(read_in.good()){                        
		getline( read_in, line );              // read in .csv line and find locs of each comma
		
		line_length = line.length();
		c1 = line.find(',');                   // c1: comma 1
		c2 = line.find(',', c1 + 1 );          // c2: comma 2
		c3 = line.find(',', c2 + 1 );          // c3: comma 3
		c4 = line.find(',', c3 + 1 );          // c4: comma 4
		c5 = line.find(',', c4 + 1 );          // c5: comma 5
       
		stringstream convert( line.substr(0, c1 ) );
		convert >> id;                         // store id of current line
	
		for(j = 0; j< max_len; j++){
			if (id == j){
				max_id = max(id, max_id);                                   // max_id is store the max number of obs from the .csv

				stringstream convert1( line.substr(c1+1,c2-1) );            // read year into array
				convert1 >> data_array[id][0];

				stringstream convert2( line.substr(c2+1,c3-1) );            // read rent into array
				convert2 >> data_array[id][1];

				stringstream convert3( line.substr(c3+1,c4-1) );            // read price into array
				convert3 >> data_array[id][2];

				stringstream convert4(line.substr(c4 + 1, c5 - 1));         // read lag_ret into array
				convert4 >> data_array[id][3];

				stringstream convert5( line.substr(c5+1, line_length-1) );  // read csf_1yr into array
				convert5 >> data_array[id][4];
			}
		}
	}

	cout << "id, year, rent, price, lagged_ret, csf_1yr"  << endl;
	
	for( j = 0; j <= max_id;  j++){                                   // move data from data_array to city data structure
	    (*city).year[j] = (double) data_array[j][0];
		(*city).rent[j] = (double) data_array[j][1];
		(*city).price[j] = (double) data_array[j][2];
		(*city).ret_lag[j] = (double) data_array[j][3];
		(*city).csf_1yr[j] = (double) data_array[j][4];               // order csf_1yr last
	
		cout << j << "," << (*city).year[j] << "," << (*city).rent[j] << ","
             << (*city).price[j] << "," << (*city).ret_lag[j] << "," << (*city).csf_1yr[j]  << endl;
	}

	cout<<"Finished load_csv "<< endl;

}
