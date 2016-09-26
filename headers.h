// headers.h
// Copyright A. Michael Sharifi, 2016 

#include <iostream>
#include <fstream>   
#include <random>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <numeric>


#include "calibration.h"
#include "vfn.h"
#include "ufnEV2.h"
#include "hdata.h"
#include <algorithm>
#include <omp.h>
#include "snodes.h"


using namespace std;

typedef struct{
    double v_opt;
	double v_opt_lagmod;
    vector<double> x_opt;
	int valid_flag;
} gen_res;

void load_simpath(void *snodes_in, double rent_in, double ph0_in,

	double ret0_in, double csf_1yr_in, int t_id, string city_init_in, int city_id_in, int age_begin_in);


void gen_VP(void *snodes_in, void *my_func_data1, void *my_func_data2 );

gen_res gen_VPw(void *snodes_in, void *vf1_in, void *vf2_in, double coh_in, vector <double> x0_in, double b_min_in, 
	double beg_equity_in, double mpmt_in );

vector<double> gen_x0(double coh, double b_min, void *vf1, void *vf2, void * evalEV21_in, vector<double> x0_in );

void store_data(void *snodes_in, void *vfn_in, string file_name, int year1_id_in, int year1_t_id_in, int T_max);

void load_csv( void *city_in, string file_name );

double ufn( double c_in, double hu_in, int pref_in);
	


