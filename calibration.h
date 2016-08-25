// calibration.h
// Includes parameters used throughout the program
// Copyright A. Michael Sharifi, 2016


/*
San Diego, San Francisco, Los Angeles, Boston, Chicago, Denver,  Miami, New York
city_id: 0,1,2,3,4,5,6,7,8
*/

const int city_begin = 0;
const int city_end = 3;
//const int city_id = 0;                    // city_id; ={0: San Diego, 1: San Francisco, 2: Los Angeles)
const int t_begin = 5;                   // begin in year 4 from .csv
const int t_end = 12;                    // = 11 to cycle through all time periods;


const double age_begin = 30.0; //45.0;  // 60
//const double age_begin = 60.0;

                                           // sd_read_in, sfr_read_in, and lax_read_in, all begin at year 4 and end at year 11
const int T_max = 5; //35; //20; //12;    // maximum number of years included in _readin.csv; 1,...,11 corresponds to years 2003,...,2014
const int w_n = 200;                      // Grid points in wealth; set = 10 for fast computation, = 2000 for precision

const double csfLev = 1.0 * ( 1.0 / 0.055 );       // Case-Shiller Index Future margin-implied leverage; (notional value contract)/(median home price)*(1/margin)
const int csfLevi = int(floor(csfLev));   // Floor for identification

const int t_n = 3;                        // possible tenure states
//const int ph_n = 9; // 5;                       // possible home price states
//const int y_n = 2;                        // possible labor income states
const int pref = 1;                       // set pref = 0 for Cobb-Douglas, = 1 for CES
const int N_control = 6;
const int N_cities = 8; // 3;                   // number of cities available

const int n_ph = 9;      // possible home price states
//const int n_ph = 5; //9;
const int n_rent = 1; // 3;
const int n_yi =  3;  // nodes in each dimension

const int n_s = n_ph * n_rent * n_yi;  // number of states

// Labor income related parameters
const double y_tax = 0.0; // 0.3;  // taxation is handled in snodes.cpp
const double y_atax = 1.0 - y_tax;
const double y_replace = 0.9388;     // From Cocco, Gomes, Maenhout (2005)
const double y_inc_city[] = { .629, .8845, .5445 };                                  // Household income by city; in $100k (Source: FFIEC)

const double w_max = 40.0; //16.05;               // maximum wealth (on grid) (100's thousands)             
const double w_min = -2.0; // 0.05;                // minimum wealth (on grid) (100's thousands) 

//const double rho = 3.0;
const double rho = 1.8;                  // Power: curvature parameter; governs risk-aversion  also = 1.0, 2.0, 4.0
const int rhoi = int(floor(rho));        // Floor for identification

const double beta = .98;   // alt: =.95               // Beta: time preferences

const double phi = 0.06;                  // moving / transaction costs in event of home sale; also =.15

// IF CES Preferences:
const double alpha_ces = -6.485; // -6.7; // .75;                // Low substitutability between C and H
const double gammad = 1.01;
const int gammai = int(floor(gammad));
const double alpha = .3;

// Housing-service related parameters
// available for San Diego, San Francisco, Los Angeles
//const double hu_med[N_cities] = { 1.5, 1.62, 1.5 }; // median square footage by city
const double hu_med[N_cities] = { 1.5, 1.62, 1.5, 1.8, 1.8, 1.59, 1.5, 1.8 }; // median square footage by city

// Home Sizes (square footage, thousands);  33-66-quintiles from AHS Data
// Assume the small house can also be rented

//const double hu_ten_store[N_cities][t_n] = { { 1.217, 1.217, 2.352 },       
//                                           { 1.206, 1.206, 2.455 },
//	       							         { 1.116, 1.116, 2.269 } };

const double hu_ten_store[N_cities][t_n] = 
{ { 0.75*1.217, 1.217, 2.352 },
{ 0.75*1.374, 1.374, 2.585 },
{ 0.75*1.206, 1.206, 2.455 },
{ 0.75*1.492, 1.492, 2.826 },
{ 0.75*1.409, 1.409, 2.781 },
{ 0.75*1.116, 1.116, 2.379 },
{ 0.75*1.18, 1.18, 2.459  },
{ 0.75*1.303, 1.303, 2.764 } };

// Housing service flow

//const double hu_ten[t_n] = { hu_ten_store[city_id][0] , hu_ten_store[city_id][1], hu_ten_store[city_id][2] };

const double hu_ten_def =  .5;

// Housing wealth weight
//const double ten_w[t_n] = { 0.0, hu_ten[1] / hu_med[city_id], hu_ten[2] / hu_med[city_id] };   //

// Housing small rental adjustment:
// this adjusts owner's equivalent rent downward for small home size

//const double rent_adj = hu_ten[0] / hu_med[city_id];  

const double h_low = 0.8;

// If Cobb-Douglas Preferences:
const double alpha_cd = 0.6;                         // Prefence weight for C: Non-durable consumption
const double calpha_cd = 0.4;                        // =1.0 - alpha_sd; Preference weight for H: Housing services

const double delta = .05;                            // Minimum down payment
const double rb = 1.0204;                             // Gross return on bonds / mortgage rate; sometimes = 1.04

//Equity approximation (Gaussian-Hermite quadrature, 2 node approximation)
const double x_mu = 0.06; // Cocco-Gomes Maenhout
const double x_std = 0.157; //0.2415; // 
const int retxn = 2;
const double retxv[] = { x_mu - x_std, x_mu + x_std }; //{ -0.10, 0.25 };
const double retxp[] = { 0.5, 0.5 };

//const double retxv[] = { -0.10, 0.25 }; //{ -0.10, 0.25 };
//const double retxp[] = { 0.5, 0.5 };

//const int retxn = 5;
//const double retxv[] = { -.6449, -.26, .12, .5, .8939 };
//const double retxp[] = { 0.0094, 0.2651, 0.5990, 0.1247, 0.0018 };

// consider adding work here on returns:
const double p_move = 0.0;                            // Probability receiving an exogenous moving shock

const double b_min_const = -20.0;                    // TODO: double-check this
const double b_motive = 1.0;                         // Strength of bequest motive

const double c_fs = .05;                             // Minimum baseline consumption (Gov Assistance)
const double coh_fs = .05;                           // Cash on hand (Gov Asssistance: Non-durable Consumption + Housing)

const double min_dpmt = .05;                         // minimum down payment
const double max_ltv = .95;                          // max loan to value

// mortgage risk criteria
const double mort_spread = .02;                       // mortgage spread above risk-free
const double pmi_dpmt = .20;                          // if down payment below this amount, pay higher APR
const double pmi_prem = 0.01;                         // 02;      // pmi premium
const double credit_prem = .18;
const double b_min_unsec = -0.4;

