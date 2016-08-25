// calibration.h
// Includes parameters used throughout the program
// Copyright A. Michael Sharifi, 2016


/*
San Diego, San Francisco, Los Angeles, Boston, Chicago, Denver,  Miami, New York
city_id: 0,1,2,3,4,5,6,7,8
*/
const int city_begin = 0;
const int city_end = 3;
const int t_begin = 5;                   // begin in year 4 from .csv
const int t_end = 12;                    // = 11 to cycle through all time periods;


const double age_begin = 30.0; //45.0;  // 60
                                           // sd_read_in, sfr_read_in, and lax_read_in, all begin at year 4 and end at year 11
const int T_max = 5; //35; //20; //12;    // maximum number of years included in _readin.csv; 1,...,11 corresponds to years 2003,...,2014
const int w_n = 200;                      // Grid points in wealth; set = 10 for fast computation, = 2000 for precision

const double csfLev = 1.0 * ( 1.0 / 0.055 );       // Case-Shiller Index Future margin-implied leverage; (notional value contract)/(median home price)*(1/margin)
const int csfLevi = int(floor(csfLev));   // Floor for identification

const int t_n = 3;                        // possible tenure states
const int pref = 1;                       // set pref = 0 for Cobb-Douglas, = 1 for CES
const int N_control = 6;
const int N_cities = 8; // 3;                   // number of cities available

const int n_ph = 9;      // possible home price states
const int n_rent = 1; // 3;  possible rent states
const int n_yi =  3;  // labor income states

const int n_s = n_ph * n_rent * n_yi;  // number of states

// Labor income related parameters
const double y_tax = 0.0;                 // 0.3;  // taxation is handled in snodes.cpp
const double y_atax = 1.0 - y_tax;
const double y_replace = 0.9388;            // From Cocco, Gomes, Maenhout (2005)
const double w_max = 40.0; //16.05;         // maximum wealth (on grid) (100's thousands)             
const double w_min = -2.0; // 0.05;          // minimum wealth (on grid) (100's thousands) 

//const double rho = 3.0;
const double rho = 1.8;                  // Power: curvature parameter; governs risk-aversion  also = 1.0, 2.0, 4.0
const int rhoi = int(floor(rho));        // Floor for identification

const double beta = .98;   // alt: =.95               // Beta: time preferences
const double phi = 0.06;                  // moving / transaction costs in event of home sale; also =.15


// Housing-service related parameters
// median square footage by city
const double hu_med[N_cities] = { 1.5, 1.62, 1.5, 1.8, 1.8, 1.59, 1.5, 1.8 }; 

// Home Sizes (square footage, thousands);  33-66-quintiles from AHS (2005) Data;
// Assume the small house can also be rented
const double hu_ten_store[N_cities][t_n] = 
{ { 1.217, 1.5, 2.352 },
{ 1.374, 1.62, 2.585 },
{ 1.206, 1.5, 2.455 },
{ 1.492, 1.8, 2.826 },
{ 1.409, 1.8, 2.781 },
{ 1.116, 1.59, 2.379 },
{ 1.18, 1.5, 2.459  },
{ 1.303, 1.8, 2.764 } };

const double hu_ten_def =  .5;  // square footage in default case

// If Cobb-Douglas Preferences:
const double alpha_cd = 0.6;                         // Prefence weight for C: Non-durable consumption
const double calpha_cd = 0.4;                        // =1.0 - alpha_sd; Preference weight for H: Housing services

const double delta = .05;                            // Minimum down payment
const double rb = 1.0204;                             // Gross return on bonds / mortgage rate; sometimes = 1.04

// IF CES Preferences:
const double alpha_ces = -6.485; // -6.7; // .75;                // Low substitutability between C and H
const double gammad = 1.01;
const int gammai = int(floor(gammad));

//Equity approximation (Gaussian-Hermite quadrature, 2 node approximation)
const double x_mu = 0.06; // Cocco-Gomes Maenhout
const double x_std = 0.157; 
const int retxn = 2;
const double retxv[] = { x_mu - x_std, x_mu + x_std }; //{ -0.10, 0.25 };
const double retxp[] = { 0.5, 0.5 };


const double p_move = 0.0;                            // Probability receiving an exogenous moving shock

const double b_min_const = -20.0;                    
const double b_motive = 1.0;                         // Strength of bequest motive

const double c_fs = .05;                             // Minimum baseline consumption (Gov Assistance)
const double coh_fs = .05;                           // Cash on hand (Gov Asssistance: Non-durable Consumption + Housing)

const double min_dpmt = .05;                         // minimum down payment
const double max_ltv = .95;                          // max loan to value

// mortgage risk criteria
const double mort_spread = .02;                       // mortgage spread above risk-free rate
const double pmi_dpmt = .20;                          // if down payment below this amount, add to mortgage spread
const double pmi_prem = 0.01;                         // pmi premium
const double credit_prem = .18;                       // unsecured credit limit
const double b_min_unsec = -0.4;                      // unsecured borrowing limit

