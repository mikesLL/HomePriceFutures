// calibration.h
// Includes parameters used throughout the program
// Copyright A. Michael Sharifi, 2016
// 
// cities: San Diego, San Francisco, Los Angeles, Boston, Chicago, Denver,  Miami, New York
// city_id: 0,1,2,3,4,5,6,7,8
// *.read_in files begin at year 5 (2007) and end at year 11 (2013)

const int city_begin = 2;
const int city_end = 2;

const int t_begin = 5;                   // begin in year 5 from .csv
const int t_end = 11;                    // = 11 to cycle through all time periods;

const int param_id = 0;  // set = 0 to define parameters here manually; set = 1, 2, 3, 4 for presets and load in main

const int age_begin_store[] = { 60, 45, 30 }; // manual age settings here
const int n_age_store[] = { 1, 1, 1, 2, 2 };
const int n_age = n_age_store[param_id]; 

const double csfLevStore[] = {1.0/0.055, 1.0/0.055, 0.0, 1.0/0.055, 0.0}; // manually set futures leverage
const double csfLev = csfLevStore[param_id];
const int w_n = 400; // Grid points in wealth; set = 200 for fast computation, = 2000 for precision

const int age_max = 65;                  // age at which household retires / annuitizes wealth  

//const double margin_store[] = { 0.0, 0.0, 0.02524, 0.032408, 0.0, 0.019866, 0.0, }; 
const double csfLev_pidxw[] = {1.0, 1.0, 0.0, 1.0, 0.0 };  // change in future notl to index weight

// csf margin requirement by city
//const double csfmarg_store[] = { 0.036444571, 0.039967956, 0.032995124, 0.033871271,
//	0.025347143, 0.023869709, 0.037148076, 0.030927638 };
const double csfmarg_store[] = { 0.036444571, 0.039967956, 0.032995124, 0.033871271,
	0.025347143, 0.023869709, 0.037148076, 0.030927638 };


//const double csfLev = 2.5 / csfLev_store[city_id];
//const double csfLev = 1.0 * ( 1.0 / 0.055 );       // Case-Shiller Index Future margin-implied leverage; (notional value contract)/(median home price)*(1/margin)
const int csfLevi = int(floor(csfLev));   // Floor for identification

const int t_n = 5; // 4;                        // possible tenure states
const int pref = 0;                       // set pref = 0 for Cobb-Douglas, = 1 for CES
const int N_control = 6;
const int N_cities = 8;                    // number of cities

const int n_ph = 9;      // possible home price states
const int n_rent = 1; // 3;  possible rent states
const int n_yi = 3; // 3;  // labor income states

const int n_s = n_ph * n_rent * n_yi;  // number of states

// Labor income related parameters
const double maint_mult = 0.98;
const double y_tax = 0.0;                 // 0.3;  // taxation is handled in snodes.cpp
const double y_atax = 1.0 - y_tax;
const double y_replace = 0.9388;            // From Cocco, Gomes, Maenhout (2005)

const double w_max = 40.0; //16.05;         // maximum wealth (on grid) (100's thousands)             
const double w_min = -2.0; // 0.0; // 0.05;          // minimum wealth (on grid) (100's thousands) 

const int w_i_zero = (int)ceil(-w_min * double(w_n) / (w_max - w_min));


//const double w_max = 40.0; //16.05;         // maximum wealth (on grid) (100's thousands)             
//const double w_min = -2.0; // 0.05;          // minimum wealth (on grid) (100's thousands) 

const double rho = 5.0;
//const double rho = 1.8;                  // Power: curvature parameter; governs risk-aversion  also = 1.0, 2.0, 4.0
const int rhoi = int(floor(rho));        // Floor for identification

const double beta = .97;   // alt: =.95               // Beta: time preferences
const double phi = 0.08;                  // moving / transaction costs in event of home sale; also =.15
const double phi_sell = 0.08; //0.10;
const double phi_buy = 0.00;  //0.02;


// Housing-service related parameters
// median square footage by city
const double hu_med[N_cities] = { 1.6, 1.7, 1.5, 1.9, 1.8, 1.6, 1.5, 1.83 };   // san fran updated

// Home Sizes (square footage, thousands);  33-66-quintiles from AHS (2005) Data;
// Assume the small house can also be rented

/*
const double hu_ten_store[N_cities][t_n] =
{ { 1.217, 1.217, 1.5, 1.9, 2.5, 2.8, 3.1, 3.6 },
{ 1.217, 1.217, 1.5, 1.9, 2.5, 2.8, 3.1, 3.6 },
{ 1.374, 1.374, 1.62,  1.0*2.585 },
{ 1.206, 1.206, 1.5,  1.0*2.455 },
{ 1.492, 1.492, 1.8,  1.0*2.826 },
{ 1.409, 1.409, 1.8,  1.0*2.781 },
{ 1.116, 1.116, 1.59,  1.0*2.379 },
{ 1.18, 1.18, 1.5,  1.0*2.459 },
{ 1.303, 1.303, 1.8,  1.0*2.764 } };
*/

/*
const double hu_ten_store[N_cities][t_n] =
{ { 1.35, 1.35, 1.7, 2.18, 3.0 }, 
{ 1.35, 1.35, 1.7, 2.18, 3.0 },    // sf updated
{ 1.2, 1.2, 1.6, 2.0, 3.0 },       // lax updated
{ 1.2, 1.2, 1.6, 2.0, 2.4 },
{ 1.2, 1.2, 1.6, 2.0, 2.4 },
{ 1.2, 1.2, 1.6, 2.2, 2.755 },    // denver updated
{ 1.2, 1.2, 1.6, 2.0, 2.4 },
{ 1.2, 1.2, 1.6, 2.0, 2.4 } };
*/

const double hu_ten_store[N_cities][t_n] = {
{1.25,	1.25,	1.60,	2.00,	2.70},
{1.35,	1.35,	1.70,	2.18,	3.00},
{1.20,	1.20,	1.50,	1.90,	2.95},
{1.50,	1.50,	1.90,	2.50,	3.945},
{1.40,	1.40,	1.80,	2.40,	3.20},
{1.20,	1.20,	1.60,	2.20,	2.755},
{1.20,	1.20,	1.50,	2.00,	3.00},
{1.40,	1.40,	1.83,	2.50,	3.60}
};


/*
const double hu_ten_store[N_cities][t_n] =
{ { 1.217, 1.217, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0 },
{ 1.374, 1.374, 1.62, 1.9, 2.2, 2.5, 2.8, 3.1 },
{ 1.206, 1.206, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, },
{ 1.492, 1.492, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3}, 
{ 1.409, 1.409, 1.8,  2.2, 2.6, 3.0, 3.4, 3.8 },
{ 1.116, 1.116, 1.59, 2.0, 2.4, 2.8, 3.2, 3.6 },
{ 1.18, 1.18, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0},
{ 1.303, 1.303, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8} };
*/


/*
const double hu_ten_store[N_cities][t_n] =
{ {1.217, 1.217, 1.5, 1.0*2.352  },
{1.374, 1.374, 1.62,  1.0*2.585 },
{ 1.206, 1.206, 1.5,  1.0*2.455 },
{1.492, 1.492, 1.8,  1.0*2.826 },
{1.409, 1.409, 1.8,  1.0*2.781 },
{1.116, 1.116, 1.59,  1.0*2.379 },
{1.18, 1.18, 1.5,  1.0*2.459 },
{1.303, 1.303, 1.8,  1.0*2.764 } };
*/

const double hu_ten_def =  .5;  // square footage in default case

// If Cobb-Douglas Preferences:
const double alpha_cd = 0.7;                         // Prefence weight for C: Non-durable consumption
const double calpha_cd = 0.3;                        // =1.0 - alpha_sd; Preference weight for H: Housing services

const double rb = 1.0304;                             // Gross return on bonds / mortgage rate; sometimes = 1.04

// IF CES Preferences:
const double alpha_ces = -6.485; // -6.7; // .75;                // Low substitutability between C and H
const double gammad = 1.01;
const int gammai = int(floor(gammad));

//Equity approximation (Gaussian-Hermite quadrature, 2 node approximation
const double x_mu = 0.04 + (rb - 1.0); // Cocco, Gomes Maenhout (2005)
const double x_std = 0.157; 
const int retxn = 2;
const double retxv[] = { x_mu - x_std, x_mu + x_std }; //{ -0.10, 0.25 };
const double retxp[] = { 0.5, 0.5 };

const double p_move = 0.0;  // also: set = 0.0        // Probability receiving an exogenous moving shock

const double b_min_const = -20.0;                    
const double b_motive = 1.0;                         // Strength of bequest motive

const double c_fs = .06;                             // Minimum baseline consumption (Gov Assistance)
const double coh_fs = .05;                           // Cash on hand (Gov Asssistance: Non-durable Consumption + Housing)

// down-payment criteria
const double delta = .20;                            // Minimum down payment
const double min_dpmt =  .20;    //0.2;                     // minimum down payment
const double max_ltv = 0.80; // .95;  0.8;                        // max loan to value

// mortgage risk criteria
const double max_lti = 0.3;
const double mort_spread = .00;                       // mortgage spread above risk-free rate
const double pmi_dpmt = 0.00;                          // if down payment below this amount, add to mortgage spread
const double pmi_prem = 0.0; //0.01;                         // pmi premium
const double credit_prem = .15;                       // unsecured credit apr
const double b_min_unsec = -2.0; //-1.0;               // unsecured borrowing limit

const double num_small = -1.0e20;

// basis risk
//int n_csf_basis = 1;// 2;

// case: no basis risk
//const int n_csf_basis = 1; 
//const double csf_basis[] = { 0.0, 0.0 };
//const double pcsf_basis[] = { 1.0, 0.0 };

//int n_csf_basis = 2;
//double csf_basis[] = { -0.045, 0.045 };

// alt:

