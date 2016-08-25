# HomePriceFutures
This project models a household's life-cycle portfolio optimization 
problem with the ability to take home price futures positions and calibrates 
the solutions to San Diego, San Francisco, and Los Angeles using data from 
2006-2013.

## Data
sdg_read_in.csv, sfr_read_in.csv, and lax_read_in.csv contain 
housing-related data for San Diego, San Francisco, and Los Angeles from 
including current rent, price, lagged home price appreciation, and futures 
price.

## Code

main.cpp is the main file in this project. main.cpp calls load_csv.cpp to load housing-related data from the .csv files. For each year, main.cpp calls load_pricepath.cpp to simulate a home price path for the next T years where T is the investor horizon. For each horizon, main.cpp then uses gen_VP.cpp to numerically solve the model using value function iteration and store_data.cpp to store the value functions and policy functions as a .csv file. main.cpp is run for one city at a time but can be run through starting years 2006,...,2013 in parallel using OpenMP framework.

calibration.h includes many of the parameters and settings for the program including city, preferences (CES vs Cobb-Douglas), and futures leverage.

The main results of interest come from running the program with and without access to home price index futures to determine the welfare benefit of hedging housing risk and modeling optimal tenure decision by year

## Compile Instructions
On Linux:
g++ *.cpp -std=c++11 -fopenmp -lpthread


## Directories
eval_perf: files for evaluating performance including evaluating welfare benefits and viewing price path

figures: stores figures generated in eval_perf for export to Latex

first_results: value function and policy data from first period only

price_results: sample simulated price paths and discretizations

vfn_results: value function and policy data from all periods

## Runtime
With 4000 gridpoints in wealth, 2 gridpoints in equity returns, 2 gridpoints in income, and 9 gridpoints in home prices, program takes around 24 hours on 1 core given a single city and starting year. 


## Acknowledgements
I would like to thank my advisors Marjorie Flavin and Johannes Wieland for their time and dedication. I hope this project meets their expectataions. This project also benefited greatly from discussions with my committee members Chris Parsons, Alexis Toda, and Ross Valkanov. My peers in the UCSD Macroeconomics, Finance, and Urban/Real Estate Groups also provided valuable discusssion and advice especially Asad Dossani, Emilien Gouin-Bonenfant, Xihan Xie, Kilian Heilmann, Irina Zhecheva, and David Stowitts.