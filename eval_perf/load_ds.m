%{
load_ds.m
This script loads all *.csv files in first_results and combines into a
single file. Scrips also retrieves file.
Copyright A. Michael Sharifi, 2016
%}

%%
function ds = load_ds( reload_data, city_str )

% path where results are stored
dir_first = '../first_results/';
dir_price = '../price_results/'; 

addpath(dir_first);
addpath(dir_price);

if ( reload_data > 0 )
    
    fileName=ls([dir_first,city_str,'*.csv']);         % load in all .csv files
    %fileName=ls([dir_first,'*.csv']);         % load in all .csv files
    N_files = size(fileName, 1);
    N_cols = 21;                               % each .csv includes 18 columns
    
    % structure holds first-year results from all .csv files
    ds_tmp{N_files} = zeros(1,N_cols);
    for id = 1:N_files
        fprintf('store file %d of %d \n', id, N_files);
        ds_tmp{id} = xlsread([dir_first, fileName(id,:)]);
    end
    
    ds = vertcat(ds_tmp{:});                   % data_struct stores data from initial horizon in all years
    save('ds_save.mat');
else
    load ds_save.mat;
end

end


