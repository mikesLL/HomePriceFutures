%{
create_ds_flat.m
This script loads all *.csv files in first_results and combines them into a
single file. 

Copyright A. Michael Sharifi, 2016
%}

%%
function ds_flat = create_ds_flat
fileName=ls('../first_results/*.csv');                      % load in all .csv files
N_files = size(fileName, 1);
N_cols = 19;                               % each .csv includes 18 columns

% create a structure to hold first-year results from all .csv files 
ds_tmp{N_files} = zeros(1,N_cols);
for id = 1:N_files
    fprintf('store file %d of %d \n', id, N_files);
    ds_tmp{id} = xlsread(fileName(id,:));    
end

ds_flat = vertcat(ds_tmp{:});                   % data_struct stores data from initial horizon in all years

save('ds_flat_save.mat','ds_flat');


end