%{
load_ds.m
This script loads all *.csv files in first_results and combines them into a
single file. The script then retrieves this file.

Copyright A. Michael Sharifi, 2016
%}
%%
function ds = load_ds( reload_data )

if ( reload_data > 0 )

fileName=ls('*.csv');                      % load in all .csv files
N_files = size(fileName, 1);
N_cols = 19;                               % each .csv includes 18 columns

% create a structure to hold first-year results from all .csv files 
ds_tmp{N_files} = zeros(1,N_cols);
for id = 1:N_files
    fprintf('store file %d of %d \n', id, N_files);
    ds_tmp{id} = xlsread(fileName(id,:));    
end
ds = vertcat(ds_tmp{:});                   % data_struct stores data from initial horizon in all years
save('ds_save.mat');
else
    load ds_save.mat;
end

end