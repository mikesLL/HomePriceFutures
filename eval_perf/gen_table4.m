%{
gen_table4.m

This file plots optimal tenure / home size and futures positions as a function
of current state variables and access to home price futures

Copyright A. Michael Sharifi
%}

function [  ten_i2, ten_i2_CSF, csf, B, B_CSF, X, X_CSF ] = ...
    gen_table4(col, age0, t_i1, city_id, ds, t_begin, t_end, p_mid  )

%ds_W = sortrows( ds(:,[10,11]) );
%idx0 = find(ds_W(:,2) >= W1, 1, 'first');
%W_i1 = ds_W(idx0,1);  % W_i1: w_i (wealth index) associated with W1

w_n = ds(1,col.w_n);
ten_i2 = zeros(w_n,t_end);
ten_i2_CSF = zeros(w_n,t_end);
csf = zeros(w_n, t_end);

B = zeros(w_n,t_end);
B_CSF = zeros(w_n,t_end);
X = zeros(w_n,t_end);
X_CSF = zeros(w_n,t_end);

idx = ( ds(:,col.csfLev) <= 0.0 );
idx_CSF = ( ds(:,col.csfLev) > 0.0 );

%%
for t = t_begin:t_end
    idx1 = all( [ ds(:,col.age0) == age0, ...   % age condition
        ds(:,col.city_id) == city_id, ...       % city_id condition
        ds(:,col.year_id) == t, ...             % year condition
        ds(:,col.hor_id) == 0, ...              % year horizon = 0
        idx, ...
        ds(:,col.t_i) == t_i1, ...              % t_i = 0
        ds(:,col.i_yi) == 1, ...                % i_yi = 1
        ds(:,col.ph_i) == p_mid ], 2);          % p_i = 2 (actual price)
    
    idx2 = all( [ ds(:,col.age0) == age0, ...   % age condition
        ds(:,col.city_id) == city_id, ...     % city_id condition
        ds(:,col.year_id) == t, ...           % year condition
        ds(:,col.hor_id) == 0, ...           % year horizon = 0
        idx_CSF, ...
        ds(:,col.t_i) == t_i1, ...               % t_i = 0
        ds(:,col.i_yi) == 1, ...                 % i_yi = 1
        ds(:,col.ph_i) == p_mid ], 2);           % p_i = 2 (actual price)
    
   
    %ten_i2(:,t) = ds(idx1,17);
    %ten_i2_CSF(:,t) = ds(idx2,17);
    ds_tmp = ds(idx1,:);
    ten_i2(:,t) = ds_tmp(1:w_n,col.t_i2);
    B(:,t) = ds_tmp(1:w_n, col.B);
    X(:,t) = ds_tmp(1:w_n, col.X);
    
    ds_tmp = ds(idx2,:);
    ten_i2_CSF(:,t) = ds_tmp(1:w_n,col.t_i2);
    B_CSF(:,t) = ds_tmp(1:w_n, col.B);
    X_CSF(:,t) = ds_tmp(1:w_n, col.X);
    
    csf(1:w_n,t) = ds_tmp(1:w_n, col.csfP) - ds_tmp(1:w_n, col.csfN);
end

end

%%
% set year and pick a wealth level
% condition: year = 2004, wealth = 40000
% want: tenure status
% idx = find(

%{
city_id, year1_id, year1t_id, rho, gamma, csfLev, w_n, t_i, ph_i, w_i,
W, C, B, X, CSFp, CSFn, t_i2, V
%}
